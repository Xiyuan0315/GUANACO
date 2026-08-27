"""Shared helpers for reading obs annotations consistently across plots."""

import hashlib
import weakref

import numpy as np
import pandas as pd

from guanaco.data.loader import obs_col


SELECTION_GROUP = "Selected / Others"
SELECTION_GROUP_LABEL = "Selection"
SELECTED_LABEL = "Selected"
OTHERS_LABEL = "Others"
SELECTION_LABELS = [SELECTED_LABEL, OTHERS_LABEL]

# id(dataset) -> (weakref to the dataset, {(col, dropna): sorted categories}).
# Keyed by id() because AnnData/MuData are unhashable, and guarded by a weakref
# so an entry is dropped the moment its dataset is collected -- this prevents a
# later object that happens to reuse the same id() from reading a stale result.
_category_cache: dict = {}


def obs_values(adata, col, overrides=None) -> pd.Series:
    """Read one obs column, optionally from a session-local synthetic column."""
    if overrides and col in overrides:
        values = overrides[col]
        if isinstance(values, pd.Series):
            return values.reindex(adata.obs_names)
        return pd.Series(values, index=adata.obs_names, name=col)
    return obs_col(adata.obs, col)


def sorted_categories(
    adata,
    col,
    *,
    dropna: bool = True,
    overrides=None,
) -> list:
    """Return an obs column's categories in one canonical order, shared by all plots.

    Every plot colours and orders categories by this list's index, so they must
    all agree -- otherwise the same label gets a different colour in the violin,
    heatmap, dotplot and stacked-bar views. Centralising the
    ``sorted(obs[col].unique())`` each callback used to hand-roll guarantees that,
    and standardises the two details that differed between them:

    - **NaN is dropped** (``dropna=True``): a missing value isn't a labelled group.
    - **Sorting matches a plain ``sorted``** -- lexical for the usual string
      categories, numeric for genuinely numeric dtypes -- but falls back to string
      order for mixed/uncomparable values so it never raises the ``TypeError`` a
      bare ``sorted`` would (this was the only reason violin used ``key=str``).

    Results for a stable dataset are memoised (the full ``adata`` is reused across
    callbacks); transient filtered *views* are a fresh object each call, so they
    miss the cache and recompute -- exactly as the inline code did before.
    """
    if overrides and col in overrides:
        values = obs_values(adata, col, overrides)
        uniques = values.dropna().unique() if dropna else values.unique()
        try:
            return sorted(uniques)
        except TypeError:
            return sorted(uniques, key=str)

    key = id(adata)
    entry = _category_cache.get(key)
    if entry is not None and entry[0]() is adata:
        hit = entry[1].get((col, dropna))
        if hit is not None:
            return list(hit)
        per_col = entry[1]
    else:
        per_col = {}
        try:
            ref = weakref.ref(adata, lambda _r, k=key: _category_cache.pop(k, None))
            _category_cache[key] = (ref, per_col)
        except TypeError:
            per_col = None  # not weak-referenceable: compute without caching

    values = obs_col(adata.obs, col)
    uniques = values.dropna().unique() if dropna else values.unique()
    try:
        ordered = sorted(uniques)
    except TypeError:
        ordered = sorted(uniques, key=str)

    if per_col is not None:
        per_col[(col, dropna)] = ordered
    return list(ordered)


def selected_cell_view(adata, selected_cells):
    """Return ``adata`` restricted to selected observation names, in data order.

    An empty selection means "all cells", matching the dashboard's lasso store
    semantics.  Matching through string values also supports non-string AnnData
    indices without changing their original dtype.

    Backed AnnData objects cannot be sliced from an existing view.  Global Data
    Filter already creates such a view, so combine the lasso positions with that
    view's positions and slice its backed parent exactly once.  This keeps the
    matrix lazy instead of materialising the filtered dataset in memory.
    """
    if selected_cells is None or len(selected_cells) == 0:
        return adata
    requested = {str(cell) for cell in selected_cells}
    mask = adata.obs_names.astype(str).isin(requested)

    if getattr(adata, "isbacked", False) and getattr(adata, "is_view", False):
        parent = getattr(adata, "_adata_ref", None)
        parent_index = getattr(adata, "_oidx", None)
        if parent is not None and parent_index is not None:
            # AnnData exposes no public mapping from view-relative positions to
            # parent positions.  ``_oidx`` is the mapping it uses internally;
            # resolving it here avoids the unsupported ``backed_view[mask]``.
            view_positions = np.arange(parent.n_obs)[parent_index]
            return parent[np.asarray(view_positions)[mask]]

    return adata[mask]


def selection_group_values(adata, selected_cells) -> pd.Series:
    """Return the session-local lasso grouping without mutating ``adata.obs``."""
    selected = {str(cell) for cell in (selected_cells or [])}
    values = np.where(
        adata.obs_names.astype(str).isin(selected),
        SELECTED_LABEL,
        OTHERS_LABEL,
    )
    return pd.Series(
        pd.Categorical(values, categories=SELECTION_LABELS, ordered=True),
        index=adata.obs_names,
        name=SELECTION_GROUP,
    )


def cell_list_signature(cells):
    """Compact, order-sensitive signature for a potentially large cell-ID list."""
    if not cells:
        return None
    payload = "\0".join(str(cell) for cell in cells).encode()
    return {"len": len(cells), "hash": hashlib.md5(payload).hexdigest()}


def selection_group_context(adata, selection_data):
    """Return the highlighting universe and its Selected/Others group values."""
    selection_data = selection_data or {}
    selected_cells = selection_data.get("selected_cells") or []
    universe_cells = selection_data.get("universe_cells")
    plot_adata = selected_cell_view(adata, universe_cells)
    return plot_adata, selection_group_values(plot_adata, selected_cells)


def selection_group_signature(selection_data):
    """Compact signature covering both the lasso and its filtered universe."""
    if not selection_data:
        return None
    return {
        "selected": cell_list_signature(selection_data.get("selected_cells")),
        "universe": cell_list_signature(selection_data.get("universe_cells")),
    }
