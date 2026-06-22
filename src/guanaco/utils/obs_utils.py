"""Shared helpers for reading obs annotations consistently across plots."""

import weakref

from guanaco.data.loader import obs_col

# id(dataset) -> (weakref to the dataset, {(col, dropna): sorted categories}).
# Keyed by id() because AnnData/MuData are unhashable, and guarded by a weakref
# so an entry is dropped the moment its dataset is collected -- this prevents a
# later object that happens to reuse the same id() from reading a stale result.
_category_cache: dict = {}


def sorted_categories(adata, col, *, dropna: bool = True) -> list:
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
