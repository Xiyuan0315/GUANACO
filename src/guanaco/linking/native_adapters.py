"""Thin linking contracts for GUANACO's existing AnnData plot builders.

The builders own transforms and visual style. This module only applies a
``ViewState`` and translates browser marks to stable cell/feature IDs.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from inspect import signature
from types import MappingProxyType
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate

from guanaco.pages.matrix.plots.dotmatrix import plot_dot_matrix
from guanaco.pages.matrix.plots.embedding import (
    _resolve_embedding_coords,
    plot_coexpression_embedding,
    plot_embedding,
)
from guanaco.pages.matrix.plots.heatmap import plot_unified_heatmap
from guanaco.pages.matrix.plots.pseudotime import plot_genes_in_pseudotime
from guanaco.pages.matrix.plots.stacked_bar import plot_stacked_bar
from guanaco.pages.matrix.plots.violin1 import (
    _build_group_sample,
    _label_masks,
    plot_ridge,
    plot_violin1,
)
from guanaco.pages.matrix.plots.violin2 import plot_violin2_new
from guanaco.pages.matrix.plots.volcano import load_volcano_payload, plot_volcano
from guanaco.utils.obs_utils import sorted_categories
from guanaco.utils.plot_style import (
    GUANACO_QUALITATIVE,
    categorical_color_map,
)

from .base import PlotAdapter, Rendered
from .data import DataSource, DataStore
from .model import MarkMembers, ViewSpec, ViewState

GRAPH_EVENTS = MappingProxyType({"click": "clickData", "select": "selectedData"})
CLICK_EVENT = MappingProxyType({"click": "clickData"})
Selector = Callable[[Mapping[str, Any]], Any]


def _source(spec: ViewSpec, store: DataStore) -> DataSource:
    source = store.source(spec.data)
    if source.kind != "anndata":
        raise TypeError(f"{spec.plot!r} requires an AnnData-like source.")
    return source


def _list(value: Any) -> list[Any]:
    if value is None:
        return []
    return [value] if isinstance(value, str) else list(value)


def _unique(values: Sequence[Any]) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            str(value) for value in values if value is not None and str(value).strip()
        )
    )


def _points(payload: Mapping[str, Any] | None) -> list[Mapping[str, Any]]:
    return [
        point
        for point in ((payload or {}).get("points") or ())
        if isinstance(point, Mapping)
    ]


def _custom(point: Mapping[str, Any]) -> list[Any]:
    value = point.get("customdata")
    if isinstance(value, np.ndarray):
        value = value.tolist()
    return list(value) if isinstance(value, (list, tuple)) else [value]


def _column(index: int) -> Selector:
    def select(point: Mapping[str, Any]) -> Any:
        row = _custom(point)
        return row[index] if len(row) > index else None

    return select


def _decode(
    payload: Mapping[str, Any] | None,
    *,
    cells: Selector | None = None,
    features: Selector | None = None,
) -> MarkMembers | None:
    points = _points(payload)
    cell_ids = _unique([cells(point) for point in points]) if cells else ()
    feature_ids = _unique([features(point) for point in points]) if features else ()
    if not cell_ids and not feature_ids:
        return None
    return MarkMembers(cells=cell_ids, features=feature_ids)


def _state(value: ViewState | None) -> ViewState:
    return ViewState.empty() if value is None else value


def _mask(adata: Any, ids: Sequence[str]) -> np.ndarray:
    wanted = set(ids)
    return np.fromiter(
        (str(cell) in wanted for cell in adata.obs_names), bool, adata.n_obs
    )


def _subset(adata: Any, state: ViewState) -> Any:
    if state.cells is None or state.is_highlighted("cell"):
        return adata
    return adata[_mask(adata, state.cells)]


def _selection(adata: Any, ids: Sequence[str]) -> pd.Series:
    wanted = set(ids)
    values = pd.Categorical(
        ["Selected" if str(cell) in wanted else "Others" for cell in adata.obs_names],
        categories=["Selected", "Others"],
        ordered=True,
    )
    return pd.Series(values, index=adata.obs_names, name="Selection")


def _selected_labels(values: pd.Series) -> list[str]:
    return [name for name in ("Selected", "Others") if (values == name).any()]


def _features(
    spec: ViewSpec, state: ViewState, option: str, *, one: bool = False
) -> list[str]:
    value: Any = state.features
    if value is None:
        value = spec.options.get(option)
    result = [str(item) for item in _list(value)]
    return result[:1] if one else result


def _valid_features(adata: Any, values: Sequence[str]) -> list[str]:
    available = set(map(str, adata.var_names))
    return [str(value) for value in values if str(value) in available]


def _feature_data(
    spec: ViewSpec,
    store: DataStore,
    value: ViewState | None,
    option: str,
    *,
    one=False,
):
    full, state = _source(spec, store).data, _state(value)
    adata = _subset(full, state)
    return (
        full,
        state,
        adata,
        _valid_features(full, _features(spec, state, option, one=one)),
    )


def _validate_features(
    spec: ViewSpec,
    store: DataStore,
    option: str,
    *,
    one: bool = False,
    groupby: bool = False,
) -> Any:
    adata = _source(spec, store).data
    genes = _features(spec, ViewState.empty(), option)
    if not genes or (one and len(genes) != 1):
        qualifier = "exactly one " if one else "at least one "
        raise ValueError(f"View {spec.id!r} requires {qualifier}valid feature.")
    missing = [gene for gene in genes if gene not in adata.var_names]
    if missing:
        raise ValueError(f"Unknown feature {missing[0]!r} in view {spec.id!r}.")
    if groupby:
        name = spec.options.get("groupby")
        if not name or name not in adata.obs.columns:
            raise ValueError(f"View {spec.id!r} requires a valid `groupby`.")
    return adata


def _builder_kwargs(
    spec: ViewSpec,
    builder: Callable[..., Any],
    *,
    aliases: Mapping[str, str] | None = None,
    **required: Any,
) -> dict[str, Any]:
    """Forward native builder options without copying their API here."""

    accepted = signature(builder).parameters
    aliases = aliases or {}
    result = {
        aliases.get(name, name): value
        for name, value in spec.options.items()
        if aliases.get(name, name) in accepted
    }
    result.update(required)
    return result


def _group_context(
    adata: Any, state: ViewState, groupby: str, labels: Any = None
) -> tuple[pd.Series, list[Any], str]:
    if state.is_highlighted("cell") and state.cells is not None:
        values = _selection(adata, state.cells)
        return values, _selected_labels(values), "Selection"
    values = adata.obs[groupby]
    return values, _list(labels) or sorted_categories(adata, groupby), groupby


def _message(spec: ViewSpec, message: str, fallback: str) -> Rendered:
    figure = go.Figure()
    figure.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font={"size": 14, "color": "#64748B"},
    )
    figure.update_xaxes(visible=False)
    figure.update_yaxes(visible=False)
    return _finish(figure, spec, fallback)


def _finish(
    figure: go.Figure,
    spec: ViewSpec,
    fallback: str,
) -> Rendered:
    figure.update_layout(title={"text": spec.title or fallback})
    return figure


def _highlight(
    figure: go.Figure,
    ids: Sequence[str] | None,
    *,
    index: Sequence[Any] | None = None,
) -> None:
    if ids is None:
        return
    wanted = set(ids)
    for trace in figure.data:
        raw = getattr(trace, "ids", None)
        if raw is None:
            raw = getattr(trace, "customdata", None)
        if raw is None or not hasattr(trace, "selectedpoints"):
            continue
        data = np.asarray(raw, dtype=object)
        values = data if data.ndim == 1 else data[:, 0]
        if index is not None:
            try:
                values = np.asarray(index, dtype=object)[values.astype(int)]
            except (IndexError, TypeError, ValueError):
                pass
        trace.selectedpoints = [
            index for index, value in enumerate(values) if str(value) in wanted
        ]
        trace.selected = {"marker": {"opacity": 1.0}}
        trace.unselected = {"marker": {"opacity": 0.14, "color": "#D1D5DB"}}


def _basis(spec: ViewSpec, source: DataSource) -> str:
    value = spec.options.get("basis") or {
        "umap": "X_umap",
        "pca": "X_pca",
        "tsne": "X_tsne",
        "spatial": "spatial",
    }.get(spec.plot, "X_umap")
    value = str(value)
    if value not in source.data.obsm and f"X_{value}" in source.data.obsm:
        return f"X_{value}"
    return value


def _embedding_xy(
    adata: Any, basis: str, spec: ViewSpec
) -> tuple[np.ndarray, np.ndarray]:
    x, y, *_ = _resolve_embedding_coords(
        adata,
        basis,
        x_axis=spec.options.get("x_axis"),
        y_axis=spec.options.get("y_axis"),
        img_key=spec.options.get("img_key"),
        library_id=spec.options.get("library_id"),
        auto_select_spatial=True,
    )
    return np.asarray(x), np.asarray(y)


def _fixed_ranges(figure: go.Figure, adata: Any, basis: str, spec: ViewSpec) -> None:
    try:
        x, y = _embedding_xy(adata, basis, spec)
    except (KeyError, TypeError, ValueError):
        return
    for values, layout_axis, update in (
        (x, figure.layout.xaxis, figure.update_xaxes),
        (y, figure.layout.yaxis, figure.update_yaxes),
    ):
        values = values[np.isfinite(values.astype(float))].astype(float)
        if values.size and not layout_axis.range:
            low, high = float(values.min()), float(values.max())
            pad = 0.04 * (high - low) or 0.5
            update(range=[low - pad, high + pad], autorange=False)


class EmbeddingAdapter(PlotAdapter):
    """Any 2D AnnData embedding, including spatial coordinates and images."""

    events = GRAPH_EVENTS
    emits = frozenset({"cell"})
    accepts = frozenset({"cell", "feature"})

    def identity_ids(self, by, spec, store):
        if by != "feature":
            return None
        adata = _source(spec, store).data
        return _unique([*adata.var_names, *adata.obs.columns])

    def validate(self, spec, store):
        source = _source(spec, store)
        if _basis(spec, source) not in source.data.obsm:
            raise ValueError(f"View {spec.id!r} cannot find its embedding in `obsm`.")
        color = spec.options.get("color")
        if not color or color not in {*source.data.obs.columns, *source.data.var_names}:
            raise ValueError(f"Embedding view {spec.id!r} requires a valid `color`.")
        if spec.options.get("color_mode", "auto") not in {
            "auto",
            "categorical",
            "continuous",
        }:
            raise ValueError("`color_mode` must be auto, categorical, or continuous.")
        if spec.options.get("render_backend", "scattergl") not in {
            "scattergl",
            "datashader",
        }:
            raise ValueError("`render_backend` must be scattergl or datashader.")
        color_range = spec.options.get("color_range")
        if color_range is not None:
            values = tuple(float(value) for value in color_range)
            if (
                len(values) != 2
                or not all(np.isfinite(values))
                or values[0] >= values[1]
            ):
                raise ValueError(
                    "`color_range` must contain two finite increasing numbers."
                )
        midpoint = spec.options.get("color_midpoint")
        if midpoint is not None and not np.isfinite(float(midpoint)):
            raise ValueError("`color_midpoint` must be finite.")

    def render(self, spec, store, state=None, *, component_id=None):
        source = _source(spec, store)
        full, state = source.data, _state(state)
        adata = _subset(full, state)
        chosen = list(state.features or ())
        color = str(chosen[0] if chosen else spec.options.get("color", ""))
        if color not in {*full.obs.columns, *full.var_names}:
            return _message(spec, f"Feature {color!r} is unavailable.", "Embedding")
        if not adata.n_obs:
            return _message(spec, "No matching cells.", color)
        palette = None
        if color in full.obs.columns:
            palette = list(
                categorical_color_map(
                    sorted_categories(full, color), spec.options.get("palette")
                ).values()
            )
        basis = _basis(spec, source)
        legend = str(spec.options.get("legend_loc", "on legend"))
        backend = (
            "scattergl"
            if state.is_highlighted("cell")
            else str(spec.options.get("render_backend", "scattergl"))
        )
        figure = plot_embedding(
            **_builder_kwargs(
                spec,
                plot_embedding,
                aliases={
                    "basis": "embedding_key",
                    "color_mode": "mode",
                    "color_map": "continuous_color_map",
                    "size": "marker_size",
                },
                adata=adata,
                adata_full=full,
                source_adata=adata,
                embedding_key=basis,
                color=color,
                discrete_color_map=palette,
                legend_show="on data" if legend == "on data" else "on legend",
                render_backend=backend,
            )
        )
        for trace in figure.data:
            marker = getattr(trace, "marker", None)
            if marker is None:
                continue
            if spec.options.get("color_range") is not None:
                marker.cmin, marker.cmax = map(float, spec.options["color_range"])
            if spec.options.get("color_midpoint") is not None:
                marker.cmid = float(spec.options["color_midpoint"])
        if state.is_highlighted("cell"):
            _highlight(figure, state.cells, index=adata.obs_names.astype(str))
        _fixed_ranges(figure, full, basis, spec)
        figure.update_layout(selectionrevision=f"{spec.id}:{color}")
        return _finish(figure, spec, color)

    def decode_event(self, event, payload, spec, store):
        del event
        first = _column(0)

        names = np.asarray(_source(spec, store).data.obs_names.astype(str))

        def cell(point):
            value = point.get("id")
            if value is not None:
                return value
            value = first(point)
            if isinstance(value, (int, np.integer)):
                return names[int(value)]
            return value

        return _decode(
            payload,
            cells=cell,
        )


class ViolinAdapter(PlotAdapter):
    """GUANACO's multi-feature stacked violin."""

    events = GRAPH_EVENTS
    emits = frozenset({"cell", "feature"})
    accepts = frozenset({"cell", "feature"})

    def validate(self, spec, store):
        _validate_features(spec, store, "keys", groupby=True)

    def render(self, spec, store, state=None, *, component_id=None):
        full, state, adata, genes = _feature_data(spec, store, state, "keys")
        if not genes or not adata.n_obs:
            message = "No matching features." if not genes else "No matching cells."
            return _message(spec, message, "Violin")
        groupby = str(spec.options["groupby"])
        values, labels, group_title = _group_context(
            adata, state, groupby, spec.options.get("labels")
        )
        try:
            figure = plot_violin1(
                **_builder_kwargs(
                    spec,
                    plot_violin1,
                    adata=adata,
                    genes=genes,
                    groupby=groupby,
                    labels=labels,
                    adata_obs=full.obs,
                    data_already_filtered=True,
                    group_values=values,
                    groupby_label_color_map=categorical_color_map(
                        labels, spec.options.get("palette")
                    ),
                )
            )
        except PreventUpdate:
            return _message(spec, "No matching cells.", "Violin")
        samples = _build_group_sample(
            _label_masks(values, labels), labels, np.random.default_rng(0)
        )
        cell_ids = {
            str(label): tuple(adata.obs_names[index].astype(str))
            for label, index in samples.items()
        }
        selectable = bool(spec.options.get("_link_source", False))
        for trace in figure.data:
            if trace.type != "violin":
                continue
            axis = str(getattr(trace, "yaxis", "y") or "y")
            try:
                number = int(axis[1:] or "1")
            except ValueError:
                number = 1
            gene = genes[min((number - 1) // 2, len(genes) - 1)]
            count = len(trace.y) if trace.y is not None else 0
            cells = cell_ids.get(str(trace.name), ())
            if len(cells) != count:
                cells = ("",) * count
            trace.customdata = np.column_stack(
                [
                    np.full(count, gene, object),
                    np.full(count, str(trace.name), object),
                    np.asarray(cells, object),
                ]
            )
            trace.ids = list(cells)
            if selectable:
                trace.update(
                    points="all", hoveron="points+violins", jitter=0.16, pointpos=0
                )
                trace.marker.update(opacity=0.5, size=4)
            trace.hovertemplate = (
                f"Feature: {gene}<br>{group_title}: %{{customdata[1]}}"
                "<br>Cell: %{customdata[2]}<br>Expression: %{y:.3g}<extra></extra>"
            )
        return _finish(figure, spec, "Expression distributions")

    def decode_event(self, event, payload, spec, store):
        del event, spec, store
        return _decode(payload, cells=_column(2), features=_column(0))


class FeatureDistributionAdapter(PlotAdapter):
    """GUANACO's single-feature ridge or comparative violin."""

    events = CLICK_EVENT
    emits = frozenset({"feature"})
    accepts = frozenset({"cell", "feature"})

    def __init__(self, kind: str) -> None:
        if kind not in {"ridge", "split_violin"}:
            raise ValueError(f"Unsupported feature distribution: {kind!r}.")
        self.kind = kind

    def validate(self, spec, store):
        adata = _validate_features(spec, store, "key", one=True, groupby=True)
        if self.kind == "split_violin":
            mode = str(spec.options.get("mode", "mode1"))
            groupby2 = spec.options.get("groupby2")
            if mode != "mode1" and not groupby2:
                raise ValueError(
                    f"Comparative violin mode {mode!r} requires `groupby2`."
                )
            if groupby2 and groupby2 not in adata.obs.columns:
                raise ValueError(f"Unknown secondary grouping column {groupby2!r}.")

    def render(self, spec, store, state=None, *, component_id=None):
        full, state, adata, genes = _feature_data(spec, store, state, "key", one=True)
        if not genes or not adata.n_obs:
            message = "No matching feature." if not genes else "No matching cells."
            return _message(spec, message, self.kind)
        gene, groupby = genes[0], str(spec.options["groupby"])
        values, labels, group_title = _group_context(
            adata, state, groupby, spec.options.get("labels")
        )
        highlighted = group_title == "Selection"
        if self.kind == "ridge":
            figure = plot_ridge(
                **_builder_kwargs(
                    spec,
                    plot_ridge,
                    adata=adata,
                    gene=gene,
                    groupby=groupby,
                    labels=labels,
                    adata_obs=full.obs,
                    data_already_filtered=True,
                    group_values=values if highlighted else None,
                    groupby_label_color_map=categorical_color_map(
                        labels, spec.options.get("palette")
                    ),
                )
            )
        else:
            mode = "mode1" if highlighted else str(spec.options.get("mode", "mode1"))
            groupby2 = None if mode == "mode1" else spec.options.get("groupby2")
            color_labels = labels + (
                _group_context(adata, state, str(groupby2))[1] if groupby2 else []
            )
            figure = plot_violin2_new(
                **_builder_kwargs(
                    spec,
                    plot_violin2_new,
                    adata=adata,
                    key=gene,
                    meta1=group_title if highlighted else groupby,
                    meta2=groupby2,
                    mode=mode,
                    labels=labels,
                    color_map=categorical_color_map(
                        color_labels, spec.options.get("palette")
                    ),
                    group_values={"Selection": values} if highlighted else None,
                )
            )
        for trace in figure.data:
            if trace.type == "violin":
                count = len(trace.x) if trace.x is not None else len(trace.y or ())
                trace.customdata = np.full((count, 1), gene, object)
        return _finish(figure, spec, f"{self.kind.replace('_', ' ').title()} · {gene}")

    def decode_event(self, event, payload, spec, store):
        del event, spec, store
        return _decode(payload, features=_column(0))


class FeatureGroupMatrixAdapter(PlotAdapter):
    """GUANACO dot/matrix plot; a mark carries a feature and its group cells."""

    emits = frozenset({"cell", "feature"})
    accepts = frozenset({"cell", "feature"})

    def __init__(self, plot_type: str) -> None:
        if plot_type not in {"dotplot", "matrixplot"}:
            raise ValueError(f"Unsupported feature matrix: {plot_type!r}.")
        self.plot_type = plot_type
        self.events = GRAPH_EVENTS if plot_type == "dotplot" else CLICK_EVENT

    def validate(self, spec, store):
        _validate_features(spec, store, "var_names", groupby=True)

    def render(self, spec, store, state=None, *, component_id=None):
        source, state = _source(spec, store), _state(state)
        adata = source.data
        genes = _valid_features(adata, _features(spec, state, "var_names"))
        empty_cells = (
            state.cells is not None
            and not state.is_highlighted("cell")
            and not state.cells
        )
        if not genes or empty_cells:
            message = "No matching features." if not genes else "No matching cells."
            return _message(spec, message, self.plot_type)
        groupby = str(spec.options["groupby"])
        values, labels, group_title = _group_context(
            adata, state, groupby, spec.options.get("labels")
        )
        selected_cells = (
            list(state.cells)
            if state.cells is not None and not state.is_highlighted("cell")
            else None
        )
        transpose = bool(spec.options.get("transpose", False))
        try:
            figure = plot_dot_matrix(
                **_builder_kwargs(
                    spec,
                    plot_dot_matrix,
                    adata=adata,
                    genes=genes,
                    groupby=groupby,
                    selected_labels=labels,
                    plot_type=self.plot_type,
                    transpose=transpose,
                    selected_cells=selected_cells,
                    group_values=values if group_title == "Selection" else None,
                )
            )
        except PreventUpdate:
            return _message(spec, "No matching cells or features.", self.plot_type)
        return _finish(figure, spec, self.plot_type.title())

    def decode_event(self, event, payload, spec, store):
        del event
        points = _points(payload)
        transpose = bool(spec.options.get("transpose", False))
        feature_key, group_key = ("y", "x") if transpose else ("x", "y")
        features = _unique([point.get(feature_key) for point in points])
        groups = _unique([point.get(group_key) for point in points])
        if not features and not groups:
            return None
        adata = _source(spec, store).data
        values = adata.obs[str(spec.options["groupby"])].astype(str)
        cells = tuple(
            dict.fromkeys(
                cell
                for group in groups
                for cell in adata.obs_names[values.eq(group)].astype(str)
            )
        )
        return MarkMembers(cells=cells, features=features)


class HeatmapAdapter(PlotAdapter):
    """GUANACO's cell-level multi-feature heatmap."""

    events = CLICK_EVENT
    emits = frozenset({"feature"})
    accepts = frozenset({"cell", "feature"})

    def validate(self, spec, store):
        _validate_features(spec, store, "var_names", groupby=True)

    def render(self, spec, store, state=None, *, component_id=None):
        full, state, adata, genes = _feature_data(spec, store, state, "var_names")
        if not genes or not adata.n_obs:
            message = "No matching features." if not genes else "No matching cells."
            return _message(spec, message, "Heatmap")
        original_group = str(spec.options["groupby"])
        values, labels, groupby = _group_context(
            full, state, original_group, spec.options.get("labels")
        )
        highlighted = groupby == "Selection"
        overrides = {groupby: values} if highlighted else None
        try:
            figure = plot_unified_heatmap(
                **_builder_kwargs(
                    spec,
                    plot_unified_heatmap,
                    aliases={"groupby": "groupby1"},
                    adata=adata,
                    genes=genes,
                    groupby1=groupby,
                    groupby2=None if highlighted else spec.options.get("groupby2"),
                    labels=labels,
                    groupby1_label_color_map=(
                        categorical_color_map(labels, spec.options.get("palette"))
                        if highlighted
                        else None
                    ),
                    adata_obs=(
                        pd.DataFrame({groupby: values}, index=full.obs_names)
                        if highlighted
                        else full.obs
                    ),
                    data_already_filtered=not highlighted,
                    color_config=list(GUANACO_QUALITATIVE),
                    obs_overrides=overrides,
                )
            )
        except PreventUpdate:
            return _message(spec, "No matching cells or features.", "Heatmap")
        return _finish(figure, spec, "Expression heatmap")

    def decode_event(self, event, payload, spec, store):
        del event, spec, store
        return _decode(payload, features=lambda point: point.get("y"))


class VolcanoAdapter(PlotAdapter):
    """GUANACO's precomputed differential-expression volcano."""

    events = GRAPH_EVENTS
    emits = frozenset({"feature"})
    accepts = frozenset({"feature"})

    def validate(self, spec, store):
        entries = load_volcano_payload(_source(spec, store).data).get("entries") or {}
        if not entries:
            raise ValueError(
                "Volcano views require precomputed differential expression."
            )
        if (
            spec.options.get("group") is not None
            and spec.options["group"] not in entries
        ):
            raise ValueError(f"Unknown volcano comparison {spec.options['group']!r}.")

    def render(self, spec, store, state=None, *, component_id=None):
        entries = load_volcano_payload(_source(spec, store).data)["entries"]
        group = spec.options.get("group") or next(iter(entries))
        figure = plot_volcano(
            **_builder_kwargs(
                spec, plot_volcano, entry_name=str(group), entry=entries[group]
            )
        )
        _highlight(figure, _state(state).features)
        return _finish(figure, spec, f"Volcano · {group}")

    def decode_event(self, event, payload, spec, store):
        del event, spec, store
        return _decode(payload, features=_column(0))


class CompositionAdapter(PlotAdapter):
    """GUANACO composition bars whose marks expand to member cells."""

    events = GRAPH_EVENTS
    emits = frozenset({"cell"})
    accepts = frozenset({"cell"})

    def validate(self, spec, store):
        adata = _source(spec, store).data
        invalid = [
            name
            for name in ("x", "color")
            if not spec.options.get(name) or spec.options[name] not in adata.obs.columns
        ]
        if invalid:
            raise ValueError(
                f"Composition view {spec.id!r} requires valid `x` and `color`."
            )

    def render(self, spec, store, state=None, *, component_id=None):
        full, state = _source(spec, store).data, _state(state)
        adata = _subset(full, state)
        if not adata.n_obs:
            return _message(spec, "No matching cells.", "Composition")
        x, color = str(spec.options["x"]), str(spec.options["color"])
        overrides = None
        if state.is_highlighted("cell") and state.cells is not None:
            color = "Selection"
            overrides = {color: _selection(adata, state.cells)}
        swap = bool(spec.options.get("swap_axes", False))
        figure = plot_stacked_bar(
            **_builder_kwargs(
                spec,
                plot_stacked_bar,
                aliases={"x": "x_meta", "color": "y_meta", "color_order": "y_order"},
                x_meta=x,
                y_meta=color,
                norm="prop"
                if spec.options.get("normalize", "proportion")
                in {"proportion", "prop", True}
                else "count",
                adata=adata,
                color_map=categorical_color_map(
                    sorted_categories(adata, color, overrides=overrides),
                    spec.options.get("palette"),
                ),
                x_order=_list(spec.options.get("x_order")) or None,
                y_order=_list(spec.options.get("color_order")) or None,
                swap_axes=swap,
                group_values=overrides,
            )
        )
        for trace in figure.data:
            if trace.type != "bar":
                continue
            groups, heights = (trace.y, trace.x) if swap else (trace.x, trace.y)
            category = str(trace.name)
            trace.customdata = np.asarray(
                [
                    [str(group), category, value]
                    for group, value in zip(groups, heights, strict=True)
                ],
                object,
            )
        return _finish(figure, spec, "Composition")

    def decode_event(self, event, payload, spec, store):
        del event
        rows = [_custom(point) for point in _points(payload)]
        rows = [row for row in rows if len(row) >= 2]
        if not rows:
            return None
        adata = _source(spec, store).data
        x, color = str(spec.options["x"]), str(spec.options["color"])
        cells = []
        for x_value, color_value, *_ in rows:
            keep = adata.obs[x].astype(str).eq(str(x_value)) & adata.obs[color].astype(
                str
            ).eq(str(color_value))
            cells.extend(adata.obs_names[keep].astype(str))
        return MarkMembers(cells=_unique(cells)) if cells else None


class PseudotimeAdapter(PlotAdapter):
    """GUANACO expression-along-pseudotime plot."""

    events = GRAPH_EVENTS
    emits = frozenset({"cell", "feature"})
    accepts = frozenset({"cell", "feature"})

    def validate(self, spec, store):
        adata = _validate_features(spec, store, "genes")
        key = str(spec.options.get("pseudotime_key", "pseudotime"))
        if key not in adata.obs.columns:
            raise ValueError(f"Pseudotime column {key!r} is unavailable.")

    def render(self, spec, store, state=None, *, component_id=None):
        _, state, adata, genes = _feature_data(spec, store, state, "genes")
        if not genes or not adata.n_obs:
            message = "No matching features." if not genes else "No matching cells."
            return _message(spec, message, "Pseudotime")
        groupby = spec.options.get("groupby")
        values = None
        if state.is_highlighted("cell") and state.cells is not None:
            groupby, values = "Selection", _selection(adata, state.cells)
        color_map = None
        if groupby:
            labels = (
                _selected_labels(values)
                if values is not None
                else sorted_categories(adata, str(groupby))
            )
            color_map = categorical_color_map(labels, spec.options.get("palette"))
        figure = plot_genes_in_pseudotime(
            **_builder_kwargs(
                spec,
                plot_genes_in_pseudotime,
                aliases={"size": "marker_size"},
                adata=adata,
                genes=genes,
                groupby=groupby,
                color_map=color_map,
                group_values=values,
            )
        )
        return _finish(figure, spec, "Expression along pseudotime")

    def decode_event(self, event, payload, spec, store):
        del event, spec, store
        return _decode(payload, cells=_column(0), features=_column(1))


class CoexpressionAdapter(EmbeddingAdapter):
    """GUANACO two-feature co-expression on any embedding."""

    emits = frozenset({"cell"})

    def identity_ids(self, by, spec, store):
        if by == "feature":
            return tuple(_source(spec, store).data.var_names.astype(str))
        return None

    def validate(self, spec, store):
        source = _source(spec, store)
        if _basis(spec, source) not in source.data.obsm:
            raise ValueError(f"View {spec.id!r} cannot find its embedding in `obsm`.")
        for name in ("gene1", "gene2"):
            if spec.options.get(name) not in source.data.var_names:
                raise ValueError(f"View {spec.id!r} requires valid `{name}`.")

    def render(self, spec, store, state=None, *, component_id=None):
        source, state = _source(spec, store), _state(state)
        full = source.data
        adata = _subset(full, state)
        chosen = list(state.features or ())
        gene1 = str(chosen[0] if chosen else spec.options.get("gene1"))
        gene2 = str(chosen[1] if len(chosen) > 1 else spec.options.get("gene2"))
        if (
            gene1 not in full.var_names
            or gene2 not in full.var_names
            or not adata.n_obs
        ):
            message = (
                "No matching feature pair." if adata.n_obs else "No matching cells."
            )
            return _message(spec, message, "Co-expression")
        basis = _basis(spec, source)
        legend = spec.options.get("legend_loc", "right")
        figure = plot_coexpression_embedding(
            **_builder_kwargs(
                spec,
                plot_coexpression_embedding,
                aliases={"basis": "embedding_key", "size": "marker_size"},
                adata=adata,
                embedding_key=basis,
                gene1=gene1,
                gene2=gene2,
                legend_show="on data" if str(legend) == "on data" else "right",
                source_adata=adata,
            )
        )
        if state.is_highlighted("cell"):
            _highlight(figure, state.cells, index=adata.obs_names.astype(str))
        _fixed_ranges(figure, full, basis, spec)
        return _finish(figure, spec, f"Co-expression · {gene1} + {gene2}")


__all__ = [
    "CoexpressionAdapter",
    "CompositionAdapter",
    "EmbeddingAdapter",
    "FeatureDistributionAdapter",
    "FeatureGroupMatrixAdapter",
    "HeatmapAdapter",
    "PseudotimeAdapter",
    "ViolinAdapter",
    "VolcanoAdapter",
]
