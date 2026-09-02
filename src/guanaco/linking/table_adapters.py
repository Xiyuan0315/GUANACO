"""Generic linked plots for indexed pandas tables.

The table index identifies atomic rows. A point represents one row; an
aggregated heatmap tile or network edge represents the rows it contains. A
row link may additionally use a shared logical key to expand one parent row
into many child rows. For mixed table/AnnData links the same key maps table
rows to native cell or feature IDs. The link's ``by`` value decides whether
those IDs mean rows, cells, or features.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import dash_cytoscape as cyto
from dash import html
from PIL import Image as PILImage

from guanaco.utils.plot_style import (
    GUANACO_QUALITATIVE,
    categorical_color_map,
)

from .base import PlotAdapter, Rendered
from .data import DataSource, DataStore
from .model import MarkMembers, ViewSpec, ViewState

PLOTLY_EVENTS = MappingProxyType({"click": "clickData", "select": "selectedData"})
NETWORK_EVENTS = MappingProxyType({"click": "tapEdgeData"})
_IDENTITIES = frozenset({"row", "cell", "feature"})
_ROW_ID = "__guanaco_row_id__"
_AGGREGATES = frozenset({"count", "sum", "mean", "median", "min", "max"})


def _state(state: ViewState | None) -> ViewState:
    return ViewState.empty() if state is None else state


def _source(spec: ViewSpec, store: DataStore) -> DataSource:
    source = store.source(spec.data)
    if source.kind != "table":
        raise TypeError(f"Table view {spec.id!r} requires a pandas DataFrame.")
    return source


def _ids_for(state: ViewState, *, highlighted: bool) -> tuple[str, ...] | None:
    """Intersect active table-index domains with the requested action."""

    domains = [
        values
        for axis, values in (
            ("row", state.rows),
            ("cell", state.cells),
            ("feature", state.features),
        )
        if values is not None and state.is_highlighted(axis) is highlighted
    ]
    if not domains:
        return None
    common = set(domains[0])
    for values in domains[1:]:
        common.intersection_update(values)
    return tuple(value for value in domains[0] if value in common)


def _frame(source: DataSource, state: ViewState) -> pd.DataFrame:
    row_ids = np.asarray(source.row_ids, dtype=object)
    selected = _ids_for(state, highlighted=False)
    if selected is None:
        frame = source.data.copy()
    else:
        mask = np.isin(row_ids, np.asarray(selected, dtype=object))
        frame = source.data.loc[mask].copy()
        row_ids = row_ids[mask]
    row_id = _ROW_ID
    while row_id in frame.columns:
        row_id = f"_{row_id}"
    frame[row_id] = row_ids
    frame.attrs["row_id"] = row_id
    return frame


def _require(frame: pd.DataFrame, columns: Sequence[Any], view_id: str) -> None:
    missing = [str(column) for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"View {view_id!r} is missing column {missing[0]!r}.")


def _points(payload: Mapping[str, Any] | None) -> list[Mapping[str, Any]]:
    if not payload:
        return []
    return [item for item in payload.get("points", ()) if isinstance(item, Mapping)]


def _custom_value(point: Mapping[str, Any]) -> str | None:
    value = point.get("customdata")
    if isinstance(value, np.ndarray):
        value = value.tolist()
    while isinstance(value, (list, tuple)):
        if not value:
            return None
        value = value[0]
    return None if value is None else str(value)


def _members(ids: Sequence[str]) -> MarkMembers | None:
    unique = tuple(dict.fromkeys(str(value) for value in ids))
    if not unique:
        return None
    return MarkMembers(rows=unique)


def _decode_points(payload: Mapping[str, Any] | None) -> MarkMembers | None:
    return _members(
        [value for point in _points(payload) if (value := _custom_value(point))]
    )


def _highlight(figure: go.Figure, ids: tuple[str, ...] | None) -> None:
    if ids is None:
        return
    selected = set(ids)
    for trace in figure.data:
        values = getattr(trace, "customdata", None)
        if values is None or not hasattr(trace, "selectedpoints"):
            continue
        array = np.asarray(values, dtype=object)
        if array.ndim > 1:
            array = array[:, 0]
        trace.selectedpoints = [
            i for i, value in enumerate(array) if str(value) in selected
        ]
        trace.selected = {"marker": {"opacity": 1.0}}
        trace.unselected = {"marker": {"opacity": 0.14, "color": "#CBD5E1"}}


def _numeric(series: pd.Series) -> bool:
    return bool(
        pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series)
    )


def _categories(series: pd.Series) -> list[Any]:
    if isinstance(series.dtype, pd.CategoricalDtype):
        present = set(series.dropna())
        return [value for value in series.cat.categories if value in present]
    return list(pd.unique(series.dropna()))


def _range(value: Any, name: str) -> list[Any]:
    if isinstance(value, (str, bytes)):
        raise ValueError(f"`{name}` must contain exactly two values.")
    try:
        result = list(value)
    except TypeError as error:
        raise ValueError(f"`{name}` must contain exactly two values.") from error
    if len(result) != 2 or any(item is None for item in result):
        raise ValueError(f"`{name}` must contain exactly two values.")
    return result


def _axis(
    full: pd.Series,
    visible: pd.Series,
    *,
    title: str,
    explicit: Any,
    fixed: bool,
    padding: float,
) -> dict[str, Any]:
    axis: dict[str, Any] = {"title": title}
    if explicit is not None:
        axis["range"] = _range(explicit, f"{full.name}_range")
    elif _numeric(full) and fixed:
        values = pd.to_numeric(full, errors="coerce").to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if len(values):
            low, high = float(values.min()), float(values.max())
            span = high - low or abs(low) or 1.0
            axis["range"] = [low - span * padding, high + span * padding]
    elif not _numeric(full):
        categories = _categories(visible)
        axis.update(type="category", categoryorder="array", categoryarray=categories)
        if categories:
            axis["range"] = [-0.5, len(categories) - 0.5]
    return axis


def _image_source(value: Any) -> Any:
    if isinstance(value, str) and value.startswith(
        ("data:image/", "http://", "https://")
    ):
        return value
    if isinstance(value, (str, Path)):
        path = Path(value)
        if not path.is_file():
            raise ValueError(f"Layout image path does not exist: {path!s}.")
        with PILImage.open(path) as image:
            return image.copy()
    if isinstance(value, np.ndarray):
        array = np.asarray(value)
        if array.ndim == 3 and array.shape[2] == 1:
            array = array[:, :, 0]
        if array.ndim not in {2, 3} or (
            array.ndim == 3 and array.shape[2] not in {3, 4}
        ):
            raise ValueError("A layout image array must be H×W, H×W×3, or H×W×4.")
        numeric = array.astype(float, copy=False)
        if not numeric.size or not np.isfinite(numeric).all():
            raise ValueError("A layout image array must contain finite values.")
        if np.issubdtype(array.dtype, np.floating) and numeric.max() <= 1:
            numeric *= 255
        return PILImage.fromarray(np.clip(np.rint(numeric), 0, 255).astype(np.uint8))
    if isinstance(value, PILImage.Image):
        return value
    raise TypeError("`layout_image.source` must be an image, array, path, or URL.")


def _layout_image(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping) or value.get("source") is None:
        raise TypeError("`layout_image` must be a mapping containing `source`.")
    result = dict(value)
    result["source"] = _image_source(result["source"])
    defaults = {
        "xref": "x",
        "yref": "y",
        "xanchor": "left",
        "yanchor": "top",
        "sizing": "stretch",
        "layer": "below",
        "opacity": 1.0,
    }
    for key, default in defaults.items():
        result.setdefault(key, default)
    return result


class _FigureAdapter(PlotAdapter):
    events = PLOTLY_EVENTS
    emits = _IDENTITIES
    accepts = _IDENTITIES

    def render(
        self,
        spec: ViewSpec,
        store: DataStore,
        state: ViewState | None = None,
        *,
        component_id: str | None = None,
    ) -> Rendered:
        del component_id
        return self._figure(spec, store, _state(state))


class TablePlotAdapter(_FigureAdapter):
    """Scatter, bar, line, and area plots with one mark per table row."""

    @staticmethod
    def _kind(spec: ViewSpec) -> str:
        return spec.plot.removeprefix("plotly.")

    def validate(self, spec: ViewSpec, store: DataStore) -> None:
        frame = _source(spec, store).data
        kind = self._kind(spec)
        if kind not in {"scatter", "bar", "line", "area"}:
            raise ValueError(f"View {spec.id!r} has unsupported plot kind {kind!r}.")
        x, y = spec.options.get("x"), spec.options.get("y")
        if not x or not y:
            raise ValueError(f"View {spec.id!r} requires `x` and `y` columns.")
        columns = [x, y] + [
            spec.options[name]
            for name in ("color", "group", "label")
            if spec.options.get(name) is not None
        ]
        size = spec.options.get("size", 8)
        if isinstance(size, str):
            columns.append(size)
        _require(frame, columns, spec.id)
        mode = str(spec.options.get("color_mode", "auto"))
        if mode not in {"auto", "categorical", "continuous"}:
            raise ValueError("`color_mode` must be auto, categorical, or continuous.")
        orientation = spec.options.get("orientation")
        if orientation is not None and orientation not in {"v", "h"}:
            raise ValueError("`orientation` must be 'v' or 'h'.")

    @staticmethod
    def _sizes(frame: pd.DataFrame, size: Any, limits: Sequence[float]) -> Any:
        if not isinstance(size, str):
            return float(size)
        values = pd.to_numeric(frame[size], errors="coerce")
        low, high = map(float, _range(limits, "size_range"))
        finite = values[np.isfinite(values)]
        if finite.empty or np.isclose(finite.min(), finite.max()):
            return pd.Series((low + high) / 2, index=frame.index)
        return low + (high - low) * (values - finite.min()) / (
            finite.max() - finite.min()
        )

    def _figure(self, spec: ViewSpec, store: DataStore, state: ViewState) -> go.Figure:
        source = _source(spec, store)
        frame = _frame(source, state)
        row_id = str(frame.attrs["row_id"])
        options = spec.options
        kind = self._kind(spec)
        x, y = str(options["x"]), str(options["y"])
        color = options.get("color")
        group = options.get("group")
        label = options.get("label")
        color_mode = str(options.get("color_mode", "auto"))
        continuous = color is not None and (
            color_mode == "continuous"
            or (color_mode == "auto" and _numeric(source.data[color]))
        )
        group_column = group or (
            color if color is not None and not continuous else None
        )
        groups = (
            [(None, frame)]
            if group_column is None
            else list(
                frame.groupby(group_column, sort=False, dropna=False, observed=True)
            )
        )
        palette: dict[Any, str] = {}
        if color is not None and not continuous:
            labels = source.data[color].astype("string").fillna("Missing")
            palette = categorical_color_map(labels.tolist(), options.get("palette"))
        sizes = self._sizes(
            frame, options.get("size", 8), options.get("size_range", (7, 22))
        )

        figure = go.Figure()
        for position, (group_value, rows) in enumerate(groups):
            if kind in {"line", "area"} and options.get("sort", True):
                rows = rows.sort_values(x, kind="stable")
            ids = rows[row_id].astype(str)
            fallback = GUANACO_QUALITATIVE[position % len(GUANACO_QUALITATIVE)]
            if continuous:
                marker_color: Any = pd.to_numeric(rows[color], errors="coerce")
            elif color is not None:
                marker_color = (
                    rows[color].astype("string").fillna("Missing").map(palette)
                )
            else:
                marker_color = fallback
            marker: dict[str, Any] = {
                "color": marker_color,
                "opacity": float(options.get("opacity", 0.82)),
            }
            if continuous:
                marker["coloraxis"] = "coloraxis"
            if kind != "bar":
                marker["size"] = (
                    sizes.loc[rows.index] if isinstance(sizes, pd.Series) else sizes
                )
            common: dict[str, Any] = {
                "x": rows[x],
                "y": rows[y],
                "ids": ids,
                "customdata": ids.to_numpy(dtype=object).reshape(-1, 1),
                "name": str(
                    group_value if group_value is not None else spec.title or spec.id
                ),
                "hovertemplate": f"{x}: %{{x}}<br>{y}: %{{y}}<extra></extra>",
            }
            if label is not None:
                common["text"] = rows[label].astype(str)
                common["hovertemplate"] = (
                    f"{label}: %{{text}}<br>" + common["hovertemplate"]
                )
            if kind == "bar":
                bar_options = {}
                if options.get("orientation") is not None:
                    bar_options["orientation"] = options["orientation"]
                figure.add_trace(go.Bar(**common, marker=marker, **bar_options))
            else:
                mode = "markers" if kind == "scatter" else "lines+markers"
                trace: dict[str, Any] = {
                    **common,
                    "marker": marker,
                    "mode": options.get("mode", mode),
                }
                if kind == "area":
                    trace.update(
                        mode=options.get("mode", "lines"),
                        fill=options.get("fill", "tozeroy"),
                        line={"color": fallback},
                    )
                figure.add_trace(go.Scatter(**trace))

        xaxis: dict[str, Any] = {"title": str(options.get("x_title", x))}
        yaxis: dict[str, Any] = {"title": str(options.get("y_title", y))}
        if kind == "scatter":
            padding = float(options.get("axis_padding", 0.05))
            if not np.isfinite(padding) or padding < 0:
                raise ValueError("`axis_padding` must be a non-negative number.")
            fixed = bool(options.get("fixed_axes", True))
            xaxis = _axis(
                source.data[x],
                frame[x],
                title=xaxis["title"],
                explicit=options.get("x_range"),
                fixed=fixed,
                padding=padding,
            )
            yaxis = _axis(
                source.data[y],
                frame[y],
                title=yaxis["title"],
                explicit=options.get("y_range"),
                fixed=fixed,
                padding=padding,
            )
            if options.get("equal_axes"):
                yaxis.update(scaleanchor="x", scaleratio=1, constrain="domain")

        layout: dict[str, Any] = {
            "barmode": options.get("barmode", "group"),
            "clickmode": "event+select",
            "dragmode": "lasso",
            "uirevision": f"{spec.id}:{kind}:{x}:{y}",
            "margin": {"b": 70, "l": 75, "r": 45, "t": 65},
            "xaxis": xaxis,
            "yaxis": yaxis,
        }
        if continuous:
            values = pd.to_numeric(source.data[color], errors="coerce")
            finite = values[np.isfinite(values)]
            color_range = options.get("color_range")
            if color_range is not None:
                cmin, cmax = map(float, _range(color_range, "color_range"))
            elif finite.empty:
                cmin = cmax = None
            else:
                cmin, cmax = float(finite.min()), float(finite.max())
            layout["coloraxis"] = {
                "colorscale": options.get("color_map", "Viridis"),
                "cmin": cmin,
                "cmax": cmax,
                "colorbar": {"title": str(options.get("colorbar_title", color))},
            }
            if options.get("color_midpoint") is not None:
                layout["coloraxis"]["cmid"] = float(options["color_midpoint"])
        figure.update_layout(**layout)
        if options.get("layout_image") is not None:
            figure.add_layout_image(_layout_image(options["layout_image"]))
        _highlight(figure, _ids_for(state, highlighted=True))
        figure.update_layout(title={"text": str(spec.title or spec.id)})
        return figure

    def decode_event(
        self,
        event: str,
        payload: Mapping[str, Any] | None,
        spec: ViewSpec,
        store: DataStore,
    ) -> MarkMembers | None:
        del event, spec, store
        return _decode_points(payload)


class TableHeatmapAdapter(_FigureAdapter):
    """Heatmap whose aggregated tiles carry their atomic table rows."""

    @staticmethod
    def _columns(spec: ViewSpec) -> tuple[str, str, str]:
        return (
            str(spec.options.get("x")),
            str(spec.options.get("y")),
            str(spec.options.get("value")),
        )

    def validate(self, spec: ViewSpec, store: DataStore) -> None:
        frame = _source(spec, store).data
        x, y, value = self._columns(spec)
        if "None" in {x, y, value}:
            raise ValueError(
                f"Heatmap view {spec.id!r} requires `x`, `y`, and `value`."
            )
        _require(frame, [x, y, value], spec.id)
        if frame[[x, y]].isna().to_numpy().any():
            raise ValueError(f"Heatmap view {spec.id!r} requires non-null coordinates.")
        aggregate = str(spec.options.get("aggregate", "mean")).lower()
        if aggregate not in _AGGREGATES:
            raise ValueError(
                f"Heatmap aggregate must be one of {', '.join(sorted(_AGGREGATES))}."
            )
        if aggregate != "count":
            numeric = pd.to_numeric(frame[value], errors="coerce")
            if (frame[value].notna() & ~np.isfinite(numeric)).any():
                raise ValueError(f"Heatmap value column {value!r} must be numeric.")

    @staticmethod
    def _axis_values(series: pd.Series, explicit: Any) -> list[Any]:
        if explicit is None:
            return _categories(series)
        result = list(explicit)
        return result + [value for value in pd.unique(series) if value not in result]

    def _tiles(
        self, spec: ViewSpec, store: DataStore, state: ViewState
    ) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
        frame = _frame(_source(spec, store), state)
        row_id = str(frame.attrs["row_id"])
        x, y, value = self._columns(spec)
        aggregate = str(spec.options.get("aggregate", "mean")).lower()
        tiles: list[dict[str, Any]] = []
        grouped = frame.groupby([x, y], sort=False, dropna=False, observed=True)
        for (x_value, y_value), rows in grouped:
            ids = tuple(rows[row_id].astype(str))
            result = (
                float(len(rows))
                if aggregate == "count"
                else float(
                    getattr(pd.to_numeric(rows[value], errors="coerce"), aggregate)()
                )
            )
            tiles.append(
                {
                    "x": x_value,
                    "y": y_value,
                    "z": result,
                    "ids": ids,
                }
            )
        return frame, tiles

    def _figure(self, spec: ViewSpec, store: DataStore, state: ViewState) -> go.Figure:
        frame, tiles = self._tiles(spec, store, state)
        x, y, value = self._columns(spec)
        xs = self._axis_values(frame[x], spec.options.get("x_order"))
        ys = self._axis_values(frame[y], spec.options.get("y_order"))
        x_at = {item: i for i, item in enumerate(xs)}
        y_at = {item: i for i, item in enumerate(ys)}
        z = np.full((len(ys), len(xs)), np.nan)
        for tile in tiles:
            row, column = y_at[tile["y"]], x_at[tile["x"]]
            z[row, column] = tile["z"]
        options: dict[str, Any] = {
            "x": xs,
            "y": ys,
            "z": z,
            "colorscale": spec.options.get("color_map", "Viridis"),
            "colorbar": {"title": str(spec.options.get("colorbar_title", value))},
            "hovertemplate": f"{x}: %{{x}}<br>{y}: %{{y}}<br>{value}: %{{z:.3g}}<extra></extra>",
        }
        if spec.options.get("color_range") is not None:
            options["zmin"], options["zmax"] = map(
                float, _range(spec.options["color_range"], "color_range")
            )
        if spec.options.get("color_midpoint") is not None:
            options["zmid"] = float(spec.options["color_midpoint"])
        figure = go.Figure(go.Heatmap(**options))
        selected = set(_ids_for(state, highlighted=True) or ())
        marked = [tile for tile in tiles if selected.intersection(tile["ids"])]
        if marked:
            figure.add_trace(
                go.Scatter(
                    x=[tile["x"] for tile in marked],
                    y=[tile["y"] for tile in marked],
                    mode="markers",
                    showlegend=False,
                    hoverinfo="skip",
                    marker={
                        "symbol": "square-open",
                        "size": spec.options.get("selection_size", 26),
                        "color": spec.options.get("selection_color", "#111827"),
                        "line": {"width": 2},
                    },
                )
            )
        figure.update_layout(
            clickmode="event+select",
            dragmode="select",
            margin={"b": 90, "l": 90, "r": 45, "t": 65},
            uirevision=f"{spec.id}:{x}:{y}:{value}",
            xaxis={
                "title": spec.options.get("x_title", x),
                "tickangle": spec.options.get("x_tickangle", -35),
            },
            yaxis={
                "title": spec.options.get("y_title", y),
                "autorange": "reversed"
                if spec.options.get("reverse_y", True)
                else True,
            },
        )
        figure.update_layout(title={"text": str(spec.title or spec.id)})
        return figure

    def decode_event(
        self,
        event: str,
        payload: Mapping[str, Any] | None,
        spec: ViewSpec,
        store: DataStore,
    ) -> MarkMembers | None:
        del event
        source = _source(spec, store)
        x, y, _value = self._columns(spec)
        ids: list[str] = []
        for point in _points(payload):
            if "x" not in point or "y" not in point:
                continue
            mask = source.data[x].eq(point["x"]) & source.data[y].eq(point["y"])
            ids.extend(np.asarray(source.row_ids, dtype=object)[mask].astype(str))
        return _members(ids)


class NetworkAdapter(PlotAdapter):
    """Directed network aggregated from source/target columns of a long table."""

    events = NETWORK_EVENTS
    emits = _IDENTITIES
    accepts = _IDENTITIES

    @staticmethod
    def _columns(spec: ViewSpec) -> tuple[str, str]:
        return (
            str(spec.options.get("source", "source")),
            str(spec.options.get("target", "target")),
        )

    def validate(self, spec: ViewSpec, store: DataStore) -> None:
        frame = _source(spec, store).data
        source, target = self._columns(spec)
        weight = spec.options.get("edge_weight")
        _require(
            frame,
            [source, target] + ([weight] if weight is not None else []),
            spec.id,
        )
        if frame[[source, target]].isna().to_numpy().any():
            raise ValueError(f"Network view {spec.id!r} requires non-null endpoints.")
        aggregate = spec.options.get(
            "aggregate", "sum" if weight is not None else "count"
        )
        if str(aggregate).lower() not in _AGGREGATES:
            raise ValueError(
                "Network aggregate must be count, sum, mean, median, min, or max."
            )
        if weight is None and str(aggregate).lower() != "count":
            raise ValueError(
                "Network aggregation other than count requires `edge_weight`."
            )

    @staticmethod
    def _weight(rows: pd.DataFrame, spec: ViewSpec) -> float:
        column = spec.options.get("edge_weight")
        aggregate = spec.options.get(
            "aggregate", "sum" if column is not None else "count"
        )
        if column is None or aggregate == "count":
            return float(len(rows))
        values = pd.to_numeric(rows[column], errors="coerce")
        result = getattr(values, str(aggregate).lower())()
        number = float(result)
        if not np.isfinite(number):
            raise ValueError("Network aggregate must return a finite number.")
        return number

    def _edges(
        self, spec: ViewSpec, store: DataStore, state: ViewState
    ) -> list[dict[str, Any]]:
        frame = _frame(_source(spec, store), state)
        row_id = str(frame.attrs["row_id"])
        source, target = self._columns(spec)
        edges = []
        grouped = frame.groupby(
            [source, target], sort=True, dropna=False, observed=True
        )
        for (sender, receiver), rows in grouped:
            ids = tuple(rows[row_id].astype(str))
            edges.append(
                {
                    "source": str(sender),
                    "target": str(receiver),
                    "ids": ids,
                    "weight": self._weight(rows, spec),
                }
            )
        return edges

    @staticmethod
    def _widths(edges: Sequence[Mapping[str, Any]]) -> list[float]:
        if not edges:
            return []
        values = np.asarray([edge["weight"] for edge in edges], dtype=float)
        if np.isclose(values.min(), values.max()):
            return [3.5] * len(edges)
        return (1.5 + 6.5 * (values - values.min()) / np.ptp(values)).tolist()

    def elements(
        self, spec: ViewSpec, store: DataStore, state: ViewState | None = None
    ) -> list[dict[str, Any]]:
        resolved = _state(state)
        edges = self._edges(spec, store, resolved)
        nodes = sorted(
            {value for edge in edges for value in (edge["source"], edge["target"])}
        )
        colors = categorical_color_map(nodes, spec.options.get("palette"))
        node_ids = {node: f"node-{position}" for position, node in enumerate(nodes)}
        elements = [
            {
                "data": {
                    "id": node_ids[node],
                    "label": node,
                    "color": colors[node],
                }
            }
            for node in nodes
        ]
        selected = set(_ids_for(resolved, highlighted=True) or ())
        for position, (edge, width) in enumerate(
            zip(edges, self._widths(edges), strict=True)
        ):
            element = {
                "data": {
                    "id": f"edge-{position}",
                    "source": node_ids[edge["source"]],
                    "target": node_ids[edge["target"]],
                    "source_value": edge["source"],
                    "target_value": edge["target"],
                    "weight": edge["weight"],
                    "width": width,
                }
            }
            if selected.intersection(edge["ids"]):
                element["classes"] = "guanaco-linked-selected"
            elements.append(element)
        return elements

    @staticmethod
    def _stylesheet() -> list[dict[str, Any]]:
        node = {
            "background-color": "data(color)",
            "border-color": "#FFFFFF",
            "border-width": 2,
            "color": "#334155",
            "font-family": "Inter, sans-serif",
            "font-size": 12,
            "height": 42,
            "label": "data(label)",
            "text-halign": "center",
            "text-margin-y": -10,
            "text-valign": "bottom",
            "width": 42,
        }
        edge = {
            "curve-style": "bezier",
            "line-color": "#94A3B8",
            "opacity": 0.76,
            "target-arrow-color": "#94A3B8",
            "target-arrow-shape": "triangle",
            "width": "data(width)",
        }
        selected = {
            "line-color": "#D55E00",
            "opacity": 1,
            "target-arrow-color": "#D55E00",
            "width": 9,
            "z-index": 9999,
        }
        return [
            {"selector": "node", "style": node},
            {"selector": "edge", "style": edge},
            {"selector": ".guanaco-linked-selected", "style": selected},
        ]

    def render(
        self,
        spec: ViewSpec,
        store: DataStore,
        state: ViewState | None = None,
        *,
        component_id: str | None = None,
    ) -> Rendered:
        raw_height = spec.options.get("height", "620px")
        height = (
            f"{raw_height}px"
            if isinstance(raw_height, (int, float))
            else str(raw_height)
        )
        layout = {
            "name": spec.options.get("layout", "circle"),
            "animate": False,
            "fit": True,
            "padding": 42,
        }
        style = {
            "backgroundColor": "white",
            "height": f"calc({height} - 54px)",
            "minHeight": "430px",
            "width": "100%",
        }
        graph = cyto.Cytoscape(
            id=component_id or spec.id,
            elements=self.elements(spec, store, state),
            layout=layout,
            stylesheet=self._stylesheet(),
            style=style,
            minZoom=0.2,
            maxZoom=4,
            responsive=True,
            boxSelectionEnabled=False,
            userPanningEnabled=True,
            userZoomingEnabled=True,
        )
        heading_style = {
            "color": "#334155",
            "fontSize": "1.2rem",
            "fontWeight": 600,
            "padding": "16px 18px 6px",
            "textAlign": "center",
        }
        heading = html.Div(
            spec.title or spec.id,
            style=heading_style,
        )
        return html.Div([heading, graph], style={"height": height, "width": "100%"})

    def decode_event(
        self,
        event: str,
        payload: Mapping[str, Any] | None,
        spec: ViewSpec,
        store: DataStore,
    ) -> MarkMembers | None:
        del event
        if (
            not payload
            or payload.get("source_value") is None
            or payload.get("target_value") is None
        ):
            return None
        source = _source(spec, store)
        sender, receiver = self._columns(spec)
        mask = source.data[sender].astype(str).eq(str(payload["source_value"]))
        mask &= source.data[receiver].astype(str).eq(str(payload["target_value"]))
        return _members(np.asarray(source.row_ids, dtype=object)[mask].astype(str))


__all__ = [
    "NETWORK_EVENTS",
    "PLOTLY_EVENTS",
    "NetworkAdapter",
    "TableHeatmapAdapter",
    "TablePlotAdapter",
]
