"""Linked circle and dot plots for precomputed ligand–receptor results."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import dash_table, html

from guanaco.data.ligand_receptor import (
    ligand_receptor_frame,
    metric_direction,
)
from guanaco.utils.cytoscape import adjacent_highlight_stylesheet

try:
    import dash_cytoscape as cyto
except Exception:
    cyto = None


MAX_DISPLAYED_INTERACTIONS = 20


def _empty(message: str, *, error: bool = False):
    class_name = "lr-empty-state lr-empty-state--error" if error else "lr-empty-state"
    return html.Div(message, className=class_name)


def empty_ligand_receptor_figure(message: str) -> go.Figure:
    figure = go.Figure()
    figure.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font={"color": "#7c838a", "size": 14},
    )
    figure.update_xaxes(visible=False)
    figure.update_yaxes(visible=False)
    figure.update_layout(
        margin={"l": 30, "r": 20, "t": 45, "b": 30},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    return figure


def build_ligand_receptor_view(
    payload: Mapping[str, Any],
    *,
    magnitude: str,
    specificity: str | None = None,
    magnitude_range: list[float] | None = None,
    specificity_range: list[float] | None = None,
    senders: list[str] | None = None,
    receivers: list[str] | None = None,
    node_colors: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Filter a canonical interaction table and aggregate circle edges."""
    frame = ligand_receptor_frame(payload)
    if magnitude not in frame.columns:
        raise ValueError("Select a numeric magnitude metric.")
    if specificity not in frame.columns:
        specificity = None

    if senders:
        frame = frame[frame["source"].isin([str(value) for value in senders])]
    if receivers:
        frame = frame[frame["target"].isin([str(value) for value in receivers])]
    frame = frame.copy()
    frame[magnitude] = pd.to_numeric(frame[magnitude], errors="coerce")
    frame = frame[frame[magnitude].notna()]
    if magnitude_range and len(magnitude_range) == 2:
        frame = frame[
            frame[magnitude].between(
                float(min(magnitude_range)),
                float(max(magnitude_range)),
                inclusive="both",
            )
        ]
    if specificity:
        frame[specificity] = pd.to_numeric(frame[specificity], errors="coerce")
        if specificity_range and len(specificity_range) == 2:
            frame = frame[
                frame[specificity].between(
                    float(min(specificity_range)),
                    float(max(specificity_range)),
                    inclusive="both",
                )
            ]
    if frame.empty:
        raise ValueError("No interactions remain after applying the current filters.")

    ascending = metric_direction(payload, magnitude) == "lower"
    frame = frame.sort_values(
        magnitude,
        ascending=ascending,
        kind="stable",
    ).reset_index(drop=True)

    pair_ids: dict[tuple[str, str], str] = {}
    edges = []
    for index, ((source, target), rows) in enumerate(
        frame.groupby(["source", "target"], sort=False)
    ):
        edge_id = f"lr-edge-{index}"
        source = str(source)
        target = str(target)
        pair_ids[(source, target)] = edge_id
        edges.append(
            {
                "id": edge_id,
                "source": source,
                "target": target,
                "pair_count": int(len(rows)),
                "magnitude": float(rows[magnitude].median()),
            }
        )

    frame["interaction_id"] = [
        f"lr-interaction-{index}" for index in range(len(frame))
    ]
    frame["edge_id"] = [
        pair_ids[(str(source), str(target))]
        for source, target in frame[["source", "target"]].itertuples(
            index=False,
            name=None,
        )
    ]
    groups = sorted(
        {
            str(value)
            for column in ("source", "target")
            for value in frame[column].dropna()
        },
        key=str.casefold,
    )
    interactions = json.loads(frame.to_json(orient="records"))
    return {
        "name": payload.get("name", "Ligand–receptor results"),
        "format": payload.get("format", "interaction table"),
        "groups": groups,
        "node_colors": dict(node_colors or {}),
        "edges": edges,
        "interactions": interactions,
        "magnitude": magnitude,
        "magnitude_direction": metric_direction(payload, magnitude),
        "specificity": specificity,
        "specificity_direction": (
            metric_direction(payload, specificity) if specificity else None
        ),
    }


def _edge_width(value: float, values: list[float]) -> float:
    if not values or min(values) == max(values):
        return 3.0
    return 1.2 + 5.0 * (value - min(values)) / (max(values) - min(values))


def _elements(view: Mapping[str, Any]) -> list[dict[str, Any]]:
    colors = view.get("node_colors", {})
    nodes = [
        {
            "data": {
                "id": group,
                "label": group,
                "color": colors.get(group, "#808080"),
            }
        }
        for group in view.get("groups", [])
    ]
    counts = [float(edge["pair_count"]) for edge in view.get("edges", [])]
    edges = []
    for edge in view.get("edges", []):
        data = dict(edge)
        data["width"] = _edge_width(float(data["pair_count"]), counts)
        edges.append({"data": data})
    return nodes + edges


def ligand_receptor_stylesheet() -> list[dict[str, Any]]:
    return [
        {
            "selector": "node",
            "style": {
                "label": "data(label)",
                "background-color": "data(color)",
                "border-width": 0,
                "color": "#343a40",
                "font-size": 11,
                "font-weight": 500,
                "height": 42,
                "text-halign": "center",
                "text-margin-y": 7,
                "text-valign": "bottom",
                "width": 42,
            },
        },
        {
            "selector": "edge",
            "style": {
                "curve-style": "bezier",
                "line-color": "#8f969d",
                "opacity": 0.72,
                "target-arrow-color": "#8f969d",
                "target-arrow-shape": "triangle",
                "arrow-scale": 1.0,
                "width": "data(width)",
            },
        },
    ]


def _selector_value(value: Any) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def ligand_receptor_highlight_stylesheet(
    view: Mapping[str, Any],
    *,
    node_id: str | None = None,
    edge_id: str | None = None,
) -> list[dict[str, Any]]:
    if node_id:
        return adjacent_highlight_stylesheet(
            ligand_receptor_stylesheet(),
            view.get("edges", []),
            node_id,
        )
    if not edge_id:
        return ligand_receptor_stylesheet()
    selected = next(
        (edge for edge in view.get("edges", []) if edge.get("id") == edge_id),
        None,
    )
    if selected is None:
        return ligand_receptor_stylesheet()
    overlays = [
        {"selector": "node", "style": {"opacity": 0.12}},
        {"selector": "edge", "style": {"opacity": 0.05}},
        {
            "selector": f'edge[id = "{_selector_value(edge_id)}"]',
            "style": {
                "line-color": "#202020",
                "opacity": 1,
                "target-arrow-color": "#202020",
                "z-index": 9999,
            },
        },
    ]
    for group in (selected["source"], selected["target"]):
        overlays.append(
            {
                "selector": f'node[id = "{_selector_value(group)}"]',
                "style": {"opacity": 1, "z-index": 9999},
            }
        )
    return ligand_receptor_stylesheet() + overlays


def build_ligand_receptor_network(
    view: Mapping[str, Any],
    *,
    component_id: str,
):
    if not view.get("edges"):
        return _empty("No precomputed interactions passed the current filters.")
    if cyto is None:
        return _empty("dash-cytoscape is required to render this network.", error=True)
    return html.Div(
        [
            html.Div(
                [
                    html.Span("Nodes = cell groups"),
                    html.Span("Arrows = sender → receiver"),
                    html.Span("Edge width = interaction count"),
                ],
                className="lr-network-legend",
            ),
            cyto.Cytoscape(
                id=component_id,
                elements=_elements(view),
                layout={
                    "name": "circle",
                    "fit": True,
                    "padding": 42,
                    "animate": False,
                },
                stylesheet=ligand_receptor_stylesheet(),
                minZoom=0.2,
                maxZoom=4,
                responsive=True,
                style={
                    "height": "100%",
                    "width": "100%",
                    "backgroundColor": "white",
                },
            ),
        ],
        className="lr-network-content",
    )


def _marker_sizes(values: np.ndarray, direction: str) -> np.ndarray:
    finite = np.isfinite(values)
    if not finite.any():
        return np.full(len(values), 11.0)
    transformed = values.astype(float, copy=True)
    if direction == "lower":
        positive = transformed[finite & (transformed > 0)]
        floor = float(positive.min() / 10) if positive.size else 1e-12
        transformed = -np.log10(np.maximum(transformed, floor))
    low = float(np.nanmin(transformed[finite]))
    high = float(np.nanmax(transformed[finite]))
    if np.isclose(low, high):
        return np.full(len(values), 13.0)
    normalized = (transformed - low) / (high - low)
    normalized[~finite] = 0
    return 8.0 + 14.0 * normalized


def metric_slider_settings(
    payload: Mapping[str, Any] | None,
    metric: str | None,
) -> dict[str, Any]:
    """Return stable linear RangeSlider settings for one numeric metric."""
    if not payload or not metric:
        return {
            "min": 0.0,
            "max": 1.0,
            "step": 0.005,
            "value": [0.0, 1.0],
            "marks": {0.0: "0", 1.0: "1"},
        }
    frame = ligand_receptor_frame(payload)
    if metric not in frame:
        return metric_slider_settings(None, None)
    values = pd.to_numeric(frame[metric], errors="coerce").to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if not finite.size:
        return metric_slider_settings(None, None)
    low = float(np.min(finite))
    high = float(np.max(finite))
    if np.isclose(low, high):
        padding = max(abs(low) * 0.05, 1.0)
        slider_low = low - padding
        slider_high = high + padding
    else:
        slider_low = low
        slider_high = high
    midpoint = (slider_low + slider_high) / 2
    span = slider_high - slider_low
    return {
        "min": slider_low,
        "max": slider_high,
        "step": span / 200,
        "value": [slider_low, slider_high],
        "marks": {
            slider_low: _format_metric(slider_low),
            midpoint: _format_metric(midpoint),
            slider_high: _format_metric(slider_high),
        },
    }


def _size_legend(
    values: np.ndarray,
    direction: str,
) -> tuple[np.ndarray, list[tuple[float, float]]]:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return _marker_sizes(values, direction), []
    legend_values = np.unique(np.quantile(finite, [0.0, 0.5, 1.0]))
    combined = np.concatenate([values, legend_values])
    combined_sizes = _marker_sizes(combined, direction)
    point_sizes = combined_sizes[: len(values)]
    legend_sizes = combined_sizes[len(values) :]
    return point_sizes, [
        (float(value), float(size))
        for value, size in zip(legend_values, legend_sizes, strict=True)
    ]


def plot_ligand_receptor_dotplot(
    view: Mapping[str, Any] | None,
    *,
    selected_node: str | None = None,
    selected_edge: str | None = None,
    colorscale: Any = "Viridis",
) -> go.Figure:
    if not view or not view.get("interactions"):
        return empty_ligand_receptor_figure(
            "Load a precomputed result to view ligand–receptor interactions."
        )
    frame = pd.DataFrame.from_records(view["interactions"])
    if selected_edge:
        frame = frame[frame["edge_id"].eq(selected_edge)]
    elif selected_node:
        frame = frame[
            frame["source"].eq(selected_node)
            | frame["target"].eq(selected_node)
        ]
    if frame.empty:
        return empty_ligand_receptor_figure(
            "No interactions match the selected cell group or edge."
        )
    interaction_count = len(frame)
    frame = frame.head(MAX_DISPLAYED_INTERACTIONS).copy()

    magnitude = str(view["magnitude"])
    specificity = view.get("specificity")
    color = pd.to_numeric(frame[magnitude], errors="coerce").to_numpy(dtype=float)
    if specificity and specificity in frame.columns:
        specificity_values = pd.to_numeric(
            frame[specificity],
            errors="coerce",
        ).to_numpy(dtype=float)
        sizes, size_legend = _size_legend(
            specificity_values,
            str(view.get("specificity_direction", "higher")),
        )
    else:
        specificity_values = np.full(len(frame), np.nan)
        sizes = np.full(len(frame), 11.0)
        size_legend = []

    customdata = np.column_stack(
        [
            frame["interaction_id"].astype(str),
            frame["edge_id"].astype(str),
            frame["source"].astype(str),
            frame["target"].astype(str),
            frame["ligand"].astype(str),
            frame["receptor"].astype(str),
            color,
            specificity_values,
        ]
    )
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=frame["source"].astype(str) + " → " + frame["target"].astype(str),
            y=frame["ligand"].astype(str) + " → " + frame["receptor"].astype(str),
            mode="markers",
            marker={
                "color": color,
                "colorscale": colorscale,
                "reversescale": view.get("magnitude_direction") == "lower",
                "size": sizes,
                "opacity": 0.88,
                "line": {"width": 0},
                "colorbar": {
                    "title": magnitude,
                    "thickness": 14,
                    "x": 1.02,
                    "xanchor": "left",
                    "y": 0.7,
                    "yanchor": "middle",
                    "len": 0.5,
                },
            },
            customdata=customdata,
            hovertemplate=(
                "%{customdata[2]} → %{customdata[3]}<br>"
                "%{customdata[4]} → %{customdata[5]}<br>"
                f"{magnitude}: %{{customdata[6]:.4g}}"
                + (
                    f"<br>{specificity}: %{{customdata[7]:.4g}}"
                    if specificity
                    else ""
                )
                + "<extra></extra>"
            ),
            showlegend=False,
        )
    )
    for value, size in size_legend:
        figure.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker={"color": "#808080", "size": size, "line": {"width": 0}},
                name=_format_metric(value),
                hoverinfo="skip",
                legendgroup="size",
                showlegend=True,
            )
        )
    figure.update_layout(
        title={
            "text": (
                f"Ligand–receptor interactions · Top "
                f"{MAX_DISPLAYED_INTERACTIONS} of {interaction_count}"
                if interaction_count > MAX_DISPLAYED_INTERACTIONS
                else "Ligand–receptor interactions"
            ),
            "x": 0.02,
        },
        clickmode="event+select",
        margin={
            "l": 160,
            "r": 185 if size_legend else 95,
            "t": 58,
            "b": 140,
        },
        legend=(
            {
                "title": {"text": str(specificity)},
                "x": 1.02,
                "xanchor": "left",
                "y": 0.36,
                "yanchor": "top",
            }
            if size_legend
            else None
        ),
        paper_bgcolor="white",
        plot_bgcolor="white",
        xaxis={"title": "Sender → receiver", "tickangle": -40},
        yaxis={"title": "Ligand → receptor", "autorange": "reversed"},
        uirevision=(
            f"lr-dot-{view.get('name')}-{magnitude}-{specificity}-"
            f"{selected_node}-{selected_edge}"
        ),
    )
    return figure


def ligand_receptor_detail(
    view: Mapping[str, Any] | None,
    *,
    node_id: str | None = None,
    edge_id: str | None = None,
    interaction_id: str | None = None,
):
    if not view or not view.get("interactions"):
        return _empty("Load a result to inspect ligand–receptor pairs.")
    rows = list(view["interactions"])
    if interaction_id:
        rows = [
            row for row in rows if row.get("interaction_id") == interaction_id
        ]
    elif edge_id:
        rows = [row for row in rows if row.get("edge_id") == edge_id]
    elif node_id:
        rows = [
            row
            for row in rows
            if str(row.get("source")) == str(node_id)
            or str(row.get("target")) == str(node_id)
        ]
    if not rows:
        return _empty("The selected interaction is no longer available.")

    magnitude = str(view["magnitude"])
    specificity = view.get("specificity")
    headers = ["Sender", "Receiver", "Ligand", "Receptor", magnitude]
    if specificity:
        headers.append(str(specificity))
    table_rows = []
    for row in rows:
        values = [
            row.get("source"),
            row.get("target"),
            row.get("ligand"),
            row.get("receptor"),
            _format_metric(row.get(magnitude)),
        ]
        if specificity:
            values.append(_format_metric(row.get(specificity)))
        table_rows.append(dict(zip(headers, values, strict=True)))
    if interaction_id or edge_id or node_id:
        title = (
            f"{len(rows):,} selected interaction"
            + ("s" if len(rows) != 1 else "")
        )
    else:
        title = f"{len(rows):,} interaction" + ("s" if len(rows) != 1 else "")
    return html.Div(
        [
            html.Div(
                title,
                className="lr-detail-title",
            ),
            html.Div(
                dash_table.DataTable(
                    columns=[{"name": header, "id": header} for header in headers],
                    data=table_rows,
                    page_action="native",
                    page_size=25,
                    sort_action="native",
                    style_table={"overflowX": "auto"},
                    style_cell={
                        "border": "0",
                        "borderBottom": "1px solid #eceeef",
                        "fontFamily": "inherit",
                        "fontSize": "0.84rem",
                        "padding": "7px 9px",
                        "textAlign": "left",
                    },
                    style_header={
                        "backgroundColor": "#f7f7f7",
                        "color": "#555b61",
                        "fontWeight": 600,
                    },
                ),
                className="lr-detail-table-wrap",
            ),
        ],
        className="lr-detail-content",
    )


def _format_metric(value: Any) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "—"
    return f"{numeric:.4g}" if np.isfinite(numeric) else "—"
