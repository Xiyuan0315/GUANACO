from __future__ import annotations

from typing import Any

from dash import html
import pandas as pd

try:
    import dash_cytoscape as cyto
except Exception:
    cyto = None


GRN_KEY = "grn"
GRN_REQUIRED_COLUMNS = {"source", "target", "regulation"}
GRN_CONTEXT_ALIASES = ("context", "cell_type", "condition", "group", "category")
GRN_WEIGHT_COLUMN = "weight"
ALL_CONTEXTS = "All"


def has_grn_data(adata: Any) -> bool:
    if adata is None or GRN_KEY not in adata.uns:
        return False
    try:
        grn_df = grn_dataframe(adata)
    except Exception:
        return False
    return not grn_df.empty


def grn_dataframe(adata: Any) -> pd.DataFrame:
    grn = adata.uns[GRN_KEY]
    if not isinstance(grn, pd.DataFrame):
        raise TypeError("adata.uns['grn'] must be a pandas DataFrame.")

    missing = sorted(GRN_REQUIRED_COLUMNS - set(grn.columns))
    if missing:
        raise ValueError(f"adata.uns['grn'] is missing required column(s): {missing}")

    df = grn.copy()
    df["source"] = df["source"].astype(str)
    df["target"] = df["target"].astype(str)
    df["regulation"] = df["regulation"].astype(str).str.strip()
    df = df[(df["source"] != "") & (df["target"] != "")]
    return df


def grn_context_column(grn_df: pd.DataFrame) -> str | None:
    for column in GRN_CONTEXT_ALIASES:
        if column in grn_df.columns:
            return column
    return None


def grn_has_weight(grn_df: pd.DataFrame) -> bool:
    if GRN_WEIGHT_COLUMN not in grn_df.columns:
        return False
    weights = pd.to_numeric(grn_df[GRN_WEIGHT_COLUMN], errors="coerce")
    return weights.notna().any()


def grn_weight_range(grn_df: pd.DataFrame) -> tuple[float | None, float | None]:
    if not grn_has_weight(grn_df):
        return None, None
    weights = pd.to_numeric(grn_df[GRN_WEIGHT_COLUMN], errors="coerce").dropna()
    if weights.empty:
        return None, None
    return float(weights.min()), float(weights.max())


def default_grn_edge_threshold(grn_df: pd.DataFrame) -> float | None:
    min_weight, max_weight = grn_weight_range(grn_df)
    if min_weight is None or max_weight is None:
        return None
    if min_weight <= 0.6 <= max_weight:
        return 0.6
    return min_weight


def grn_weight_step(grn_df: pd.DataFrame) -> float:
    min_weight, max_weight = grn_weight_range(grn_df)
    if min_weight is None or max_weight is None:
        return 0.05
    weight_span = max_weight - min_weight
    if weight_span <= 1:
        return 0.05
    if weight_span <= 10:
        return 0.1
    return 1.0


def grn_context_options(adata: Any) -> tuple[list[dict[str, str]], str | None, bool]:
    if not has_grn_data(adata):
        return [{"label": ALL_CONTEXTS, "value": ALL_CONTEXTS}], None, False

    grn_df = grn_dataframe(adata)
    context_column = grn_context_column(grn_df)
    if context_column is None:
        options = [{"label": ALL_CONTEXTS, "value": ALL_CONTEXTS}]
    else:
        contexts = sorted(str(value) for value in grn_df[context_column].dropna().unique())
        options = [{"label": ALL_CONTEXTS, "value": ALL_CONTEXTS}] + [
            {"label": value, "value": value} for value in contexts
        ]
    return options, context_column, grn_has_weight(grn_df)


def _regulation_sign(value: Any) -> str:
    text = str(value).strip().lower()
    if text in {"+", "activation", "activate", "positive", "up"}:
        return "activation"
    if text in {"-", "repression", "repress", "negative", "down"}:
        return "repression"
    return "unknown"


def _legend():
    item_style = {"display": "inline-flex", "alignItems": "center", "marginRight": "14px", "marginBottom": "6px"}
    swatch_style = {
        "display": "inline-block",
        "width": "12px",
        "height": "12px",
        "borderRadius": "50%",
        "marginRight": "6px",
        "border": "1px solid rgba(0,0,0,0.15)",
    }
    return html.Div(
        [
            html.Div(
                [
                    html.Div([html.Span(style={**swatch_style, "backgroundColor": "#E07A5F"}), html.Span("Source")], style=item_style),
                    html.Div([html.Span(style={**swatch_style, "backgroundColor": "#3D5A80"}), html.Span("Target")], style=item_style),
                    html.Div([html.Span(style={**swatch_style, "backgroundColor": "#2A9D8F"}), html.Span("Activation (+)")], style=item_style),
                    html.Div([html.Span(style={**swatch_style, "backgroundColor": "#C0392B"}), html.Span("Repression (-)")], style=item_style),
                ],
                style={"display": "flex", "flexWrap": "wrap"},
            ),
        ],
        style={"padding": "8px 10px 0 10px"},
    )


def _empty_grn_component(message: str):
    return html.Div(
        message,
        style={
            "padding": "16px",
            "color": "#7A1F1F",
            "backgroundColor": "#FDECEC",
            "border": "1px solid #F5C2C7",
            "borderRadius": "6px",
            "margin": "8px",
        },
    )


def _filter_grn_edges(
    grn_df: pd.DataFrame,
    selected_context: str,
    context_column: str | None,
    edge_threshold: float | None,
) -> pd.DataFrame:
    filtered = grn_df
    if context_column is not None and selected_context and selected_context != ALL_CONTEXTS:
        filtered = filtered[filtered[context_column].astype(str) == str(selected_context)]

    if edge_threshold is not None and grn_has_weight(filtered):
        weights = pd.to_numeric(filtered[GRN_WEIGHT_COLUMN], errors="coerce")
        filtered = filtered[weights >= float(edge_threshold)]
    return filtered


def _edge_widths(filtered: pd.DataFrame) -> pd.Series:
    if not grn_has_weight(filtered):
        return pd.Series([3.0] * len(filtered), index=filtered.index)

    weights = pd.to_numeric(filtered[GRN_WEIGHT_COLUMN], errors="coerce")
    finite_weights = weights.dropna()
    if finite_weights.empty:
        return pd.Series([3.0] * len(filtered), index=filtered.index)

    min_weight = float(finite_weights.min())
    max_weight = float(finite_weights.max())
    if min_weight == max_weight:
        return pd.Series([4.5 if pd.notna(value) else 3.0 for value in weights], index=filtered.index)

    scaled = 1.5 + ((weights - min_weight) / (max_weight - min_weight)) * 5.0
    return scaled.fillna(3.0)


def grn_graph(
    adata: Any,
    *,
    selected_context=ALL_CONTEXTS,
    edge_threshold: float | None = None,
    layout_name="cose",
    node_font_size=12,
):
    """Backend-agnostic GRN graph for cytoscape.js.

    Returns ``{"elements", "stylesheet", "layout", "legend", "directed"}`` so it
    can drive either the Dash ``cyto.Cytoscape`` component or the ipycytoscape
    notebook widget (``gc.pl.grn``). Raises ``ValueError`` when there is nothing
    to render.
    """
    grn_df = grn_dataframe(adata)
    context_column = grn_context_column(grn_df)
    filtered = _filter_grn_edges(grn_df, selected_context, context_column, edge_threshold)
    if filtered.empty:
        raise ValueError("No GRN edges match the current filters.")

    source_nodes = set(filtered["source"].astype(str))
    nodes = {}
    for source in filtered["source"].astype(str):
        nodes[source] = {"id": source, "label": source, "node_type": "source"}
    for target in filtered["target"].astype(str):
        node_type = "source" if target in source_nodes else "target"
        nodes.setdefault(target, {"id": target, "label": target, "node_type": node_type})

    widths = _edge_widths(filtered)
    edge_elements = []
    for row_index, row in filtered.iterrows():
        source = str(row["source"])
        target = str(row["target"])
        sign = _regulation_sign(row["regulation"])
        weight = row.get(GRN_WEIGHT_COLUMN)
        edge_elements.append(
            {
                "data": {
                    "id": f"{source}-{target}-{row_index}",
                    "source": source,
                    "target": target,
                    "weight": None if pd.isna(weight) else float(weight) if GRN_WEIGHT_COLUMN in filtered.columns else None,
                    "sign": sign,
                    "width": float(widths.loc[row_index]),
                }
            }
        )

    node_elements = [{"data": node_data} for node_data in nodes.values()]

    stylesheet = [
        {
            "selector": "node",
            "style": {
                "label": "data(label)",
                "font-size": node_font_size,
                "font-weight": "bold",
                "color": "#22333B",
                "text-valign": "bottom",
                "text-halign": "center",
                "text-margin-y": 8,
                "border-width": 1,
                "border-color": "white",
            },
        },
        {
            "selector": 'node[node_type = "source"]',
            "style": {
                "shape": "round-rectangle",
                "width": 54,
                "height": 34,
                "background-color": "#E07A5F",
            },
        },
        {
            "selector": 'node[node_type = "target"]',
            "style": {
                "shape": "ellipse",
                "width": 34,
                "height": 34,
                "background-color": "#3D5A80",
                "color": "#22333B",
            },
        },
        {
            "selector": "edge",
            "style": {
                "curve-style": "bezier",
                "target-arrow-shape": "triangle",
                "arrow-scale": 1.2,
                "width": "data(width)",
                "opacity": 0.78,
            },
        },
        {
            "selector": 'edge[sign = "activation"]',
            "style": {
                "line-color": "#2A9D8F",
                "target-arrow-color": "#2A9D8F",
            },
        },
        {
            "selector": 'edge[sign = "repression"]',
            "style": {
                "line-color": "#C0392B",
                "target-arrow-color": "#C0392B",
                "target-arrow-shape": "tee",
            },
        },
        {
            "selector": 'edge[sign = "unknown"]',
            "style": {
                "line-color": "#7A7A7A",
                "target-arrow-color": "#7A7A7A",
            },
        },
    ]

    return {
        "elements": node_elements + edge_elements,
        "stylesheet": stylesheet,
        "layout": {"name": layout_name, "fit": True, "padding": 40},
        "legend": {
            "kind": "categorical",
            "title": "Gene regulatory network",
            "entries": [
                {"label": "Source", "color": "#E07A5F"},
                {"label": "Target", "color": "#3D5A80"},
                {"label": "Activation (+)", "color": "#2A9D8F"},
                {"label": "Repression (−)", "color": "#C0392B"},
            ],
        },
        "directed": True,
    }


def build_grn_cytoscape(
    adata: Any,
    *,
    component_id="grn-cytoscape-view",
    selected_context=ALL_CONTEXTS,
    edge_threshold: float | None = None,
    layout_name="cose",
):
    """Dash GRN component (wraps :func:`grn_graph`)."""
    if cyto is None:
        return _empty_grn_component("dash-cytoscape is not installed. Add dash-cytoscape to render GRN plots.")
    try:
        graph = grn_graph(
            adata,
            selected_context=selected_context,
            edge_threshold=edge_threshold,
            layout_name=layout_name,
        )
    except Exception as exc:
        return _empty_grn_component(str(exc))

    cytoscape_component = cyto.Cytoscape(
        id=component_id,
        elements=graph["elements"],
        layout=graph["layout"],
        stylesheet=graph["stylesheet"],
        minZoom=0.2,
        maxZoom=3.0,
        responsive=True,
        style={"width": "100%", "height": "100%", "flexGrow": "1", "backgroundColor": "white"},
    )

    return html.Div(
        [
            _legend(),
            cytoscape_component,
        ],
        style={"height": "100%", "width": "100%", "display": "flex", "flexDirection": "column"},
    )
