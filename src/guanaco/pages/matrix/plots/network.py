"""Cytoscape rendering for database-derived exploratory networks."""

from __future__ import annotations

from typing import Any

from dash import html
from plotly.colors import sample_colorscale

from guanaco.utils.colors import resolve_continuous_colorscale
from guanaco.utils.cytoscape import adjacent_highlight_stylesheet

try:
    import dash_cytoscape as cyto
except Exception:
    cyto = None


NETWORK_LABELS = {
    "ppi": "PPI",
    "tf-gene": "TF–gene · TRRUST",
    "metabolite": "Metabolite–gene · KEGG",
    "mirna": "miRNA–mRNA · miRTarBase",
}

NETWORK_LAYOUTS = {"cose", "breadthfirst", "circle", "concentric"}


def _empty_component(message: str, *, error: bool = False):
    return html.Div(
        message,
        className="network-empty-state network-empty-state--error"
        if error
        else "network-empty-state",
    )


def network_stylesheet() -> list[dict[str, Any]]:
    mirna_color = "#7b4f9d"
    stylesheet = [
        {
            "selector": "node",
            "style": {
                "label": "data(label)",
                "font-size": 11,
                "font-weight": 500,
                "min-zoomed-font-size": 24,
                "color": "#343a40",
                "text-valign": "bottom",
                "text-halign": "center",
                "text-margin-y": 7,
                "width": 10,
                "height": 10,
                "background-color": "#c7c9cc",
                "border-width": 0,
            },
        },
        {
            "selector": "node.network-measured-node",
            "style": {
                "background-color": "data(expression_color)",
            },
        },
        {
            "selector": "node.network-query-node",
            "style": {
                "border-width": 3,
                "border-color": "#4a4a4a",
                "min-zoomed-font-size": 0,
            },
        },
        {
            "selector": "node.network-hub-node",
            "style": {
                "font-weight": 600,
                "min-zoomed-font-size": 14,
            },
        },
        {
            "selector": 'node[entity_type = "tf"]',
            "style": {"shape": "round-rectangle", "width": 48, "height": 32},
        },
        {
            "selector": 'node[entity_type = "metabolite"]',
            "style": {
                "shape": "diamond",
                "background-color": "#d8b365",
                "width": 38,
                "height": 38,
            },
        },
        {
            "selector": 'node[entity_type = "mirna"], node[entity_type = "microrna"]',
            "style": {
                "shape": "triangle",
                "background-color": mirna_color,
            },
        },
        {
            "selector": "edge",
            "style": {
                "curve-style": "bezier",
                "line-color": "#9a9a9a",
                "target-arrow-color": "#9a9a9a",
                "target-arrow-shape": "none",
                "width": "data(width)",
                "opacity": 0.7,
            },
        },
        {
            "selector": "edge[directed = true]",
            "style": {"target-arrow-shape": "triangle", "arrow-scale": 1.05},
        },
        {
            "selector": 'edge[effect = "activation"]',
            "style": {"line-color": "#2a9d8f", "target-arrow-color": "#2a9d8f"},
        },
        {
            "selector": 'edge[effect = "inhibition"]',
            "style": {
                "line-color": "#c0392b",
                "target-arrow-color": "#c0392b",
                "target-arrow-shape": "tee",
            },
        },
        {
            "selector": "edge.network-undirected-edge",
            "style": {"target-arrow-shape": "none"},
        },
    ]
    stylesheet.append(
        {
            "selector": "node[degree >= 0]",
            "style": {
                "width": "data(node_size)",
                "height": "data(node_size)",
            },
        }
    )
    return stylesheet


def network_highlight_stylesheet(
    graph: dict[str, Any],
    selected_node: str,
) -> list[dict[str, Any]]:
    return adjacent_highlight_stylesheet(
        network_stylesheet(), graph.get("edges", []), selected_node
    )


def _edge_width(
    score: Any,
    scores: list[float],
    *,
    network_type: str | None = None,
) -> float:
    try:
        value = float(score)
    except (TypeError, ValueError):
        return 2.0
    if network_type == "ppi":
        # STRING confidence is an absolute 0–1 scale. Guanaco requests edges
        # from 0.4 upward, so keep their widths stable across queries, views,
        # and slider thresholds instead of renormalizing each displayed graph.
        normalized = (min(1.0, max(0.4, value)) - 0.4) / 0.6
        return 1.2 + 4.0 * normalized
    if not scores:
        return 2.0
    minimum, maximum = min(scores), max(scores)
    if minimum == maximum:
        return 3.0
    return 1.2 + 4.0 * (value - minimum) / (maximum - minimum)


def _expression_color(value, minimum, maximum, colorscale):
    if value is None:
        return "#c7c9cc"
    normalized = (
        0.5 if minimum == maximum else (float(value) - minimum) / (maximum - minimum)
    )
    return sample_colorscale(colorscale, [float(max(0.0, min(1.0, normalized)))])[0]


def _hub_node_ids(nodes: list[dict[str, Any]]) -> set[str]:
    if not nodes:
        return set()
    hub_count = min(8, max(2, round(len(nodes) ** 0.5)))
    candidates = [node for node in nodes if not node.get("is_query")]
    candidates.sort(
        key=lambda node: (
            -int(node.get("degree") or 0),
            str(node.get("id", "")).casefold(),
        )
    )
    return {
        str(node["id"])
        for node in candidates[:hub_count]
        if int(node.get("degree") or 0) > 0
    }


def _elements(
    graph: dict[str, Any],
    continuous_colormap: str = "viridis",
) -> list[dict[str, Any]]:
    colorscale = resolve_continuous_colorscale(continuous_colormap or "viridis")
    minimum = float(graph.get("expression_min", 0.0))
    maximum = float(graph.get("expression_max", 1.0))
    hub_ids = _hub_node_ids(graph.get("nodes", []))
    nodes = []
    for node in graph.get("nodes", []):
        data = dict(node)
        data["expression_color"] = _expression_color(
            data.get("expression"), minimum, maximum, colorscale
        )
        element = {"data": data}
        classes = []
        if data.get("measured") is True:
            classes.append("network-measured-node")
        if data.get("is_query") is True:
            classes.append("network-query-node")
        if str(data.get("id")) in hub_ids:
            classes.append("network-hub-node")
        if classes:
            element["classes"] = " ".join(classes)
        nodes.append(element)
    scores = [
        float(edge["score"])
        for edge in graph.get("edges", [])
        if edge.get("score") is not None
    ]
    edges = []
    for edge in graph.get("edges", []):
        data = dict(edge)
        data["width"] = _edge_width(
            data.get("score"),
            scores,
            network_type=graph.get("network_type"),
        )
        element = {"data": data}
        if not data.get("directed"):
            element["classes"] = "network-undirected-edge"
        edges.append(element)
    return nodes + edges


def _legend(graph: dict[str, Any], continuous_colormap: str = "viridis"):
    network_type = graph.get("network_type")
    network_label = NETWORK_LABELS.get(network_type, "Network")
    if network_type == "ppi":
        network_label = (
            f"{network_label} · STRING {graph.get('ppi_network_type', 'functional')}"
        )
    type_items = []
    node_types = {node.get("entity_type") for node in graph.get("nodes", [])}
    if "metabolite" in node_types:
        type_items.append(html.Span("◆ Metabolite", className="network-legend-item"))
    type_items.append(
        html.Span(
            "Node size = network degree",
            className="network-legend-item",
        )
    )
    colorscale = resolve_continuous_colorscale(continuous_colormap or "viridis")
    gradient_colors = sample_colorscale(colorscale, [index / 8 for index in range(9)])
    gradient_style = {
        "background": f"linear-gradient(90deg, {', '.join(gradient_colors)})"
    }
    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        f"{network_label} · mean expression",
                        className="network-legend-title",
                    ),
                    html.Span(
                        className="network-expression-gradient",
                        style=gradient_style,
                    ),
                    html.Span(f"{graph['expression_min']:.3g}"),
                    html.Span("→"),
                    html.Span(f"{graph['expression_max']:.3g}"),
                ],
                className="network-expression-legend",
            ),
            html.Div(
                [
                    html.Span("Outlined = query gene", className="network-legend-item"),
                    html.Span("Grey = not measured", className="network-legend-item"),
                    *type_items,
                ],
                className="network-type-legend",
            ),
        ],
        className="network-legend",
    )


def network_layout(
    layout_name: str | None,
    network_type: str | None = None,
) -> dict[str, Any]:
    resolved = layout_name if layout_name in NETWORK_LAYOUTS else "cose"
    layout = {"name": resolved, "fit": True, "padding": 36, "animate": False}
    if resolved == "cose":
        if network_type == "tf-gene":
            layout.update(
                {
                    "nodeRepulsion": 18000,
                    "idealEdgeLength": 130,
                    "gravity": 0.08,
                    "componentSpacing": 120,
                    "numIter": 2500,
                    "randomize": True,
                }
            )
        else:
            layout.update({"nodeRepulsion": 6500, "idealEdgeLength": 90})
    elif resolved == "breadthfirst":
        layout.update({"directed": True, "spacingFactor": 1.25})
    return layout


def build_key_regulator_component(enrichment: dict[str, Any], *, max_rows: int = 25):
    method = html.Div(
        (
            "One-sided Fisher exact test · BH-FDR · "
            f"background: {enrichment['background_size']:,} measurable RNA genes · "
            f"query: {enrichment['valid_query_count']}/{enrichment['input_query_count']} measurable"
        ),
        className="network-enrichment-method",
    )
    results = enrichment.get("results", [])
    if not results:
        return html.Div(
            [
                method,
                html.Div(enrichment.get("reason", "No enriched regulator found.")),
            ],
            className="network-enrichment-empty",
        )

    rows = []
    for result in results[:max_rows]:
        rows.append(
            html.Tr(
                [
                    html.Td(result["tf"], className="network-enrichment-tf"),
                    html.Td(str(result["overlap_count"])),
                    html.Td(str(result["target_count"])),
                    html.Td(f"{result['p_value']:.3g}"),
                    html.Td(f"{result['fdr']:.3g}"),
                    html.Td(", ".join(result["overlap_genes"])),
                ],
                className=(
                    "network-enrichment-significant" if result["fdr"] < 0.05 else ""
                ),
            )
        )
    footer = None
    if len(results) > max_rows:
        footer = html.Div(
            f"Showing the top {max_rows} of {len(results)} tested key regulators.",
            className="network-enrichment-footnote",
        )
    return html.Div(
        [
            method,
            html.Div(
                html.Table(
                    [
                        html.Thead(
                            html.Tr(
                                [
                                    html.Th("TF"),
                                    html.Th("Overlap"),
                                    html.Th("Targets"),
                                    html.Th("P value"),
                                    html.Th("FDR"),
                                    html.Th("Overlapping genes"),
                                ]
                            )
                        ),
                        html.Tbody(rows),
                    ],
                    className="network-enrichment-table",
                ),
                className="network-enrichment-table-wrap",
            ),
            footer,
        ],
        className="network-enrichment-results",
    )


def build_network_cytoscape(
    graph: dict[str, Any],
    *,
    component_id: str,
    layout_name: str = "cose",
    continuous_colormap: str = "viridis",
):
    if cyto is None:
        return _empty_component(
            "dash-cytoscape is required to render Network.", error=True
        )
    if not graph.get("nodes"):
        return _empty_component("No database interactions were found for these genes.")
    return html.Div(
        [
            _legend(graph, continuous_colormap),
            cyto.Cytoscape(
                id=component_id,
                elements=_elements(graph, continuous_colormap),
                layout=network_layout(layout_name, graph.get("network_type")),
                stylesheet=network_stylesheet(),
                minZoom=0.2,
                maxZoom=4.0,
                responsive=True,
                style={"width": "100%", "height": "100%", "backgroundColor": "white"},
            ),
        ],
        className="network-graph-content",
    )
