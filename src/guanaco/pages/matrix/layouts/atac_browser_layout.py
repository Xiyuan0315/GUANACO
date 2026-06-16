from dash import dcc, html
import dash_bootstrap_components as dbc

from guanaco.pages.matrix.plots.atac_browser import (
    default_region,
    format_locus,
    has_genomic_peak_features,
)
from guanaco.utils.ui_helpers import graph_flex_container, responsive_graph_grid


def generate_atac_browser_layout(adata, prefix, gene_annotation_path=None):
    if not has_genomic_peak_features(adata):
        return dbc.Alert(
            "Peak Browser needs peak variables named like chr1:1000-1500 or var columns chrom/start/end.",
            color="secondary",
            style={"marginTop": "10px"},
        )

    region = default_region(adata)
    locus = format_locus(str(region["chrom"]), int(region["start"]), int(region["end"]))
    region = {**region, "nav": locus}
    return html.Div(
        [
            dcc.Store(id=f"{prefix}-atac-browser-region-store", data=region),
            html.Div(
                [
                    html.Label("Gene / locus", style={"fontWeight": "bold", "marginBottom": "5px", "display": "block"}),
                    html.Div(
                        [
                            dcc.Input(
                                id=f"{prefix}-atac-browser-locus",
                                value=locus,
                                debounce=True,
                                placeholder="chr1:100,000-200,000 or a gene name (e.g. IL7R)",
                                style={"width": "320px", "height": "38px"},
                            ),
                            dbc.Button(
                                "Update plot",
                                id=f"{prefix}-atac-browser-go",
                                color="primary",
                                size="sm",
                                style={"height": "38px", "whiteSpace": "nowrap"},
                            ),
                            dbc.Tooltip(
                                "You can type a location (e.g. chr1:100,000-200,000) or a gene name here.",
                                target=f"{prefix}-atac-browser-locus",
                                placement="bottom",
                            ),
                            html.Span(
                                id=f"{prefix}-atac-browser-message",
                                style={"fontSize": "13px", "color": "#555", "marginLeft": "6px"},
                            ),
                        ],
                        style={"display": "flex", "gap": "8px", "alignItems": "center", "flexWrap": "wrap"},
                    ),
                ],
                style={"marginBottom": "12px"},
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Expression metric", style={"fontWeight": "bold", "marginBottom": "5px", "display": "block"}),
                            dbc.RadioItems(
                                id=f"{prefix}-atac-browser-metric",
                                options=[
                                    {"label": "Mean expression", "value": "mean"},
                                    {"label": "% detected", "value": "detection"},
                                ],
                                value="mean",
                                inline=True,
                                style={"fontSize": "13px"},
                            ),
                        ],
                        style={"minWidth": "240px"},
                    ),
                    html.Div(
                        [
                            html.Label("Y-axis scale", style={"fontWeight": "bold", "marginBottom": "5px", "display": "block"}),
                            dbc.RadioItems(
                                id=f"{prefix}-atac-browser-yscale",
                                options=[
                                    {"label": "Auto", "value": "auto"},
                                    {"label": "Shared", "value": "shared"},
                                ],
                                value="shared",
                                inline=True,
                                style={"fontSize": "13px"},
                            ),
                        ],
                        style={"minWidth": "180px"},
                    ),
                ],
                style={"display": "flex", "gap": "32px", "flexWrap": "wrap", "marginBottom": "12px"},
            ),
            responsive_graph_grid(
                f"{prefix}-atac-browser-grid",
                f"{prefix}-atac-browser-grid-item",
                graph_flex_container(
                    f"{prefix}-atac-browser-graph",
                    config={
                        # Scroll/trackpad zoom disabled -- it jittered. Navigation
                        # is left/right drag (pan); use the modebar zoom buttons /
                        # box-zoom to zoom in.
                        "scrollZoom": False,
                        "displaylogo": False,
                        "modeBarButtonsToRemove": ["select2d", "lasso2d"],
                    },
                    container_id=f"{prefix}-atac-browser-grid-item",
                ),
            ),
        ],
        style={"padding": "20px"},
    )
