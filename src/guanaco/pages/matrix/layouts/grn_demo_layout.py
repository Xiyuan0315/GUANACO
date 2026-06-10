from dash import dcc, html

from guanaco.pages.matrix.plots.grn_demo import (
    ALL_CONTEXTS,
    default_grn_edge_threshold,
    grn_context_options,
    grn_dataframe,
    grn_weight_range,
    grn_weight_step,
)
from guanaco.utils.ui_helpers import component_flex_container, labeled_dropdown


def generate_grn_demo_layout(adata, prefix):
    context_options, context_column, has_weight = grn_context_options(adata)
    if has_weight:
        grn_df = grn_dataframe(adata)
        min_weight, max_weight = grn_weight_range(grn_df)
        threshold_value = default_grn_edge_threshold(grn_df)
        threshold_step = grn_weight_step(grn_df)
    else:
        min_weight, max_weight = None, None
        threshold_value = None
        threshold_step = 0.05
    layout_options = [
        {"label": "cose", "value": "cose"},
        {"label": "breadthfirst", "value": "breadthfirst"},
        {"label": "circle", "value": "circle"},
        {"label": "concentric", "value": "concentric"},
    ]
    context_label = f"{context_column.replace('_', ' ').title()}:" if context_column else "Context:"

    controls = html.Div(
        [
            labeled_dropdown(
                context_label,
                f"{prefix}-grn-demo-context",
                context_options,
                value=ALL_CONTEXTS,
                clearable=False,
                dropdown_style={"width": "220px"},
                wrapper_style={"marginBottom": "10px"},
            ),
            labeled_dropdown(
                "Layout:",
                f"{prefix}-grn-demo-layout",
                layout_options,
                value="cose",
                clearable=False,
                dropdown_style={"width": "220px"},
                wrapper_style={"marginBottom": "10px"},
            ),
            html.Div(
                [
                    html.Label("Edge Threshold:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                    dcc.Input(
                        id=f"{prefix}-grn-demo-threshold",
                        type="number",
                        min=min_weight,
                        max=max_weight,
                        step=threshold_step,
                        value=threshold_value,
                        style={"width": "140px"},
                    ),
                ],
                id=f"{prefix}-grn-demo-threshold-wrapper",
                style={"marginBottom": "10px", "display": "block" if has_weight else "none"},
            ),
            html.Button(
                "⬇ Download SVG",
                id=f"{prefix}-grn-demo-download-svg",
                n_clicks=0,
                style={"border": "1px solid #ccc", "borderRadius": "5px",
                       "padding": "5px 10px", "backgroundColor": "white", "cursor": "pointer"},
            ),
        ]
    )

    # The GRN layout is circular, so use a square viewport (instead of a wide
    # rectangle) -- a circle in a wide-but-short box gets its top/bottom clipped
    # on zoom. Responsive square, capped at the viewport height and centered.
    graph_container = html.Div(
        component_flex_container(f"{prefix}-grn-demo"),
        style={
            "width": "min(100%, 80vh)",
            "aspectRatio": "1 / 1",
            "minHeight": "480px",
            "margin": "0 auto",
        },
    )

    return html.Div(
        [dcc.Store(id=f"{prefix}-grn-demo-rendered-key"), controls, graph_container],
        style={"padding": "20px", "marginBottom": "15px"},
    )
