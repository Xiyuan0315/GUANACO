from dash import dcc, html
import pandas as pd

from guanaco.utils.ui_helpers import component_flex_container, labeled_dropdown, labeled_radioitems


def generate_paga_layout(adata, prefix):
    categorical_obs_columns = [
        col for col in adata.obs.columns if not pd.api.types.is_numeric_dtype(adata.obs[col])
    ]
    obs_options = [{"label": col, "value": col} for col in categorical_obs_columns]

    controls = html.Div(
        [
            labeled_radioitems(
                "Color By:",
                f"{prefix}-paga-color-mode",
                [
                    {"label": "obs", "value": "obs"},
                    {"label": "gene", "value": "gene"},
                ],
                value="obs",
                inline=True,
                wrapper_style={"marginBottom": "10px"},
            ),
            html.Div(
                labeled_dropdown(
                    "obs Column:",
                    f"{prefix}-paga-obs-dropdown",
                    obs_options,
                    value=categorical_obs_columns[0] if categorical_obs_columns else None,
                    clearable=False,
                    dropdown_style={"width": "240px"},
                ),
                id=f"{prefix}-paga-obs-wrapper",
                style={"marginBottom": "10px"},
            ),
            html.Div(
                labeled_dropdown(
                    "Gene:",
                    f"{prefix}-paga-gene-dropdown",
                    [{"label": adata.var_names[0], "value": adata.var_names[0]}] if adata.n_vars > 0 else [],
                    value=None,
                    placeholder="Search gene",
                    clearable=True,
                    dropdown_style={"width": "240px"},
                ),
                id=f"{prefix}-paga-gene-wrapper",
                style={"display": "none", "marginBottom": "10px"},
            ),
            html.Div(
                [
                    html.Label("Connectivity Threshold:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                    dcc.Input(
                        id=f"{prefix}-paga-threshold",
                        type="number",
                        min=0,
                        max=1,
                        step=0.01,
                        value=0.03,
                        style={"width": "140px"},
                    ),
                ],
                style={"marginBottom": "10px"},
            ),
            html.Button(
                "⬇ Download SVG",
                id=f"{prefix}-paga-download-svg",
                n_clicks=0,
                style={"border": "1px solid #ccc", "borderRadius": "5px",
                       "padding": "5px 10px", "backgroundColor": "white", "cursor": "pointer"},
            ),
        ]
    )

    # Cytoscape handles its own zoom/pan, so the graph just fills a
    # screen-adaptive box instead of a drag-to-resize container.
    graph_container = html.Div(
        component_flex_container(f"{prefix}-paga"),
        style={"width": "100%", "height": "75vh", "minHeight": "480px"},
    )

    return html.Div(
        [
            dcc.Store(id=f"{prefix}-paga-rendered-key"),
            controls,
            graph_container,
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )
