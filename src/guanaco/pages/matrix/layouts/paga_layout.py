import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.data.loader import get_discrete_labels
from guanaco.utils.ui_helpers import (
    component_flex_container,
    labeled_dropdown,
    labeled_radioitems,
)


def generate_paga_layout(adata, prefix):
    categorical_obs_columns = get_discrete_labels(adata)
    obs_options = [{"label": col, "value": col} for col in categorical_obs_columns]

    color_mode_control = labeled_radioitems(
        "Color By:",
        f"{prefix}-paga-color-mode",
        [
            {"label": "obs", "value": "obs"},
            {"label": "gene", "value": "gene"},
        ],
        value="obs",
        inline=True,
        wrapper_style={"marginBottom": "10px"},
    )
    obs_selection = html.Div(
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
    )
    gene_selection = html.Div(
        labeled_dropdown(
            "Gene:",
            f"{prefix}-paga-gene-dropdown",
            (
                [{"label": adata.var_names[0], "value": adata.var_names[0]}]
                if adata.n_vars > 0
                else []
            ),
            value=None,
            placeholder="Search gene",
            clearable=True,
            dropdown_style={"width": "240px"},
        ),
        id=f"{prefix}-paga-gene-wrapper",
        style={"display": "none", "marginBottom": "10px"},
    )
    threshold_control = html.Div(
        [
            html.Label(
                "Connectivity Threshold:",
                style={"fontWeight": "bold", "marginBottom": "5px"},
            ),
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
    )
    download_button = html.Button(
        "⬇ Download SVG",
        id=f"{prefix}-paga-download-svg",
        n_clicks=0,
        style={
            "border": "1px solid #ccc",
            "borderRadius": "5px",
            "padding": "5px 10px",
            "backgroundColor": "white",
            "cursor": "pointer",
        },
    )

    controls = html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(color_mode_control, xs=12, lg=4),
                    dbc.Col([obs_selection, gene_selection], xs=12, lg=4),
                    dbc.Col(threshold_control, xs=12, lg=4),
                ],
                style={"marginBottom": "5px"},
            ),
            dbc.Row(
                dbc.Col(download_button, xs=12),
                style={"marginBottom": "15px"},
            ),
        ],
        style={
            "marginBottom": "15px",
            "paddingBottom": "10px",
            "borderBottom": "1px solid #eee",
        },
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
