import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.utils.ui_helpers import labeled_dropdown, labeled_radioitems


def generate_pseudotime_layout(prefix):
    return html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.Div(
                                [
                                    html.Label("Minimum Expression Threshold:", className="control-label"),
                                    dcc.Slider(
                                        id=f"{prefix}-pseudotime-min-expr-slider",
                                        min=0,
                                        max=5,
                                        step=0.1,
                                        value=0.5,
                                        marks=None,
                                        tooltip={"placement": "bottom", "always_visible": True},
                                        className="dbc-slider",
                                    ),
                                ],
                                style={"marginBottom": "20px"},
                            ),
                            labeled_radioitems(
                                "Transformation:",
                                f"{prefix}-pseudotime-transformation",
                                [
                                    {"label": "None", "value": "none"},
                                    {"label": "Log", "value": "log"},
                                    {"label": "Z-score", "value": "z_score"},
                                ],
                                value="none",
                                inline=True,
                                radio_style={"fontSize": "14px"},
                                wrapper_style={"marginBottom": "20px"},
                            ),
                            labeled_dropdown(
                                "Pseudotime Column:",
                                f"{prefix}-pseudotime-key-dropdown",
                                [],
                                value=None,
                                placeholder="Select pseudotime column",
                                clearable=False,
                                dropdown_style={"marginBottom": "15px"},
                                wrapper_style={"marginBottom": "20px"},
                            ),
                        ],
                        width=12,
                    )
                ],
                style={"marginBottom": "20px", "padding": "15px", "backgroundColor": "#f8f9fa", "borderRadius": "5px"},
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dcc.Loading(
                                id=f"{prefix}-pseudotime-loading",
                                type="default",
                                children=[
                                    dcc.Graph(
                                        id=f"{prefix}-pseudotime-plot",
                                        style={"height": "auto", "minHeight": "600px"},
                                        config={
                                            "displayModeBar": True,
                                            "displaylogo": False,
                                            "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                                            "toImageButtonOptions": {
                                                "format": "png",
                                                "filename": "pseudotime_plot",
                                                "height": 800,
                                                "width": 1200,
                                                "scale": 2,
                                            },
                                        },
                                    )
                                ],
                            )
                        ],
                        width=12,
                    )
                ]
            ),
        ]
    )
