import dash_bootstrap_components as dbc
import dash_draggable
from dash import dcc, html

from guanaco.utils.ui_helpers import graph_flex_container, labeled_dropdown, labeled_radioitems


def generate_pseudotime_layout(prefix):
    return html.Div(
        [
            # Tracks the cache key of the figure currently shown, so the callback
            # can skip recomputing/redrawing on tab switches when nothing changed.
            dcc.Store(id=f"{prefix}-pseudotime-rendered-key"),
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
                                "Continuous Variable:",
                                f"{prefix}-pseudotime-key-dropdown",
                                [],
                                value=None,
                                placeholder="Select continuous obs variable",
                                clearable=False,
                                dropdown_style={"marginBottom": "15px"},
                                wrapper_style={"marginBottom": "20px"},
                            ),
                        ],
                        width=12,
                    )
                ],
                style={"marginBottom": "20px", "padding": "15px", "backgroundColor": "white"},
            ),
            dcc.Loading(
                id=f"{prefix}-pseudotime-loading",
                type="circle",
                children=[
                    dash_draggable.GridLayout(
                        id=f"{prefix}-draggable-pseudotime",
                        className="grid-layout-no-border",
                        children=[graph_flex_container(f"{prefix}-pseudotime-plot")],
                        isResizable=True,
                        isDraggable=True,
                        height=30,
                        gridCols=12,
                        style={"backgroundColor": "transparent", "padding": "0px", "border": "none", "boxShadow": "none"},
                    )
                ],
            ),
        ]
    )
