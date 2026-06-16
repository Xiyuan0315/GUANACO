import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.utils.ui_helpers import (
    LOADING_OVERLAY_STYLE,
    graph_flex_container,
    labeled_dropdown,
    responsive_graph_grid,
)


def generate_pseudotime_layout(prefix):
    return html.Div(
        [
            # Tracks the cache key of the figure currently shown, so the callback
            # can skip recomputing/redrawing on tab switches when nothing changed.
            dcc.Store(id=f"{prefix}-pseudotime-rendered-key"),
            dbc.Row(
                [
                    dbc.Col(
                        labeled_dropdown(
                            "Continuous Variable:",
                            f"{prefix}-pseudotime-key-dropdown",
                            [],
                            value=None,
                            placeholder="Select continuous obs variable",
                            clearable=False,
                        ),
                        width=6,
                    ),
                    dbc.Col(
                        html.Div(
                            [
                                html.Label("Minimum Expression Threshold:", style={"fontWeight": "bold"}),
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
                        ),
                        width=6,
                    ),
                ],
                style={"marginBottom": "20px", "padding": "15px", "backgroundColor": "white"},
            ),
            dcc.Loading(
                id=f"{prefix}-pseudotime-loading",
                type="circle",
                overlay_style=LOADING_OVERLAY_STYLE,
                children=[
                    responsive_graph_grid(
                        f"{prefix}-pseudotime-grid",
                        f"{prefix}-pseudotime-grid-item",
                        graph_flex_container(
                            f"{prefix}-pseudotime-plot",
                            container_id=f"{prefix}-pseudotime-grid-item",
                        ),
                    )
                ],
            ),
        ]
    )
