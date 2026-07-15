"""Layout for the IGV exploratory plot."""

import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.utils.ui_helpers import LOADING_OVERLAY_STYLE


def build_igv_layout(prefix, genome_tracks):
    """Build the modality-scoped IGV session and motif controls."""
    return html.Div(
        style={"display": "flex", "flexDirection": "column"},
        children=[
            dbc.Row(
                className="g-0",
                children=[
                    dbc.Col(
                        xs=12,
                        lg=8,
                        style={"padding": "18px"},
                        children=[
                            html.P("Select the IGV session to display below:"),
                            dcc.Dropdown(
                                id=f"{prefix}-igv-genome-select",
                                options=[
                                    {"label": session, "value": session}
                                    for session in genome_tracks
                                ],
                                value=None,
                                placeholder="Select an IGV session...",
                            ),
                            html.Hr(),
                            dcc.Loading(
                                id=f"{prefix}-igv-container",
                                overlay_style=LOADING_OVERLAY_STYLE,
                            ),
                        ],
                    ),
                    dbc.Col(
                        xs=12,
                        lg=4,
                        style={
                            "padding": "18px",
                            "borderLeft": "1px solid #ccc",
                        },
                        children=[
                            html.H4("Motif Search Box"),
                            html.P("Search with motif id, from JASPAR database"),
                            html.Div(
                                style={
                                    "display": "flex",
                                    "alignItems": "center",
                                    "gap": "10px",
                                },
                                children=[
                                    dcc.Input(
                                        id=f"{prefix}-search-input",
                                        type="text",
                                        placeholder="Enter a motif id (e.g., MA1972.1)",
                                        style={
                                            "flex": "1",
                                            "marginBottom": "10px",
                                        },
                                    ),
                                    html.Button(
                                        "Search",
                                        id=f"{prefix}-search-button",
                                        n_clicks=0,
                                        style={
                                            "padding": "10px 20px",
                                            "whiteSpace": "nowrap",
                                        },
                                    ),
                                ],
                            ),
                            html.Div(
                                id=f"{prefix}-search-results",
                                style={"marginTop": "20px"},
                                children=[],
                            ),
                        ],
                    ),
                ],
            )
        ],
    )

