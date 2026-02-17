import dash_bootstrap_components as dbc
import dash_draggable
import plotly.express as px
from dash import dcc, html

from guanaco.config import common_config


def _continuous_colormap_options():
    options = [{"label": scale, "value": scale} for scale in px.colors.named_colorscales()]
    try:
        import colorcet as cc

        cc_names = sorted(
            name for name in cc.palette.keys()
            if "glasbey" not in name.lower()
        )
        options.extend(
            {"label": f"colorcet/{name}", "value": f"cc:{name}"}
            for name in cc_names
        )
    except Exception:
        # colorcet is optional; keep Plotly maps if unavailable.
        pass
    return options


def generate_dotplot_layout(prefix):
    colorscales = _continuous_colormap_options()

    row1 = dbc.Row(
        [
            html.Div(
                [
                    html.Label("Transformation: ", style={"fontWeight": "bold", "marginBottom": "5px"}),
                    dbc.RadioItems(
                        id=f"{prefix}-dotplot-log-or-zscore",
                        options=[{"label": "None", "value": "None"}, {"label": "Log", "value": "log"}],
                        value="None",
                        inline=True,
                        style={"marginRight": "15px"},
                    ),
                ]
            ),
            html.Div(
                [
                    html.Label("Standardization:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                    dbc.RadioItems(
                        id=f"{prefix}-dotplot-standardization",
                        options=[
                            {"label": "None", "value": "None"},
                            {"label": "By variable", "value": "var"},
                            {"label": "By label", "value": "group"},
                        ],
                        value="None",
                        inline=True,
                    ),
                ]
            ),
        ],
        style={"marginBottom": "15px", "paddingBottom": "10px", "borderBottom": "1px solid #eee"},
    )

    row2 = dbc.Row(
        [
            dbc.Col(
                [
                    html.Label("Clustering:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                    dbc.RadioItems(
                        id=f"{prefix}-dotplot-cluster-mode",
                        options=[
                            {"label": "None", "value": "none"},
                            {"label": "Rows", "value": "row"},
                            {"label": "Cols", "value": "col"},
                            {"label": "Both", "value": "both"},
                        ],
                        value="none",
                        inline=True,
                        style={"marginBottom": "5px"},
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label("Linkage Method", style={"fontSize": "0.9em", "marginRight": "5px"}),
                                    dcc.Dropdown(
                                        id=f"{prefix}-dotplot-cluster-method",
                                        options=[{"label": m.capitalize(), "value": m} for m in ["average", "ward", "complete", "single"]],
                                        value="average",
                                        clearable=False,
                                        style={"width": "100px"},
                                    ),
                                ],
                                style={"display": "inline-block", "marginRight": "15px"},
                            ),
                            html.Div(
                                [
                                    html.Label("Distance Metric", style={"fontSize": "0.9em", "marginRight": "5px"}),
                                    dcc.Dropdown(
                                        id=f"{prefix}-dotplot-cluster-metric",
                                        options=[
                                            {"label": "Correlation", "value": "correlation"},
                                            {"label": "Euclidean", "value": "euclidean"},
                                            {"label": "Cosine", "value": "cosine"},
                                        ],
                                        value="correlation",
                                        clearable=False,
                                        style={"width": "120px"},
                                    ),
                                ],
                                style={"display": "inline-block"},
                            ),
                        ]
                    ),
                ],
                width=6,
                style={"borderRight": "1px solid #eee", "paddingRight": "15px"},
            ),
            dbc.Col(
                [
                    html.Label("Appearance:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                    html.Div(
                        [
                            dbc.Checklist(
                                options=[{"label": "Matrix Plot", "value": "matrixplot"}],
                                value=[],
                                id=f"{prefix}-plot-type-switch",
                                inline=True,
                                switch=True,
                                style={"display": "inline-block", "marginRight": "15px"},
                            ),
                            dbc.Checklist(
                                options=[{"label": "Swap Axes", "value": "swap"}],
                                value=[],
                                id=f"{prefix}-dotplot-transpose",
                                inline=True,
                                switch=True,
                                style={"display": "inline-block", "marginRight": "15px"},
                            ),
                        ],
                        style={"marginBottom": "5px"},
                    ),
                    html.Div(
                        [
                            html.Label("Continuous ColorMap:", style={"fontWeight": "bold", "marginBottom": "5px", "fontSize": "0.9em"}),
                            dcc.Dropdown(
                                id=f"{prefix}-dotmatrix-color-map-dropdown",
                                options=colorscales,
                                value="viridis",
                                clearable=False,
                                style={"width": "200px", "marginBottom": "10px"},
                            ),
                        ]
                    ),
                ],
                width=6,
                style={"paddingLeft": "15px"},
            ),
        ],
        style={"marginBottom": "10px", "backgroundColor": "#f8f9fa", "padding": "10px", "borderRadius": "5px"},
    )

    draggable_container = dash_draggable.GridLayout(
        id=f"{prefix}-draggable-dotplot",
        className="grid-layout-no-border",
        children=[
            html.Div(
                [
                    dcc.Graph(
                        id=f"{prefix}-dotplot",
                        config=common_config,
                        responsive=True,
                        style={"min-height": "0", "flex-grow": "1", "backgroundColor": "transparent"},
                    )
                ],
                style={
                    "backgroundColor": "transparent",
                    "border": "none",
                    "boxShadow": "none",
                    "height": "100%",
                    "width": "100%",
                    "display": "flex",
                    "flex-direction": "column",
                    "flex-grow": "0",
                },
            )
        ],
        isResizable=True,
        isDraggable=True,
        height=30,
        gridCols=12,
        style={"backgroundColor": "transparent", "padding": "0px", "border": "none", "boxShadow": "none"},
    )

    return html.Div([row1, row2, draggable_container], style={"padding": "20px", "marginBottom": "15px"})
