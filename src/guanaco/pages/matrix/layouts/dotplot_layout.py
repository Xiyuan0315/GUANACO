import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.utils.plot_config import dotplot_config
from guanaco.utils.render_guard import rendered_key_store
from guanaco.utils.ui_helpers import responsive_graph_grid


def generate_dotplot_layout(prefix):
    row1 = dbc.Row(
        [
            html.Div(
                [
                    html.Label("Standardization:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                    dbc.RadioItems(
                        id=f"{prefix}-dotplot-standardization",
                        options=[
                            {"label": "None", "value": "None"},
                            {"label": "Min–max", "value": "minmax"},
                            {"label": "Z-score", "value": "zscore"},
                        ],
                        value="minmax",
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
                ],
                width=6,
                style={"paddingLeft": "15px"},
            ),
        ],
        style={"marginBottom": "10px", "padding": "10px"},
    )

    dotplot_advanced_toggle = dbc.Button(
        "▸ More options",
        id=f"{prefix}-dotplot-options-toggle",
        color="link",
        size="sm",
        style={"padding": "2px 0", "textDecoration": "none", "fontWeight": "bold", "marginBottom": "10px"},
    )

    dotplot_advanced_panel = dbc.Collapse(
        row2,
        id=f"{prefix}-dotplot-options-collapse",
        is_open=False,
    )

    # Width-responsive grid (scales with the screen instead of capping at ~1200px).
    # Default opens at 3/4 width x ~510px; drag-resize is bounded by the helper's
    # min/max. The grid item carries the id so the per-item min/max reach
    # react-grid-layout. Fresh grid id (was "-draggable-dotplot") so a stale saved
    # layout without these constraints can't override them.
    grid_item_id = f"{prefix}-dotplot-grid-item"
    dotplot_graph_item = html.Div(
        [
            dcc.Graph(
                id=f"{prefix}-dotplot",
                config=dotplot_config,
                responsive=True,
                style={"min-height": "0", "flex-grow": "1", "backgroundColor": "transparent"},
            )
        ],
        id=grid_item_id,
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
    draggable_container = responsive_graph_grid(
        f"{prefix}-dotplot-grid", grid_item_id, dotplot_graph_item
    )

    return html.Div(
        [rendered_key_store(prefix, "dotplot"), row1, dotplot_advanced_toggle, dotplot_advanced_panel, draggable_container],
        style={"padding": "20px", "marginBottom": "15px"},
    )
