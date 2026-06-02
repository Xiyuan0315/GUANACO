import dash_bootstrap_components as dbc
import dash_draggable
from dash import dcc, html

from guanaco.pages.matrix.plots.volcano import (
    DEFAULT_PADJ_THRESHOLD,
    DEFAULT_TOP_N,
    DEFAULT_X_THRESHOLD,
    volcano_entry_options,
    x_axis_options,
)
from guanaco.plot_config import common_config
from guanaco.utils.ui_helpers import labeled_dropdown


def _number_input(input_id, label, value, *, min_value=None, max_value=None, step=None, min_width="160px"):
    return html.Div(
        [
            html.Label(label, style={"fontWeight": "bold", "marginBottom": "5px"}),
            dcc.Input(
                id=input_id,
                type="number",
                value=value,
                min=min_value,
                max=max_value,
                step=step,
                style={"width": "100%", "minWidth": min_width},
            ),
        ],
        style={"marginBottom": "10px"},
    )


def generate_volcano_layout(adata, prefix):
    entry_options, default_entry = volcano_entry_options(adata)

    controls = html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        labeled_dropdown(
                            "Group / contrast:",
                            f"{prefix}-volcano-entry-dropdown",
                            entry_options,
                            value=default_entry,
                            placeholder="Select a DE result",
                            clearable=False,
                            dropdown_style={"minWidth": "220px"},
                            wrapper_style={"marginBottom": "10px"},
                        ),
                        width=12,
                        lg=4,
                    ),
                    dbc.Col(
                        labeled_dropdown(
                            "X axis:",
                            f"{prefix}-volcano-x-axis-dropdown",
                            x_axis_options(),
                            value="logfoldchange",
                            clearable=False,
                            dropdown_style={"minWidth": "180px"},
                            wrapper_style={"marginBottom": "10px"},
                        ),
                        width=12,
                        lg=3,
                    ),
                    dbc.Col(
                        _number_input(
                            f"{prefix}-volcano-top-n",
                            "Top labels:",
                            DEFAULT_TOP_N,
                            min_value=0,
                            step=1,
                            min_width="160px",
                        ),
                        width=12,
                        lg=2,
                    ),
                ],
                style={"marginBottom": "8px"},
            ),
            dbc.Row(
                [
                    dbc.Col(
                        _number_input(
                            f"{prefix}-volcano-padj-threshold",
                            "P-value threshold:",
                            DEFAULT_PADJ_THRESHOLD,
                            min_value=0,
                            max_value=1,
                            step=0.01,
                        ),
                        width=12,
                        lg=3,
                    ),
                    dbc.Col(
                        _number_input(
                            f"{prefix}-volcano-x-threshold",
                            "X-axis threshold:",
                            DEFAULT_X_THRESHOLD,
                            min_value=0,
                            step=0.25,
                        ),
                        width=12,
                        lg=3,
                    ),
                ]
            ),
        ],
        style={"marginBottom": "15px", "paddingBottom": "10px", "borderBottom": "1px solid #eee"},
    )

    volcano_panel = html.Div(
        [
            html.Div(
                [
                    dcc.Graph(
                        id=f"{prefix}-volcano-plot",
                        config=common_config,
                        responsive=True,
                        style={"min-height": "0", "flex-grow": "1", "backgroundColor": "transparent"},
                    )
                ],
                style={
                    "flex": "3 1 620px",
                    "minWidth": "0",
                    "height": "100%",
                    "display": "flex",
                    "flexDirection": "column",
                },
            ),
            html.Div(
                [
                    html.Div(
                        id=f"{prefix}-volcano-deg-summary",
                        style={"marginBottom": "12px"},
                    ),
                    dbc.Button(
                        "Download DEGs CSV",
                        id=f"{prefix}-volcano-download-button",
                        n_clicks=0,
                        color="primary",
                        size="sm",
                        style={"width": "100%"},
                    ),
                    dcc.Download(id=f"{prefix}-volcano-degs-download"),
                ],
                style={
                    "flex": "1 1 260px",
                    "minWidth": "240px",
                    "padding": "14px",
                    "border": "1px solid rgba(31, 41, 51, 0.12)",
                    "borderRadius": "8px",
                    "backgroundColor": "#f8f9fa",
                    "alignSelf": "flex-start",
                },
            ),
        ],
        style={
            "height": "100%",
            "width": "100%",
            "display": "flex",
            "gap": "16px",
            "flexWrap": "wrap",
        },
    )

    draggable_container = dash_draggable.GridLayout(
        id=f"{prefix}-draggable-volcano",
        className="grid-layout-no-border",
        children=[
            html.Div(
                volcano_panel,
                style={
                    "height": "100%",
                    "width": "100%",
                    "backgroundColor": "transparent",
                    "border": "none",
                    "boxShadow": "none",
                },
            )
        ],
        isResizable=True,
        isDraggable=True,
        height=30,
        gridCols=12,
        style={"backgroundColor": "transparent", "padding": "0px", "border": "none", "boxShadow": "none"},
    )

    return html.Div([controls, draggable_container], style={"padding": "20px", "marginBottom": "15px"})
