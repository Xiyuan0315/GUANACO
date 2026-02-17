import dash_ag_grid as dag
import dash_draggable
from dash import dcc, html

from guanaco.config import common_config
from guanaco.utils.ui_helpers import labeled_dropdown, labeled_radioitems


def generate_stacked_bar_layout(discrete_label_list, prefix):
    x_axis_dropdown = labeled_dropdown(
        "Cell info (x-axis):",
        f"{prefix}-stacked-bar-x-axis",
        [{"label": meta, "value": meta} for meta in discrete_label_list],
        value=discrete_label_list[0] if discrete_label_list else None,
        clearable=False,
        placeholder="Select metadata for x-axis",
        dropdown_style={"marginBottom": "15px"},
        wrapper_style={"marginBottom": "15px"},
    )

    norm_box = labeled_radioitems(
        "Plot value:",
        f"{prefix}-norm-box",
        [
            {"label": "Proportion", "value": "prop"},
            {"label": "Count", "value": "count"},
        ],
        value="prop",
        inline=True,
        wrapper_style={"marginBottom": "15px"},
    )

    draggable_bar = dash_draggable.GridLayout(
        id=f"{prefix}-draggable",
        className="grid-layout-no-border",
        children=[
            html.Div(
                children=[
                    dcc.Graph(
                        id=f"{prefix}-stacked-bar-plot",
                        config=common_config,
                        responsive=True,
                        style={"flex-grow": "1"},
                    )
                ],
                style={"height": "100%", "width": "100%", "display": "flex", "flex-direction": "column", "flex-grow": "0"},
            ),
        ],
        isResizable=True,
        isDraggable=True,
        height=30,
        gridCols=12,
    )

    info_text = html.Div(
        [
            html.P(
                [
                    html.I(className="fas fa-info-circle", style={"marginRight": "5px"}),
                    "The stacked layers come from 'Select Annotation' and 'Select Labels' in the left control panel",
                ],
                style={"color": "#6c757d", "fontSize": "14px", "marginBottom": "15px"},
            )
        ]
    )

    x_axis_order_component = html.Div(
        [
            html.Label("X-axis group order:", style={"fontWeight": "bold", "marginBottom": "10px"}),
            dag.AgGrid(
                id=f"{prefix}-stacked-bar-x-order-grid",
                rowData=[],
                columnDefs=[],
                defaultColDef={
                    "sortable": False,
                    "filter": False,
                    "resizable": True,
                    "suppressMenu": True,
                    "headerHeight": 40,
                    "minWidth": 120,
                    "width": 150,
                    "headerClass": "ag-header-cell-center",
                },
                dashGridOptions={
                    "headerHeight": 40,
                    "rowHeight": 0,
                    "suppressRowClickSelection": True,
                    "suppressCellSelection": True,
                    "suppressMovableColumns": False,
                    "animateRows": False,
                    "suppressHorizontalScroll": False,
                    "onColumnMoved": True,
                    "suppressLoadingOverlay": True,
                    "suppressNoRowsOverlay": True,
                    "suppressDisplayTotal": True,
                },
                style={"height": "40px", "marginBottom": "10px", "overflow": "hidden"},
                className="ag-theme-alpine",
            ),
            html.P(
                "Drag column headers to reorder x-axis groups.",
                style={"fontSize": "12px", "color": "#6c757d", "marginTop": "5px", "marginBottom": "15px"},
            ),
        ]
    )

    column_order_store = dcc.Store(id=f"{prefix}-x-axis-column-order-store", data=[])

    return html.Div(
        [
            column_order_store,
            info_text,
            x_axis_dropdown,
            norm_box,
            draggable_bar,
            x_axis_order_component,
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )
