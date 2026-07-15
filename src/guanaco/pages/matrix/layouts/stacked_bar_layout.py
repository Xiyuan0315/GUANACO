import dash_ag_grid as dag
import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.utils.obs_utils import sorted_categories
from guanaco.utils.plot_config import common_config
from guanaco.utils.ui_helpers import (
    labeled_dropdown,
    labeled_radioitems,
    responsive_graph_grid,
)


def _reference_categories(adata, annotation):
    if not annotation:
        return [], None
    categories = [str(value) for value in sorted_categories(adata, annotation)]
    return [{"label": value, "value": value} for value in categories], None


def generate_stacked_bar_layout(adata, discrete_label_list, prefix):
    default_x = discrete_label_list[0] if discrete_label_list else None
    default_stack_by = (
        discrete_label_list[1]
        if len(discrete_label_list) > 1
        else default_x
    )
    sample_unit_list = [str(column) for column in adata.obs.columns]
    reference_options, _ = _reference_categories(adata, default_stack_by)
    view_mode = labeled_radioitems(
        "View:",
        f"{prefix}-composition-view",
        [
            {"label": "Stacked bars", "value": "bars"},
            {"label": "Hierarchy", "value": "hierarchy"},
        ],
        value="bars",
        inline=True,
        wrapper_style={"marginBottom": "15px"},
    )
    x_group_dropdown = labeled_dropdown(
        "Group bars by:",
        f"{prefix}-stacked-bar-x-group",
        [{"label": meta, "value": meta} for meta in discrete_label_list],
        value=default_x,
        clearable=False,
        placeholder="Select the x-axis grouping",
        label_id=f"{prefix}-composition-primary-label",
        wrapper_style={"marginBottom": "15px"},
    )
    stack_by_label_id = f"{prefix}-stack-by-label"
    stack_by_dropdown = labeled_dropdown(
        "Stack bars by:",
        f"{prefix}-stacked-bar-stack-by",
        [{"label": meta, "value": meta} for meta in discrete_label_list],
        value=default_stack_by,
        clearable=False,
        placeholder="Select metadata for the stacked color layers",
        label_id=stack_by_label_id,
        dropdown_style={"marginBottom": "15px"},
        wrapper_style={"marginBottom": "15px"},
    )
    stack_by_tooltip = dbc.Tooltip(
        "X-axis groups come from 'Group bars by'. Colors show the 'Stack bars by' variable.",
        target=stack_by_label_id,
        id=f"{prefix}-composition-selector-tooltip",
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

    sample_unit_dropdown = labeled_dropdown(
        "Sample unit:",
        f"{prefix}-composition-sample-unit",
        [{"label": meta, "value": meta} for meta in sample_unit_list],
        value=None,
        clearable=True,
        placeholder="Select donor, patient, or replicate",
        wrapper_style={"marginBottom": "10px"},
    )
    reference_dropdown = labeled_dropdown(
        "ALR reference population:",
        f"{prefix}-composition-alr-reference",
        reference_options,
        value=None,
        clearable=False,
        placeholder="Select a reference population",
        wrapper_style={"marginBottom": "10px"},
    )

    draggable_bar = responsive_graph_grid(
        f"{prefix}-stacked-bar-grid",
        f"{prefix}-stacked-bar-grid-item",
        html.Div(
            id=f"{prefix}-stacked-bar-grid-item",
            children=[
                dcc.Graph(
                    id=f"{prefix}-stacked-bar-plot",
                    config=common_config,
                    responsive=True,
                    style={"flex-grow": "1"},
                )
            ],
            style={
                "height": "100%",
                "width": "100%",
                "display": "flex",
                "flex-direction": "column",
                "flex-grow": "0",
            },
        ),
    )

    x_axis_order_title_id = f"{prefix}-x-axis-order-title"
    x_axis_order_component = html.Div(
        [
            html.Label(
                "X-axis group order:",
                id=x_axis_order_title_id,
                style={"fontWeight": "bold", "marginBottom": "10px"},
            ),
            dbc.Tooltip(
                "Drag column headers to reorder x-axis groups.",
                target=x_axis_order_title_id,
            ),
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
        ],
        id=f"{prefix}-stacked-bar-x-order-container",
    )

    column_order_store = dcc.Store(id=f"{prefix}-x-axis-column-order-store", data=[])

    controls_row = dbc.Row(
        [
            dbc.Col(view_mode, xs=12, lg=True),
            dbc.Col(x_group_dropdown, xs=12, lg=True),
            dbc.Col(stack_by_dropdown, xs=12, lg=True),
            dbc.Col(
                norm_box,
                id=f"{prefix}-composition-value-control",
                xs=12,
                lg=True,
            ),
        ],
        style={"marginBottom": "15px"},
    )
    differential_abundance_controls = dbc.Row(
        [
            dbc.Col(sample_unit_dropdown, xs=12, lg=4),
            dbc.Col(reference_dropdown, xs=12, lg=4),
            dbc.Col(
                html.Div(
                    "Tests use sample-level cell counts. Every pair of x-axis "
                    "groups is compared and FDR-corrected together.",
                    style={"color": "#6c757d", "fontSize": "13px"},
                ),
                xs=12,
                lg=4,
                style={"display": "flex", "alignItems": "center"},
            ),
        ],
        id=f"{prefix}-composition-da-controls",
        style={"marginBottom": "10px"},
    )
    differential_abundance_panel = html.Div(
        [
            dbc.Button(
                "▸ Differential abundance",
                id=f"{prefix}-composition-da-toggle",
                n_clicks=0,
                color="link",
                style={
                    "fontSize": "18px",
                    "fontWeight": "600",
                    "padding": "8px 0",
                    "textDecoration": "none",
                    "color": "#5f6368",
                },
            ),
            html.Div(
                html.Div(
                    [
                        differential_abundance_controls,
                        dcc.Store(id=f"{prefix}-composition-da-rendered-key"),
                        html.Div(
                            dcc.Graph(
                                id=f"{prefix}-composition-da-plot",
                                config=common_config,
                                responsive=True,
                                style={"height": "100%", "width": "100%"},
                            ),
                            id=f"{prefix}-composition-da-results",
                            className="composition-da-results",
                            style={"display": "none"},
                        ),
                    ],
                    style={"width": "100%"},
                ),
                id=f"{prefix}-composition-da-collapse",
                style={"display": "none", "width": "100%"},
            ),
        ],
        id=f"{prefix}-composition-da-panel",
        style={"marginTop": "24px", "width": "100%"},
    )

    return html.Div(
        [
            column_order_store,
            dcc.Store(id=f"{prefix}-stacked-bar-rendered-key"),
            controls_row,
            stack_by_tooltip,
            draggable_bar,
            x_axis_order_component,
            differential_abundance_panel,
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )
