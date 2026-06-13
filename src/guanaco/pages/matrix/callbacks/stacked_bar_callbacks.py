import plotly.graph_objects as go
from dash import Input, Output, State, no_update

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import sorted_categories
from guanaco.utils.render_guard import signature


def register_stacked_bar_callbacks(
    app,
    adata,
    prefix,
    *,
    filter_data,
    plot_stacked_bar,
    palette_json,
    color_config,
):
    # The x-axis is "Select Annotation"; the x-axis groups are its "Select Labels".
    # This grid lets the user reorder those groups by dragging the column headers.
    @app.callback(
        [Output(f"{prefix}-stacked-bar-x-order-grid", "columnDefs"), Output(f"{prefix}-stacked-bar-x-order-grid", "rowData")],
        [
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-tabs", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
        ],
    )
    def update_x_axis_order_grid(selected_labels, annotation, active_tab, selected_cells):
        if active_tab != "stacked-bar-tab" or not annotation:
            return [], []

        if selected_labels:
            x_values = [str(v) for v in selected_labels]
        else:
            src = adata[selected_cells] if selected_cells else adata
            x_values = [str(v) for v in sorted_categories(src, annotation)]

        column_defs = [
            {
                "field": val,
                "headerName": val,
                "width": 150,
                "minWidth": 120,
                "suppressMovable": False,
                "headerClass": "ag-header-cell-center",
                "resizable": True,
            }
            for val in x_values
        ]
        return column_defs, []

    @app.callback(
        Output(f"{prefix}-x-axis-column-order-store", "data"),
        Input(f"{prefix}-stacked-bar-x-order-grid", "columnState"),
        prevent_initial_call=True,
    )
    def update_column_order(column_state):
        if not column_state:
            return []
        return [col["colId"] for col in column_state if "colId" in col]

    @app.callback(
        [
            Output(f"{prefix}-stacked-bar-plot", "figure"),
            Output(f"{prefix}-stacked-bar-rendered-key", "data"),
        ],
        [
            Input(f"{prefix}-norm-box", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-stacked-bar-stack-by", "value"),
            Input(f"{prefix}-x-axis-column-order-store", "data"),
        ],
        [
            State(f"{prefix}-stacked-bar-plot", "figure"),
            State(f"{prefix}-stacked-bar-rendered-key", "data"),
            State(f"{prefix}-selected-cells-store", "data"),
        ],
    )
    def update_stacked_bar(norm, discrete_color_map, cells_hash, active_tab, annotation, selected_labels, stack_by, x_axis_order, current_figure, rendered_key, selected_cells):
        if active_tab != "stacked-bar-tab":
            return no_update, no_update

        # x-axis = "Select Annotation"; stacked color layers = "Stack bars by".
        if not annotation or not stack_by:
            fig = go.Figure()
            fig.add_annotation(
                text="Select an annotation (x-axis) in the left panel and a 'Stack bars by' variable",
                xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                font=dict(size=14), xanchor="center", yanchor="middle",
            )
            fig.update_layout(plot_bgcolor="white", paper_bgcolor="white", xaxis=dict(visible=False), yaxis=dict(visible=False))
            return fig, None

        cache_key = signature(
            "stacked-bar", norm, discrete_color_map, cells_hash,
            annotation, selected_labels, stack_by, x_axis_order,
        )
        if cache_key == rendered_key and current_figure:
            return no_update, no_update

        # Keep only the cells in the x-axis groups to display (Select Labels of the
        # x-axis annotation). The stack variable keeps all of its layers.
        filtered_adata = filter_data(adata, annotation, selected_labels, selected_cells)

        # X-axis group order: dragged order from the grid, else the selected labels,
        # else every category of the x-axis annotation.
        if x_axis_order:
            final_x_order = x_axis_order
        elif selected_labels:
            final_x_order = [str(v) for v in selected_labels]
        else:
            final_x_order = [str(v) for v in sorted_categories(adata, annotation)]

        # Color the stacked layers (the "Stack bars by" variable). Resolve the
        # palette like the scatter/embedding so the same category gets the same
        # color everywhere; str() keys match plot_stacked_bar's astype(str).
        stack_categories = sorted_categories(adata, stack_by)
        discrete_palette = resolve_discrete_palette(
            discrete_color_map, len(stack_categories), default=color_config
        )
        fixed_color_map = {
            str(cat): discrete_palette[i % len(discrete_palette)]
            for i, cat in enumerate(stack_categories)
        }

        fig = plot_stacked_bar(
            x_meta=annotation,
            y_meta=stack_by,
            norm=norm,
            adata=filtered_adata,
            color_map=fixed_color_map,
            y_order=None,
            x_order=final_x_order,
        )
        return fig, cache_key
