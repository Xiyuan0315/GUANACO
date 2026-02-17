import plotly.graph_objects as go
from dash import Input, Output, State, ALL, html
import dash_bootstrap_components as dbc


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
    @app.callback(
        [Output(f"{prefix}-x-axis-draggable-grid", "children"), Output(f"{prefix}-x-axis-groups-state", "data")],
        [Input(f"{prefix}-x-meta-dropdown", "value"), Input(f"{prefix}-single-cell-label-selection", "value"), Input(f"{prefix}-single-cell-tabs", "value")],
        [State(f"{prefix}-single-cell-annotation-dropdown", "value"), State(f"{prefix}-x-axis-groups-state", "data")],
    )
    def update_x_axis_groups_grid(x_meta, selected_labels, active_tab, y_meta, current_state):
        if active_tab != "stacked-bar-tab" or not x_meta:
            return [], {}

        x_values = sorted(adata.obs[x_meta].unique())
        if selected_labels and len(selected_labels) < len(x_values):
            if set(selected_labels).issubset(set(x_values)):
                x_values = selected_labels

        x_values = [str(val) for val in x_values]
        if not current_state:
            current_state = {}

        new_state = {val: current_state.get(val, True) for val in x_values}

        items = []
        for value in x_values:
            is_enabled = new_state.get(value, True)
            items.append(
                html.Div(
                    [
                        html.Div(
                            value,
                            style={"fontWeight": "bold", "marginBottom": "5px", "color": "#000" if is_enabled else "#999"},
                        ),
                        dbc.Switch(id={"type": f"{prefix}-x-group-switch", "index": value}, value=is_enabled, style={"marginTop": "5px"}),
                    ],
                    key=f"x-group-{value}",
                    style={
                        "padding": "10px",
                        "backgroundColor": "#fff" if is_enabled else "#f5f5f5",
                        "border": "2px solid #007bff" if is_enabled else "1px solid #ddd",
                        "borderRadius": "5px",
                        "textAlign": "center",
                        "cursor": "move",
                        "height": "100%",
                    },
                )
            )

        if not items:
            items = [html.Div("No groups available", key="empty-msg", style={"color": "#6c757d", "fontStyle": "italic", "padding": "20px"})]

        return items, new_state

    @app.callback(
        Output(f"{prefix}-x-axis-groups-state", "data", allow_duplicate=True),
        [Input({"type": f"{prefix}-x-group-switch", "index": ALL}, "value")],
        [State({"type": f"{prefix}-x-group-switch", "index": ALL}, "id"), State(f"{prefix}-x-axis-groups-state", "data")],
        prevent_initial_call=True,
    )
    def update_group_state(switch_values, switch_ids, current_state):
        if not switch_ids:
            return current_state
        new_state = current_state.copy() if current_state else {}
        for i, switch_id in enumerate(switch_ids):
            if "index" in switch_id:
                group_name = switch_id["index"]
                new_state[group_name] = switch_values[i]
        return new_state

    @app.callback(
        [Output(f"{prefix}-stacked-bar-x-order-grid", "columnDefs"), Output(f"{prefix}-stacked-bar-x-order-grid", "rowData")],
        [Input(f"{prefix}-stacked-bar-x-axis", "value"), Input(f"{prefix}-single-cell-tabs", "value"), Input(f"{prefix}-selected-cells-store", "data")],
    )
    def update_x_axis_order_grid(x_axis_meta, active_tab, selected_cells):
        if active_tab != "stacked-bar-tab" or not x_axis_meta:
            return [], []

        if selected_cells:
            filtered_adata = adata[selected_cells]
            x_values = sorted(filtered_adata.obs[x_axis_meta].unique())
        else:
            x_values = sorted(adata.obs[x_axis_meta].unique())
        x_values_str = [str(val) for val in x_values]

        column_defs = []
        for val in x_values_str:
            column_defs.append(
                {
                    "field": val,
                    "headerName": val,
                    "width": 150,
                    "minWidth": 120,
                    "suppressMovable": False,
                    "headerClass": "ag-header-cell-center",
                    "resizable": True,
                }
            )
        return column_defs, []

    @app.callback(
        Output(f"{prefix}-x-axis-column-order-store", "data"),
        Input(f"{prefix}-stacked-bar-x-order-grid", "columnState"),
        prevent_initial_call=True,
    )
    def update_column_order(column_state):
        if not column_state:
            return []
        column_order = []
        for col in column_state:
            if "colId" in col:
                column_order.append(col["colId"])
        return column_order

    @app.callback(
        Output(f"{prefix}-stacked-bar-plot", "figure"),
        [
            Input(f"{prefix}-norm-box", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-stacked-bar-x-axis", "value"),
            Input(f"{prefix}-x-axis-column-order-store", "data"),
        ],
        [State(f"{prefix}-stacked-bar-plot", "figure")],
    )
    def update_stacked_bar(norm, discrete_color_map, selected_cells, active_tab, annotation, selected_labels, x_axis_meta, x_axis_order, current_figure):
        if active_tab != "stacked-bar-tab":
            return current_figure if current_figure else go.Figure()

        if not annotation or not x_axis_meta:
            fig = go.Figure()
            fig.add_annotation(
                text="Please select both annotation (for stacking) and x-axis metadata",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
                font=dict(size=14),
                xanchor="center",
                yanchor="middle",
            )
            fig.update_layout(plot_bgcolor="white", paper_bgcolor="white", xaxis=dict(visible=False), yaxis=dict(visible=False))
            return fig

        filtered_adata = filter_data(adata, annotation, selected_labels, selected_cells)

        if x_axis_order and len(x_axis_order) > 0:
            final_x_order = x_axis_order
        elif x_axis_meta:
            x_values = sorted(adata.obs[x_axis_meta].unique())
            final_x_order = [str(val) for val in x_values]
        else:
            final_x_order = None

        all_categories = sorted(adata.obs[annotation].unique())
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            fixed_color_map = {cat: discrete_palette[i % len(discrete_palette)] for i, cat in enumerate(all_categories)}
        else:
            fixed_color_map = {cat: color_config[i % len(color_config)] for i, cat in enumerate(all_categories)}

        return plot_stacked_bar(
            x_meta=x_axis_meta,
            y_meta=annotation,
            norm=norm,
            adata=filtered_adata,
            color_map=fixed_color_map,
            y_order=selected_labels,
            x_order=final_x_order,
        )
