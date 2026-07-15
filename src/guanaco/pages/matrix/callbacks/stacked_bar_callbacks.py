import plotly.graph_objects as go
from dash import Input, Output, State, no_update

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import sorted_categories
from guanaco.utils.render_guard import signature


_STACKED_BAR_TAB = "stacked-bar-tab"
_EXPLORATORY_WORKSPACE = "dataset-exploration"
_MISSING_STACKED_BAR_SELECTION = (
    "Select a 'Group bars by' annotation and a 'Stack bars by' annotation"
)
_BAR_SELECTOR_TOOLTIP = (
    "X-axis groups come from 'Group bars by'. Colors show the 'Stack bars by' "
    "variable."
)
_HIERARCHY_SELECTOR_TOOLTIP = (
    "Parent level defines the first hierarchy layer. Child level defines the "
    "categories nested beneath each parent."
)


def _composition_control_state(view_mode):
    if view_mode == "hierarchy":
        return (
            "Parent level:",
            "Child level:",
            {"display": "none"},
            _HIERARCHY_SELECTOR_TOOLTIP,
            {"display": "none"},
            {"display": "none"},
        )
    return (
        "Group bars by:",
        "Stack bars by:",
        {},
        _BAR_SELECTOR_TOOLTIP,
        {"marginBottom": "10px"},
        {"marginTop": "24px", "width": "100%"},
    )


def _empty_stacked_bar_figure(message):
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        showarrow=False,
        font=dict(size=14),
        xanchor="center",
        yanchor="middle",
    )
    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    return fig


def _differential_abundance_results_style(fig=None):
    height = getattr(getattr(fig, "layout", None), "height", None) or 360
    return {
        "display": "block",
        "height": f"{int(height)}px",
        "minHeight": f"{int(height)}px",
        "width": "100%",
    }


def _x_axis_values(adata, annotation):
    return [str(value) for value in sorted_categories(adata, annotation)]


def _column_defs(x_values):
    return [
        {
            "field": value,
            "headerName": value,
            "width": 150,
            "minWidth": 120,
            "suppressMovable": False,
            "headerClass": "ag-header-cell-center",
            "resizable": True,
        }
        for value in x_values
    ]


def _stacked_bar_x_order(adata, annotation, x_axis_order):
    observed = _x_axis_values(adata, annotation)
    if not x_axis_order:
        return observed
    requested = [str(value) for value in x_axis_order]
    return [value for value in requested if value in observed] + [
        value for value in observed if value not in requested
    ]


def _stacked_bar_color_map(adata, stack_by, discrete_color_map, color_config):
    stack_categories = sorted_categories(adata, stack_by)
    discrete_palette = resolve_discrete_palette(
        discrete_color_map, len(stack_categories), default=color_config
    )
    return {
        str(category): discrete_palette[i % len(discrete_palette)]
        for i, category in enumerate(stack_categories)
    }


def _reference_population_options(adata, annotation, current_value=None):
    if not annotation:
        return [], None

    categories = [str(value) for value in sorted_categories(adata, annotation)]
    options = [{"label": value, "value": value} for value in categories]
    if current_value is not None and str(current_value) in categories:
        return options, str(current_value)

    return options, None


def register_stacked_bar_callbacks(
    app,
    adata,
    prefix,
    *,
    calculate_alr_welch,
    plot_composition_hierarchy,
    plot_composition_differential_abundance,
    plot_stacked_bar,
    palette_json,
    color_config,
    resolve_plot_adata_from_filter,
    hash_list_signature,
):
    @app.callback(
        [
            Output(f"{prefix}-composition-da-collapse", "style"),
            Output(f"{prefix}-composition-da-toggle", "children"),
        ],
        Input(f"{prefix}-composition-da-toggle", "n_clicks"),
        State(f"{prefix}-composition-da-collapse", "style"),
        prevent_initial_call=True,
    )
    def toggle_differential_abundance(_n_clicks, current_style):
        is_open = (current_style or {}).get("display") != "none"
        if is_open:
            return {"display": "none", "width": "100%"}, "▸ Differential abundance"
        return {"display": "block", "width": "100%"}, "▾ Differential abundance"

    @app.callback(
        [
            Output(f"{prefix}-composition-primary-label", "children"),
            Output(f"{prefix}-stack-by-label", "children"),
            Output(f"{prefix}-composition-value-control", "style"),
            Output(f"{prefix}-composition-selector-tooltip", "children"),
            Output(f"{prefix}-composition-da-controls", "style"),
            Output(f"{prefix}-composition-da-panel", "style"),
        ],
        Input(f"{prefix}-composition-view", "value"),
    )
    def update_composition_controls(view_mode):
        return _composition_control_state(view_mode)

    @app.callback(
        [
            Output(f"{prefix}-composition-alr-reference", "options"),
            Output(f"{prefix}-composition-alr-reference", "value"),
        ],
        [
            Input(f"{prefix}-stacked-bar-stack-by", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-visualization-workspace-tabs", "value"),
            Input(f"{prefix}-exploratory-tabs", "value"),
        ],
        State(f"{prefix}-composition-alr-reference", "value"),
    )
    def update_reference_population(
        annotation,
        filtered_data,
        active_workspace,
        active_tab,
        current_value,
    ):
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
        ):
            return no_update, no_update
        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        return _reference_population_options(
            plot_adata, annotation, current_value
        )

    # This grid lets the user reorder all groups in the selected x-axis annotation.
    @app.callback(
        [
            Output(f"{prefix}-stacked-bar-x-order-grid", "columnDefs"),
            Output(f"{prefix}-stacked-bar-x-order-grid", "rowData"),
            Output(f"{prefix}-stacked-bar-x-order-container", "style"),
        ],
        [
            Input(f"{prefix}-stacked-bar-x-group", "value"),
            Input(f"{prefix}-composition-view", "value"),
            Input(f"{prefix}-exploratory-tabs", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-visualization-workspace-tabs", "value"),
        ],
    )
    def update_x_axis_order_grid(
        annotation,
        view_mode,
        active_tab,
        filtered_data,
        active_workspace,
    ):
        order_style = {"display": "none"} if view_mode == "hierarchy" else {}
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
            or not annotation
        ):
            return [], [], order_style
        if view_mode == "hierarchy":
            return no_update, no_update, order_style

        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        return _column_defs(_x_axis_values(plot_adata, annotation)), [], order_style

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
            Input(f"{prefix}-composition-view", "value"),
            Input(f"{prefix}-norm-box", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-visualization-workspace-tabs", "value"),
            Input(f"{prefix}-exploratory-tabs", "value"),
            Input(f"{prefix}-stacked-bar-x-group", "value"),
            Input(f"{prefix}-stacked-bar-stack-by", "value"),
            Input(f"{prefix}-x-axis-column-order-store", "data"),
            Input(f"{prefix}-global-filtered-data", "data"),
        ],
        [
            State(f"{prefix}-stacked-bar-plot", "figure"),
            State(f"{prefix}-stacked-bar-rendered-key", "data"),
        ],
    )
    def update_stacked_bar(
        view_mode,
        norm,
        discrete_color_map,
        active_workspace,
        active_tab,
        annotation,
        stack_by,
        x_axis_order,
        filtered_data,
        current_figure,
        rendered_key,
    ):
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
        ):
            return no_update, no_update

        # x-axis = "Select Annotation"; stacked color layers = "Stack bars by".
        if not annotation or not stack_by:
            return _empty_stacked_bar_figure(_MISSING_STACKED_BAR_SELECTION), None

        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        if plot_adata.n_obs == 0:
            return _empty_stacked_bar_figure(
                "No cells match the Global Data Filter."
            ), None

        filter_key = hash_list_signature(
            (filtered_data or {}).get("cell_indices")
        )

        cache_key = signature(
            "composition",
            view_mode,
            norm if view_mode == "bars" else None,
            discrete_color_map,
            annotation,
            stack_by,
            x_axis_order if view_mode == "bars" else None,
            filter_key,
        )
        if cache_key == rendered_key and current_figure:
            return no_update, no_update

        # X-axis group order: dragged order from the grid, otherwise every
        # category of the selected x-axis annotation.
        final_x_order = _stacked_bar_x_order(
            plot_adata, annotation, x_axis_order
        )

        # Color the stacked layers (the "Stack bars by" variable). Resolve the
        # palette like the scatter/embedding so the same category gets the same
        # color everywhere; str() keys match plot_stacked_bar's astype(str).
        child_color_map = _stacked_bar_color_map(
            adata, stack_by, discrete_color_map, color_config
        )
        if view_mode == "hierarchy":
            parent_color_map = _stacked_bar_color_map(
                adata, annotation, discrete_color_map, color_config
            )
            fig = plot_composition_hierarchy(
                parent_meta=annotation,
                child_meta=stack_by,
                adata=plot_adata,
                parent_color_map=parent_color_map,
                child_color_map=child_color_map,
                parent_order=final_x_order,
            )
        else:
            fig = plot_stacked_bar(
                x_meta=annotation,
                y_meta=stack_by,
                norm=norm,
                adata=plot_adata,
                color_map=child_color_map,
                y_order=None,
                x_order=final_x_order,
            )
        return fig, cache_key

    @app.callback(
        [
            Output(f"{prefix}-composition-da-plot", "figure"),
            Output(f"{prefix}-composition-da-rendered-key", "data"),
            Output(f"{prefix}-composition-da-results", "style"),
        ],
        [
            Input(f"{prefix}-composition-view", "value"),
            Input(f"{prefix}-visualization-workspace-tabs", "value"),
            Input(f"{prefix}-exploratory-tabs", "value"),
            Input(f"{prefix}-stacked-bar-x-group", "value"),
            Input(f"{prefix}-stacked-bar-stack-by", "value"),
            Input(f"{prefix}-composition-sample-unit", "value"),
            Input(f"{prefix}-composition-alr-reference", "value"),
            Input(f"{prefix}-x-axis-column-order-store", "data"),
            Input(f"{prefix}-global-filtered-data", "data"),
        ],
        [
            State(f"{prefix}-composition-da-plot", "figure"),
            State(f"{prefix}-composition-da-rendered-key", "data"),
        ],
    )
    def update_differential_abundance(
        view_mode,
        active_workspace,
        active_tab,
        annotation,
        stack_by,
        sample_unit,
        reference_population,
        x_axis_order,
        filtered_data,
        current_figure,
        rendered_key,
    ):
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
            or view_mode != "bars"
        ):
            return no_update, no_update, no_update

        if not annotation or not stack_by or not sample_unit or not reference_population:
            return no_update, None, {"display": "none"}

        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        final_x_order = _stacked_bar_x_order(
            plot_adata, annotation, x_axis_order
        )
        filter_key = hash_list_signature(
            (filtered_data or {}).get("cell_indices")
        )
        cache_key = signature(
            "composition-da",
            annotation,
            stack_by,
            sample_unit,
            reference_population,
            final_x_order,
            filter_key,
        )
        if cache_key == rendered_key and current_figure:
            return no_update, no_update, no_update

        try:
            results = calculate_alr_welch(
                plot_adata,
                group_key=annotation,
                population_key=stack_by,
                sample_key=sample_unit,
                reference_population=reference_population,
                group_order=final_x_order,
            )
        except ValueError as error:
            fig = _empty_stacked_bar_figure(str(error))
            return fig, cache_key, _differential_abundance_results_style(fig)

        fig = plot_composition_differential_abundance(results)
        return fig, cache_key, _differential_abundance_results_style(fig)
