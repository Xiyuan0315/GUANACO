import plotly.graph_objects as go
from dash import Input, Output, State, no_update

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import (
    SELECTION_GROUP,
    SELECTION_GROUP_LABEL,
    selected_cell_view,
    selection_group_context,
    sorted_categories,
)
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


def _option_enabled(selection, option):
    if isinstance(selection, (list, tuple, set)):
        return option in selection
    return selection == option


def _composition_view_mode(selection):
    return "hierarchy" if _option_enabled(selection, "hierarchy") else "bars"


def _composition_control_state(view_selection):
    if _composition_view_mode(view_selection) == "hierarchy":
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


def _x_axis_values(adata, annotation, group_values=None):
    return [
        str(value)
        for value in sorted_categories(
            adata,
            annotation,
            overrides=group_values,
        )
    ]


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


def _stacked_bar_x_order(
    adata,
    annotation,
    x_axis_order,
    group_values=None,
):
    observed = _x_axis_values(adata, annotation, group_values)
    if not x_axis_order:
        return observed
    requested = [str(value) for value in x_axis_order]
    return [value for value in requested if value in observed] + [
        value for value in observed if value not in requested
    ]


def _stacked_bar_color_map(
    adata,
    stack_by,
    discrete_color_map,
    color_config,
    group_values=None,
):
    stack_categories = sorted_categories(
        adata,
        stack_by,
        overrides=group_values,
    )
    discrete_palette = resolve_discrete_palette(
        discrete_color_map, len(stack_categories), default=color_config
    )
    return {
        str(category): discrete_palette[i % len(discrete_palette)]
        for i, category in enumerate(stack_categories)
    }


def _reference_population_options(
    adata,
    annotation,
    current_value=None,
    group_values=None,
):
    if not annotation:
        return [], None

    categories = [
        str(value)
        for value in sorted_categories(
            adata,
            annotation,
            overrides=group_values,
        )
    ]
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
    def resolve_selection_context(
        filtered_data,
        selected_cells,
        highlighted_cells,
    ):
        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        if selected_cells:
            return selected_cell_view(plot_adata, selected_cells), None
        if highlighted_cells:
            plot_adata, selection_values = selection_group_context(
                plot_adata,
                highlighted_cells,
            )
            return plot_adata, {SELECTION_GROUP: selection_values}
        return plot_adata, None

    @app.callback(
        [
            Output(f"{prefix}-composition-da-collapse", "style"),
            Output(f"{prefix}-composition-da-toggle-label", "children"),
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
    def update_composition_controls(view_selection):
        return _composition_control_state(view_selection)

    @app.callback(
        Output(f"{prefix}-stacked-bar-x-group", "options"),
        Output(f"{prefix}-stacked-bar-x-group", "value"),
        Output(f"{prefix}-stacked-bar-stack-by", "options"),
        Output(f"{prefix}-stacked-bar-stack-by", "value"),
        Input(f"{prefix}-selection-group-hash", "data"),
        State(f"{prefix}-stacked-bar-x-group", "options"),
        State(f"{prefix}-stacked-bar-x-group", "value"),
        State(f"{prefix}-stacked-bar-stack-by", "options"),
        State(f"{prefix}-stacked-bar-stack-by", "value"),
    )
    def sync_lasso_grouping(
        selection_group_hash,
        x_options,
        x_value,
        stack_options,
        stack_value,
    ):
        selection_option = {
            "label": SELECTION_GROUP_LABEL,
            "value": SELECTION_GROUP,
        }
        base_x = [
            option
            for option in (x_options or [])
            if option.get("value") != SELECTION_GROUP
        ]
        base_stack = [
            option
            for option in (stack_options or [])
            if option.get("value") != SELECTION_GROUP
        ]
        if selection_group_hash:
            return (
                [selection_option, *base_x],
                x_value,
                [selection_option, *base_stack],
                SELECTION_GROUP,
            )

        next_x = (
            base_x[0]["value"]
            if x_value == SELECTION_GROUP and base_x
            else x_value
        )
        if stack_value == SELECTION_GROUP:
            next_stack = (
                base_stack[1]["value"]
                if len(base_stack) > 1
                else (base_stack[0]["value"] if base_stack else None)
            )
        else:
            next_stack = stack_value
        return base_x, next_x, base_stack, next_stack

    @app.callback(
        [
            Output(f"{prefix}-composition-alr-reference", "options"),
            Output(f"{prefix}-composition-alr-reference", "value"),
        ],
        [
            Input(f"{prefix}-stacked-bar-stack-by", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
            Input(f"{prefix}-visualization-workspace-tabs", "value"),
            Input(f"{prefix}-exploratory-tabs", "value"),
        ],
        [
            State(f"{prefix}-selected-cells-store", "data"),
            State(f"{prefix}-selection-group-store", "data"),
            State(f"{prefix}-composition-alr-reference", "value"),
        ],
    )
    def update_reference_population(
        annotation,
        filtered_data,
        _selected_cells_hash,
        _selection_group_hash,
        active_workspace,
        active_tab,
        selected_cells,
        highlighted_cells,
        current_value,
    ):
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
        ):
            return no_update, no_update
        plot_adata, group_values = resolve_selection_context(
            filtered_data,
            selected_cells,
            highlighted_cells,
        )
        return _reference_population_options(
            plot_adata,
            annotation,
            current_value,
            group_values,
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
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
            Input(f"{prefix}-visualization-workspace-tabs", "value"),
        ],
        [
            State(f"{prefix}-selected-cells-store", "data"),
            State(f"{prefix}-selection-group-store", "data"),
        ],
    )
    def update_x_axis_order_grid(
        annotation,
        view_selection,
        active_tab,
        filtered_data,
        _selected_cells_hash,
        _selection_group_hash,
        active_workspace,
        selected_cells,
        highlighted_cells,
    ):
        view_mode = _composition_view_mode(view_selection)
        order_style = {"display": "none"} if view_mode == "hierarchy" else {}
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
            or not annotation
        ):
            return [], [], order_style
        if view_mode == "hierarchy":
            return no_update, no_update, order_style

        plot_adata, group_values = resolve_selection_context(
            filtered_data,
            selected_cells,
            highlighted_cells,
        )
        return (
            _column_defs(
                _x_axis_values(plot_adata, annotation, group_values)
            ),
            [],
            order_style,
        )

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
            Input(f"{prefix}-composition-swap-axes", "value"),
            Input(f"{prefix}-norm-box", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-visualization-workspace-tabs", "value"),
            Input(f"{prefix}-exploratory-tabs", "value"),
            Input(f"{prefix}-stacked-bar-x-group", "value"),
            Input(f"{prefix}-stacked-bar-stack-by", "value"),
            Input(f"{prefix}-x-axis-column-order-store", "data"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
        ],
        [
            State(f"{prefix}-selected-cells-store", "data"),
            State(f"{prefix}-selection-group-store", "data"),
            State(f"{prefix}-stacked-bar-plot", "figure"),
            State(f"{prefix}-stacked-bar-rendered-key", "data"),
        ],
    )
    def update_stacked_bar(
        view_selection,
        swap_axes_selection,
        norm,
        discrete_color_map,
        active_workspace,
        active_tab,
        annotation,
        stack_by,
        x_axis_order,
        filtered_data,
        selected_cells_hash,
        selection_group_hash,
        selected_cells,
        highlighted_cells,
        current_figure,
        rendered_key,
    ):
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
        ):
            return no_update, no_update

        view_mode = _composition_view_mode(view_selection)
        swap_axes = _option_enabled(swap_axes_selection, "swap")

        # x-axis = "Select Annotation"; stacked color layers = "Stack bars by".
        if not annotation or not stack_by:
            return _empty_stacked_bar_figure(_MISSING_STACKED_BAR_SELECTION), None

        plot_adata, group_values = resolve_selection_context(
            filtered_data,
            selected_cells,
            highlighted_cells,
        )
        if plot_adata.n_obs == 0:
            return _empty_stacked_bar_figure(
                "No cells match the active filters."
            ), None

        filter_key = hash_list_signature(
            (filtered_data or {}).get("cell_indices")
        )

        cache_key = signature(
            "composition",
            view_mode,
            swap_axes if view_mode == "bars" else None,
            norm if view_mode == "bars" else None,
            discrete_color_map,
            annotation,
            stack_by,
            x_axis_order if view_mode == "bars" else None,
            filter_key,
            selected_cells_hash,
            selection_group_hash,
        )
        if cache_key == rendered_key and current_figure:
            return no_update, no_update

        # X-axis group order: dragged order from the grid, otherwise every
        # category of the selected x-axis annotation.
        final_x_order = _stacked_bar_x_order(
            plot_adata,
            annotation,
            x_axis_order,
            group_values,
        )

        # Color the stacked layers (the "Stack bars by" variable). Resolve the
        # palette like the scatter/embedding so the same category gets the same
        # color everywhere; str() keys match plot_stacked_bar's astype(str).
        child_color_map = _stacked_bar_color_map(
            plot_adata,
            stack_by,
            discrete_color_map,
            color_config,
            group_values,
        )
        if view_mode == "hierarchy":
            parent_color_map = _stacked_bar_color_map(
                plot_adata,
                annotation,
                discrete_color_map,
                color_config,
                group_values,
            )
            fig = plot_composition_hierarchy(
                parent_meta=annotation,
                child_meta=stack_by,
                adata=plot_adata,
                parent_color_map=parent_color_map,
                child_color_map=child_color_map,
                parent_order=final_x_order,
                group_values=group_values,
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
                swap_axes=swap_axes,
                group_values=group_values,
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
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
        ],
        [
            State(f"{prefix}-selected-cells-store", "data"),
            State(f"{prefix}-selection-group-store", "data"),
            State(f"{prefix}-composition-da-plot", "figure"),
            State(f"{prefix}-composition-da-rendered-key", "data"),
        ],
    )
    def update_differential_abundance(
        view_selection,
        active_workspace,
        active_tab,
        annotation,
        stack_by,
        sample_unit,
        reference_population,
        x_axis_order,
        filtered_data,
        selected_cells_hash,
        selection_group_hash,
        selected_cells,
        highlighted_cells,
        current_figure,
        rendered_key,
    ):
        view_mode = _composition_view_mode(view_selection)
        if (
            active_workspace != _EXPLORATORY_WORKSPACE
            or active_tab != _STACKED_BAR_TAB
            or view_mode != "bars"
        ):
            return no_update, no_update, no_update

        if not annotation or not stack_by or not sample_unit or not reference_population:
            return no_update, None, {"display": "none"}

        plot_adata, group_values = resolve_selection_context(
            filtered_data,
            selected_cells,
            highlighted_cells,
        )
        final_x_order = _stacked_bar_x_order(
            plot_adata,
            annotation,
            x_axis_order,
            group_values,
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
            selected_cells_hash,
            selection_group_hash,
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
                group_values=group_values,
            )
        except ValueError as error:
            fig = _empty_stacked_bar_figure(str(error))
            return fig, cache_key, _differential_abundance_results_style(fig)

        fig = plot_composition_differential_abundance(results)
        return fig, cache_key, _differential_abundance_results_style(fig)
