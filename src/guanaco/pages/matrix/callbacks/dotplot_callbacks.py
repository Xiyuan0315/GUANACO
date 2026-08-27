from dash import Input, Output, State, no_update

from guanaco.utils.obs_utils import SELECTION_GROUP, selection_group_context


_DOTPLOT_TAB = "dotplot-tab"
_DEFAULT_MATRIX_LAYER = "X"
_DEFAULT_CLUSTER_MODE = "none"
_DEFAULT_CLUSTER_METHOD = "average"
_DEFAULT_CLUSTER_METRIC = "correlation"


def _resolve_layer(data_layer):
    return data_layer if data_layer and data_layer != _DEFAULT_MATRIX_LAYER else None


def _option_enabled(values, option):
    return bool(values and option in values)


def _plot_type(plot_type_selection):
    return (
        "matrixplot"
        if _option_enabled(plot_type_selection, "matrixplot")
        else "dotplot"
    )


def _cluster_settings(cluster_mode, cluster_method, cluster_metric):
    return (
        cluster_mode or _DEFAULT_CLUSTER_MODE,
        cluster_method or _DEFAULT_CLUSTER_METHOD,
        cluster_metric or _DEFAULT_CLUSTER_METRIC,
    )


def register_dotplot_callbacks(
    app,
    adata,
    prefix,
    *,
    filter_data,
    plot_dot_matrix,
    make_cache_key,
    hash_list_signature,
    cached_figure_get,
    cached_figure_set,
    multiomics_source=None,
):
    @app.callback(
        [
            Output(f"{prefix}-dotplot", "figure"),
            Output(f"{prefix}-dotplot-rendered-key", "data"),
        ],
        [
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-plot-type-switch", "value"),
            Input(f"{prefix}-data-layer", "data"),
            Input(f"{prefix}-dotplot-standardization", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-dotplot-cluster-mode", "value"),
            Input(f"{prefix}-dotplot-cluster-method", "value"),
            Input(f"{prefix}-dotplot-cluster-metric", "value"),
            Input(f"{prefix}-dotplot-transpose", "value"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
            Input(f"{prefix}-marker-tabs", "value"),
        ],
        [
            State(f"{prefix}-dotplot", "figure"),
            State(f"{prefix}-dotplot-rendered-key", "data"),
            State(f"{prefix}-selected-cells-store", "data"),
            State(f"{prefix}-selection-group-store", "data"),
        ],
    )
    def update_dotplot(
        selected_genes,
        selected_annotation,
        selected_labels,
        plot_type_selection,
        data_layer,
        standardization,
        color_map,
        cluster_mode,
        cluster_method,
        cluster_metric,
        transpose_selection,
        cells_hash,
        selection_group_hash,
        active_tab,
        current_figure,
        rendered_key,
        selected_cells,
        highlighted_cells,
    ):
        if active_tab != _DOTPLOT_TAB:
            return no_update, no_update

        transpose = _option_enabled(transpose_selection, "swap")
        plot_type_str = _plot_type(plot_type_selection)
        cluster_mode, cluster_method, cluster_metric = _cluster_settings(
            cluster_mode, cluster_method, cluster_metric
        )
        layer = _resolve_layer(data_layer)
        cache_key = make_cache_key(
            "dotplot",
            adata,
            selected_genes=hash_list_signature(selected_genes),
            selected_annotation=selected_annotation,
            selected_labels=hash_list_signature(selected_labels),
            plot_type=plot_type_str,
            data_layer=data_layer,
            standardization=standardization,
            color_map=color_map,
            cluster_mode=cluster_mode,
            cluster_method=cluster_method,
            cluster_metric=cluster_metric,
            transpose=transpose,
            selected_cells=cells_hash,
            selection_group=selection_group_hash,
            is_backed=bool(hasattr(adata, "isbacked") and adata.isbacked),
            n_obs=adata.n_obs,
        )
        if rendered_key == cache_key and current_figure:
            return no_update, no_update

        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return cached_fig, cache_key

        plot_adata = (
            multiomics_source.materialize(selected_genes)
            if multiomics_source is not None
            else adata
        )
        group_values = None
        if selected_annotation == SELECTION_GROUP and highlighted_cells:
            plot_adata, group_values = selection_group_context(
                plot_adata, highlighted_cells
            )
        fig = plot_dot_matrix(
            plot_adata,
            selected_genes,
            selected_annotation,
            selected_labels,
            layer=layer,
            standardization=standardization,
            color_map=color_map,
            plot_type=plot_type_str,
            cluster=cluster_mode,
            method=cluster_method,
            metric=cluster_metric,
            transpose=transpose,
            selected_cells=selected_cells,
            group_values=group_values,
        )
        cached_figure_set(cache_key, fig)
        return fig, cache_key

    @app.callback(
        Output(f"{prefix}-dotplot-options-collapse", "is_open"),
        Output(f"{prefix}-dotplot-options-toggle", "children"),
        Input(f"{prefix}-dotplot-options-toggle", "n_clicks"),
        State(f"{prefix}-dotplot-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_dotplot_options(_n_clicks, is_open):
        now_open = not is_open
        label = "▾ More options" if now_open else "▸ More options"
        return now_open, label
