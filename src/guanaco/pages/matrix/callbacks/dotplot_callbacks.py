from dash import Input, Output, State, no_update



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
            Input(f"{prefix}-data-layer", "value"),
            Input(f"{prefix}-dotplot-standardization", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-dotplot-cluster-mode", "value"),
            Input(f"{prefix}-dotplot-cluster-method", "value"),
            Input(f"{prefix}-dotplot-cluster-metric", "value"),
            Input(f"{prefix}-dotplot-transpose", "value"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
        ],
        [
            State(f"{prefix}-dotplot", "figure"),
            State(f"{prefix}-dotplot-rendered-key", "data"),
            State(f"{prefix}-selected-cells-store", "data"),
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
        active_tab,
        current_figure,
        rendered_key,
        selected_cells,
    ):
        if active_tab != "dotplot-tab":
            return no_update, no_update

        transpose = "swap" in transpose_selection if transpose_selection else False
        plot_type_str = "matrixplot" if plot_type_selection and "matrixplot" in plot_type_selection else "dotplot"
        layer = data_layer if data_layer and data_layer != "X" else None
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
            cluster_mode=cluster_mode or "none",
            cluster_method=cluster_method or "average",
            cluster_metric=cluster_metric or "correlation",
            transpose=transpose,
            selected_cells=cells_hash,
            is_backed=bool(hasattr(adata, "isbacked") and adata.isbacked),
            n_obs=adata.n_obs,
        )
        if rendered_key == cache_key and current_figure:
            return no_update, no_update

        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return cached_fig, cache_key

        fig = plot_dot_matrix(
            adata,
            selected_genes,
            selected_annotation,
            selected_labels,
            layer=layer,
            standardization=standardization,
            color_map=color_map,
            plot_type=plot_type_str,
            cluster=cluster_mode or "none",
            method=cluster_method or "average",
            metric=cluster_metric or "correlation",
            transpose=transpose,
            selected_cells=selected_cells,
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
    def toggle_dotplot_options(n_clicks, is_open):
        now_open = not is_open
        label = "▾ More options" if now_open else "▸ More options"
        return now_open, label
