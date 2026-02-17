from dash import Input, Output, State
import plotly.graph_objects as go


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
        Output(f"{prefix}-dotplot", "figure"),
        [
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-plot-type-switch", "value"),
            Input(f"{prefix}-dotplot-log-or-zscore", "value"),
            Input(f"{prefix}-dotplot-standardization", "value"),
            Input(f"{prefix}-dotmatrix-color-map-dropdown", "value"),
            Input(f"{prefix}-dotplot-cluster-mode", "value"),
            Input(f"{prefix}-dotplot-cluster-method", "value"),
            Input(f"{prefix}-dotplot-cluster-metric", "value"),
            Input(f"{prefix}-dotplot-transpose", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
        ],
        [State(f"{prefix}-dotplot", "figure")],
    )
    def update_dotplot(
        selected_genes,
        selected_annotation,
        selected_labels,
        plot_type_selection,
        transformation,
        standardization,
        color_map,
        cluster_mode,
        cluster_method,
        cluster_metric,
        transpose_selection,
        selected_cells,
        active_tab,
        current_figure,
    ):
        if active_tab != "dotplot-tab":
            return current_figure if current_figure else go.Figure()

        transpose = "swap" in transpose_selection if transpose_selection else False
        plot_type_str = "matrixplot" if plot_type_selection and "matrixplot" in plot_type_selection else "dotplot"
        cache_key = make_cache_key(
            "dotplot",
            adata,
            selected_genes=hash_list_signature(selected_genes),
            selected_annotation=selected_annotation,
            selected_labels=hash_list_signature(selected_labels),
            plot_type=plot_type_str,
            transformation=transformation,
            standardization=standardization,
            color_map=color_map,
            cluster_mode=cluster_mode or "none",
            cluster_method=cluster_method or "average",
            cluster_metric=cluster_metric or "correlation",
            transpose=transpose,
            selected_cells=hash_list_signature(selected_cells),
            is_backed=bool(hasattr(adata, "isbacked") and adata.isbacked),
            n_obs=adata.n_obs,
        )
        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return cached_fig

        if adata.n_obs > 10000 or (hasattr(adata, "isbacked") and adata.isbacked):
            fig = plot_dot_matrix(
                adata,
                selected_genes,
                selected_annotation,
                selected_labels,
                transformation=transformation,
                standardization=standardization,
                color_map=color_map,
                plot_type=plot_type_str,
                cluster=cluster_mode or "none",
                method=cluster_method or "average",
                metric=cluster_metric or "correlation",
                transpose=transpose,
            )
        else:
            filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
            fig = plot_dot_matrix(
                filtered_adata,
                selected_genes,
                selected_annotation,
                selected_labels,
                transformation=transformation,
                standardization=standardization,
                color_map=color_map,
                plot_type=plot_type_str,
                cluster=cluster_mode or "none",
                method=cluster_method or "average",
                metric=cluster_metric or "correlation",
                transpose=transpose,
            )
        cached_figure_set(cache_key, fig)
        return fig
