from dash import Input, Output, State
import plotly.graph_objects as go


def register_heatmap_callbacks(
    app,
    adata,
    prefix,
    *,
    filter_data,
    plot_unified_heatmap,
    palette_json,
    make_cache_key,
    hash_list_signature,
    cached_figure_get,
    cached_figure_set,
):
    @app.callback(
        Output(f"{prefix}-heatmap-secondary-colormap-wrapper", "style"),
        Input(f"{prefix}-heatmap-label-dropdown", "value"),
    )
    def toggle_secondary_colormap_dropdown(secondary_annotation):
        if secondary_annotation and secondary_annotation != "None":
            return {"display": "block"}
        return {"display": "none"}

    @app.callback(
        Output(f"{prefix}-heatmap", "figure"),
        [
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-heatmap-transformation", "value"),
            Input(f"{prefix}-heatmap-colorscale-dropdown", "value"),
            Input(f"{prefix}-heatmap-label-dropdown", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-heatmap-secondary-colormap-dropdown", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
        ],
        [State(f"{prefix}-heatmap", "figure")],
    )
    def update_heatmap(
        selected_genes,
        selected_annotation,
        selected_labels,
        transformation,
        heatmap_color,
        secondary_annotation,
        discrete_color_map,
        secondary_colormap,
        selected_cells,
        active_tab,
        current_figure,
    ):
        if active_tab != "heatmap-tab":
            return current_figure if current_figure else go.Figure()

        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
        groupby1_label_color_map = None
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            unique_labels1 = sorted(adata.obs[selected_annotation].unique())
            groupby1_label_color_map = {
                label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels1)
            }
        groupby2_label_color_map = None
        if secondary_annotation and secondary_annotation != "None" and secondary_annotation != selected_annotation:
            unique_labels2 = sorted(adata.obs[secondary_annotation].unique())
            if secondary_colormap:
                secondary_palette = palette_json["color_palettes"][secondary_colormap]
                groupby2_label_color_map = {
                    label: secondary_palette[i % len(secondary_palette)] for i, label in enumerate(unique_labels2)
                }

        labels_need_post_filter = bool(selected_cells) and bool(selected_labels)
        common_kwargs = dict(
            adata=filtered_adata,
            genes=selected_genes,
            groupby1=selected_annotation,
            groupby2=secondary_annotation
            if secondary_annotation and secondary_annotation != "None" and secondary_annotation != selected_annotation
            else None,
            labels=selected_labels,
            log=(transformation == "log"),
            z_score=(transformation == "z_score"),
            boundary=1,
            color_map=heatmap_color,
            groupby1_label_color_map=groupby1_label_color_map,
            groupby2_label_color_map=groupby2_label_color_map,
            transformation=transformation,
            adata_obs=adata.obs,
            data_already_filtered=not labels_need_post_filter,
        )
        return plot_unified_heatmap(**common_kwargs)
