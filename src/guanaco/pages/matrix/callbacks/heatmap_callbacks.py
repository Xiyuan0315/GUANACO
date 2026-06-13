from dash import Input, Output, State, no_update

from guanaco.utils.colors import resolve_discrete_palette


# Clientside double-click reset for the heatmap. Plotly's own reset is unreliable
# here: the axes use constrain='domain' and are linked via shared_xaxes, so the
# default 'reset+autosize' can leave the view partly zoomed. On a double-click
# (relayout emits *.autorange:true) we force every axis back to the figure's
# configured range -- the full heatmap extent -- which always restores cleanly.
# Output is a throwaway store; the real work is the Plotly.relayout side effect.
_HEATMAP_RESET_JS = """
function(relayout, figure) {
    const noUpdate = window.dash_clientside.no_update;
    if (!relayout || !figure) return noUpdate;
    // Double-click reset emits *.autorange:true. Ignore zoom/pan (emit ranges) and
    // responsive resize (emits 'autosize'), so this only fires on an actual reset.
    const isReset = relayout['xaxis.autorange'] === true
                 || relayout['yaxis.autorange'] === true;
    if (!isReset) return noUpdate;
    // Debounce so the relayout we trigger below doesn't re-enter this handler.
    const now = Date.now();
    if (window.__guanacoHeatmapReset && (now - window.__guanacoHeatmapReset) < 350) return noUpdate;
    window.__guanacoHeatmapReset = now;

    const wrap = document.getElementById('__HEATMAP_ID__');
    if (!wrap || !window.Plotly) return noUpdate;
    const gd = wrap.classList.contains('js-plotly-plot') ? wrap : wrap.querySelector('.js-plotly-plot');
    if (!gd) return noUpdate;

    // figure.layout holds the server-set (original) ranges -- the full view.
    const lay = (figure && figure.layout) || {};
    const upd = {};
    ['xaxis','xaxis2','xaxis3','yaxis','yaxis2','yaxis3'].forEach(function(ax) {
        if (lay[ax] && lay[ax].range) {
            upd[ax + '.range'] = lay[ax].range.slice();
            upd[ax + '.autorange'] = false;
        }
    });
    if (Object.keys(upd).length > 0) window.Plotly.relayout(gd, upd);
    return noUpdate;
}
"""


def register_heatmap_callbacks(
    app,
    adata,
    prefix,
    *,
    filter_data,
    plot_unified_heatmap,
    palette_json,
    color_config=None,
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
        [
            Output(f"{prefix}-heatmap", "figure"),
            Output(f"{prefix}-heatmap-rendered-key", "data"),
        ],
        [
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-heatmap-standardization", "value"),
            Input(f"{prefix}-data-layer", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-heatmap-label-dropdown", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-heatmap-secondary-colormap-dropdown", "value"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
        ],
        [
            State(f"{prefix}-heatmap", "figure"),
            State(f"{prefix}-heatmap-rendered-key", "data"),
            State(f"{prefix}-selected-cells-store", "data"),
        ],
    )
    def update_heatmap(
        selected_genes,
        selected_annotation,
        selected_labels,
        standardization,
        data_layer,
        heatmap_color,
        secondary_annotation,
        discrete_color_map,
        secondary_colormap,
        cells_hash,
        active_tab,
        current_figure,
        rendered_key,
        selected_cells,
    ):
        if active_tab != "heatmap-tab":
            return no_update, no_update

        layer = data_layer if data_layer and data_layer != "X" else None
        cache_key = make_cache_key(
            "heatmap",
            adata,
            selected_genes=hash_list_signature(selected_genes),
            selected_annotation=selected_annotation,
            selected_labels=hash_list_signature(selected_labels),
            standardization=standardization,
            data_layer=data_layer,
            heatmap_color=heatmap_color,
            secondary_annotation=secondary_annotation,
            discrete_color_map=discrete_color_map,
            secondary_colormap=secondary_colormap,
            selected_cells=cells_hash,
            is_backed=bool(hasattr(adata, "isbacked") and adata.isbacked),
            n_obs=adata.n_obs,
        )
        if rendered_key == cache_key and current_figure:
            return no_update, no_update

        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return cached_fig, cache_key

        # Gene-independent signature for the per-gene binned cache: stays stable as
        # the user adds/removes genes (so unchanged genes stay cached) but changes
        # with the cell selection / grouping / standardization / layer, even though
        # a fresh adata[selected_cells] view would otherwise change identity each render.
        plan_sig = make_cache_key(
            "heatmap-plan",
            adata,
            selected_annotation=selected_annotation,
            selected_labels=hash_list_signature(selected_labels),
            standardization=standardization,
            data_layer=data_layer,
            selected_cells=cells_hash,
            n_obs=adata.n_obs,
        )

        # Keep the source AnnData stable when possible so gene-vector cache hits are maximized.
        # filter_data returns the same cached view for the same selected_cells, avoiding
        # repeated obs/var DataFrame slicing for purely cosmetic parameter changes.
        plot_adata = filter_data(adata, None, None, selected_cells)

        groupby1_label_color_map = None
        if discrete_color_map:
            unique_labels1 = sorted(adata.obs[selected_annotation].unique())
            discrete_palette = resolve_discrete_palette(discrete_color_map, len(unique_labels1))
            groupby1_label_color_map = {
                label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels1)
            }
        groupby2_label_color_map = None
        if secondary_annotation and secondary_annotation != "None" and secondary_annotation != selected_annotation:
            unique_labels2 = sorted(adata.obs[secondary_annotation].unique())
            if secondary_colormap:
                secondary_palette = resolve_discrete_palette(secondary_colormap, len(unique_labels2))
                groupby2_label_color_map = {
                    label: secondary_palette[i % len(secondary_palette)] for i, label in enumerate(unique_labels2)
                }

        labels_need_post_filter = bool(selected_labels)
        common_kwargs = dict(
            adata=plot_adata,
            genes=selected_genes,
            groupby1=selected_annotation,
            groupby2=secondary_annotation
            if secondary_annotation and secondary_annotation != "None" and secondary_annotation != selected_annotation
            else None,
            labels=selected_labels,
            standardization=standardization,
            layer=layer,
            boundary=1,
            color_map=heatmap_color,
            groupby1_label_color_map=groupby1_label_color_map,
            groupby2_label_color_map=groupby2_label_color_map,
            cache_sig=plan_sig,
            adata_obs=adata.obs,
            data_already_filtered=not labels_need_post_filter,
            color_config=color_config,
        )
        fig = plot_unified_heatmap(**common_kwargs)
        cached_figure_set(cache_key, fig)
        return fig, cache_key

    # Clientside double-click reset (see _HEATMAP_RESET_JS).
    app.clientside_callback(
        _HEATMAP_RESET_JS.replace("__HEATMAP_ID__", f"{prefix}-heatmap"),
        Output(f"{prefix}-heatmap-reset-link", "data"),
        Input(f"{prefix}-heatmap", "relayoutData"),
        State(f"{prefix}-heatmap", "figure"),
        prevent_initial_call=True,
    )
