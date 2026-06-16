from dataclasses import dataclass

from dash import Input, Output, State, no_update

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import sorted_categories


_HEATMAP_TAB = "heatmap-tab"
_DEFAULT_MATRIX_LAYER = "X"
_NO_SECONDARY_ANNOTATION = "None"
_HEATMAP_BOUNDARY = 1


@dataclass(frozen=True)
class _HeatmapRequest:
    selected_genes: object
    selected_annotation: object
    selected_labels: object
    standardization: object
    data_layer: object
    heatmap_color: object
    secondary_annotation: object
    discrete_color_map: object
    secondary_colormap: object
    cells_hash: object
    active_tab: object
    selected_cells: object


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


def _resolve_layer(data_layer):
    return data_layer if data_layer and data_layer != _DEFAULT_MATRIX_LAYER else None


def _is_backed(adata):
    return bool(hasattr(adata, "isbacked") and adata.isbacked)


def _uses_secondary_annotation(primary_annotation, secondary_annotation):
    return (
        bool(secondary_annotation)
        and secondary_annotation != _NO_SECONDARY_ANNOTATION
        and secondary_annotation != primary_annotation
    )


def _label_color_map(adata, annotation, palette_name):
    if not annotation or not palette_name:
        return None
    unique_labels = sorted_categories(adata, annotation)
    palette = resolve_discrete_palette(palette_name, len(unique_labels))
    if not palette:
        return None
    return {label: palette[i % len(palette)] for i, label in enumerate(unique_labels)}


def _heatmap_cache_key(
    make_cache_key,
    hash_list_signature,
    adata,
    request,
):
    return make_cache_key(
        "heatmap",
        adata,
        selected_genes=hash_list_signature(request.selected_genes),
        selected_annotation=request.selected_annotation,
        selected_labels=hash_list_signature(request.selected_labels),
        standardization=request.standardization,
        data_layer=request.data_layer,
        heatmap_color=request.heatmap_color,
        secondary_annotation=request.secondary_annotation,
        discrete_color_map=request.discrete_color_map,
        secondary_colormap=request.secondary_colormap,
        selected_cells=request.cells_hash,
        is_backed=_is_backed(adata),
        n_obs=adata.n_obs,
    )


def _heatmap_plan_signature(
    make_cache_key,
    hash_list_signature,
    adata,
    request,
):
    # Gene-independent so unchanged per-gene work stays cacheable as the gene list changes.
    return make_cache_key(
        "heatmap-plan",
        adata,
        selected_annotation=request.selected_annotation,
        selected_labels=hash_list_signature(request.selected_labels),
        standardization=request.standardization,
        data_layer=request.data_layer,
        selected_cells=request.cells_hash,
        n_obs=adata.n_obs,
    )


def _heatmap_label_color_maps(adata, request):
    primary_color_map = _label_color_map(
        adata,
        request.selected_annotation,
        request.discrete_color_map,
    )
    secondary_color_map = None
    if _uses_secondary_annotation(request.selected_annotation, request.secondary_annotation):
        secondary_color_map = _label_color_map(
            adata,
            request.secondary_annotation,
            request.secondary_colormap,
        )
    return primary_color_map, secondary_color_map


def _heatmap_kwargs(
    request,
    *,
    plot_adata,
    layer,
    label_color_maps,
    plan_sig,
    adata_obs,
    color_config,
):
    groupby1_label_color_map, groupby2_label_color_map = label_color_maps
    return dict(
        adata=plot_adata,
        genes=request.selected_genes,
        groupby1=request.selected_annotation,
        groupby2=request.secondary_annotation
        if _uses_secondary_annotation(request.selected_annotation, request.secondary_annotation)
        else None,
        labels=request.selected_labels,
        standardization=request.standardization,
        layer=layer,
        boundary=_HEATMAP_BOUNDARY,
        color_map=request.heatmap_color,
        groupby1_label_color_map=groupby1_label_color_map,
        groupby2_label_color_map=groupby2_label_color_map,
        cache_sig=plan_sig,
        adata_obs=adata_obs,
        data_already_filtered=not bool(request.selected_labels),
        color_config=color_config,
    )


def _build_heatmap_figure(
    request,
    adata,
    *,
    filter_data,
    plot_unified_heatmap,
    make_cache_key,
    hash_list_signature,
    color_config,
):
    plan_sig = _heatmap_plan_signature(
        make_cache_key,
        hash_list_signature,
        adata,
        request,
    )
    # Keep the source AnnData stable when possible so gene-vector cache hits are maximized.
    # filter_data returns the same cached view for the same selected_cells, avoiding
    # repeated obs/var DataFrame slicing for purely cosmetic parameter changes.
    plot_adata = filter_data(adata, None, None, request.selected_cells)
    common_kwargs = _heatmap_kwargs(
        request,
        plot_adata=plot_adata,
        layer=_resolve_layer(request.data_layer),
        label_color_maps=_heatmap_label_color_maps(adata, request),
        plan_sig=plan_sig,
        adata_obs=adata.obs,
        color_config=color_config,
    )
    return plot_unified_heatmap(**common_kwargs)


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
        if secondary_annotation and secondary_annotation != _NO_SECONDARY_ANNOTATION:
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
        request = _HeatmapRequest(
            selected_genes, selected_annotation, selected_labels, standardization,
            data_layer, heatmap_color, secondary_annotation, discrete_color_map,
            secondary_colormap, cells_hash, active_tab, selected_cells,
        )
        if request.active_tab != _HEATMAP_TAB:
            return no_update, no_update

        cache_key = _heatmap_cache_key(
            make_cache_key,
            hash_list_signature,
            adata,
            request,
        )
        if rendered_key == cache_key and current_figure:
            return no_update, no_update

        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return cached_fig, cache_key

        fig = _build_heatmap_figure(
            request,
            adata,
            filter_data=filter_data,
            plot_unified_heatmap=plot_unified_heatmap,
            make_cache_key=make_cache_key,
            hash_list_signature=hash_list_signature,
            color_config=color_config,
        )
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
