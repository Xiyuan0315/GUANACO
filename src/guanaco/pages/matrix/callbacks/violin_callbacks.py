from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import sorted_categories


_CURRENT_CACHE_KEY = "current_key"
_MAX_VIOLIN_CACHE_ENTRIES = 10

_MODE_EXPLANATIONS = {
    "mode1": "Compare expression across groups in obs1 only. Obs2 will be ignored.",
    "mode2": "Within each obs1 group on the x-axis, split/group the violin by obs2 and compare.",
    "mode3": "Linear model treating obs2 as a confounder: expression ~ obs1 + obs2",
    "mode4": "Mixed model treating obs2 as a random effect (e.g. donor, batch, replicate) to account for non-independent samples.",
}

_BASE_TEST_OPTIONS = (("Auto", "auto"), ("None", "none"))
_TWO_LEVEL_TEST_OPTIONS = (("Mann-Whitney U", "mwu-test"), ("T-test", "ttest"))
_MULTI_LEVEL_TEST_OPTIONS = (("Kruskal-Wallis", "kw-test"), ("ANOVA", "anova"))
_MODEL_TEST_OPTIONS = {
    "mode3": (
        ("Linear Model", "linear-model"),
        ("Linear Model with Interaction", "linear-model-interaction"),
    ),
    "mode4": (("Mixed Model", "mixed-model"),),
}


def _dropdown_options(option_pairs):
    return [{"label": label, "value": value} for label, value in option_pairs]


def _resolve_layer(data_layer):
    return data_layer if data_layer and data_layer != "X" else None


def _has_checklist_value(values, value):
    return value in values if values else False


def _annotation_level_count(adata, annotation):
    if not annotation or annotation not in adata.obs:
        return 0
    return adata.obs[annotation].nunique()


def _test_options_for_mode(adata, mode, meta1, meta2):
    options = _dropdown_options(_BASE_TEST_OPTIONS)
    if mode == "mode1":
        n_levels = _annotation_level_count(adata, meta1)
        if n_levels == 0:
            return options
        comparison_options = _TWO_LEVEL_TEST_OPTIONS if n_levels == 2 else _MULTI_LEVEL_TEST_OPTIONS
        return options + _dropdown_options(comparison_options)
    if mode == "mode2":
        n_levels = _annotation_level_count(adata, meta2)
        if n_levels == 0:
            return options
        comparison_options = _TWO_LEVEL_TEST_OPTIONS if n_levels == 2 else _MULTI_LEVEL_TEST_OPTIONS
        return options + _dropdown_options(comparison_options)
    return options + _dropdown_options(_MODEL_TEST_OPTIONS.get(mode, ()))


def _resolve_violin1_color_map(adata, annotation, palette_name):
    if not palette_name or not annotation:
        return None
    unique_labels = sorted_categories(adata, annotation)
    discrete_palette = resolve_discrete_palette(palette_name, len(unique_labels))
    return {label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels)}


def _resolve_violin2_meta2(mode, meta2):
    if mode == "mode1" or meta2 == "none":
        return None
    return meta2


def _resolve_violin2_palette(filtered_adata, meta1, meta2, palette_name, color_config):
    n_colors = max(
        (_annotation_level_count(filtered_adata, annotation) for annotation in (meta1, meta2) if annotation),
        default=0,
    )
    return resolve_discrete_palette(palette_name, n_colors, default=color_config)


def _size_violin1_figure(fig, selected_genes, selected_labels):
    num_genes = len(selected_genes) if selected_genes else 0
    num_categories = len(selected_labels) if selected_labels else 0
    fig.update_layout(
        height=min(1000, max(300, 80 * num_genes)),
        width=min(500, max(200, 110 * num_categories)),
        margin=dict(l=130, r=10, t=30, b=30),
    )


def _prune_violin_cache(cache_data):
    figure_keys = [key for key in cache_data if key != _CURRENT_CACHE_KEY]
    max_figures = max(_MAX_VIOLIN_CACHE_ENTRIES - 1, 0)
    for key in figure_keys[:-max_figures]:
        cache_data.pop(key, None)


def _store_current_violin_figure(cache_data, cache_key, fig):
    cache_data[cache_key] = fig.to_dict()
    cache_data[_CURRENT_CACHE_KEY] = cache_key
    _prune_violin_cache(cache_data)
    return cache_data


def register_violin_callbacks(
    app,
    adata,
    prefix,
    *,
    filter_data,
    plot_violin1,
    plot_violin2_new,
    palette_json,
    var_names,
    var_names_lower,
    color_config=None,
):
    @app.callback(
        Output(f"{prefix}-violin-plot-cache-store", "data"),
        [
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-data-layer", "value"),
            Input(f"{prefix}-show-box1", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
        ],
        [
            State(f"{prefix}-violin-plot-cache-store", "data"),
            State(f"{prefix}-selected-cells-store", "data"),
        ],
    )
    def update_violin_cache(
        selected_genes,
        selected_annotation,
        selected_labels,
        data_layer,
        show_box_plot,
        discrete_color_map,
        cells_hash,
        active_tab,
        current_cache,
        selected_cells,
    ):
        # Lazy: only build the violin figure when its tab is active, so the default
        # (e.g. dot plot) view doesn't also pay to compute violins it shares inputs
        # with. The tab itself is an Input, so switching to violin triggers the build.
        if active_tab != "violin-tab":
            return no_update
        layer = _resolve_layer(data_layer)
        cache_key = f"{selected_genes}_{selected_annotation}_{selected_labels}_{data_layer}_{show_box_plot}_{discrete_color_map}_{cells_hash}"

        if current_cache is None:
            current_cache = {}

        if current_cache.get(_CURRENT_CACHE_KEY) == cache_key:
            return current_cache

        color_map = _resolve_violin1_color_map(adata, selected_annotation, discrete_color_map)
        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)

        fig = plot_violin1(
            filtered_adata,
            selected_genes,
            selected_annotation,
            labels=selected_labels,
            layer=layer,
            show_box=_has_checklist_value(show_box_plot, "show"),
            groupby_label_color_map=color_map,
            adata_obs=adata.obs,
            data_already_filtered=True,
        )

        _size_violin1_figure(fig, selected_genes, selected_labels)
        return _store_current_violin_figure(current_cache, cache_key, fig)

    @app.callback(
        [
            Output(f"{prefix}-violin-plot1", "figure"),
            Output(f"{prefix}-violin1-rendered-key", "data"),
        ],
        [
            Input(f"{prefix}-violin-plot-cache-store", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
        ],
        [
            State(f"{prefix}-violin-plot1", "figure"),
            State(f"{prefix}-violin1-rendered-key", "data"),
        ],
    )
    def display_violin1(cache_data, active_tab, current_figure, rendered_key):
        if active_tab != "violin-tab":
            return no_update, no_update

        current_key = cache_data.get(_CURRENT_CACHE_KEY) if cache_data else None
        # Already showing the figure for this key: don't redraw on a tab switch.
        if current_key and current_key == rendered_key and current_figure:
            return no_update, no_update
        if current_key and current_key in cache_data:
            # Return the cached figure dict as-is: Dash accepts a plain dict for a
            # Graph 'figure', which skips the expensive go.Figure(...) re-validation of
            # the violin figure on every switch to the tab.
            return cache_data[current_key], current_key
        return no_update, no_update

    @app.callback(Output(f"{prefix}-mode-explanation", "children"), Input(f"{prefix}-mode-selection", "value"))
    def update_mode_explanation(mode):
        return _MODE_EXPLANATIONS.get(mode, "")

    @app.callback(
        Output(f"{prefix}-test-method-selection", "options"),
        Output(f"{prefix}-test-method-selection", "value"),
        [Input(f"{prefix}-mode-selection", "value"), Input(f"{prefix}-meta1-selection", "value"), Input(f"{prefix}-meta2-selection", "value")],
    )
    def update_test_methods(mode, meta1, meta2):
        return _test_options_for_mode(adata, mode, meta1, meta2), "auto"

    @app.callback(
        Output(f"{prefix}-violin2-gene-selection", "options"),
        Input(f"{prefix}-violin2-gene-selection", "search_value"),
    )
    def update_violin_genes_dropdown(search_value):
        if not search_value:
            raise PreventUpdate
        q = search_value.lower()
        matching_labels = [label for label, label_l in zip(var_names, var_names_lower) if q in label_l]
        return [{"label": label, "value": label} for label in matching_labels[:10]]

    @app.callback(
        [Output(f"{prefix}-meta2-selection", "disabled"), Output(f"{prefix}-meta2-selection", "value")],
        Input(f"{prefix}-mode-selection", "value"),
    )
    def toggle_meta2_dropdown(mode):
        if mode == "mode1":
            return True, "none"
        return False, no_update

    @app.callback(
        Output(f"{prefix}-split-violin-options-collapse", "is_open"),
        Output(f"{prefix}-split-violin-options-toggle", "children"),
        Input(f"{prefix}-split-violin-options-toggle", "n_clicks"),
        State(f"{prefix}-split-violin-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_split_violin_options(_n_clicks, is_open):
        now_open = not is_open
        label = "▾ More options" if now_open else "▸ More options"
        return now_open, label

    @app.callback(
        Output(f"{prefix}-violin-plot2", "figure"),
        [
            Input(f"{prefix}-violin2-gene-selection", "value"),
            Input(f"{prefix}-meta1-selection", "value"),
            Input(f"{prefix}-meta2-selection", "value"),
            Input(f"{prefix}-mode-selection", "value"),
            Input(f"{prefix}-test-method-selection", "value"),
            Input(f"{prefix}-show-box2", "value"),
            Input(f"{prefix}-data-layer", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
        ],
    )
    def update_violin2(gene_selection, meta1, meta2, mode, test_method, show_box2, data_layer, selected_palette_name, selected_cells):
        # Not tab-gated: all of this plot's controls live on its own tab, so it only
        # recomputes in response to its own inputs -- it was never part of the
        # eager-load issue that the (left-panel-driven) violin1 cache had.
        layer = _resolve_layer(data_layer)
        if selected_cells:
            filtered_adata = filter_data(adata, None, None, selected_cells)
        else:
            filtered_adata = adata

        meta2 = _resolve_violin2_meta2(mode, meta2)

        if mode in ["mode2", "mode3", "mode4"] and meta2 is None:
            raise PreventUpdate

        # Resolve the categorical palette the same way every other plot does -- the
        # selected discrete colormap, falling back to the dataset color_config -- so
        # the violins share the app's default colors instead of a private palette.
        palette = _resolve_violin2_palette(filtered_adata, meta1, meta2, selected_palette_name, color_config)

        fig = plot_violin2_new(
            filtered_adata,
            key=gene_selection,
            meta1=meta1,
            meta2=meta2,
            mode=mode,
            layer=layer,
            show_box=_has_checklist_value(show_box2, "show"),
            test_method=test_method,
            labels=None,
            color_map=None,
            palette=palette,
        )
        return fig
