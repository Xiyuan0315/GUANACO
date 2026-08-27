from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import (
    SELECTION_GROUP,
    SELECTION_GROUP_LABEL,
    SELECTION_LABELS,
    selected_cell_view,
    selection_group_context,
    sorted_categories,
)
from guanaco.utils.search import ranked_substring_matches
from guanaco.data.loader import obs_col


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


def _annotation_level_count(adata, annotation, group_values=None):
    if not annotation:
        return 0
    if group_values and annotation in group_values:
        return group_values[annotation].nunique()
    if annotation not in adata.obs:
        return 0
    return obs_col(adata.obs, annotation).nunique()


def _test_options_for_mode(adata, mode, meta1, meta2):
    options = _dropdown_options(_BASE_TEST_OPTIONS)
    if mode == "mode1":
        n_levels = _annotation_level_count(adata, meta1)
        if n_levels == 0:
            return options
        comparison_options = (
            _TWO_LEVEL_TEST_OPTIONS if n_levels == 2 else _MULTI_LEVEL_TEST_OPTIONS
        )
        return options + _dropdown_options(comparison_options)
    if mode == "mode2":
        n_levels = _annotation_level_count(adata, meta2)
        if n_levels == 0:
            return options
        comparison_options = (
            _TWO_LEVEL_TEST_OPTIONS if n_levels == 2 else _MULTI_LEVEL_TEST_OPTIONS
        )
        return options + _dropdown_options(comparison_options)
    return options + _dropdown_options(_MODEL_TEST_OPTIONS.get(mode, ()))


def _resolve_violin1_color_map(adata, annotation, palette_name):
    if not palette_name or not annotation:
        return None
    unique_labels = (
        SELECTION_LABELS
        if annotation == SELECTION_GROUP
        else sorted_categories(adata, annotation)
    )
    discrete_palette = resolve_discrete_palette(palette_name, len(unique_labels))
    return {
        label: discrete_palette[i % len(discrete_palette)]
        for i, label in enumerate(unique_labels)
    }


def _resolve_violin2_meta2(mode, meta2):
    if mode == "mode1" or meta2 == "none":
        return None
    return meta2


def _resolve_violin2_palette(
    filtered_adata,
    meta1,
    meta2,
    palette_name,
    color_config,
    group_values=None,
):
    n_colors = max(
        (
            _annotation_level_count(filtered_adata, annotation, group_values)
            for annotation in (meta1, meta2)
            if annotation
        ),
        default=0,
    )
    return resolve_discrete_palette(palette_name, n_colors, default=color_config)


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


def _violin1_graph_style(figure):
    if isinstance(figure, dict):
        height = figure.get("layout", {}).get("height")
    else:
        height = getattr(getattr(figure, "layout", None), "height", None)
    height = int(height or 400)
    return {
        "width": "100%",
        "height": f"{height}px",
        "minHeight": "400px",
    }


def register_marker_violin_callbacks(
    app,
    adata,
    prefix,
    *,
    filter_data,
    plot_violin1,
    multiomics_source=None,
):
    @app.callback(
        Output(f"{prefix}-violin-plot-cache-store", "data"),
        [
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-data-layer", "data"),
            Input(f"{prefix}-show-box1", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
            Input(f"{prefix}-marker-tabs", "value"),
        ],
        [
            State(f"{prefix}-violin-plot-cache-store", "data"),
            State(f"{prefix}-selected-cells-store", "data"),
            State(f"{prefix}-selection-group-store", "data"),
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
        selection_group_hash,
        active_tab,
        current_cache,
        selected_cells,
        highlighted_cells,
    ):
        # Lazy: only build the violin figure when its tab is active, so the default
        # (e.g. dot plot) view doesn't also pay to compute violins it shares inputs
        # with. The tab itself is an Input, so switching to violin triggers the build.
        if active_tab != "violin-tab":
            return no_update
        layer = _resolve_layer(data_layer)
        cache_key = (
            f"{selected_genes}_{selected_annotation}_{selected_labels}_{data_layer}_"
            f"{show_box_plot}_{discrete_color_map}_{cells_hash}_{selection_group_hash}"
        )

        if current_cache is None:
            current_cache = {}

        if current_cache.get(_CURRENT_CACHE_KEY) == cache_key:
            return current_cache

        source_adata = (
            multiomics_source.materialize(selected_genes)
            if multiomics_source is not None
            else adata
        )
        color_map = _resolve_violin1_color_map(
            source_adata, selected_annotation, discrete_color_map
        )
        filtering_by_selection_group = selected_annotation == SELECTION_GROUP
        group_values = None
        if filtering_by_selection_group and highlighted_cells:
            source_adata, group_values = selection_group_context(
                source_adata, highlighted_cells
            )
        filtered_adata = filter_data(
            source_adata,
            None if filtering_by_selection_group else selected_annotation,
            None if filtering_by_selection_group else selected_labels,
            selected_cells,
        )

        fig = plot_violin1(
            filtered_adata,
            selected_genes,
            selected_annotation,
            labels=selected_labels,
            layer=layer,
            show_box=_has_checklist_value(show_box_plot, "show"),
            groupby_label_color_map=color_map,
            adata_obs=source_adata.obs,
            data_already_filtered=True,
            group_values=group_values,
        )

        return _store_current_violin_figure(current_cache, cache_key, fig)

    @app.callback(
        [
            Output(f"{prefix}-violin-plot1", "figure"),
            Output(f"{prefix}-violin1-rendered-key", "data"),
            Output(f"{prefix}-violin-plot1", "style"),
        ],
        [
            Input(f"{prefix}-violin-plot-cache-store", "data"),
            Input(f"{prefix}-marker-tabs", "value"),
        ],
        [
            State(f"{prefix}-violin-plot1", "figure"),
            State(f"{prefix}-violin1-rendered-key", "data"),
        ],
    )
    def display_violin1(cache_data, active_tab, current_figure, rendered_key):
        if active_tab != "violin-tab":
            return no_update, no_update, no_update

        current_key = cache_data.get(_CURRENT_CACHE_KEY) if cache_data else None
        # Already showing the figure for this key: don't redraw on a tab switch.
        if current_key and current_key == rendered_key and current_figure:
            return no_update, no_update, no_update
        if current_key and current_key in cache_data:
            # Return the cached figure dict as-is: Dash accepts a plain dict for a
            # Graph 'figure', which skips the expensive go.Figure(...) re-validation of
            # the violin figure on every switch to the tab.
            figure = cache_data[current_key]
            return figure, current_key, _violin1_graph_style(figure)
        return no_update, no_update, no_update


def register_comparative_violin_callbacks(
    app,
    adata,
    prefix,
    *,
    plot_violin2_new,
    var_names,
    var_names_lower,
    color_config=None,
    resolve_plot_adata_from_filter=None,
):
    if resolve_plot_adata_from_filter is None:
        def resolve_plot_adata_from_filter(_filtered):
            return adata

    @app.callback(
        Output(f"{prefix}-mode-explanation", "children"),
        Input(f"{prefix}-mode-selection", "value"),
    )
    def update_mode_explanation(mode):
        return _MODE_EXPLANATIONS.get(mode, "")

    @app.callback(
        Output(f"{prefix}-test-method-selection", "options"),
        Output(f"{prefix}-test-method-selection", "value"),
        [
            Input(f"{prefix}-mode-selection", "value"),
            Input(f"{prefix}-meta1-selection", "value"),
            Input(f"{prefix}-meta2-selection", "value"),
        ],
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
        matching_labels = ranked_substring_matches(
            var_names,
            search_value,
            limit=10,
            match_values=var_names_lower,
        )
        return [{"label": label, "value": label} for label in matching_labels]

    @app.callback(
        Output(f"{prefix}-meta2-selection", "disabled"),
        Input(f"{prefix}-mode-selection", "value"),
    )
    def toggle_meta2_dropdown(mode):
        return mode == "mode1"

    @app.callback(
        Output(f"{prefix}-meta1-selection", "options"),
        Output(f"{prefix}-meta1-selection", "value"),
        Output(f"{prefix}-meta2-selection", "options"),
        Output(f"{prefix}-meta2-selection", "value"),
        Input(f"{prefix}-selection-group-hash", "data"),
        State(f"{prefix}-meta1-selection", "options"),
        State(f"{prefix}-meta1-selection", "value"),
        State(f"{prefix}-meta2-selection", "options"),
        State(f"{prefix}-meta2-selection", "value"),
    )
    def sync_lasso_grouping(
        selection_group_hash,
        meta1_options,
        meta1_value,
        meta2_options,
        meta2_value,
    ):
        selection_option = {
            "label": SELECTION_GROUP_LABEL,
            "value": SELECTION_GROUP,
        }
        base_meta1 = [
            option
            for option in (meta1_options or [])
            if option.get("value") != SELECTION_GROUP
        ]
        base_meta2 = [
            option
            for option in (meta2_options or [])
            if option.get("value") != SELECTION_GROUP
        ]
        if selection_group_hash:
            return (
                [selection_option, *base_meta1],
                SELECTION_GROUP,
                [*base_meta2[:1], selection_option, *base_meta2[1:]],
                meta2_value,
            )

        next_meta1 = (
            base_meta1[0]["value"]
            if meta1_value == SELECTION_GROUP and base_meta1
            else meta1_value
        )
        next_meta2 = "none" if meta2_value == SELECTION_GROUP else meta2_value
        return base_meta1, next_meta1, base_meta2, next_meta2

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
            Input(f"{prefix}-violin2-data-layer", "data"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
            Input(f"{prefix}-exploratory-tabs", "value"),
        ],
        [
            State(f"{prefix}-selected-cells-store", "data"),
            State(f"{prefix}-selection-group-store", "data"),
        ],
    )
    def update_violin2(
        gene_selection,
        meta1,
        meta2,
        mode,
        test_method,
        show_box2,
        data_layer,
        selected_palette_name,
        filtered_data,
        _selected_cells_hash,
        _selection_group_hash,
        active_tab,
        selected_cells,
        highlighted_cells,
    ):
        if active_tab != "split-violin-tab":
            raise PreventUpdate
        layer = _resolve_layer(data_layer)
        filtered_adata = resolve_plot_adata_from_filter(filtered_data)
        group_values = None
        if selected_cells:
            filtered_adata = selected_cell_view(filtered_adata, selected_cells)
        elif highlighted_cells:
            filtered_adata, selection_values = selection_group_context(
                filtered_adata,
                highlighted_cells,
            )
            group_values = {SELECTION_GROUP: selection_values}

        meta2 = _resolve_violin2_meta2(mode, meta2)

        if mode in ["mode2", "mode3", "mode4"] and meta2 is None:
            raise PreventUpdate

        # Resolve the categorical palette the same way every other plot does -- the
        # selected discrete colormap, falling back to the dataset color_config -- so
        # the violins share the app's default colors instead of a private palette.
        palette = _resolve_violin2_palette(
            filtered_adata,
            meta1,
            meta2,
            selected_palette_name,
            color_config,
            group_values,
        )

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
            group_values=group_values,
        )
        return fig


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
    resolve_plot_adata_from_filter=None,
):
    """Compatibility wrapper registering both violin control models."""
    if resolve_plot_adata_from_filter is None:
        def resolve_plot_adata_from_filter(_filtered):
            return adata

    register_marker_violin_callbacks(
        app,
        adata,
        prefix,
        filter_data=filter_data,
        plot_violin1=plot_violin1,
    )
    register_comparative_violin_callbacks(
        app,
        adata,
        prefix,
        plot_violin2_new=plot_violin2_new,
        var_names=var_names,
        var_names_lower=var_names_lower,
        color_config=color_config,
        resolve_plot_adata_from_filter=resolve_plot_adata_from_filter,
    )
