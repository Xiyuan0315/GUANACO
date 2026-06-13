import dash
from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import sorted_categories



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
        Output(f"{prefix}-violin2-group-selection", "options"),
        Output(f"{prefix}-violin2-group-selection", "value"),
        [Input(f"{prefix}-meta1-selection", "value"), Input(f"{prefix}-selected-cells-store", "data")],
    )
    def update_group_labels(selected_column, selected_cells):
        src = adata[selected_cells] if selected_cells else adata
        unique_labels = sorted_categories(src, selected_column)
        options = [{"label": str(label), "value": str(label)} for label in unique_labels]
        values = [str(label) for label in unique_labels]
        return options, values

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
        layer = data_layer if data_layer and data_layer != "X" else None
        cache_key = f"{selected_genes}_{selected_annotation}_{selected_labels}_{data_layer}_{show_box_plot}_{discrete_color_map}_{cells_hash}"

        if current_cache is None:
            current_cache = {}

        if "current_key" in current_cache and current_cache["current_key"] == cache_key:
            return current_cache

        color_map = None
        if discrete_color_map:
            unique_labels = sorted_categories(adata, selected_annotation)
            discrete_palette = resolve_discrete_palette(discrete_color_map, len(unique_labels))
            color_map = {label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels)}

        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)

        fig = plot_violin1(
            filtered_adata,
            selected_genes,
            selected_annotation,
            labels=selected_labels,
            layer=layer,
            show_box="show" in show_box_plot if show_box_plot else False,
            groupby_label_color_map=color_map,
            adata_obs=adata.obs,
            data_already_filtered=True,
        )

        num_genes = len(selected_genes) if selected_genes else 0
        num_categories = len(selected_labels) if selected_labels else 0
        fig.update_layout(
            height=min(1000, max(300, 80 * num_genes)),
            width=min(500, max(200, 110 * num_categories)),
            margin=dict(l=130, r=10, t=30, b=30),
        )

        if len(current_cache) > 10:
            keys_to_remove = list(current_cache.keys())[:-9]
            for key in keys_to_remove:
                if key != "current_key":
                    current_cache.pop(key, None)

        current_cache[cache_key] = fig.to_dict()
        current_cache["current_key"] = cache_key
        return current_cache

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

        current_key = cache_data.get("current_key") if cache_data else None
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
        explanations = {
            "mode1": "Compare expression across groups in obs1 only. Obs2 will be ignored.",
            "mode2": "Within each obs1 group on the x-axis, split/group the violin by obs2 and compare.",
            "mode3": "Linear model treating obs2 as a confounder: expression ~ obs1 + obs2",
            "mode4": "Mixed model treating obs2 as a random effect (e.g. donor, batch, replicate) to account for non-independent samples.",
        }
        return explanations.get(mode, "")

    @app.callback(
        Output(f"{prefix}-test-method-selection", "options"),
        Output(f"{prefix}-test-method-selection", "value"),
        [Input(f"{prefix}-mode-selection", "value"), Input(f"{prefix}-meta1-selection", "value"), Input(f"{prefix}-meta2-selection", "value")],
    )
    def update_test_methods(mode, meta1, meta2):
        # Auto is always available and is the default; None disables the test.
        base_options = [
            {"label": "Auto", "value": "auto"},
            {"label": "None", "value": "none"},
        ]

        if mode == "mode1":
            n_levels = len(adata.obs[meta1].unique()) if meta1 else 0
            if n_levels == 2:
                options = base_options + [{"label": "Mann-Whitney U", "value": "mwu-test"}, {"label": "T-test", "value": "ttest"}]
            else:
                options = base_options + [{"label": "Kruskal-Wallis", "value": "kw-test"}, {"label": "ANOVA", "value": "anova"}]
        elif mode == "mode2":
            if meta2 and meta2 != "none":
                n_levels = len(adata.obs[meta2].unique())
                if n_levels == 2:
                    options = base_options + [{"label": "Mann-Whitney U", "value": "mwu-test"}, {"label": "T-test", "value": "ttest"}]
                else:
                    options = base_options + [{"label": "Kruskal-Wallis", "value": "kw-test"}, {"label": "ANOVA", "value": "anova"}]
            else:
                options = base_options
        elif mode == "mode3":
            options = base_options + [
                {"label": "Linear Model", "value": "linear-model"},
                {"label": "Linear Model with Interaction", "value": "linear-model-interaction"},
            ]
        elif mode == "mode4":
            options = base_options + [{"label": "Mixed Model", "value": "mixed-model"}]
        else:
            options = base_options

        return options, "auto"

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
        return False, dash.no_update

    @app.callback(
        Output(f"{prefix}-split-violin-options-collapse", "is_open"),
        Output(f"{prefix}-split-violin-options-toggle", "children"),
        Input(f"{prefix}-split-violin-options-toggle", "n_clicks"),
        State(f"{prefix}-split-violin-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_split_violin_options(n_clicks, is_open):
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
        layer = data_layer if data_layer and data_layer != "X" else None
        if selected_cells:
            filtered_adata = filter_data(adata, None, None, selected_cells)
        else:
            filtered_adata = adata

        if mode == "mode1":
            meta2 = None
        elif meta2 == "none":
            meta2 = None

        if mode in ["mode2", "mode3", "mode4"] and meta2 is None:
            raise PreventUpdate

        # Resolve the categorical palette the same way every other plot does -- the
        # selected discrete colormap, falling back to the dataset color_config -- so
        # the violins share the app's default colors instead of a private palette.
        n_colors = 0
        if meta1:
            n_colors = max(n_colors, filtered_adata.obs[meta1].nunique())
        if meta2:
            n_colors = max(n_colors, filtered_adata.obs[meta2].nunique())
        palette = resolve_discrete_palette(selected_palette_name, n_colors, default=color_config)

        fig = plot_violin2_new(
            filtered_adata,
            key=gene_selection,
            meta1=meta1,
            meta2=meta2,
            mode=mode,
            layer=layer,
            show_box="show" in show_box2 if show_box2 else False,
            test_method=test_method,
            labels=None,
            color_map=None,
            palette=palette,
        )
        return fig
