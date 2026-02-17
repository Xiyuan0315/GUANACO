import dash
import plotly.graph_objects as go
from dash import Input, Output, State
from dash.exceptions import PreventUpdate


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
):
    @app.callback(
        Output(f"{prefix}-violin2-group-selection", "options"),
        Output(f"{prefix}-violin2-group-selection", "value"),
        [Input(f"{prefix}-meta1-selection", "value"), Input(f"{prefix}-selected-cells-store", "data")],
    )
    def update_group_labels(selected_column, selected_cells):
        if selected_cells:
            filtered_adata = adata[selected_cells]
            unique_labels = sorted(filtered_adata.obs[selected_column].dropna().unique(), key=str)
        else:
            unique_labels = sorted(adata.obs[selected_column].dropna().unique(), key=str)
        options = [{"label": str(label), "value": str(label)} for label in unique_labels]
        values = [str(label) for label in unique_labels]
        return options, values

    @app.callback(
        Output(f"{prefix}-violin-plot-cache-store", "data"),
        [
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-violin-log-or-zscore", "value"),
            Input(f"{prefix}-show-box1", "value"),
            Input(f"{prefix}-show-scatter1", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
        ],
        [State(f"{prefix}-violin-plot-cache-store", "data")],
    )
    def update_violin_cache(
        selected_genes,
        selected_annotation,
        selected_labels,
        transformation,
        show_box_plot,
        show_scatter1,
        discrete_color_map,
        selected_cells,
        current_cache,
    ):
        cache_key = f"{selected_genes}_{selected_annotation}_{selected_labels}_{transformation}_{show_box_plot}_{show_scatter1}_{discrete_color_map}_{selected_cells}"

        if current_cache is None:
            current_cache = {}

        if "current_key" in current_cache and current_cache["current_key"] == cache_key:
            return current_cache

        color_map = None
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            unique_labels = sorted(adata.obs[selected_annotation].unique())
            color_map = {label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels)}

        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)

        fig = plot_violin1(
            filtered_adata,
            selected_genes,
            selected_annotation,
            labels=selected_labels,
            transformation=transformation,
            show_box="show" in show_box_plot if show_box_plot else False,
            show_points="show" in show_scatter1 if show_scatter1 else False,
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
        Output(f"{prefix}-violin-plot1", "figure"),
        [Input(f"{prefix}-violin-plot-cache-store", "data"), Input(f"{prefix}-single-cell-tabs", "value")],
        [State(f"{prefix}-violin-plot1", "figure")],
    )
    def display_violin1(cache_data, active_tab, current_figure):
        if active_tab != "violin-tab":
            return current_figure if current_figure else go.Figure()
        if cache_data and "current_key" in cache_data:
            current_key = cache_data["current_key"]
            if current_key in cache_data:
                return go.Figure(cache_data[current_key])
        return current_figure if current_figure else go.Figure()

    @app.callback(Output(f"{prefix}-mode-explanation", "children"), Input(f"{prefix}-mode-selection", "value"))
    def update_mode_explanation(mode):
        explanations = {
            "mode1": "Compare expression across groups in obs1 only. Obs2 will be ignored.",
            "mode2": "Create facets by obs1, compare obs2 groups within each facet.",
            "mode3": "Linear model treating obs2 as a confounder: expression ~ obs1 + obs2",
            "mode4": "Mixed model treating obs2 as random effect: expression ~ meta1 + (1|obs2)",
        }
        return explanations.get(mode, "")

    @app.callback(
        Output(f"{prefix}-test-method-selection", "options"),
        Output(f"{prefix}-test-method-selection", "value"),
        [Input(f"{prefix}-mode-selection", "value"), Input(f"{prefix}-meta1-selection", "value"), Input(f"{prefix}-meta2-selection", "value")],
    )
    def update_test_methods(mode, meta1, meta2):
        base_options = [{"label": "None", "value": "none"}]

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
        Output(f"{prefix}-violin-plot2", "figure"),
        [
            Input(f"{prefix}-violin2-gene-selection", "value"),
            Input(f"{prefix}-meta1-selection", "value"),
            Input(f"{prefix}-meta2-selection", "value"),
            Input(f"{prefix}-mode-selection", "value"),
            Input(f"{prefix}-test-method-selection", "value"),
            Input(f"{prefix}-show-box2", "value"),
            Input(f"{prefix}-show-scatter2", "value"),
            Input(f"{prefix}-violin2-log-or-zscore", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
        ],
    )
    def update_violin2(gene_selection, meta1, meta2, mode, test_method, show_box2, show_points, transformation, selected_palette_name, selected_cells):
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

        palette = None
        if selected_palette_name:
            palette = palette_json["color_palettes"].get(selected_palette_name)

        return plot_violin2_new(
            filtered_adata,
            key=gene_selection,
            meta1=meta1,
            meta2=meta2,
            mode=mode,
            transformation=transformation,
            show_box="show" in show_box2 if show_box2 else False,
            show_points="show" in show_points if show_points else False,
            test_method=test_method,
            labels=None,
            color_map=None,
            palette=palette,
        )
