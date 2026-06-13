from dash import Input, Output, State, no_update
import plotly.graph_objects as go

from guanaco.utils.colors import resolve_discrete_palette



def register_pseudotime_callbacks(
    app,
    adata,
    prefix,
    *,
    filter_data,
    plot_genes_in_pseudotime,
    palette_json,
    color_config,
    make_cache_key,
    hash_list_signature,
    cached_figure_get,
    cached_figure_set,
):
    @app.callback(
        [
            Output(f"{prefix}-pseudotime-plot", "figure"),
            Output(f"{prefix}-pseudotime-rendered-key", "data"),
        ],
        [
            Input(f"{prefix}-single-cell-tabs", "value"),
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-pseudotime-min-expr-slider", "value"),
            Input(f"{prefix}-data-layer", "value"),
            Input(f"{prefix}-pseudotime-key-dropdown", "value"),
            Input(f"{prefix}-marker-size-slider", "value"),
            Input(f"{prefix}-opacity-slider", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-hash", "data"),
        ],
        [
            State(f"{prefix}-pseudotime-plot", "figure"),
            State(f"{prefix}-pseudotime-rendered-key", "data"),
            State(f"{prefix}-selected-cells-store", "data"),
        ],
    )
    def update_pseudotime_plot(
        selected_tab,
        selected_genes,
        selected_annotation,
        selected_labels,
        min_expr,
        data_layer,
        pseudotime_key,
        marker_size,
        opacity,
        discrete_color_map,
        cells_hash,
        current_figure,
        rendered_key,
        selected_cells,
    ):
        if selected_tab != "pseudotime-tab":
            # Not the active tab: leave whatever is there untouched.
            return no_update, no_update

        if not selected_genes:
            fig = go.Figure()
            fig.add_annotation(
                text="<b>Please select genes to plot</b><br><br>Use the 'Select Variables' dropdown on the left to choose genes",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
                font=dict(size=16),
                align="center",
            )
            fig.update_layout(
                plot_bgcolor="#f8f9fa",
                paper_bgcolor="white",
                height=400,
                margin=dict(t=50, b=50, l=50, r=50),
            )
            return fig, None

        cache_key = make_cache_key(
            "pseudotime",
            adata,
            selected_genes=hash_list_signature(selected_genes),
            selected_annotation=selected_annotation,
            selected_labels=hash_list_signature(selected_labels),
            min_expr=min_expr,
            data_layer=data_layer,
            pseudotime_key=pseudotime_key,
            marker_size=marker_size,
            opacity=opacity,
            discrete_color_map=discrete_color_map,
            selected_cells=cells_hash,
        )
        # Already showing the figure for these exact parameters: do nothing,
        # so a plain tab switch neither recomputes nor redraws.
        if rendered_key == cache_key and current_figure:
            return no_update, no_update

        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return cached_fig, cache_key

        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)

        all_categories = sorted(adata.obs[selected_annotation].unique())
        if discrete_color_map:
            discrete_palette = resolve_discrete_palette(discrete_color_map, len(all_categories))
            color_map = {cat: discrete_palette[i % len(discrete_palette)] for i, cat in enumerate(all_categories)}
        else:
            color_map = {cat: color_config[i % len(color_config)] for i, cat in enumerate(all_categories)}

        if not pseudotime_key:
            numeric_cols = []
            for col in filtered_adata.obs.columns:
                if filtered_adata.obs[col].dtype in ["float32", "float64", "int32", "int64"]:
                    if filtered_adata.obs[col].nunique() > 50:
                        numeric_cols.append(col)

            pseudotime_cols = [col for col in numeric_cols if "pseudotime" in col.lower() or "dpt" in col.lower()]

            if pseudotime_cols:
                pseudotime_key = pseudotime_cols[0]
            elif numeric_cols:
                pseudotime_key = numeric_cols[0]
            else:
                fig = go.Figure()
                fig.add_annotation(
                    text="<b>No continuous variable found</b><br><br>Please ensure your data has a numeric obs column",
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    font=dict(size=16),
                    align="center",
                )
                fig.update_layout(
                    plot_bgcolor="#f8f9fa",
                    paper_bgcolor="white",
                    height=400,
                    margin=dict(t=50, b=50, l=50, r=50),
                )
                return fig, None

        fig = plot_genes_in_pseudotime(
            filtered_adata,
            selected_genes,
            pseudotime_key=pseudotime_key,
            groupby=selected_annotation,
            min_expr=min_expr,
            transformation="none",
            layer=data_layer if data_layer and data_layer != "X" else None,
            color_map=color_map,
            marker_size=marker_size,
            opacity=opacity,
        )
        cached_figure_set(cache_key, fig)
        return fig, cache_key

    @app.callback(
        [Output(f"{prefix}-pseudotime-key-dropdown", "options"), Output(f"{prefix}-pseudotime-key-dropdown", "value")],
        Input(f"{prefix}-single-cell-tabs", "value"),
    )
    def update_pseudotime_key_options(active_tab):
        if active_tab == "pseudotime-tab":
            numeric_cols = []
            for col in adata.obs.columns:
                if adata.obs[col].dtype in ["float32", "float64", "int32", "int64"]:
                    if adata.obs[col].nunique() > 50:
                        numeric_cols.append(col)
            pseudotime_cols = [col for col in numeric_cols if "pseudotime" in col.lower() or "dpt" in col.lower()]
            other_cols = [col for col in numeric_cols if col not in pseudotime_cols]
            all_cols = pseudotime_cols + other_cols
            options = [{"label": col, "value": col} for col in all_cols]
            if not options:
                options = [{"label": "No continuous variables found", "value": None}]
                return options, None
            default_value = all_cols[0] if all_cols else None
            return options, default_value
        return [], None
