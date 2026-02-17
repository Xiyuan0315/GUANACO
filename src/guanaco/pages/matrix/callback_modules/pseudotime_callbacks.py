from dash import Input, Output, State
import plotly.graph_objects as go


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
        Output(f"{prefix}-pseudotime-plot", "figure"),
        [
            Input(f"{prefix}-single-cell-tabs", "value"),
            Input(f"{prefix}-single-cell-genes-selection", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-pseudotime-min-expr-slider", "value"),
            Input(f"{prefix}-pseudotime-transformation", "value"),
            Input(f"{prefix}-pseudotime-key-dropdown", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
        ],
        [State(f"{prefix}-pseudotime-plot", "figure")],
    )
    def update_pseudotime_plot(
        selected_tab,
        selected_genes,
        selected_annotation,
        selected_labels,
        min_expr,
        transformation,
        pseudotime_key,
        discrete_color_map,
        selected_cells,
        current_figure,
    ):
        if selected_tab != "pseudotime-tab":
            return current_figure if current_figure else go.Figure()

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
            return fig

        cache_key = make_cache_key(
            "pseudotime",
            adata,
            selected_genes=hash_list_signature(selected_genes),
            selected_annotation=selected_annotation,
            selected_labels=hash_list_signature(selected_labels),
            min_expr=min_expr,
            transformation=transformation,
            pseudotime_key=pseudotime_key,
            discrete_color_map=discrete_color_map,
            selected_cells=hash_list_signature(selected_cells),
        )
        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return cached_fig

        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)

        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            all_categories = sorted(filtered_adata.obs[selected_annotation].unique())
            color_map = {cat: discrete_palette[i % len(discrete_palette)] for i, cat in enumerate(all_categories)}
        else:
            all_categories = sorted(filtered_adata.obs[selected_annotation].unique())
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
                    text="<b>No pseudotime column found</b><br><br>Please ensure your data has a numeric column with pseudotime values",
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
                return fig

        fig = plot_genes_in_pseudotime(
            filtered_adata,
            selected_genes,
            pseudotime_key=pseudotime_key,
            groupby=selected_annotation,
            min_expr=min_expr,
            transformation=transformation,
            color_map=color_map,
        )
        cached_figure_set(cache_key, fig)
        return fig

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
                options = [{"label": "No pseudotime columns found", "value": None}]
                return options, None
            default_value = all_cols[0] if all_cols else None
            return options, default_value
        return [], None
