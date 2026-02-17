import numpy as np
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html, callback_context, exceptions
from dash.exceptions import PreventUpdate


def register_scatter_callbacks(
    app,
    adata,
    prefix,
    *,
    embedding_render_backend,
    initialize_scatter_components,
    apply_relayout,
    is_continuous_annotation,
    resolve_plot_adata_from_filter,
    search_combined,
    obs_columns,
    obs_columns_lower,
    var_names,
    var_names_lower,
    palette_json,
    color_config,
    plot_embedding,
    plot_coexpression_embedding,
    make_cache_key,
    filtered_data_signature,
    cached_figure_get,
    cached_figure_set,
):
    @app.callback(
        Output(f"{prefix}-controls-container", "style"),
        Output(f"{prefix}-toggle-button", "children"),
        Input(f"{prefix}-toggle-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def toggle_controls(n_clicks):
        if n_clicks % 2 == 1:
            return {"display": "block"}, "Hide controls"
        return {"display": "none"}, "More controls"

    @app.callback(
        [
            Output(f"{prefix}-coordinates-dropdowns", "children"),
            Output(f"{prefix}-x-axis", "value"),
            Output(f"{prefix}-y-axis", "value"),
        ],
        Input(f"{prefix}-clustering-dropdown", "value"),
    )
    def update_coordinates_dropdowns(selected_clustering):
        _, embedding_columns, _, _ = initialize_scatter_components(adata)
        selected_columns = embedding_columns[selected_clustering]
        options = [{"label": col, "value": col} for col in selected_columns]
        x_value = selected_columns[0]
        y_value = selected_columns[1] if len(selected_columns) > 1 else selected_columns[0]
        return (
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("X-axis:"),
                            dcc.Dropdown(
                                id=f"{prefix}-x-axis",
                                options=options,
                                value=x_value,
                                clearable=False,
                                style={"fontSize": "14px"},
                            ),
                        ],
                        style={"flex": "1", "paddingRight": "10px"},
                    ),
                    html.Div(
                        [
                            html.Label("Y-axis:"),
                            dcc.Dropdown(
                                id=f"{prefix}-y-axis",
                                options=options,
                                value=y_value,
                                clearable=False,
                                style={"fontSize": "14px"},
                            ),
                        ],
                        style={"flex": "1", "paddingLeft": "10px"},
                    ),
                ],
                style={"display": "flex", "marginBottom": "15px"},
            ),
            x_value,
            y_value,
        )

    @app.callback(
        [
            Output(f"{prefix}-spatial-imgkey-container", "style"),
            Output(f"{prefix}-spatial-imgkey-dropdown", "options"),
            Output(f"{prefix}-spatial-imgkey-dropdown", "value"),
        ],
        Input(f"{prefix}-clustering-dropdown", "value"),
        State(f"{prefix}-spatial-imgkey-dropdown", "value"),
    )
    def update_spatial_imgkey_dropdown(selected_clustering, current_value):
        if selected_clustering != "spatial" or "spatial" not in adata.uns:
            return {"display": "none", "marginBottom": "15px"}, [], None

        spatial = adata.uns.get("spatial", {})
        image_keys = set()
        for lib_data in spatial.values():
            image_keys.update(lib_data.get("images", {}).keys())

        image_keys = sorted(image_keys)
        options = [{"label": key, "value": key} for key in image_keys]

        if not image_keys:
            return {"display": "none", "marginBottom": "15px"}, [], None

        if current_value in image_keys:
            selected_value = current_value
        elif "hires" in image_keys:
            selected_value = "hires"
        elif "highres" in image_keys:
            selected_value = "highres"
        else:
            selected_value = image_keys[0]

        return {"display": "block", "marginBottom": "15px"}, options, selected_value

    @app.callback(
        Output(f"{prefix}-annotation-dropdown", "options"),
        Input(f"{prefix}-annotation-dropdown", "search_value"),
    )
    def update_annotation_dropdown(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        all_matches = search_combined(
            obs_columns, obs_columns_lower, var_names, var_names_lower, search_value, limit=10
        )
        return [{"label": item, "value": item} for item in all_matches]

    @app.callback(
        Output(f"{prefix}-scatter-gene-selection", "options"),
        Input(f"{prefix}-scatter-gene-selection", "search_value"),
    )
    def update_scatter_gene_selection(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        all_matches = search_combined(
            obs_columns, obs_columns_lower, var_names, var_names_lower, search_value, limit=10
        )
        return [{"label": item, "value": item} for item in all_matches]

    @app.callback(
        Output(f"{prefix}-scatter-gene2-selection", "options"),
        Input(f"{prefix}-scatter-gene2-selection", "search_value"),
    )
    def update_scatter_gene2_selection(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        q = search_value.lower()
        matching_genes = [gene for gene, gene_l in zip(var_names, var_names_lower) if q in gene_l]
        return [{"label": gene, "value": gene} for gene in matching_genes[:20]]

    @app.callback(
        [Output(f"{prefix}-gene2-container", "style"), Output(f"{prefix}-threshold-container", "style")],
        Input(f"{prefix}-coexpression-toggle", "value"),
    )
    def toggle_coexpression_controls(mode):
        if mode == "coexpression":
            return {"display": "block"}, {"display": "block"}
        return {"display": "none"}, {"display": "none"}

    @app.callback(
        Output(f"{prefix}-annotation-scatter", "figure"),
        [
            Input(f"{prefix}-clustering-dropdown", "value"),
            Input(f"{prefix}-x-axis", "value"),
            Input(f"{prefix}-y-axis", "value"),
            Input(f"{prefix}-annotation-dropdown", "value"),
            Input(f"{prefix}-marker-size-slider", "value"),
            Input(f"{prefix}-opacity-slider", "value"),
            Input(f"{prefix}-scatter-legend-toggle", "value"),
            Input(f"{prefix}-axis-toggle", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-scatter-log-or-zscore", "value"),
            Input(f"{prefix}-plot-order", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-spatial-imgkey-dropdown", "value"),
            Input(f"{prefix}-gene-scatter", "relayoutData"),
        ],
        State(f"{prefix}-annotation-scatter", "relayoutData"),
    )
    def update_annotation_scatter(
        clustering_method,
        x_axis,
        y_axis,
        annotation,
        marker_size,
        opacity,
        legend_show,
        axis_show,
        discrete_color_map,
        transformation,
        order,
        continuous_color_map,
        filtered_data,
        spatial_img_key,
        gene_relayout,
        annotation_relayout,
    ):
        if not annotation:
            raise exceptions.PreventUpdate

        triggered_prop = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        # Keep left/right plots synced when the right plot was zoomed/reset.
        cross_relayout = (
            gene_relayout
            if triggered_prop == f"{prefix}-gene-scatter.relayoutData"
            else None
        )
        # Preserve current view on style-only updates; reset on geometry/data changes.
        if triggered_prop in {
            f"{prefix}-clustering-dropdown.value",
            f"{prefix}-x-axis.value",
            f"{prefix}-y-axis.value",
            f"{prefix}-global-filtered-data.data",
            f"{prefix}-spatial-imgkey-dropdown.value",
        }:
            self_relayout = None
        else:
            self_relayout = annotation_relayout
        effective_relayout = cross_relayout if cross_relayout is not None else self_relayout

        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        filtered_sig = filtered_data_signature(filtered_data)

        is_continuous = annotation in adata.var_names or is_continuous_annotation(adata, annotation)
        render_backend = embedding_render_backend
        discrete_palette = color_config if discrete_color_map is None else palette_json["color_palettes"][discrete_color_map]
        use_transformation = transformation if annotation in adata.var_names else None
        cache_key = make_cache_key(
            "annotation_scatter",
            adata,
            clustering_method=clustering_method,
            x_axis=x_axis,
            y_axis=y_axis,
            annotation=annotation,
            marker_size=marker_size,
            opacity=opacity,
            render_backend=render_backend,
            legend_show=legend_show,
            axis_show=axis_show,
            discrete_color_map=discrete_color_map,
            transformation=use_transformation if is_continuous else None,
            order=order if is_continuous else None,
            continuous_color_map=continuous_color_map or "Viridis",
            filtered_data=filtered_sig,
            spatial_img_key=spatial_img_key,
            mode="continuous" if is_continuous else "categorical",
        )
        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return apply_relayout(cached_fig, effective_relayout)

        fig = plot_embedding(
            adata=plot_adata,
            adata_full=adata,
            embedding_key=clustering_method,
            color=annotation,
            x_axis=x_axis,
            y_axis=y_axis,
            mode="continuous" if is_continuous else "categorical",
            transformation=use_transformation if is_continuous else None,
            order=order if is_continuous else None,
            continuous_color_map=continuous_color_map or "Viridis",
            discrete_color_map=discrete_palette,
            marker_size=marker_size,
            opacity=opacity,
            render_backend=render_backend,
            legend_show=legend_show,
            axis_show=axis_show,
            img_key=spatial_img_key,
        )
        cached_figure_set(cache_key, fig)
        return apply_relayout(fig, effective_relayout)

    @app.callback(
        Output(f"{prefix}-gene-scatter", "figure"),
        [
            Input(f"{prefix}-scatter-gene-selection", "value"),
            Input(f"{prefix}-clustering-dropdown", "value"),
            Input(f"{prefix}-x-axis", "value"),
            Input(f"{prefix}-y-axis", "value"),
            Input(f"{prefix}-scatter-log-or-zscore", "value"),
            Input(f"{prefix}-plot-order", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-marker-size-slider", "value"),
            Input(f"{prefix}-opacity-slider", "value"),
            Input(f"{prefix}-annotation-scatter", "relayoutData"),
            Input(f"{prefix}-axis-toggle", "value"),
            Input(f"{prefix}-coexpression-toggle", "value"),
            Input(f"{prefix}-scatter-gene2-selection", "value"),
            Input(f"{prefix}-gene1-threshold-slider", "value"),
            Input(f"{prefix}-gene2-threshold-slider", "value"),
            Input(f"{prefix}-scatter-legend-toggle", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-spatial-imgkey-dropdown", "value"),
        ],
        State(f"{prefix}-gene-scatter", "relayoutData"),
    )
    def update_gene_scatter(
        gene_name,
        clustering,
        x_axis,
        y_axis,
        transformation,
        order,
        color_map,
        marker_size,
        opacity,
        annotation_relayout,
        axis_show,
        coexpression_mode,
        gene2_name,
        threshold1,
        threshold2,
        legend_show,
        discrete_color_map,
        filtered_data,
        spatial_img_key,
        gene_relayout,
    ):
        if not gene_name:
            raise exceptions.PreventUpdate

        triggered_prop = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        # Keep left/right plots synced when the left plot was zoomed/reset.
        cross_relayout = (
            annotation_relayout
            if triggered_prop == f"{prefix}-annotation-scatter.relayoutData"
            else None
        )
        # Preserve current view on style-only updates; reset on geometry/data changes.
        if triggered_prop in {
            f"{prefix}-clustering-dropdown.value",
            f"{prefix}-x-axis.value",
            f"{prefix}-y-axis.value",
            f"{prefix}-global-filtered-data.data",
            f"{prefix}-spatial-imgkey-dropdown.value",
        }:
            self_relayout = None
        else:
            self_relayout = gene_relayout
        effective_relayout = cross_relayout if cross_relayout is not None else self_relayout

        render_backend = embedding_render_backend
        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        cache_key = make_cache_key(
            "gene_scatter",
            adata,
            gene_name=gene_name,
            clustering=clustering,
            x_axis=x_axis,
            y_axis=y_axis,
            transformation=transformation,
            order=order,
            color_map=color_map or "Viridis",
            marker_size=marker_size,
            opacity=opacity,
            render_backend=render_backend,
            axis_show=axis_show,
            coexpression_mode=coexpression_mode,
            gene2_name=gene2_name,
            threshold1=threshold1,
            threshold2=threshold2,
            legend_show=legend_show,
            discrete_color_map=discrete_color_map,
            filtered_data=filtered_data_signature(filtered_data),
            spatial_img_key=spatial_img_key,
        )
        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return apply_relayout(cached_fig, effective_relayout)

        if gene_name in adata.var_names:
            if coexpression_mode == "coexpression" and gene2_name:
                fig = plot_coexpression_embedding(
                    adata=plot_adata,
                    embedding_key=clustering,
                    gene1=gene_name,
                    gene2=gene2_name,
                    x_axis=x_axis,
                    y_axis=y_axis,
                    threshold1=threshold1,
                    threshold2=threshold2,
                    transformation=transformation,
                    color_map=None,
                    marker_size=marker_size,
                    opacity=opacity,
                    legend_show=legend_show,
                    axis_show=axis_show,
                )
            else:
                fig = plot_embedding(
                    adata=plot_adata,
                    embedding_key=clustering,
                    color=gene_name,
                    x_axis=x_axis,
                    y_axis=y_axis,
                    mode="continuous",
                    transformation=transformation,
                    order=order,
                    continuous_color_map=color_map or "Viridis",
                    marker_size=marker_size,
                    opacity=opacity,
                    render_backend=render_backend,
                    annotation=None,
                    axis_show=axis_show,
                    img_key=spatial_img_key,
                )
        elif is_continuous_annotation(adata, gene_name):
            fig = plot_embedding(
                adata=plot_adata,
                embedding_key=clustering,
                color=gene_name,
                x_axis=x_axis,
                y_axis=y_axis,
                mode="continuous",
                transformation=None,
                order=order,
                continuous_color_map=color_map or "Viridis",
                marker_size=marker_size,
                opacity=opacity,
                render_backend=render_backend,
                axis_show=axis_show,
                img_key=spatial_img_key,
            )
        else:
            if discrete_color_map is None:
                discrete_color_map_value = color_config
            else:
                discrete_color_map_value = palette_json["color_palettes"][discrete_color_map]
            fig = plot_embedding(
                adata=plot_adata,
                adata_full=adata,
                embedding_key=clustering,
                color=gene_name,
                x_axis=x_axis,
                y_axis=y_axis,
                mode="categorical",
                discrete_color_map=discrete_color_map_value,
                marker_size=marker_size,
                opacity=opacity,
                render_backend=render_backend,
                legend_show=legend_show,
                axis_show=axis_show,
                img_key=spatial_img_key,
            )
        cached_figure_set(cache_key, fig)
        return apply_relayout(fig, effective_relayout)

    @app.callback(
        [
            Output(f"{prefix}-gene1-threshold-slider", "min"),
            Output(f"{prefix}-gene1-threshold-slider", "max"),
            Output(f"{prefix}-gene1-threshold-slider", "value"),
            Output(f"{prefix}-gene2-threshold-slider", "min"),
            Output(f"{prefix}-gene2-threshold-slider", "max"),
            Output(f"{prefix}-gene2-threshold-slider", "value"),
        ],
        [
            Input(f"{prefix}-scatter-gene-selection", "value"),
            Input(f"{prefix}-scatter-gene2-selection", "value"),
            Input(f"{prefix}-coexpression-toggle", "value"),
            Input(f"{prefix}-scatter-log-or-zscore", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
        ],
    )
    def update_threshold_ranges(gene1, gene2, coexpression_mode, transformation, filtered_data):
        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        default_min, default_max, default_value = 0, 1, 0.5
        from guanaco.utils.gene_extraction_utils import extract_gene_expression, apply_transformation

        if gene1 and gene1 in plot_adata.var_names:
            gene1_expr = extract_gene_expression(plot_adata, gene1)
            if transformation:
                gene1_expr = apply_transformation(gene1_expr, transformation, copy=True)
            if gene1_expr.max() > gene1_expr.min():
                gene1_min = float(gene1_expr.min())
                gene1_max = float(gene1_expr.max())
                gene1_value = (gene1_min + gene1_max) / 2
            else:
                gene1_min, gene1_max, gene1_value = default_min, default_max, default_value
        else:
            gene1_min, gene1_max, gene1_value = default_min, default_max, default_value

        if coexpression_mode == "coexpression" and gene2 and gene2 in plot_adata.var_names:
            gene2_expr = extract_gene_expression(plot_adata, gene2)
            if transformation:
                gene2_expr = apply_transformation(gene2_expr, transformation, copy=True)
            if gene2_expr.max() > gene2_expr.min():
                gene2_min = float(gene2_expr.min())
                gene2_max = float(gene2_expr.max())
                gene2_value = (gene2_min + gene2_max) / 2
            else:
                gene2_min, gene2_max, gene2_value = default_min, default_max, default_value
        else:
            gene2_min, gene2_max, gene2_value = default_min, default_max, default_value

        return gene1_min, gene1_max, gene1_value, gene2_min, gene2_max, gene2_value

    @app.callback(
        Output(f"{prefix}-selected-cells-store", "data"),
        Output(f"{prefix}-selection-status", "children"),
        [Input(f"{prefix}-update-plots-button", "n_clicks")],
        [
            State(f"{prefix}-annotation-scatter", "selectedData"),
            State(f"{prefix}-annotation-dropdown", "value"),
            State(f"{prefix}-global-filtered-data", "data"),
        ],
        prevent_initial_call=True,
    )
    def store_selected_cells(n_clicks, selected_data, current_annotation, filtered_data):
        if n_clicks == 0:
            return None, ""

        plot_adata = resolve_plot_adata_from_filter(filtered_data)

        if not selected_data or not selected_data.get("points"):
            all_indices = plot_adata.obs.index.tolist()
            n_cells = len(all_indices)
            status_msg = dbc.Alert(
                f"✓ All {n_cells} cells from scatter plot selected. Other plots updated.",
                color="info",
                dismissable=True,
                duration=4000,
            )
            return all_indices, status_msg

        selected_points = selected_data["points"]
        selected_indices = []

        if current_annotation in adata.var_names or is_continuous_annotation(plot_adata, current_annotation):
            for point in selected_points:
                if "customdata" in point:
                    customdata = point["customdata"]
                    if isinstance(customdata, (list, tuple)) and len(customdata) > 1:
                        cell_idx = int(customdata[1])
                    else:
                        cell_idx = int(customdata)
                    selected_indices.append(plot_adata.obs.index[cell_idx])
                else:
                    point_number = point.get("pointNumber", 0)
                    selected_indices.append(plot_adata.obs.index[point_number])
        else:
            category_to_indices = {
                cat: plot_adata.obs.index[idx].tolist() for cat, idx in plot_adata.obs.groupby(current_annotation).indices.items()
            }
            unique_categories = sorted(category_to_indices.keys())
            for point in selected_points:
                curve_number = point.get("curveNumber", 0)
                point_number = point.get("pointNumber", 0)
                if curve_number == 0:
                    continue

                if "customdata" in point:
                    customdata = point["customdata"]
                    if isinstance(customdata, (list, tuple)):
                        category = customdata[0]
                    else:
                        category = customdata
                    category_indices = category_to_indices.get(category, [])
                    if point_number < len(category_indices):
                        selected_indices.append(category_indices[point_number])
                else:
                    category_index = curve_number - 1
                    if category_index >= 0 and category_index < len(unique_categories):
                        selected_category = unique_categories[category_index]
                        category_indices = category_to_indices.get(selected_category, [])
                        if point_number < len(category_indices):
                            selected_indices.append(category_indices[point_number])

        if selected_indices:
            n_selected = len(selected_indices)
            status_msg = dbc.Alert(
                f"✓ {n_selected} cells selected from {current_annotation}. Other plots updated.",
                color="success",
                dismissable=True,
                duration=4000,
            )
            return selected_indices, status_msg
        return None, ""

    @app.callback(Output(f"{prefix}-download-menu", "disabled"), [Input(f"{prefix}-selected-cells-store", "data")])
    def toggle_download_menu(selected_cells):
        return not bool(selected_cells)

    @app.callback(
        Output(f"{prefix}-download-cells-data", "data"),
        [Input(f"{prefix}-download-cellids", "n_clicks")],
        [State(f"{prefix}-selected-cells-store", "data")],
        prevent_initial_call=True,
    )
    def download_selected_cells(n_clicks_txt, selected_cells):
        ctx = callback_context
        if not ctx.triggered or not selected_cells:
            raise PreventUpdate
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if f"{prefix}-download-cellids" in button_id:
            content = "\n".join(selected_cells)
            return dict(content=content, filename="selected_cells.txt")
        raise PreventUpdate
