import numpy as np
import hashlib
import json
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, callback_context, exceptions
from dash.exceptions import PreventUpdate

from guanaco.utils.colors import resolve_discrete_palette



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
    def _filtered_cell_indices(filtered_data):
        if (
            filtered_data
            and filtered_data.get("cell_indices") is not None
            and filtered_data.get("n_cells", adata.n_obs) < adata.n_obs
        ):
            return np.asarray(filtered_data["cell_indices"], dtype=np.int64)
        return None

    def _list_signature(values):
        if values is None:
            return None
        if not isinstance(values, (list, tuple)):
            return values
        if len(values) == 0:
            return {"len": 0, "hash": None}
        payload = json.dumps(list(values), sort_keys=False, default=str, separators=(",", ":"))
        return {"len": len(values), "hash": hashlib.md5(payload.encode("utf-8")).hexdigest()}

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

    def _coordinate_dropdown_children(id_prefix, selected_clustering):
        _, embedding_columns, _, _ = initialize_scatter_components(adata)
        selected_columns = embedding_columns[selected_clustering]
        options = [{"label": col, "value": col} for col in selected_columns]
        x_value = selected_columns[0]
        y_value = selected_columns[1] if len(selected_columns) > 1 else selected_columns[0]
        display_style = {"display": "flex", "marginBottom": "15px"} if selected_clustering == "X_pca" else {"display": "none"}
        return (
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("X-axis:"),
                            dcc.Dropdown(
                                id=f"{id_prefix}-x-axis",
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
                                id=f"{id_prefix}-y-axis",
                                options=options,
                                value=y_value,
                                clearable=False,
                                style={"fontSize": "14px"},
                            ),
                        ],
                        style={"flex": "1", "paddingLeft": "10px"},
                    ),
                ],
                style=display_style,
            ),
            x_value,
            y_value,
        )

    @app.callback(
        [
            Output(f"{prefix}-coordinates-dropdowns", "children"),
            Output(f"{prefix}-x-axis", "value"),
            Output(f"{prefix}-y-axis", "value"),
        ],
        Input(f"{prefix}-clustering-dropdown", "value"),
    )
    def update_coordinates_dropdowns(selected_clustering):
        return _coordinate_dropdown_children(prefix, selected_clustering)

    @app.callback(
        [
            Output(f"{prefix}-right-coordinates-dropdowns", "children"),
            Output(f"{prefix}-right-x-axis", "value"),
            Output(f"{prefix}-right-y-axis", "value"),
        ],
        Input(f"{prefix}-right-clustering-dropdown", "value"),
    )
    def update_right_coordinates_dropdowns(selected_clustering):
        return _coordinate_dropdown_children(f"{prefix}-right", selected_clustering)

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
        [
            State(f"{prefix}-annotation-scatter", "relayoutData"),
            State(f"{prefix}-right-clustering-dropdown", "value"),
            State(f"{prefix}-right-x-axis", "value"),
            State(f"{prefix}-right-y-axis", "value"),
        ],
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
        right_clustering,
        right_x_axis,
        right_y_axis,
    ):
        if not annotation:
            raise exceptions.PreventUpdate

        triggered_prop = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        # Keep left/right plots synced when the right plot was zoomed/reset.
        same_embedding_view = (
            clustering_method == right_clustering
            and x_axis == right_x_axis
            and y_axis == right_y_axis
        )
        cross_relayout = gene_relayout if (
            triggered_prop == f"{prefix}-gene-scatter.relayoutData" and same_embedding_view
        ) else None
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
        filtered_cell_idx = _filtered_cell_indices(filtered_data)
        filtered_sig = filtered_data_signature(filtered_data)

        is_continuous = annotation in adata.var_names or is_continuous_annotation(adata, annotation)
        render_backend = embedding_render_backend
        n_annotation_categories = adata.obs[annotation].nunique() if annotation in adata.obs.columns else 0
        discrete_palette = resolve_discrete_palette(
            discrete_color_map,
            n_annotation_categories,
            default=color_config,
        )
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
            source_adata=adata,
            cell_indices=filtered_cell_idx,
        )
        cached_figure_set(cache_key, fig)
        return apply_relayout(fig, effective_relayout)

    @app.callback(
        Output(f"{prefix}-gene-scatter", "figure"),
        [
            Input(f"{prefix}-scatter-gene-selection", "value"),
            Input(f"{prefix}-right-clustering-dropdown", "value"),
            Input(f"{prefix}-right-x-axis", "value"),
            Input(f"{prefix}-right-y-axis", "value"),
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
            Input(f"{prefix}-left-highlighted-cells-store", "data"),
        ],
        [
            State(f"{prefix}-gene-scatter", "relayoutData"),
            State(f"{prefix}-clustering-dropdown", "value"),
            State(f"{prefix}-x-axis", "value"),
            State(f"{prefix}-y-axis", "value"),
        ],
    )
    def update_gene_scatter(
        gene_name,
        right_clustering,
        right_x_axis,
        right_y_axis,
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
        highlighted_cells,
        gene_relayout,
        left_clustering,
        left_x_axis,
        left_y_axis,
    ):
        if not gene_name:
            raise exceptions.PreventUpdate

        triggered_prop = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        # Keep left/right plots synced when the left plot was zoomed/reset.
        same_embedding_view = (
            right_clustering == left_clustering
            and right_x_axis == left_x_axis
            and right_y_axis == left_y_axis
        )
        cross_relayout = annotation_relayout if (
            triggered_prop == f"{prefix}-annotation-scatter.relayoutData" and same_embedding_view
        ) else None
        # Preserve current view on style-only updates; reset on geometry/data changes.
        if triggered_prop in {
            f"{prefix}-right-clustering-dropdown.value",
            f"{prefix}-right-x-axis.value",
            f"{prefix}-right-y-axis.value",
            f"{prefix}-global-filtered-data.data",
            f"{prefix}-spatial-imgkey-dropdown.value",
        }:
            self_relayout = None
        else:
            self_relayout = gene_relayout
        effective_relayout = cross_relayout if cross_relayout is not None else self_relayout

        render_backend = embedding_render_backend
        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        filtered_cell_idx = _filtered_cell_indices(filtered_data)
        highlighted_cell_ids = (highlighted_cells or {}).get("cell_ids")
        cache_key = make_cache_key(
            "gene_scatter",
            adata,
            gene_name=gene_name,
            clustering=right_clustering,
            x_axis=right_x_axis,
            y_axis=right_y_axis,
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
            highlighted_cells=_list_signature(highlighted_cell_ids),
        )
        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return apply_relayout(cached_fig, effective_relayout)

        if gene_name in adata.var_names:
            if coexpression_mode == "coexpression" and gene2_name:
                fig = plot_coexpression_embedding(
                    adata=plot_adata,
                    embedding_key=right_clustering,
                    gene1=gene_name,
                    gene2=gene2_name,
                    x_axis=right_x_axis,
                    y_axis=right_y_axis,
                    threshold1=threshold1,
                    threshold2=threshold2,
                    transformation=transformation,
                    color_map=None,
                    marker_size=marker_size,
                    opacity=opacity,
                    legend_show=legend_show,
                    axis_show=axis_show,
                    source_adata=adata,
                    cell_indices=filtered_cell_idx,
                    highlighted_cell_ids=highlighted_cell_ids,
                )
            else:
                fig = plot_embedding(
                    adata=plot_adata,
                    embedding_key=right_clustering,
                    color=gene_name,
                    x_axis=right_x_axis,
                    y_axis=right_y_axis,
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
                    source_adata=adata,
                    cell_indices=filtered_cell_idx,
                    highlighted_cell_ids=highlighted_cell_ids,
                )
        elif is_continuous_annotation(adata, gene_name):
            fig = plot_embedding(
                adata=plot_adata,
                embedding_key=right_clustering,
                color=gene_name,
                x_axis=right_x_axis,
                y_axis=right_y_axis,
                mode="continuous",
                transformation=None,
                order=order,
                continuous_color_map=color_map or "Viridis",
                marker_size=marker_size,
                opacity=opacity,
                render_backend=render_backend,
                axis_show=axis_show,
                img_key=spatial_img_key,
                source_adata=adata,
                cell_indices=filtered_cell_idx,
                highlighted_cell_ids=highlighted_cell_ids,
            )
        else:
            n_gene_categories = adata.obs[gene_name].nunique() if gene_name in adata.obs.columns else 0
            discrete_color_map_value = resolve_discrete_palette(
                discrete_color_map,
                n_gene_categories,
                default=color_config,
            )
            fig = plot_embedding(
                adata=plot_adata,
                adata_full=adata,
                embedding_key=right_clustering,
                color=gene_name,
                x_axis=right_x_axis,
                y_axis=right_y_axis,
                mode="categorical",
                discrete_color_map=discrete_color_map_value,
                marker_size=marker_size,
                opacity=opacity,
                render_backend=render_backend,
                legend_show=legend_show,
                axis_show=axis_show,
                img_key=spatial_img_key,
                source_adata=adata,
                cell_indices=filtered_cell_idx,
                highlighted_cell_ids=highlighted_cell_ids,
            )
        cached_figure_set(cache_key, fig)
        return apply_relayout(fig, effective_relayout)

    def _extract_cell_ids_from_customdata(customdata, plot_adata=None):
        if customdata is None:
            return None
        if isinstance(customdata, (list, tuple)):
            if customdata and str(customdata[0]) in adata.obs_names:
                return str(customdata[0])
            if customdata and plot_adata is not None:
                try:
                    return str(plot_adata.obs.index[int(customdata[0])])
                except Exception:
                    return None
        try:
            if str(customdata) in adata.obs_names:
                return str(customdata)
            if plot_adata is not None:
                return str(plot_adata.obs.index[int(customdata)])
        except Exception:
            return None
        return None

    def _unique_cell_ids(cell_ids):
        seen = set()
        unique = []
        for cell_id in cell_ids:
            if cell_id is not None and cell_id not in seen:
                unique.append(cell_id)
                seen.add(cell_id)
        return unique

    @app.callback(
        Output(f"{prefix}-left-highlighted-cells-store", "data"),
        [
            Input(f"{prefix}-annotation-scatter", "selectedData"),
            Input(f"{prefix}-annotation-scatter", "restyleData"),
        ],
        [
            State(f"{prefix}-annotation-scatter", "figure"),
            State(f"{prefix}-annotation-dropdown", "value"),
            State(f"{prefix}-global-filtered-data", "data"),
            State(f"{prefix}-left-highlighted-cells-store", "data"),
        ],
    )
    def update_left_highlighted_cells(selected_data, restyle_data, current_figure, current_annotation, filtered_data, current_store):
        triggered_prop = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        plot_adata = resolve_plot_adata_from_filter(filtered_data)

        if triggered_prop == f"{prefix}-annotation-scatter.selectedData":
            if not selected_data or not selected_data.get("points"):
                return None
            cell_ids = _unique_cell_ids(
                _extract_cell_ids_from_customdata(point.get("customdata"), plot_adata)
                for point in selected_data.get("points", [])
            )
            return {"cell_ids": cell_ids, "source": "lasso", "hidden_traces": (current_store or {}).get("hidden_traces", [])}

        if triggered_prop != f"{prefix}-annotation-scatter.restyleData" or not restyle_data or not current_figure:
            return current_store

        hidden_traces = set((current_store or {}).get("hidden_traces", []))
        update = restyle_data[0] if len(restyle_data) > 0 else {}
        trace_indices = restyle_data[1] if len(restyle_data) > 1 else []
        visible_values = update.get("visible")
        if not isinstance(visible_values, list):
            visible_values = [visible_values] * len(trace_indices)

        for trace_index, visible in zip(trace_indices, visible_values):
            trace_index = int(trace_index)
            if visible in ("legendonly", False):
                hidden_traces.add(trace_index)
            else:
                hidden_traces.discard(trace_index)

        traces = current_figure.get("data", [])
        selectable_trace_indices = [
            idx for idx, trace in enumerate(traces)
            if trace.get("name") != "Background" and trace.get("showlegend", False)
        ]
        if not hidden_traces:
            return None

        visible_trace_indices = [idx for idx in selectable_trace_indices if idx not in hidden_traces]
        cell_ids = []
        for idx in visible_trace_indices:
            for customdata in traces[idx].get("customdata", []) or []:
                cell_ids.append(_extract_cell_ids_from_customdata(customdata, plot_adata))
        cell_ids = _unique_cell_ids(cell_ids)
        return {"cell_ids": cell_ids, "source": "legend", "hidden_traces": sorted(hidden_traces)}

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
        filtered_cell_idx = _filtered_cell_indices(filtered_data)
        default_min, default_max, default_value = 0, 1, 0.5
        from guanaco.utils.gene_extraction_utils import extract_gene_expression, apply_transformation

        if gene1 and gene1 in adata.var_names:
            gene1_expr = extract_gene_expression(adata, gene1)
            if filtered_cell_idx is not None:
                gene1_expr = gene1_expr[filtered_cell_idx]
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

        if coexpression_mode == "coexpression" and gene2 and gene2 in adata.var_names:
            gene2_expr = extract_gene_expression(adata, gene2)
            if filtered_cell_idx is not None:
                gene2_expr = gene2_expr[filtered_cell_idx]
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
                        if customdata and str(customdata[0]) in plot_adata.obs.index:
                            selected_indices.append(str(customdata[0]))
                            continue
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
