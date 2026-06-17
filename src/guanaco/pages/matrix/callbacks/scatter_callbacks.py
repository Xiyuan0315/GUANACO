import numpy as np
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, callback_context, exceptions, no_update

from guanaco.utils.colors import resolve_discrete_palette

# Auto-dismiss delay (ms) for the cell-selection status toasts.
_SELECTION_ALERT_DURATION_MS = 4000


def _is_reset_relayout(relayout):
    """True for a double-click 'reset axes' (emits autorange) vs. zoom/pan (ranges)."""
    return bool(relayout) and (
        "xaxis.autorange" in relayout or "yaxis.autorange" in relayout
    )


# Clientside reset-link: when either scatter is reset (double-click -> autorange),
# reset it in the browser, and reset the other panel too when both share the same
# dimension reduction. Done client-side so it resets the clicked plot (the server
# cross-link only ever reset the *other* one) and never gates the initial render.
# For a raster (categorical datashader = a layout image with no autorange-able
# data) it restores the image's data-coord extent instead of autoranging to
# nothing. __LEFT_ID__/__RIGHT_ID__ are substituted with the two graph ids.
_AXIS_RESET_LINK_JS = """
function(leftRelayout, rightRelayout, leftClu, rightClu, leftX, rightX, leftY, rightY) {
    const noUpdate = window.dash_clientside.no_update;
    const ctx = window.dash_clientside.callback_context;
    if (!ctx || !ctx.triggered || ctx.triggered.length === 0) return noUpdate;
    const trig = ctx.triggered[0];
    const rl = trig.value;
    if (!rl) return noUpdate;
    // Only react to a reset (double-click emits autorange); ignore zoom/pan.
    if (!('xaxis.autorange' in rl) && !('yaxis.autorange' in rl)) return noUpdate;
    // Debounce so the relayout we trigger below doesn't re-enter this callback.
    const now = Date.now();
    if (window.__guanacoAxisReset && (now - window.__guanacoAxisReset) < 350) return noUpdate;
    window.__guanacoAxisReset = now;

    const LEFT = '__LEFT_ID__', RIGHT = '__RIGHT_ID__';
    const prop = trig.prop_id || '';
    // Same dimension reduction => the two plots are linked.
    const sameEmbedding = (leftClu === rightClu) && (leftX === rightX) && (leftY === rightY);
    // Always reset the plot that was double-clicked; reset the other only if linked.
    const targets = [];
    if (prop.indexOf(LEFT + '.') === 0) {
        targets.push(LEFT);
        if (sameEmbedding) targets.push(RIGHT);
    } else if (prop.indexOf(RIGHT + '.') === 0) {
        targets.push(RIGHT);
        if (sameEmbedding) targets.push(LEFT);
    } else {
        return noUpdate;
    }

    function resetGraph(id) {
        const wrap = document.getElementById(id);
        if (!wrap) return;
        const gd = wrap.classList.contains('js-plotly-plot') ? wrap : wrap.querySelector('.js-plotly-plot');
        if (!gd || !window.Plotly) return;
        const lay = gd.layout || {};
        const imgs = lay.images || [];
        if (imgs.length > 0 && imgs[0].sizex != null && imgs[0].sizey != null) {
            // Raster: restore the image extent (x from left edge, y from top edge).
            const im = imgs[0];
            window.Plotly.relayout(gd, {
                'xaxis.range': [im.x, im.x + im.sizex],
                'yaxis.range': [im.y - im.sizey, im.y],
                'xaxis.autorange': false,
                'yaxis.autorange': false
            });
        } else {
            window.Plotly.relayout(gd, {'xaxis.autorange': true, 'yaxis.autorange': true});
        }
    }
    targets.forEach(resetGraph);
    return noUpdate;
}
"""


# Client-side cross-highlight: when the left plot's selection/legend changes, grey
# out the de-selected cells on the *right* plot without rebuilding it. The store
# carries the visible cells' row positions; each right-plot trace carries those
# positions (continuous/coexpression: customdata is the position; categorical:
# customdata's last column). We translate positions -> per-trace point indices and
# set them via Plotly.restyle(selectedpoints), so unselected points dim to the
# trace's `unselected.marker.opacity` -- a pure re-style, no server round-trip and
# no figure rebuild. Re-runs on a right-plot rebuild (figure Input) to re-apply the
# current highlight. __RIGHT_ID__ is substituted with the right graph's id.
_RIGHT_HIGHLIGHT_JS = """
function(highlightData, _rightFigure) {
    const noUpdate = window.dash_clientside.no_update;
    const RIGHT = '__RIGHT_ID__';
    const wrap = document.getElementById(RIGHT);
    if (!wrap) return noUpdate;
    const gd = wrap.classList.contains('js-plotly-plot') ? wrap : wrap.querySelector('.js-plotly-plot');
    if (!gd || !window.Plotly || !gd.data) return noUpdate;

    const positions = (highlightData && highlightData.positions) ? highlightData.positions : null;
    const posSet = positions ? new Set(positions) : null;

    const traceIdx = [];
    const selected = [];
    for (let t = 0; t < gd.data.length; t++) {
        const tr = gd.data[t];
        if (!tr.customdata) continue;
        traceIdx.push(t);
        if (posSet === null) {
            selected.push(null);  // no highlight -> clear selection, all points normal
            continue;
        }
        const cd = tr.customdata;
        const picked = [];
        for (let i = 0; i < cd.length; i++) {
            const ci = cd[i];
            const pos = Array.isArray(ci) ? ci[ci.length - 1] : ci;
            if (posSet.has(pos)) picked.push(i);
        }
        selected.push(picked);
    }
    if (traceIdx.length === 0) return noUpdate;
    window.Plotly.restyle(gd, {selectedpoints: selected}, traceIdx);
    return noUpdate;
}
"""


# Client-side debounce for legend select/deselect. Toggling legend entries fires a
# `restyleData` event per click; recomputing the cross-highlight server-side on every
# one makes rapid clicking (or quickly hiding several labels in a row) stutter. This
# collapses a burst of clicks into a single trailing update: each event (re)arms a
# 250 ms timer, and only when clicking settles do we read the *live* trace
# visibilities off the left graph and push the full set of hidden label names to the
# debounce store via set_props. The server then recomputes the highlight once.
# __LEFT_ID__ / __HIDDEN_STORE_ID__ are substituted with the graph id and store id.
_LEGEND_DEBOUNCE_JS = """
function(restyleData) {
    const noUpdate = window.dash_clientside.no_update;
    const LEFT = '__LEFT_ID__';
    const STORE = '__HIDDEN_STORE_ID__';
    const wrap = document.getElementById(LEFT);
    if (!wrap) return noUpdate;
    window.__guanacoLegendTimers = window.__guanacoLegendTimers || {};
    if (window.__guanacoLegendTimers[LEFT]) {
        clearTimeout(window.__guanacoLegendTimers[LEFT]);
    }
    window.__guanacoLegendTimers[LEFT] = setTimeout(function() {
        const gd = wrap.classList.contains('js-plotly-plot') ? wrap : wrap.querySelector('.js-plotly-plot');
        if (!gd || !gd.data) return;
        const hidden = [];
        for (let t = 0; t < gd.data.length; t++) {
            const tr = gd.data[t];
            if (!tr.name) continue;
            if (tr.visible === 'legendonly' || tr.visible === false) hidden.push(tr.name);
        }
        if (window.dash_clientside.set_props) {
            window.dash_clientside.set_props(STORE, {data: {labels: hidden, t: Date.now()}});
        }
    }, 250);
    return noUpdate;
}
"""


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
):
    def _filtered_cell_indices(filtered_data):
        if (
            filtered_data
            and filtered_data.get("cell_indices") is not None
            and filtered_data.get("n_cells", adata.n_obs) < adata.n_obs
        ):
            return np.asarray(filtered_data["cell_indices"], dtype=np.int64)
        return None

    def _resolve_layer(data_layer):
        # "X" is the sentinel for the default matrix, i.e. no explicit layer.
        return data_layer if data_layer and data_layer != "X" else None

    def _triggered_prop():
        return callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""

    def _combined_search_options(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        all_matches = search_combined(
            obs_columns, obs_columns_lower, var_names, var_names_lower, search_value, limit=10
        )
        return [{"label": item, "value": item} for item in all_matches]

    def _resolve_discrete_palette_for(annotation, discrete_color_map):
        n_categories = adata.obs[annotation].nunique() if annotation in adata.obs.columns else 0
        return resolve_discrete_palette(discrete_color_map, n_categories, default=color_config)

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
        return _combined_search_options(search_value)

    @app.callback(
        Output(f"{prefix}-scatter-gene-selection", "options"),
        Input(f"{prefix}-scatter-gene-selection", "search_value"),
    )
    def update_scatter_gene_selection(search_value):
        return _combined_search_options(search_value)

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
            Input(f"{prefix}-data-layer", "value"),
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
        data_layer,
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

        layer = _resolve_layer(data_layer)

        triggered_prop = _triggered_prop()

        # Continuous renders ignore the discrete colormap and legend toggle, so a
        # change to either control produces an identical figure -- skip the rebuild
        # those shared controls would otherwise force when this panel is continuous.
        is_continuous = annotation in adata.var_names or is_continuous_annotation(adata, annotation)
        if is_continuous and triggered_prop in {
            f"{prefix}-discrete-color-map-dropdown.value",
            f"{prefix}-scatter-legend-toggle.value",
        }:
            return no_update
        # Conversely, a categorical render ignores the continuous colormap, so don't
        # rebuild the discrete embedding when that dropdown changes.
        if not is_continuous and triggered_prop == f"{prefix}-scatter-color-map-dropdown.value":
            return no_update

        if triggered_prop == f"{prefix}-gene-scatter.relayoutData" and _is_reset_relayout(gene_relayout):
            return no_update
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

        render_backend = embedding_render_backend
        discrete_palette = _resolve_discrete_palette_for(annotation, discrete_color_map)
        fig = plot_embedding(
            adata=plot_adata,
            adata_full=adata,
            embedding_key=clustering_method,
            color=annotation,
            x_axis=x_axis,
            y_axis=y_axis,
            mode="continuous" if is_continuous else "categorical",
            layer=layer,
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
        return apply_relayout(fig, effective_relayout)

    @app.callback(
        Output(f"{prefix}-gene-scatter", "figure"),
        [
            Input(f"{prefix}-scatter-gene-selection", "value"),
            Input(f"{prefix}-right-clustering-dropdown", "value"),
            Input(f"{prefix}-right-x-axis", "value"),
            Input(f"{prefix}-right-y-axis", "value"),
            Input(f"{prefix}-data-layer", "value"),
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
        data_layer,
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
        left_clustering,
        left_x_axis,
        left_y_axis,
    ):
        if not gene_name:
            raise exceptions.PreventUpdate

        layer = _resolve_layer(data_layer)

        triggered_prop = _triggered_prop()

        # Classify this panel's render mode up-front. The discrete colormap only feeds
        # the categorical render; the legend toggle only feeds categorical/coexpression.
        # When the right plot is showing a continuous gene/annotation, a change to
        # either control yields an identical figure -- skip the full rebuild it would
        # otherwise force (these are shared controls the left plot legitimately uses).
        if gene_name in adata.var_names:
            right_mode = "coexpression" if (coexpression_mode == "coexpression" and gene2_name) else "continuous"
        elif is_continuous_annotation(adata, gene_name):
            right_mode = "continuous"
        else:
            right_mode = "categorical"
        if triggered_prop == f"{prefix}-discrete-color-map-dropdown.value" and right_mode != "categorical":
            return no_update
        if triggered_prop == f"{prefix}-scatter-legend-toggle.value" and right_mode == "continuous":
            return no_update
        # The continuous colormap only feeds the continuous render (coexpression and
        # categorical ignore it), so don't rebuild those modes when it changes.
        if triggered_prop == f"{prefix}-scatter-color-map-dropdown.value" and right_mode != "continuous":
            return no_update

        # A reset of the other plot is handled entirely by the clientside reset-link
        # (resets both panels in the browser). Don't rebuild/sync it here, so the
        # server can't clobber the raster reset with an autorange.
        if triggered_prop == f"{prefix}-annotation-scatter.relayoutData" and _is_reset_relayout(annotation_relayout):
            return no_update
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
        # The left plot's selection no longer rebuilds this figure: the cross-highlight
        # (grey-out of deselected cells) is applied client-side via selectedpoints, so
        # the right plot always renders the full cell set and keeps every legend entry.

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
                    layer=layer,
                    color_map=None,
                    marker_size=marker_size,
                    opacity=opacity,
                    legend_show=legend_show,
                    axis_show=axis_show,
                    # Pass the chosen spatial image so co-expression uses the same
                    # hires/lowres tissue image (and its scalefactor) as the other
                    # panels -- without it the coexpression view auto-selected hires
                    # regardless of the dropdown, so the two panels' coordinate ranges
                    # diverged and zoom no longer lined up.
                    img_key=spatial_img_key,
                    source_adata=adata,
                    cell_indices=filtered_cell_idx,
                )
            else:
                fig = plot_embedding(
                    adata=plot_adata,
                    embedding_key=right_clustering,
                    color=gene_name,
                    x_axis=right_x_axis,
                    y_axis=right_y_axis,
                    mode="continuous",
                    layer=layer,
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
                )
        elif is_continuous_annotation(adata, gene_name):
            fig = plot_embedding(
                adata=plot_adata,
                embedding_key=right_clustering,
                color=gene_name,
                x_axis=right_x_axis,
                y_axis=right_y_axis,
                mode="continuous",
                order=order,
                continuous_color_map=color_map or "Viridis",
                marker_size=marker_size,
                opacity=opacity,
                render_backend=render_backend,
                axis_show=axis_show,
                img_key=spatial_img_key,
                source_adata=adata,
                cell_indices=filtered_cell_idx,
            )
        else:
            discrete_color_map_value = _resolve_discrete_palette_for(gene_name, discrete_color_map)
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
            )
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
            Input(f"{prefix}-legend-hidden-store", "data"),
            Input(f"{prefix}-global-filtered-data", "data"),
        ],
        [
            State(f"{prefix}-annotation-dropdown", "value"),
            State(f"{prefix}-left-highlighted-cells-store", "data"),
        ],
    )
    def update_left_highlighted_cells(selected_data, legend_hidden, filtered_data, current_annotation, current_store):
        triggered_prop = _triggered_prop()
        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        src = plot_adata if plot_adata is not None else adata

        # The highlight is expressed as row positions in the current plotted data, which
        # the global filter re-bases. Clear it on a filter change so the client-side
        # cross-highlight never applies stale positions to the rebuilt right plot.
        if triggered_prop == f"{prefix}-global-filtered-data.data":
            return None

        if triggered_prop == f"{prefix}-annotation-scatter.selectedData":
            if not selected_data or not selected_data.get("points"):
                return None
            cell_ids = _unique_cell_ids(
                _extract_cell_ids_from_customdata(point.get("customdata"), plot_adata)
                for point in selected_data.get("points", [])
            )
            # Map lassoed cell ids to row positions in the plotted data: the right
            # plot's traces carry these positions, so it highlights them client-side.
            positions = src.obs.index.get_indexer(cell_ids)
            return {
                "positions": positions[positions >= 0].tolist(),
                "source": "lasso",
                "hidden_labels": (current_store or {}).get("hidden_labels", []),
            }

        # Legend select/deselect (debounced clientside): the store carries the settled
        # set of hidden label *names* read off the live legend, so we recompute the
        # visible-cell positions once per click-burst instead of on every toggle.
        if triggered_prop != f"{prefix}-legend-hidden-store.data" or not legend_hidden:
            return current_store

        hidden_labels = set(legend_hidden.get("labels", []))
        # Nothing hidden, or the panel is showing a continuous gene/obs (no categorical
        # legend, so the column isn't in obs) -> clear the highlight.
        if not hidden_labels or not current_annotation or current_annotation not in src.obs.columns:
            return None

        # Visible cells = obs rows whose annotation value is NOT in hidden_labels.
        # Emit their row positions (single C-level pass) for the client-side
        # cross-highlight; positions index the right plot's trace customdata.
        visible_mask = ~src.obs[current_annotation].astype(str).isin(hidden_labels)
        positions = np.flatnonzero(visible_mask.to_numpy())
        return {"positions": positions.tolist(), "source": "legend", "hidden_labels": sorted(hidden_labels)}

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
            Input(f"{prefix}-data-layer", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
        ],
    )
    def update_threshold_ranges(gene1, gene2, coexpression_mode, data_layer, filtered_data):
        filtered_cell_idx = _filtered_cell_indices(filtered_data)
        default_min, default_max, default_value = 0, 1, 0.5
        layer = _resolve_layer(data_layer)
        from guanaco.utils.gene_extraction_utils import extract_gene_expression

        if gene1 and gene1 in adata.var_names:
            gene1_expr = extract_gene_expression(adata, gene1, layer=layer)
            if filtered_cell_idx is not None:
                gene1_expr = gene1_expr[filtered_cell_idx]
            if gene1_expr.max() > gene1_expr.min():
                gene1_min = float(gene1_expr.min())
                gene1_max = float(gene1_expr.max())
                gene1_value = (gene1_min + gene1_max) / 2
            else:
                gene1_min, gene1_max, gene1_value = default_min, default_max, default_value
        else:
            gene1_min, gene1_max, gene1_value = default_min, default_max, default_value

        if coexpression_mode == "coexpression" and gene2 and gene2 in adata.var_names:
            gene2_expr = extract_gene_expression(adata, gene2, layer=layer)
            if filtered_cell_idx is not None:
                gene2_expr = gene2_expr[filtered_cell_idx]
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
                duration=_SELECTION_ALERT_DURATION_MS,
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
            for point in selected_points:
                if "customdata" in point:
                    try:
                        row_idx = int(point["customdata"])
                        selected_indices.append(plot_adata.obs.index[row_idx])
                    except (IndexError, ValueError, TypeError):
                        pass

        if selected_indices:
            n_selected = len(selected_indices)
            status_msg = dbc.Alert(
                f"✓ {n_selected} cells selected from {current_annotation}. Other plots updated.",
                color="success",
                dismissable=True,
                duration=_SELECTION_ALERT_DURATION_MS,
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
            raise exceptions.PreventUpdate
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if f"{prefix}-download-cellids" in button_id:
            content = "\n".join(selected_cells)
            return dict(content=content, filename="selected_cells.txt")
        raise exceptions.PreventUpdate

    # Reset-link (see _AXIS_RESET_LINK_JS): client-side double-click reset, linked
    # to the other panel when both use the same dimension reduction.
    app.clientside_callback(
        _AXIS_RESET_LINK_JS
        .replace("__LEFT_ID__", f"{prefix}-annotation-scatter")
        .replace("__RIGHT_ID__", f"{prefix}-gene-scatter"),
        Output(f"{prefix}-axis-reset-link", "data"),
        Input(f"{prefix}-annotation-scatter", "relayoutData"),
        Input(f"{prefix}-gene-scatter", "relayoutData"),
        State(f"{prefix}-clustering-dropdown", "value"),
        State(f"{prefix}-right-clustering-dropdown", "value"),
        State(f"{prefix}-x-axis", "value"),
        State(f"{prefix}-right-x-axis", "value"),
        State(f"{prefix}-y-axis", "value"),
        State(f"{prefix}-right-y-axis", "value"),
        prevent_initial_call=True,
    )

    # Cross-highlight (see _RIGHT_HIGHLIGHT_JS): grey out the left plot's de-selected
    # cells on the right plot client-side. Fires when the left selection changes and
    # when the right plot is rebuilt (re-applies the current highlight to the new figure).
    app.clientside_callback(
        _RIGHT_HIGHLIGHT_JS.replace("__RIGHT_ID__", f"{prefix}-gene-scatter"),
        Output(f"{prefix}-right-highlight-link", "data"),
        Input(f"{prefix}-left-highlighted-cells-store", "data"),
        Input(f"{prefix}-gene-scatter", "figure"),
        prevent_initial_call=True,
    )

    # Legend debounce (see _LEGEND_DEBOUNCE_JS): collapse a burst of legend clicks
    # into one trailing update written to the hidden-store, which the server callback
    # above turns into the cross-highlight -- so rapid select/deselect no longer fires
    # a server round-trip per click.
    app.clientside_callback(
        _LEGEND_DEBOUNCE_JS
        .replace("__LEFT_ID__", f"{prefix}-annotation-scatter")
        .replace("__HIDDEN_STORE_ID__", f"{prefix}-legend-hidden-store"),
        Output(f"{prefix}-legend-debounce-dummy", "data"),
        Input(f"{prefix}-annotation-scatter", "restyleData"),
        prevent_initial_call=True,
    )
