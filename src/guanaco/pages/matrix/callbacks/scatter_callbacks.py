import numpy as np
import dash_bootstrap_components as dbc
from dash import Input, Output, Patch, State, dcc, html, callback_context, exceptions, no_update

from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.search import ranked_substring_matches
from guanaco.utils.obs_utils import obs_col

# Auto-dismiss delay (ms) for the cell-selection status toasts.
_SELECTION_ALERT_DURATION_MS = 4000


# Coordinate-only interactions never need a new figure. Comparing the live axis
# values prevents feedback loops without dropping quick, successive zoom events.
_AXIS_RESET_LINK_JS = """
function(leftRelayout, rightRelayout, leftClu, rightClu, leftX, rightX, leftY, rightY, _filtered, _image, current) {
    const noUpdate = window.dash_clientside.no_update;
    const ctx = window.dash_clientside.callback_context;
    if (!ctx || !ctx.triggered || ctx.triggered.length === 0) return noUpdate;
    const viewports = {...(current || {})};
    const changed = ctx.triggered.map(item => item.prop_id);
    const prefix = '__PREFIX__';
    const geometry = ['clustering-dropdown', 'x-axis', 'y-axis'];
    const clearBoth = changed.includes(prefix + '-global-filtered-data.data') ||
        changed.includes(prefix + '-spatial-imgkey-dropdown.value');
    const clearLeft = clearBoth || geometry.some(name => changed.includes(prefix + '-' + name + '.value'));
    const clearRight = clearBoth || geometry.some(name => changed.includes(prefix + '-right-' + name + '.value'));
    if (clearLeft) delete viewports.left;
    if (clearRight) delete viewports.right;
    if (clearLeft || clearRight) return viewports;
    const trig = ctx.triggered[0];
    const rl = trig.value;
    if (!rl) return noUpdate;
    const LEFT = '__LEFT_ID__', RIGHT = '__RIGHT_ID__';
    const prop = trig.prop_id || '';
    const sameEmbedding = (leftClu === rightClu) && (leftX === rightX) && (leftY === rightY);
    const source = prop === LEFT + '.relayoutData' ? LEFT :
        prop === RIGHT + '.relayoutData' ? RIGHT : null;
    if (!source) return noUpdate;
    const reset = rl['xaxis.autorange'] === true || rl['yaxis.autorange'] === true;
    if (!reset && !['xaxis', 'yaxis'].some(axis =>
        rl[axis + '.range'] || (rl[axis + '.range[0]'] != null && rl[axis + '.range[1]'] != null)
    )) return noUpdate;
    const targets = reset ? [source] : [];
    if (sameEmbedding) targets.push(source === LEFT ? RIGHT : LEFT);
    for (const id of targets) {
        const wrap = document.getElementById(id);
        if (!wrap) continue;
        const gd = wrap.classList.contains('js-plotly-plot') ? wrap : wrap.querySelector('.js-plotly-plot');
        if (!gd || !window.Plotly) continue;
        const lay = gd.layout || {};
        const im = (lay.images || [])[0];
        const raster = reset && im && im.sizex != null && im.sizey != null;
        const update = {};
        for (const axis of ['xaxis', 'yaxis']) {
            const current = lay[axis] || {};
            let range = rl[axis + '.range'];
            if (!range && rl[axis + '.range[0]'] != null && rl[axis + '.range[1]'] != null) {
                range = [rl[axis + '.range[0]'], rl[axis + '.range[1]']];
            }
            if (raster) range = axis === 'xaxis' ? [im.x, im.x + im.sizex] : [im.y - im.sizey, im.y];
            if (reset && !raster) {
                if (current.autorange !== true) update[axis + '.autorange'] = true;
            } else if (range && (!current.range || current.range[0] !== range[0] ||
                       current.range[1] !== range[1] || current.autorange !== false)) {
                update[axis + '.range'] = range;
                update[axis + '.autorange'] = false;
            }
        }
        if (Object.keys(update).length) window.Plotly.relayout(gd, update);
    }
    // Keep ranges separate from relayoutData: a later font/size relayout replaces
    // that event payload, but must not lose the viewport on the next gene load.
    for (const id of new Set([source, ...targets])) {
        const wrap = document.getElementById(id);
        const gd = wrap && (wrap.classList.contains('js-plotly-plot') ? wrap : wrap.querySelector('.js-plotly-plot'));
        if (!gd || !gd.layout) continue;
        const ranges = {};
        for (const axis of ['xaxis', 'yaxis']) {
            const value = gd.layout[axis] || {};
            if (value.autorange === true) ranges[axis + '.autorange'] = true;
            else if (value.range) {
                ranges[axis + '.range'] = value.range.slice();
                ranges[axis + '.autorange'] = false;
            }
        }
        viewports[id === LEFT ? 'left' : 'right'] = ranges;
    }
    return viewports;
}
"""


_SCATTER_STYLE_JS = """
function(size, opacity, axisShow, leftFigure, rightFigure) {
    const changed = (window.dash_clientside.callback_context.triggered || []).map(item => item.prop_id);
    const redraw = !changed.length || changed.some(prop => prop.endsWith('.figure'));
    const resize = redraw || changed.some(prop => prop.endsWith('-marker-size-slider.value'));
    const fade = redraw || changed.some(prop => prop.endsWith('-opacity-slider.value'));
    const axes = redraw || changed.some(prop => prop.endsWith('-axis-toggle.value'));
    for (const id of ['__LEFT_ID__', '__RIGHT_ID__', '__PSEUDOTIME_ID__']) {
        const wrap = document.getElementById(id);
        if (!wrap || !window.Plotly) continue;
        const gd = wrap.classList.contains('js-plotly-plot') ? wrap : wrap.querySelector('.js-plotly-plot');
        if (!gd || !gd.data) continue;
        const indices = [], selectedOpacity = [];
        gd.data.forEach((trace, i) => {
            if (trace.type !== 'scattergl' || !trace.mode || !trace.mode.includes('markers')) return;
            indices.push(i);
            // Categorical plots follow the opacity slider even when selected;
            // continuous/coexpression plots intentionally keep selections opaque.
            selectedOpacity.push(trace.legendgroup != null ? opacity : 1);
        });
        const style = {};
        if (resize) style['marker.size'] = size;
        if (fade) {
            style['marker.opacity'] = opacity;
            if (id !== '__PSEUDOTIME_ID__') style['selected.marker.opacity'] = selectedOpacity;
        }
        if (indices.length && (resize || fade)) window.Plotly.restyle(gd, style, indices);
        if (axes && id !== '__PSEUDOTIME_ID__') window.Plotly.relayout(gd, {
            'xaxis.tickfont.color': axisShow ? 'black' : 'rgba(0,0,0,0)',
            'yaxis.tickfont.color': axisShow ? 'black' : 'rgba(0,0,0,0)'
        });
    }
    return window.dash_clientside.no_update;
}
"""


def register_scatter_display_callbacks(app, prefix):
    """Shared browser-only styles and controls for paired/unpaired panels."""
    app.clientside_callback(
        _SCATTER_STYLE_JS
        .replace("__LEFT_ID__", f"{prefix}-annotation-scatter")
        .replace("__RIGHT_ID__", f"{prefix}-gene-scatter")
        .replace("__PSEUDOTIME_ID__", f"{prefix}-pseudotime-plot"),
        Output(f"{prefix}-scatter-style-link", "data"),
        Input(f"{prefix}-marker-size-slider", "value"),
        Input(f"{prefix}-opacity-slider", "value"),
        Input(f"{prefix}-axis-toggle", "value"),
        Input(f"{prefix}-annotation-scatter", "figure"),
        Input(f"{prefix}-gene-scatter", "figure"),
        Input(f"{prefix}-pseudotime-plot", "figure", allow_optional=True),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        """function(n) {
            return n % 2 ? [{display: 'block'}, 'Hide controls'] :
                           [{display: 'none'}, 'More controls'];
        }""",
        Output(f"{prefix}-controls-container", "style"),
        Output(f"{prefix}-toggle-button", "children"),
        Input(f"{prefix}-toggle-button", "n_clicks"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        """function(mode) {
            const display = mode === 'coexpression' ? 'block' : 'none';
            return [{display: display}, {display: display}];
        }""",
        Output(f"{prefix}-gene2-container", "style"),
        Output(f"{prefix}-threshold-container", "style"),
        Input(f"{prefix}-coexpression-toggle", "value"),
    )


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


# Client-side render-metadata extractor for the right (gene) scatter. Reads the
# figure that's already in the browser and stores only two booleans/counts --
# whether a tissue image is present and how many traces there are. The gene-switch
# server callback reads this tiny store (instead of State(figure)) to decide whether
# it can patch just the point trace, so the multi-MB base64 image never round-trips
# to the server on a spatial gene switch.
_GENE_SCATTER_META_JS = """
function(figure) {
    if (!figure || !figure.layout) return {hasImage: false, nTraces: 0};
    const imgs = figure.layout.images || [];
    const data = figure.data || [];
    return {hasImage: imgs.length > 0, nTraces: data.length};
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
    if (!restyleData || !restyleData[0] || !('visible' in restyleData[0])) return noUpdate;
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
    multiomics_source=None,
):
    # Raster marker size is baked into the image; only that backend needs a
    # server render on style changes. States preserve styles on later data loads.
    marker_dependency = Input if embedding_render_backend == "datashader" else State
    register_scatter_display_callbacks(app, prefix)

    def _is_feature(value):
        if multiomics_source is not None:
            return multiomics_source.is_feature(value)
        return value in adata.var_names if value else False

    def _materialize(features, embedding):
        if multiomics_source is None:
            return adata
        return multiomics_source.materialize(features, embedding=embedding)

    def _filtered_adata(source_adata, filtered_data):
        indices = _filtered_cell_indices(filtered_data)
        return source_adata[indices] if indices is not None else source_adata

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
        n_categories = obs_col(adata.obs, annotation).nunique() if annotation in adata.obs.columns else 0
        return resolve_discrete_palette(discrete_color_map, n_categories, default=color_config)


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
            # Skip non-library scalar flags (e.g. squidpy's "is_single" bool) that
            # can sit alongside the real per-library dicts in adata.uns["spatial"].
            if not isinstance(lib_data, dict):
                continue
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
        matching_genes = ranked_substring_matches(
            var_names,
            search_value,
            limit=20,
            match_values=var_names_lower,
        )
        return [{"label": gene, "value": gene} for gene in matching_genes]


    @app.callback(
        Output(f"{prefix}-annotation-scatter", "figure"),
        [
            Input(f"{prefix}-clustering-dropdown", "value"),
            Input(f"{prefix}-x-axis", "value"),
            Input(f"{prefix}-y-axis", "value"),
            Input(f"{prefix}-annotation-dropdown", "value"),
            marker_dependency(f"{prefix}-marker-size-slider", "value"),
            marker_dependency(f"{prefix}-opacity-slider", "value"),
            Input(f"{prefix}-scatter-legend-toggle", "value"),
            State(f"{prefix}-axis-toggle", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-data-layer", "data"),
            Input(f"{prefix}-plot-order", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-spatial-imgkey-dropdown", "value"),
        ],
        [
            State(f"{prefix}-axis-reset-link", "data"),
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
        viewports,
    ):
        if not annotation:
            raise exceptions.PreventUpdate

        layer = _resolve_layer(data_layer)

        triggered_prop = _triggered_prop()

        # Continuous renders ignore the discrete colormap and legend toggle, so a
        # change to either control produces an identical figure -- skip the rebuild
        # those shared controls would otherwise force when this panel is continuous.
        is_continuous = _is_feature(annotation) or is_continuous_annotation(adata, annotation)
        if is_continuous and triggered_prop in {
            f"{prefix}-discrete-color-map-dropdown.value",
            f"{prefix}-scatter-legend-toggle.value",
        }:
            return no_update
        # Conversely, a categorical render ignores the continuous colormap, so don't
        # rebuild the discrete embedding when that dropdown changes.
        if not is_continuous and triggered_prop == f"{prefix}-scatter-color-map-dropdown.value":
            return no_update

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
            self_relayout = (viewports or {}).get("left")

        source_adata = _materialize(
            [annotation] if _is_feature(annotation) else [], clustering_method
        )
        plot_adata = _filtered_adata(source_adata, filtered_data)
        filtered_cell_idx = _filtered_cell_indices(filtered_data)

        render_backend = embedding_render_backend
        discrete_palette = _resolve_discrete_palette_for(annotation, discrete_color_map)
        fig = plot_embedding(
            adata=plot_adata,
            adata_full=source_adata,
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
            source_adata=source_adata,
            cell_indices=filtered_cell_idx,
        )
        return apply_relayout(fig, self_relayout)

    @app.callback(
        Output(f"{prefix}-gene-scatter", "figure"),
        [
            Input(f"{prefix}-scatter-gene-selection", "value"),
            Input(f"{prefix}-right-clustering-dropdown", "value"),
            Input(f"{prefix}-right-x-axis", "value"),
            Input(f"{prefix}-right-y-axis", "value"),
            Input(f"{prefix}-data-layer", "data"),
            Input(f"{prefix}-plot-order", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            marker_dependency(f"{prefix}-marker-size-slider", "value"),
            marker_dependency(f"{prefix}-opacity-slider", "value"),
            State(f"{prefix}-axis-toggle", "value"),
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
            State(f"{prefix}-axis-reset-link", "data"),
            # Lightweight render metadata ({hasImage, nTraces}) maintained client-side
            # from the figure. Used instead of State(..., "figure") so a spatial gene
            # switch doesn't ship the multi-MB base64 tissue image back to the server.
            State(f"{prefix}-gene-scatter-meta", "data"),
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
        axis_show,
        coexpression_mode,
        gene2_name,
        threshold1,
        threshold2,
        legend_show,
        discrete_color_map,
        filtered_data,
        spatial_img_key,
        viewports,
        current_meta,
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
        if _is_feature(gene_name):
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
        if (
            triggered_prop in {
                f"{prefix}-scatter-gene2-selection.value",
                f"{prefix}-gene1-threshold-slider.value",
                f"{prefix}-gene2-threshold-slider.value",
            }
            and right_mode != "coexpression"
        ):
            return no_update

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
            self_relayout = (viewports or {}).get("right")

        render_backend = embedding_render_backend
        requested_features = [gene_name] if _is_feature(gene_name) else []
        if coexpression_mode == "coexpression" and _is_feature(gene2_name):
            requested_features.append(gene2_name)
        source_adata = _materialize(requested_features, right_clustering)
        plot_adata = _filtered_adata(source_adata, filtered_data)
        filtered_cell_idx = _filtered_cell_indices(filtered_data)
        # The left plot's selection no longer rebuilds this figure: the cross-highlight
        # (grey-out of deselected cells) is applied client-side via selectedpoints, so
        # the right plot always renders the full cell set and keeps every legend entry.

        # A spatial single-gene switch changes only the WebGL point trace; the tissue
        # image underneath is identical. When the right plot is already a spatial
        # scatter with the image in place (per the lightweight meta store -- so we
        # never ship the multi-MB figure back to the server), patch just the point
        # trace instead of rebuilding/re-encoding/re-sending the static background.
        patch_spatial_points = (
            triggered_prop == f"{prefix}-scatter-gene-selection.value"
            and right_clustering == "spatial"
            and bool(current_meta)
            and current_meta.get("hasImage")
            and current_meta.get("nTraces") == 1
        )

        if _is_feature(gene_name):
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
                    source_adata=source_adata,
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
                    # Spatial points are always a WebGL trace over the tissue image;
                    # never datashader. Forcing scattergl matters when the background
                    # is suppressed for the patch below, because passing
                    # spatial_image=None would otherwise let the datashader path turn
                    # the points into a raster image (and drop the real trace).
                    render_backend="scattergl" if right_clustering == "spatial" else render_backend,
                    annotation=None,
                    axis_show=axis_show,
                    img_key=spatial_img_key,
                    source_adata=source_adata,
                    cell_indices=filtered_cell_idx,
                    include_spatial_background=not patch_spatial_points,
                )
                # A spatial gene switch changes only the WebGL point trace. Patch
                # that trace in place so Dash does not resend/redecode the static
                # tissue bitmap or recreate the Plotly image layer.
                if patch_spatial_points:
                    patch = Patch()
                    patch["data"][0] = fig.data[0].to_plotly_json()
                    patch["layout"]["title"] = fig.layout.title.to_plotly_json()
                    return patch
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
                source_adata=source_adata,
                cell_indices=filtered_cell_idx,
            )
        else:
            discrete_color_map_value = _resolve_discrete_palette_for(gene_name, discrete_color_map)
            fig = plot_embedding(
                adata=plot_adata,
                adata_full=source_adata,
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
                source_adata=source_adata,
                cell_indices=filtered_cell_idx,
            )
        return apply_relayout(fig, self_relayout)

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
        visible_mask = ~obs_col(src.obs, current_annotation).astype(str).isin(hidden_labels)
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
            Input(f"{prefix}-data-layer", "data"),
            Input(f"{prefix}-global-filtered-data", "data"),
        ],
    )
    def update_threshold_ranges(gene1, gene2, coexpression_mode, data_layer, filtered_data):
        if coexpression_mode != "coexpression":
            raise exceptions.PreventUpdate

        filtered_cell_idx = _filtered_cell_indices(filtered_data)
        default_min, default_max, default_value = 0, 1, 0.5
        layer = _resolve_layer(data_layer)
        from guanaco.utils.gene_extraction_utils import extract_gene_expression

        if gene1 and _is_feature(gene1):
            threshold_adata = _materialize([gene1], None)
            gene1_expr = extract_gene_expression(threshold_adata, gene1, layer=layer)
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

        if coexpression_mode == "coexpression" and gene2 and _is_feature(gene2):
            threshold_adata = _materialize([gene2], None)
            gene2_expr = extract_gene_expression(threshold_adata, gene2, layer=layer)
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
        Output(f"{prefix}-selection-group-store", "data"),
        Output(f"{prefix}-selection-status", "children"),
        [
            Input(f"{prefix}-highlight-plots-button", "n_clicks"),
            Input(f"{prefix}-filter-plots-button", "n_clicks"),
        ],
        [
            State(f"{prefix}-annotation-scatter", "selectedData"),
            State(f"{prefix}-annotation-dropdown", "value"),
            State(f"{prefix}-global-filtered-data", "data"),
        ],
        prevent_initial_call=True,
    )
    def store_selected_cells(
        _highlight_clicks,
        _filter_clicks,
        selected_data,
        current_annotation,
        filtered_data,
    ):
        triggered_id = callback_context.triggered_id
        highlighting = triggered_id == f"{prefix}-highlight-plots-button"
        if not highlighting and triggered_id != f"{prefix}-filter-plots-button":
            return no_update, no_update, no_update

        plot_adata = resolve_plot_adata_from_filter(filtered_data)

        if not selected_data or not selected_data.get("points"):
            selected_indices = plot_adata.obs.index.tolist()
        else:
            selected_points = selected_data["points"]
            selected_indices = []

            if _is_feature(current_annotation) or is_continuous_annotation(plot_adata, current_annotation):
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
            if highlighting:
                n_others = max(plot_adata.n_obs - n_selected, 0)
                message = (
                    f"✓ Highlighting applied: {n_selected} Selected and "
                    f"{n_others} Others."
                )
            else:
                message = f"✓ Filtering applied to {n_selected} cells."
            status_msg = dbc.Alert(
                message,
                color="success",
                dismissable=True,
                duration=_SELECTION_ALERT_DURATION_MS,
            )
            if highlighting:
                universe_cells = None
                if filtered_data and filtered_data.get("cell_indices") is not None:
                    universe_cells = plot_adata.obs.index.tolist()
                selection_data = {
                    "selected_cells": selected_indices,
                    "universe_cells": universe_cells,
                }
                return None, selection_data, status_msg
            return selected_indices, None, status_msg
        return None, None, ""

    app.clientside_callback(
        """function(selected, grouped) {
            return grouped ? [false, true] : selected ? [true, false] : [true, true];
        }""",
        [
            Output(f"{prefix}-highlight-plots-button", "outline"),
            Output(f"{prefix}-filter-plots-button", "outline"),
        ],
        [
            Input(f"{prefix}-selected-cells-hash", "data"),
            Input(f"{prefix}-selection-group-hash", "data"),
        ],
    )

    # Sync zoom, pan and reset in the browser; retain only the small viewport
    # snapshot for future server-side data changes.
    app.clientside_callback(
        _AXIS_RESET_LINK_JS
        .replace("__PREFIX__", prefix)
        .replace("__LEFT_ID__", f"{prefix}-annotation-scatter")
        .replace("__RIGHT_ID__", f"{prefix}-gene-scatter"),
        Output(f"{prefix}-axis-reset-link", "data"),
        Input(f"{prefix}-annotation-scatter", "relayoutData"),
        Input(f"{prefix}-gene-scatter", "relayoutData"),
        Input(f"{prefix}-clustering-dropdown", "value"),
        Input(f"{prefix}-right-clustering-dropdown", "value"),
        Input(f"{prefix}-x-axis", "value"),
        Input(f"{prefix}-right-x-axis", "value"),
        Input(f"{prefix}-y-axis", "value"),
        Input(f"{prefix}-right-y-axis", "value"),
        Input(f"{prefix}-global-filtered-data", "data"),
        Input(f"{prefix}-spatial-imgkey-dropdown", "value"),
        State(f"{prefix}-axis-reset-link", "data"),
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

    # Render-metadata store (see _GENE_SCATTER_META_JS): keep a tiny {hasImage,
    # nTraces} summary of the right-plot figure in sync, so the gene-switch callback
    # can read it as State without shipping the full figure (and its base64 image).
    app.clientside_callback(
        _GENE_SCATTER_META_JS,
        Output(f"{prefix}-gene-scatter-meta", "data"),
        Input(f"{prefix}-gene-scatter", "figure"),
        prevent_initial_call=False,
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
