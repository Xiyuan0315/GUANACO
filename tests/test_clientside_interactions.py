"""Execute the browser callbacks as JavaScript, and check Dash routing boundaries.

Node tests exercise event/state semantics, not browser rendering performance.
"""

import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from dash import Dash

from guanaco.linking.clientside import CELL_HIGHLIGHT_JS
from guanaco.linking.model import link, view
from guanaco.linking.runtime import LinkedView
from guanaco.pages.matrix.callbacks.scatter_callbacks import (
    _AXIS_RESET_LINK_JS,
    _SCATTER_STYLE_JS,
    register_scatter_callbacks,
)
from guanaco.pages.matrix.callbacks.unpaired_multiomics_callbacks import (
    register_unpaired_multiomics_callbacks,
)


def run_js(function, body):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is required to execute clientside callbacks")
    script = """
    const assert = require('node:assert/strict');
    const graphs = {}, calls = [];
    global.window = {
        dash_clientside: {no_update: 'NO_UPDATE', callback_context: {triggered: []}},
        Plotly: {
            relayout: (gd, update) => {
                calls.push({kind: 'relayout', id: gd.id, update});
                for (const [key, value] of Object.entries(update)) {
                    const parts = key.split('.');
                    let dest = gd.layout;
                    for (const part of parts.slice(0, -1)) dest = dest[part] ||= {};
                    dest[parts.at(-1)] = value;
                }
                return Promise.resolve();
            },
            restyle: (gd, update, indices) => {
                calls.push({kind: 'restyle', id: gd.id, update, indices});
                return Promise.resolve();
            }
        }
    };
    global.document = {getElementById: id => graphs[id]};
    function graph(id, data = []) {
        return graphs[id] = {id, data, layout: {
            xaxis: {range: [0, 10], autorange: false},
            yaxis: {range: [0, 10], autorange: false}
        }, classList: {contains: () => true}};
    }
    function trigger(prop_id, value) {
        window.dash_clientside.callback_context.triggered = [{prop_id, value}];
    }
    """ + f"\nconst callback = ({function});\n" + body
    result = subprocess.run([node, "-e", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_zoom_sync_ignores_unrelated_views_and_does_not_echo():
    run_js(_AXIS_RESET_LINK_JS, """
        graph('__LEFT_ID__'); graph('__RIGHT_ID__');
        const range = {'xaxis.range[0]': 2, 'xaxis.range[1]': 4};
        trigger('__LEFT_ID__.relayoutData', range);
        callback(range, null, 'umap', 'pca', 0, 0, 1, 1);
        assert.equal(calls.length, 0);
        callback(range, null, 'umap', 'umap', 0, 0, 1, 1);
        assert.deepEqual(graphs.__RIGHT_ID__.layout.xaxis.range, [2, 4]);
        callback(range, null, 'umap', 'umap', 0, 0, 1, 1);
        assert.equal(calls.length, 1);
        trigger('__LEFT_ID__.relayoutData', {autosize: true});
        callback(null, null, 'umap', 'umap', 0, 0, 1, 1);
        assert.equal(calls.length, 1);
        trigger('__RIGHT_ID__.relayoutData', {'yaxis.range': [8, 3]});
        callback(null, null, 'umap', 'umap', 0, 0, 1, 1);
        assert.deepEqual(graphs.__LEFT_ID__.layout.yaxis.range, [8, 3]);
    """)


def test_reset_restores_raster_extent_and_normal_autorange():
    run_js(_AXIS_RESET_LINK_JS, """
        graph('__LEFT_ID__').layout.images = [{x: 1, y: 9, sizex: 8, sizey: 6}];
        graph('__RIGHT_ID__');
        trigger('__LEFT_ID__.relayoutData', {'xaxis.autorange': true});
        callback(null, null, 'spatial', 'umap', 0, 0, 1, 1);
        assert.deepEqual(graphs.__LEFT_ID__.layout.xaxis.range, [1, 9]);
        assert.deepEqual(graphs.__LEFT_ID__.layout.yaxis.range, [3, 9]);
        assert.equal(calls.length, 1);
        callback(null, null, 'spatial', 'spatial', 0, 0, 1, 1);
        assert.equal(graphs.__RIGHT_ID__.layout.xaxis.autorange, true);
        assert.equal(calls.length, 2);
    """)


def test_viewport_snapshot_survives_style_events_but_resets_on_geometry_changes():
    run_js(_AXIS_RESET_LINK_JS, """
        graph('__LEFT_ID__').layout.xaxis.range = [2, 4];
        graph('__RIGHT_ID__');
        trigger('__LEFT_ID__.relayoutData', {'xaxis.range': [2, 4]});
        const saved = callback(null, null, 'umap', 'umap', 0, 0, 1, 1);
        assert.deepEqual(saved.left['xaxis.range'], [2, 4]);
        assert.deepEqual(saved.right['xaxis.range'], [2, 4]);
        trigger('__RIGHT_ID__.relayoutData', {'xaxis.tickfont.color': 'black'});
        assert.equal(callback(null, null, 'umap', 'umap', 0, 0, 1, 1, null, null, saved), 'NO_UPDATE');
        trigger('__PREFIX__-right-clustering-dropdown.value', 'pca');
        const cleared = callback(null, null, 'umap', 'pca', 0, 0, 1, 1, null, null, saved);
        assert.deepEqual(cleared.left, saved.left);
        assert.equal(cleared.right, undefined);
        trigger('__PREFIX__-global-filtered-data.data', {});
        assert.deepEqual(callback(null, null, 'umap', 'pca', 0, 0, 1, 1, null, null, saved), {});
    """)


def test_scatter_style_updates_only_point_traces_and_preserves_data_and_images():
    run_js(_SCATTER_STYLE_JS, """
        const left = graph('__LEFT_ID__', [
            {type: 'scattergl', mode: 'markers', legendgroup: 'A', x: [1,2]},
            {type: 'scattergl', mode: 'text', x: [1]},
            {type: 'heatmap', z: [[1]]}
        ]);
        left.layout.images = [{source: 'tissue-image'}];
        const right = graph('__RIGHT_ID__', [{type: 'scattergl', mode: 'markers', x: [3]}]);
        const pseudotime = graph('__PSEUDOTIME_ID__', [
            {type: 'scattergl', mode: 'markers', x: [1]},
            {type: 'scatter', mode: 'lines', x: [1,2]}
        ]);
        const before = JSON.stringify(left.data);
        callback(9, 0.4, false, null, null);
        const styles = calls.filter(c => c.kind === 'restyle');
        assert.deepEqual(styles[0].indices, [0]);
        assert.equal(styles[0].update['marker.size'], 9);
        assert.deepEqual(styles[0].update['selected.marker.opacity'], [0.4]);
        assert.deepEqual(styles[1].update['selected.marker.opacity'], [1]);
        assert.deepEqual(styles[2].indices, [0]);
        assert.equal(styles[2].update['selected.marker.opacity'], undefined);
        assert.equal(pseudotime.layout.xaxis.tickfont, undefined);
        assert.equal(left.layout.xaxis.tickfont.color, 'rgba(0,0,0,0)');
        assert.equal(JSON.stringify(left.data), before);
        assert.equal(left.layout.images[0].source, 'tissue-image');
        assert.deepEqual(left.layout.xaxis.range, [0, 10]);
        calls.length = 0;
        trigger('p-axis-toggle.value', true);
        callback(9, 0.4, true, null, null);
        assert.equal(calls.filter(c => c.kind === 'restyle').length, 0);
        assert.equal(left.layout.xaxis.tickfont.color, 'black');
    """)


def test_cell_highlight_handles_click_lasso_clear_and_nonmatching_ids():
    run_js(CELL_HIGHLIGHT_JS, """
        graph('__TARGET_ID__', [
            {type: 'scattergl', ids: ['c3', 'c1']},
            {type: 'scattergl', ids: ['c2']}
        ]);
        trigger('__SOURCE_ID__.clickData', {points: [{id: 'c1'}]});
        let ids = callback(null, null, null, null);
        assert.deepEqual(ids, ['c1']);
        assert.deepEqual(calls.at(-1).update.selectedpoints, [[1], []]);
        window.dash_clientside.callback_context.triggered = [
            {prop_id: '__SOURCE_ID__.clickData', value: null},
            {prop_id: '__SOURCE_ID__.selectedData', value: {points: [{id: 'c2'}, {id: 'c3'}, {id: 'c2'}]}}
        ];
        ids = callback(null, null, null, ids);
        assert.deepEqual(ids, ['c2', 'c3']);
        assert.deepEqual(calls.at(-1).update.selectedpoints, [[0], [0]]);
        trigger('__TARGET_ID__.figure', {});
        assert.deepEqual(callback(null, null, {}, ids), ids);
        trigger('__SOURCE_ID__.selectedData', {points: [{id: 'missing'}]});
        callback(null, null, null, ids);
        assert.deepEqual(calls.at(-1).update.selectedpoints, [[], []]);
        trigger('__SOURCE_ID__.selectedData', null);
        assert.equal(callback(null, null, null, ids), null);
        assert.deepEqual(calls.at(-1).update.selectedpoints, [null, null]);
    """)


def scatter_app(backend="scattergl", unpaired=False):
    app = Dash(__name__)
    if unpaired:
        register_unpaired_multiomics_callbacks(app, None, "p", embedding_render_backend=backend)
    else:
        register_scatter_callbacks(
            app, None, "p", embedding_render_backend=backend,
            initialize_scatter_components=None, apply_relayout=None,
            is_continuous_annotation=None, resolve_plot_adata_from_filter=None,
            search_combined=None, obs_columns=[], obs_columns_lower=[],
            var_names=[], var_names_lower=[], palette_json={}, color_config=[],
            plot_embedding=None, plot_coexpression_embedding=None,
        )
    return app


@pytest.mark.parametrize("unpaired", [False, True])
@pytest.mark.parametrize("backend", ["scattergl", "datashader"])
def test_display_changes_do_not_trigger_server_figures_except_raster_styles(backend, unpaired):
    app = scatter_app(backend, unpaired)
    for plot in ["annotation-scatter", "gene-scatter"]:
        entry = app.callback_map[f"p-{plot}.figure"]
        inputs = {(item["id"], item["property"]) for item in entry["inputs"]}
        states = {(item["id"], item["property"]) for item in entry["state"]}
        assert not any(prop == "relayoutData" for _, prop in inputs)
        assert ("p-axis-toggle", "value") in states
        for control in ["marker-size-slider", "opacity-slider"]:
            assert (f"p-{control}", "value") in (inputs if backend == "datashader" else states)
        assert not any(prop == "figure" for _, prop in states)
    toggles = [entry for key, entry in app.callback_map.items()
               if "controls-container.style" in key or "gene2-container.style" in key]
    assert len(toggles) == 2
    assert all("callback" not in entry for entry in toggles)


def linked_embeddings(*, action=None, backend="scattergl", extra_link=False):
    data = AnnData(np.array([[1.], [3.], [2.]]),
                   obs=pd.DataFrame(index=["c1", "c2", "c3"]),
                   var=pd.DataFrame(index=["G1"]))
    data.obsm["X_umap"] = np.array([[0., 1.], [1., 0.], [2., 2.]])
    views = [view("umap", "source", data="rna", color="G1"),
             view("umap", "target", data="protein", color="G1", render_backend=backend)]
    links = [link("source", "target", action=action)]
    if extra_link:
        views.append(view("umap", "other", data="rna", color="G1"))
        links.append(link("other", "target"))
    return LinkedView({"rna": data, "protein": data[[2, 0, 1]].copy()},
                      views=views, links=links, prefix="demo")


def test_native_highlight_uses_stable_ids_without_a_server_redraw():
    demo = linked_embeddings()
    app = demo.create_app()
    assert "demo-view-target.figure" not in app.callback_map
    assert "callback" not in app.callback_map["demo-view-target-highlight.data"]
    assert "callback" in app.callback_map["demo-state.data"]  # notebook state remains available
    for spec in demo.views:
        fig = demo._render(spec)
        source = demo.store.source(spec.data).data
        trace = fig.data[0]
        assert trace.customdata is None
        assert set(trace.ids) == set(source.obs_names)
        for id, x in zip(trace.ids, trace.x, strict=True):
            assert x == source.obsm["X_umap"][source.obs_names.get_loc(id), 0]


@pytest.mark.parametrize("kwargs", [{"action": "filter"}, {"backend": "datashader"}, {"extra_link": True}])
def test_complex_or_raster_links_retain_general_server_registration(kwargs):
    demo = linked_embeddings(**kwargs)
    # Registration can classify the native target without rendering a raster.
    demo._component_kinds["target"] = "figure"
    app = Dash(__name__)
    demo.register(app)
    assert "callback" in app.callback_map["demo-view-target.figure"]
    assert "demo-view-target-highlight.data" not in app.callback_map
