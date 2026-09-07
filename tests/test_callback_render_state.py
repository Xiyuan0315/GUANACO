"""Render keys must gate updates without uploading the previously drawn figure."""

import inspect
from unittest.mock import Mock

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest
from dash import Dash, no_update

from guanaco.pages.matrix.callbacks.dotplot_callbacks import register_dotplot_callbacks
from guanaco.pages.matrix.callbacks.heatmap_callbacks import register_heatmap_callbacks
from guanaco.pages.matrix.callbacks.pseudotime_callbacks import register_pseudotime_callbacks
from guanaco.pages.matrix.callbacks.stacked_bar_callbacks import register_stacked_bar_callbacks
from guanaco.pages.matrix.callbacks.violin_callbacks import (
    register_marker_violin_callbacks,
    register_ridge_callbacks,
)
from guanaco.pages.matrix.callbacks.volcano_callbacks import register_volcano_callbacks


def _register(registrar):
    app = Dash(__name__)
    source = ad.AnnData(
        X=np.ones((2, 2), dtype=np.float32),
        obs=pd.DataFrame({"group": ["A", "B"]}, index=["a", "b"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    plot = Mock(side_effect=lambda *args, **kwargs: go.Figure())
    dependencies = {
        "filter_data": lambda source, *args: source,
        "palette_json": {},
        "color_config": ["red", "blue"],
        "make_cache_key": lambda kind, source, **kwargs: repr((kind, kwargs)),
        "hash_list_signature": lambda values: tuple(values or []),
        "cached_figure_get": lambda key: None,
        "cached_figure_set": lambda key, figure: None,
        "calculate_alr_welch": Mock(),
        "resolve_plot_adata_from_filter": lambda filtered: source,
    }
    kwargs = {}
    for name, parameter in inspect.signature(registrar).parameters.items():
        if parameter.kind != inspect.Parameter.KEYWORD_ONLY:
            continue
        if name.startswith("plot_"):
            kwargs[name] = plot
        elif name in dependencies:
            kwargs[name] = dependencies[name]
    registrar(app, source, "p", **kwargs)
    return app, plot


@pytest.mark.parametrize("registrar, outputs", [
    (register_dotplot_callbacks, ["p-dotplot.figure"]),
    (register_heatmap_callbacks, ["p-heatmap.figure"]),
    (register_pseudotime_callbacks, ["p-pseudotime-plot.figure"]),
    (register_marker_violin_callbacks, ["p-violin-plot1.figure"]),
    (register_ridge_callbacks, ["p-ridge-plot.figure"]),
    (register_stacked_bar_callbacks, [
        "p-stacked-bar-plot.figure", "p-composition-da-plot.figure",
    ]),
    (register_volcano_callbacks, ["p-volcano-plot.figure"]),
])
def test_render_callbacks_do_not_request_full_figure(registrar, outputs):
    app, _ = _register(registrar)
    for output in outputs:
        entry = next(v for k, v in app.callback_map.items() if output in k)
        assert all(s["property"] != "figure" for s in entry["state"])
        assert any("rendered-key" in s["id"] for s in entry["state"])


def test_render_key_preserves_initial_render_tab_switch_and_changed_request():
    app, plot = _register(register_dotplot_callbacks)
    entry = next(v for k, v in app.callback_map.items() if "p-dotplot.figure" in k)
    callback = entry["callback"].__wrapped__
    args = dict.fromkeys(inspect.signature(callback).parameters)
    args.update(
        selected_genes=["G1"], selected_annotation="group", active_tab="dotplot-tab",
    )
    figure, key = callback(**args)
    assert isinstance(figure, go.Figure)
    assert key
    assert plot.call_count == 1

    args["rendered_key"] = key
    assert callback(**args) == (no_update, no_update)
    args["active_tab"] = "heatmap-tab"
    assert callback(**args) == (no_update, no_update)
    args["active_tab"] = "dotplot-tab"
    assert callback(**args) == (no_update, no_update)
    assert plot.call_count == 1

    args["selected_genes"] = ["G2"]
    _, next_key = callback(**args)
    assert next_key != key
    assert plot.call_count == 2


def test_marker_violin_renders_once_without_browser_figure_cache():
    app, plot = _register(register_marker_violin_callbacks)
    assert not any("plot-cache-store" in str(entry) for entry in app.callback_map.values())
    entry = next(v for k, v in app.callback_map.items() if "p-violin-plot1.figure" in k)
    callback = entry["callback"].__wrapped__
    args = dict.fromkeys(inspect.signature(callback).parameters)
    args.update(selected_genes=["G1"], selected_annotation="group", active_tab="violin-tab")
    figure, key, style = callback(**args)
    assert isinstance(figure, go.Figure)
    assert style["height"] == "400px"
    assert plot.call_count == 1
    args["rendered_key"] = key
    assert callback(**args) == (no_update,) * 3
    args["active_tab"] = "dotplot-tab"
    assert callback(**args) == (no_update,) * 3
    args.update(active_tab="violin-tab", selected_genes=["G2"])
    assert callback(**args)[1] != key
    assert plot.call_count == 2


def test_pseudotime_shared_style_controls_do_not_recompute_data():
    app, _ = _register(register_pseudotime_callbacks)
    entry = next(v for k, v in app.callback_map.items() if "p-pseudotime-plot.figure" in k)
    inputs = {item["id"] for item in entry["inputs"]}
    states = {item["id"] for item in entry["state"]}
    for name in ["p-marker-size-slider", "p-opacity-slider"]:
        assert name not in inputs
        assert name in states
