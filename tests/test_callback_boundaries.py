import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from dash import Dash, html
from dash.exceptions import PreventUpdate

import guanaco.pages.visualizations.callbacks as visualization_callbacks
import guanaco.pages.visualizations.plots.igv.callbacks as igv_callbacks
from guanaco.pages.matrix.callbacks.atac_browser_callbacks import (
    register_atac_browser_callbacks,
)
from guanaco.pages.matrix.callbacks.paga_callbacks import register_paga_callbacks
from guanaco.pages.matrix.callbacks.stacked_bar_callbacks import (
    register_stacked_bar_callbacks,
)
from guanaco.pages.matrix.callbacks.violin_callbacks import register_violin_callbacks
from guanaco.pages.visualizations.plots.igv.callbacks import register_igv_callbacks
from guanaco.utils.obs_utils import SELECTION_GROUP


FORBIDDEN_EXPLORATORY_INPUTS = {
    "p-single-cell-genes-selection",
    "p-single-cell-annotation-dropdown",
    "p-single-cell-label-selection",
    "p-selected-cells-store",
}


def _adata(*, peaks=False):
    var_names = (
        ["chr1:1000-1500", "chr1:2000-2500", "chr1:3000-3500"]
        if peaks
        else ["G1", "G2", "G3"]
    )
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["A", "B", "A", "B"]),
            "condition": pd.Categorical(["control", "control", "case", "case"]),
        },
        index=[f"cell-{index}" for index in range(4)],
    )
    return AnnData(
        X=np.ones((4, 3), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=var_names),
    )


def _inputs_for(app, output_fragment):
    callback = next(
        value for key, value in app.callback_map.items() if output_fragment in key
    )
    return {item["id"] for item in callback["inputs"]}


def _assert_exploratory_boundary(inputs):
    assert not inputs & FORBIDDEN_EXPLORATORY_INPUTS
    assert "p-exploratory-tabs" in inputs


def test_stacked_bar_owns_its_grouping_inputs():
    app = Dash(__name__)
    source = _adata()
    rendered = {}

    def resolve_filtered(filtered_data):
        indices = (filtered_data or {}).get("cell_indices")
        return source if indices is None else source[indices]

    def calculate_alr(source, *_args, **_kwargs):
        rendered["analysis_n_obs"] = source.n_obs

    def plot_stacked_bar(**kwargs):
        rendered["bar_n_obs"] = kwargs["adata"].n_obs
        rendered["bar_group_values"] = kwargs.get("group_values")

    register_stacked_bar_callbacks(
        app,
        source,
        "p",
        calculate_alr_welch=calculate_alr,
        plot_composition_differential_abundance=lambda *_args, **_kwargs: None,
        plot_composition_hierarchy=lambda **_kwargs: None,
        plot_stacked_bar=plot_stacked_bar,
        palette_json={},
        color_config=["red", "blue"],
        resolve_plot_adata_from_filter=resolve_filtered,
        hash_list_signature=lambda values: tuple(values or []),
    )

    inputs = _inputs_for(app, "p-stacked-bar-plot.figure")
    _assert_exploratory_boundary(inputs)
    assert "p-stacked-bar-x-group" in inputs
    assert "p-composition-view" in inputs
    assert "p-composition-swap-axes" in inputs
    assert "p-visualization-workspace-tabs" in inputs
    assert "p-global-filtered-data" in inputs
    assert {"p-selected-cells-hash", "p-selection-group-hash"} <= inputs
    assert "p-stacked-bar-x-labels" not in inputs

    analysis_inputs = _inputs_for(app, "p-composition-da-plot.figure")
    _assert_exploratory_boundary(analysis_inputs)
    assert {
        "p-visualization-workspace-tabs",
        "p-stacked-bar-x-group",
        "p-stacked-bar-stack-by",
        "p-composition-sample-unit",
        "p-composition-alr-reference",
        "p-global-filtered-data",
    } <= analysis_inputs
    assert "p-composition-da-panel" not in analysis_inputs

    collapse_callback = next(
        value
        for key, value in app.callback_map.items()
        if "p-composition-da-collapse.style" in key
    )["callback"].__wrapped__
    collapse_style, toggle_label = collapse_callback(
        1, {"display": "none", "width": "100%"}
    )
    assert collapse_style == {"display": "block", "width": "100%"}
    assert toggle_label == "▾ Differential abundance"

    analysis_callback = next(
        value
        for key, value in app.callback_map.items()
        if "p-composition-da-plot.figure" in key
    )["callback"].__wrapped__
    figure, rendered_key, results_style = analysis_callback(
        "bars",
        "dataset-exploration",
        "stacked-bar-tab",
        "condition",
        "cell_type",
        "sample",
        "A",
        [],
        {"cell_indices": [0, 1], "n_cells": 2},
        None,
        None,
        None,
        None,
        None,
        None,
    )
    assert figure is None
    assert rendered_key
    assert results_style == {
        "display": "block",
        "height": "360px",
        "minHeight": "360px",
        "width": "100%",
    }
    assert rendered["analysis_n_obs"] == 2

    bar_callback = next(
        value
        for key, value in app.callback_map.items()
        if "p-stacked-bar-plot.figure" in key
    )["callback"].__wrapped__
    _figure, bar_key = bar_callback(
        "bars",
        [],
        "prop",
        "Plotly",
        "dataset-exploration",
        "stacked-bar-tab",
        "condition",
        "cell_type",
        [],
        {"cell_indices": [0, 1], "n_cells": 2},
        None,
        None,
        None,
        None,
        None,
        None,
    )
    assert bar_key
    assert rendered["bar_n_obs"] == 2
    assert "p-discrete-color-map-dropdown" not in analysis_inputs

    _figure, lasso_key = bar_callback(
        "bars",
        [],
        "prop",
        "Plotly",
        "dataset-exploration",
        "stacked-bar-tab",
        "condition",
        "cell_type",
        [],
        {"cell_indices": None, "n_cells": 4},
        {"len": 1, "hash": "filtered"},
        None,
        ["cell-1"],
        None,
        None,
        None,
    )
    assert lasso_key
    assert rendered["bar_n_obs"] == 1
    assert rendered["bar_group_values"] is None

    _figure, highlight_key = bar_callback(
        "bars",
        [],
        "prop",
        "Plotly",
        "dataset-exploration",
        "stacked-bar-tab",
        "condition",
        SELECTION_GROUP,
        [],
        {"cell_indices": None, "n_cells": 4},
        None,
        {"selected": {"len": 1, "hash": "highlighted"}},
        None,
        {
            "selected_cells": ["cell-1"],
            "universe_cells": None,
        },
        None,
        None,
    )
    assert highlight_key
    assert rendered["bar_n_obs"] == 4
    assert rendered["bar_group_values"][SELECTION_GROUP].tolist() == [
        "Others",
        "Selected",
        "Others",
        "Others",
    ]


def test_peak_browser_owns_its_grouping_inputs():
    app = Dash(__name__)
    register_atac_browser_callbacks(app, _adata(peaks=True), "p")

    inputs = _inputs_for(app, "p-atac-browser-graph.figure")
    _assert_exploratory_boundary(inputs)
    assert {"p-atac-browser-groupby", "p-atac-browser-labels"} <= inputs


def test_paga_does_not_inherit_scatter_selection():
    app = Dash(__name__)
    register_paga_callbacks(
        app,
        _adata(),
        "p",
        build_paga_cytoscape=lambda *_args, **_kwargs: None,
        color_config=["red", "blue"],
    )

    _assert_exploratory_boundary(_inputs_for(app, "p-paga.children"))


def test_comparative_violin_uses_local_layer_and_lasso_selection():
    adata = _adata()
    app = Dash(__name__)
    rendered = {}

    def plot_violin2(source, **kwargs):
        rendered["n_obs"] = source.n_obs
        rendered["meta1"] = kwargs["meta1"]
        rendered["group_values"] = kwargs["group_values"]

    register_violin_callbacks(
        app,
        adata,
        "p",
        filter_data=lambda source, *_args: source,
        plot_violin1=lambda *_args, **_kwargs: None,
        plot_violin2_new=plot_violin2,
        palette_json={},
        var_names=list(adata.var_names),
        var_names_lower=[name.lower() for name in adata.var_names],
        color_config=["red", "blue"],
        resolve_plot_adata_from_filter=lambda filtered: (
            adata
            if not (filtered or {}).get("cell_indices")
            else adata[(filtered or {})["cell_indices"]]
        ),
    )

    inputs = _inputs_for(app, "p-violin-plot2.figure")
    _assert_exploratory_boundary(inputs)
    assert "p-violin2-data-layer" in inputs
    assert "p-data-layer" not in inputs
    assert {
        "p-global-filtered-data",
        "p-selected-cells-hash",
        "p-selection-group-hash",
    } <= inputs

    callback = next(
        value
        for key, value in app.callback_map.items()
        if key == "p-violin-plot2.figure"
    )["callback"].__wrapped__
    callback(
        "G1",
        "cell_type",
        "none",
        "mode1",
        "none",
        [],
        "X",
        "Plotly",
        {"cell_indices": None, "n_cells": 4},
        {"len": 1, "hash": "filtered"},
        None,
        "split-violin-tab",
        ["cell-1"],
        None,
    )
    assert rendered["n_obs"] == 1
    assert rendered["group_values"] is None

    callback(
        "G1",
        SELECTION_GROUP,
        "none",
        "mode1",
        "none",
        [],
        "X",
        "Plotly",
        {"cell_indices": None, "n_cells": 4},
        None,
        {"selected": {"len": 1, "hash": "highlighted"}},
        "split-violin-tab",
        None,
        {
            "selected_cells": ["cell-1"],
            "universe_cells": None,
        },
    )
    assert rendered["n_obs"] == 4
    assert rendered["meta1"] == SELECTION_GROUP
    assert rendered["group_values"][SELECTION_GROUP].tolist() == [
        "Others",
        "Selected",
        "Others",
        "Others",
    ]


def test_igv_uses_the_common_exploratory_workspace_guard():
    app = Dash(__name__)
    rendered = {}

    def igv_factory(**kwargs):
        rendered.update(kwargs)
        return html.Div("IGV")

    register_igv_callbacks(
        app,
        {"sample": [{"name": "accessibility"}]},
        {"label": "hg38"},
        "p",
        discrete_color_prefix="p",
        igv_factory=igv_factory,
    )

    inputs = _inputs_for(app, "p-igv-container.children")
    _assert_exploratory_boundary(inputs)
    assert "p-visualization-workspace-tabs" in inputs

    callback = next(
        value
        for key, value in app.callback_map.items()
        if "p-igv-container.children" in key
    )["callback"].__wrapped__
    with pytest.raises(PreventUpdate):
        callback("sample", "feature-analysis", "igv-tab", "plotly/Plotly")

    result = callback(
        "sample",
        "dataset-exploration",
        "igv-tab",
        "plotly/Plotly",
    )
    assert result.children.children == "IGV"
    assert rendered["genome"] == "hg38"
    assert rendered["tracks"][0]["name"] == "accessibility"


def test_one_visualization_registrar_routes_matrix_and_igv_callbacks(monkeypatch):
    registered = []

    def matrix_registrar(_app, _adata, prefix, **kwargs):
        registered.append(("matrix", prefix, set(kwargs["enabled_components"])))

    def igv_registrar(_app, _tracks, _reference, prefix, **_kwargs):
        registered.append(("igv", prefix))

    monkeypatch.setattr(
        visualization_callbacks,
        "matrix_callbacks",
        matrix_registrar,
    )
    monkeypatch.setattr(igv_callbacks, "register_igv_callbacks", igv_registrar)

    enabled = visualization_callbacks.register_visualization_callbacks(
        object(),
        _adata(),
        "p",
        optional_plot_components=["dotplot"],
        genome_tracks={"sample": [{"name": "accessibility"}]},
        ref_track={"label": "hg38"},
    )

    assert enabled == ("dotplot", "igv")
    assert registered == [
        ("matrix", "p", {"dotplot"}),
        ("igv", "p"),
    ]
