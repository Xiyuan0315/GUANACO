from __future__ import annotations

import builtins
import inspect

import anndata as ad
import guanaco as gc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from guanaco.linking.model import LinkSpec, ViewSpec, view


@pytest.fixture
def adata() -> ad.AnnData:
    result = ad.AnnData(
        X=np.asarray([[1.0, 0.0], [0.0, 2.0], [3.0, 1.0]]),
        obs=pd.DataFrame(
            {"cell_type": pd.Categorical(["T", "B", "T"])},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=["CD4", "MS4A1"]),
    )
    result.obsm["X_umap"] = np.asarray([[0, 0], [1, 2], [2, 0]])
    result.obsm["X_fa"] = np.asarray([[0, 0], [3, 4], [5, 0]])
    return result


def test_scanpy_style_plot_helpers_build_linked_view_specs() -> None:
    specs = [
        gc.pl.umap(id="cells", color="cell_type", height="360px"),
        gc.pl.embedding(id="factor", basis="fa", color="cell_type"),
        gc.pl.spatial(id="tissue", color="cell_type"),
        gc.pl.violin(id="violin", keys=["CD4"], groupby="cell_type"),
        gc.pl.heatmap(id="heat", var_names=["CD4"], groupby="cell_type"),
        gc.pl.dotplot(id="dot", var_names=["CD4"], groupby="cell_type"),
        gc.pl.matrixplot(id="matrix", var_names=["CD4"], groupby="cell_type"),
        gc.pl.volcano(id="de", group="comparison"),
    ]

    assert all(isinstance(spec, ViewSpec) for spec in specs)
    assert specs[0] == view("umap", id="cells", color="cell_type", height="360px")
    assert specs[1].options["basis"] == "fa"
    assert specs[2].plot == "spatial"
    assert specs[3].options["keys"] == ["CD4"]


def test_umap_standalone_and_linked_forms_use_the_same_plot_contract(adata) -> None:
    standalone = gc.pl.umap(adata, color="cell_type", return_fig=True)
    workspace = gc.pl.linked_view(
        adata,
        views=[gc.pl.umap(id="cells", color="cell_type")],
        links=[],
    )
    linked = workspace._render(workspace.views[0])

    assert isinstance(standalone, go.Figure)
    assert isinstance(linked, go.Figure)
    assert [trace.type for trace in standalone.data] == [
        trace.type for trace in linked.data
    ]
    assert all(
        np.array_equal(left.customdata, right.customdata)
        for left, right in zip(standalone.data, linked.data, strict=True)
    )


def test_named_dot_and_matrix_plots_keep_their_native_mark_types(adata) -> None:
    dot = gc.pl.dotplot(adata, ["CD4", "MS4A1"], "cell_type", return_fig=True)
    matrix = gc.pl.matrixplot(adata, ["CD4", "MS4A1"], "cell_type", return_fig=True)

    assert {trace.type for trace in dot.data} == {"scatter"}
    assert [trace.type for trace in matrix.data] == ["heatmap"]


def test_standalone_plotting_does_not_import_linking(monkeypatch, adata) -> None:
    original_import = builtins.__import__

    def reject_linking(name, *args, **kwargs):
        if name.startswith("guanaco.linking"):
            raise AssertionError(f"standalone plotting imported {name}")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", reject_linking)

    figure = gc.pl.embedding(adata, basis="fa", color="cell_type", return_fig=True)

    assert isinstance(figure, go.Figure)
    assert max(float(value) for trace in figure.data for value in trace.x) == 5.0


def test_public_linking_functions_have_small_explicit_signatures() -> None:
    assert tuple(inspect.signature(gc.pl.view).parameters) == (
        "plot",
        "id",
        "data",
        "title",
        "options",
    )
    assert tuple(inspect.signature(gc.pl.link).parameters) == (
        "source",
        "target",
        "by",
        "action",
        "key",
    )
    assert tuple(inspect.signature(gc.pl.linked_view).parameters) == (
        "data",
        "views",
        "links",
        "layout",
        "registry",
        "prefix",
        "title",
    )
    assert gc.pl.link("cells", "expression") == LinkSpec("cells", "expression")


def test_registry_contract_describes_identity_capabilities() -> None:
    contract = gc.pl.linked_plot_types("umap")

    assert contract["emits"] == ("cell",)
    assert contract["accepts"] == ("cell", "feature")
    assert "actions" not in contract
    assert "target_parameters" not in contract
