import numpy as np
import pandas as pd
from anndata import AnnData
from dash import Dash

from guanaco.pages.matrix.callbacks.spatial_relationships_callbacks import (
    register_spatial_relationships_callbacks,
)
from guanaco.pages.matrix.layouts.spatial_relationships_layout import (
    generate_spatial_relationships_layout,
)
from guanaco.pages.matrix.plots.spatial_relationships import (
    default_spatial_pair,
    load_spatial_relationship_payload,
    plot_co_occurrence_group,
    plot_neighborhood_enrichment,
    spatial_pair_from_click,
    spatial_relationship_keys,
)
from guanaco.pages.visualizations.registry import resolve_plot_components


def _adata():
    categories = ["A", "B", "C"]
    adata = AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(
                    ["A", "A", "B", "B", "C", "C"],
                    categories=categories,
                )
            },
            index=[f"spot-{index}" for index in range(6)],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    zscore = np.array(
        [
            [8.0, 4.0, -2.0],
            [4.0, 7.0, 1.0],
            [-2.0, 1.0, 6.0],
        ]
    )
    adata.uns["cluster_nhood_enrichment"] = {
        "zscore": zscore,
        "count": np.arange(9).reshape(3, 3),
    }
    occurrence = np.zeros((3, 3, 3), dtype=float)
    for conditional in range(3):
        for observed in range(3):
            occurrence[conditional, observed] = [
                1 + conditional,
                1 + observed,
                1 + conditional + observed,
            ]
    adata.uns["cluster_co_occurrence"] = {
        "occ": occurrence,
        "interval": np.array([0.0, 10.0, 20.0, 30.0]),
    }
    return adata


def _walk(component):
    if component is None:
        return
    if isinstance(component, (list, tuple)):
        for child in component:
            yield from _walk(child)
        return
    yield component
    children = getattr(component, "children", None)
    if children is not None:
        yield from _walk(children)


def _by_id(component, component_id):
    return next(
        item for item in _walk(component) if getattr(item, "id", None) == component_id
    )


def test_spatial_relationship_capability_requires_both_compatible_results():
    adata = _adata()
    assert spatial_relationship_keys(adata) == ["cluster"]
    assert resolve_plot_components(
        adata,
        ["spatial-relationships"],
    ) == ("spatial-relationships",)

    del adata.uns["cluster_co_occurrence"]
    assert spatial_relationship_keys(adata) == []
    assert resolve_plot_components(adata, ["spatial-relationships"]) == ()


def test_linked_spatial_figures_show_all_curves_and_emphasize_heatmap_pair():
    payload = load_spatial_relationship_payload(_adata(), "cluster")
    assert default_spatial_pair(payload) == ("A", "B")
    assert spatial_pair_from_click(
        {"points": [{"x": "C", "y": "B"}]},
        payload,
    ) == ("B", "C")

    heatmap = plot_neighborhood_enrichment(payload)
    clustered_columns = list(heatmap.data[0].x)
    assert set(clustered_columns) == {"A", "B", "C"}
    assert clustered_columns != ["A", "B", "C"]
    assert heatmap.data[0].z[0][clustered_columns.index("B")] == 4
    assert any(trace.xaxis == "x2" for trace in heatmap.data[1:])

    original_order = plot_neighborhood_enrichment(
        payload,
        cluster_columns=False,
    )
    assert list(original_order.data[0].x) == ["A", "B", "C"]
    assert len(original_order.data) == 1
    assert original_order.layout.yaxis.domain == (0, 1)
    assert "xaxis2" not in original_order.to_dict()["layout"]

    curve = plot_co_occurrence_group(payload, "B", emphasized="C")
    assert [trace.name for trace in curve.data] == ["A", "B", "C"]
    assert list(curve.data[0].x) == [10.0, 20.0, 30.0]
    assert list(curve.data[2].y) == [2.0, 3.0, 4.0]
    assert curve.data[2].line.width > curve.data[0].line.width
    assert all(trace.hoverinfo == "skip" for trace in curve.data)
    assert curve.layout.hovermode is False
    assert "group | B" in curve.layout.title.text


def test_spatial_relationship_layout_and_callbacks_keep_pair_link_local():
    adata = _adata()
    layout = generate_spatial_relationships_layout(adata, "p")
    assert _by_id(layout, "p-spatial-relationships-groupby").value == "cluster"
    assert _by_id(layout, "p-spatial-cluster-columns").value == ["cluster"]
    assert _by_id(layout, "p-spatial-nhood-heatmap") is not None
    assert _by_id(layout, "p-spatial-co-occurrence-curve") is not None

    app = Dash(__name__)
    register_spatial_relationships_callbacks(app, adata, "p")
    heatmap_entry = next(
        value
        for key, value in app.callback_map.items()
        if key.startswith("..p-spatial-nhood-heatmap.figure")
    )
    heatmap_inputs = {item["id"] for item in heatmap_entry["inputs"]}
    assert heatmap_inputs == {
        "p-spatial-relationships-groupby",
        "p-spatial-cluster-columns",
        "p-visualization-workspace-tabs",
        "p-exploratory-tabs",
    }
    heatmap, cleared_click = heatmap_entry["callback"].__wrapped__(
        "cluster",
        ["cluster"],
        "dataset-exploration",
        "spatial-relationships-tab",
    )
    assert heatmap.data[0].type == "heatmap"
    assert cleared_click is None

    curve_entry = next(
        value
        for key, value in app.callback_map.items()
        if key == "p-spatial-co-occurrence-curve.figure"
    )
    curve_inputs = {item["id"] for item in curve_entry["inputs"]}
    assert "p-spatial-nhood-heatmap" in curve_inputs
    assert "p-discrete-color-map-dropdown" in curve_inputs
    assert "p-global-filtered-data" not in curve_inputs
    assert "p-selected-cells-store" not in curve_inputs
    curve = curve_entry["callback"].__wrapped__(
        {"points": [{"x": "C", "y": "B"}]},
        "cluster",
        "Plotly",
        "dataset-exploration",
        "spatial-relationships-tab",
    )
    assert [trace.name for trace in curve.data] == ["A", "B", "C"]
    assert curve.data[2].line.width > curve.data[0].line.width
