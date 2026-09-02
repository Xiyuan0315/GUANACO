from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from guanaco.linking.data import DataStore
from guanaco.linking.model import MarkMembers, ViewState, view
from guanaco.linking.native_adapters import (
    CoexpressionAdapter,
    CompositionAdapter,
    EmbeddingAdapter,
    FeatureDistributionAdapter,
    FeatureGroupMatrixAdapter,
    HeatmapAdapter,
    PseudotimeAdapter,
    ViolinAdapter,
)


@pytest.fixture
def adata() -> ad.AnnData:
    result = ad.AnnData(
        X=np.asarray(
            [
                [0.0, 1.0, 2.0],
                [1.0, 0.0, 2.0],
                [2.0, 3.0, 0.0],
                [0.0, 4.0, 1.0],
                [3.0, 0.0, 1.0],
                [0.0, 2.0, 0.0],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(["T", "T", "B", "B", "Mono", "Mono"]),
                "condition": pd.Categorical(["ctrl", "stim"] * 3),
                "pseudotime": np.linspace(0.0, 1.0, 6),
            },
            index=[f"cell-{index}" for index in range(6)],
        ),
        var=pd.DataFrame(index=["CD4", "MS4A1", "LYZ"]),
    )
    result.obsm["X_umap"] = np.asarray(
        [[0, 0], [1, 0], [4, 4], [5, 4], [8, 0], [9, 0]], dtype=np.float32
    )
    result.obsm["X_fa"] = result.obsm["X_umap"] * 2
    result.obsm["spatial"] = result.obsm["X_umap"] * 10
    return result


@pytest.fixture
def store(adata: ad.AnnData) -> DataStore:
    return DataStore.from_data(adata)


def test_embedding_maps_cell_positions_for_highlight_and_filter(store: DataStore):
    spec = view("umap", id="cells", color="cell_type")
    adapter = EmbeddingAdapter()
    adapter.validate(spec, store)

    base = adapter.render(spec, store, component_id="ignored-by-plotly-runtime")
    assert isinstance(base, go.Figure)
    assert sum(len(trace.x) for trace in base.data) == 6

    highlighted = adapter.render(
        spec,
        store,
        ViewState(cells=("cell-1", "cell-4"), highlight_axes=("cell",)),
    )
    selected_positions = {
        int(trace.customdata[index])
        for trace in highlighted.data
        for index in (trace.selectedpoints or ())
    }
    assert selected_positions == {1, 4}
    assert highlighted.layout.xaxis.range == base.layout.xaxis.range
    assert highlighted.layout.yaxis.range == base.layout.yaxis.range

    filtered = adapter.render(
        spec,
        store,
        ViewState(cells=("cell-1", "cell-4")),
    )
    assert sum(len(trace.x) for trace in filtered.data) == 2
    assert filtered.layout.xaxis.range == base.layout.xaxis.range
    assert filtered.layout.yaxis.range == base.layout.yaxis.range


def test_embedding_accepts_any_basis_and_feature_context(store: DataStore):
    spec = view("umap", id="cells", color="cell_type", basis="fa")
    adapter = EmbeddingAdapter()
    adapter.validate(spec, store)
    figure = adapter.render(spec, store, ViewState(features=("CD4",)))
    assert figure.layout.title.text == "CD4"
    assert max(float(value) for trace in figure.data for value in trace.x) == 18.0


def test_embedding_decodes_stable_point_ids_to_cells(store: DataStore):
    spec = view("umap", id="cells", color="cell_type")
    adapter = EmbeddingAdapter()
    direct = adapter.decode_event(
        "select",
        {"points": [{"customdata": 2}, {"id": "cell-3"}]},
        spec,
        store,
    )
    assert direct == MarkMembers(cells=("cell-2", "cell-3"))


def test_datashader_target_uses_points_when_cells_must_be_highlighted(store):
    spec = view("umap", id="cells", color="cell_type", render_backend="datashader")
    figure = EmbeddingAdapter().render(
        spec,
        store,
        ViewState(cells=("cell-1",), highlight_axes=("cell",)),
    )

    assert {trace.type for trace in figure.data} == {"scattergl"}
    assert sum(len(trace.selectedpoints or ()) for trace in figure.data) == 1


def test_coexpression_feature_domain_and_cell_highlight_are_exact(store):
    spec = view("coexpression", id="coexpression", gene1="CD4", gene2="MS4A1")
    adapter = CoexpressionAdapter()
    features = adapter.identity_ids("feature", spec, store)
    figure = adapter.render(
        spec,
        store,
        ViewState(cells=("cell-4",), highlight_axes=("cell",)),
    )

    assert set(features) == {"CD4", "MS4A1", "LYZ"}
    assert sum(len(trace.selectedpoints or ()) for trace in figure.data) == 1


def test_violin_highlight_is_selected_vs_others_with_shared_bandwidth(store: DataStore):
    spec = view("violin", id="expression", keys=["CD4"], groupby="cell_type")
    adapter = ViolinAdapter()
    figure = adapter.render(
        spec,
        store,
        ViewState(
            cells=("cell-0", "cell-2", "cell-4"),
            highlight_axes=("cell",),
        ),
    )
    violins = [trace for trace in figure.data if trace.type == "violin"]
    assert {trace.name for trace in violins} == {"Selected", "Others"}
    bandwidths = {trace.bandwidth for trace in violins}
    assert len(bandwidths) == 1
    assert next(iter(bandwidths)) > 0
    assert adapter.decode_event(
        "click", {"points": [{"customdata": ["CD4", "Selected"]}]}, spec, store
    ) == MarkMembers(features=("CD4",))


def test_link_source_violin_exposes_sampled_cells_for_lasso(store: DataStore):
    spec = view(
        "violin",
        id="expression",
        keys=["CD4"],
        groupby="cell_type",
        _link_source=True,
    )
    adapter = ViolinAdapter()

    figure = adapter.render(spec, store)
    violins = [trace for trace in figure.data if trace.type == "violin"]
    points = [
        list(row)
        for trace in violins
        for row in np.asarray(trace.customdata, dtype=object)
    ]
    selected = points[:2]
    members = adapter.decode_event(
        "select",
        {"points": [{"customdata": row} for row in selected]},
        spec,
        store,
    )

    assert all(trace.points == "all" for trace in violins)
    assert all(len(row) == 3 and row[2].startswith("cell-") for row in points)
    assert members is not None
    assert members.cells == tuple(dict.fromkeys(row[2] for row in selected))
    assert members.features == ("CD4",)


def test_ridge_accepts_key_and_feature_context(store: DataStore):
    spec = view("ridge", id="ridge", key="CD4", groupby="cell_type")
    adapter = FeatureDistributionAdapter("ridge")
    adapter.validate(spec, store)
    figure = adapter.render(spec, store, ViewState(features=("MS4A1",)))
    assert "MS4A1" in figure.layout.title.text


def test_dotplot_filter_recomputes_real_fraction_and_mark_members(store: DataStore):
    spec = view(
        "dotplot",
        id="dot",
        var_names=["CD4", "MS4A1"],
        groupby="cell_type",
    )
    adapter = FeatureGroupMatrixAdapter("dotplot")
    figure = adapter.render(
        spec,
        store,
        ViewState(cells=("cell-0", "cell-1")),
    )
    trace = figure.data[0]
    index = next(
        index
        for index, (feature, group) in enumerate(zip(trace.x, trace.y, strict=True))
        if feature == "CD4" and group == "T"
    )
    assert float(np.asarray(trace.customdata)[index, 1]) == pytest.approx(0.5)

    members = adapter.decode_event(
        "click",
        {"points": [{"x": "CD4", "y": "T"}]},
        spec,
        store,
    )
    assert members == MarkMembers(cells=("cell-0", "cell-1"), features=("CD4",))


def test_dotplot_highlight_compares_selected_and_others(store: DataStore):
    spec = view("dotplot", id="dot", var_names=["CD4"], groupby="cell_type")
    figure = FeatureGroupMatrixAdapter("dotplot").render(
        spec,
        store,
        ViewState(cells=("cell-0", "cell-1"), highlight_axes=("cell",)),
    )
    trace = figure.data[0]
    assert set(trace.y) == {"Selected", "Others"}
    selected = list(trace.y).index("Selected")
    assert float(np.asarray(trace.customdata)[selected, 1]) == pytest.approx(0.5)


def test_matrixplot_click_uses_heatmap_axes_without_customdata(store: DataStore):
    spec = view(
        "matrixplot",
        id="matrix",
        var_names=["CD4", "MS4A1"],
        groupby="cell_type",
    )
    adapter = FeatureGroupMatrixAdapter("matrixplot")
    adapter.render(spec, store)

    members = adapter.decode_event(
        "click",
        {"points": [{"x": "MS4A1", "y": "B", "z": 3.5}]},
        spec,
        store,
    )

    assert members == MarkMembers(
        cells=("cell-2", "cell-3"),
        features=("MS4A1",),
    )


def test_heatmap_feature_context_and_click_identity(store: DataStore):
    spec = view("heatmap", id="heat", var_names=["CD4"], groupby="cell_type")
    adapter = HeatmapAdapter()
    figure = adapter.render(spec, store, ViewState(features=("MS4A1",)))
    assert set(figure.data[0].y) == {"MS4A1"}
    assert adapter.decode_event(
        "click",
        {"points": [{"x": 2, "y": "MS4A1"}]},
        spec,
        store,
    ) == MarkMembers(features=("MS4A1",))


def test_composition_aggregate_expands_to_atomic_cells(store: DataStore):
    spec = view(
        "composition",
        id="composition",
        x="condition",
        color="cell_type",
    )
    adapter = CompositionAdapter()
    members = adapter.decode_event(
        "click",
        {"points": [{"customdata": ["ctrl", "T", 1]}]},
        spec,
        store,
    )
    assert members == MarkMembers(cells=("cell-0",))


def test_pseudotime_points_carry_both_cell_and_feature(store: DataStore):
    spec = view(
        "pseudotime",
        id="trajectory",
        genes=["CD4"],
        pseudotime_key="pseudotime",
        groupby="cell_type",
        min_expr=0,
    )
    adapter = PseudotimeAdapter()
    figure = adapter.render(spec, store)
    point = next(
        row
        for trace in figure.data
        if trace.customdata is not None
        for row in np.asarray(trace.customdata, dtype=object)
    )
    members = adapter.decode_event(
        "click", {"points": [{"customdata": point.tolist()}]}, spec, store
    )
    assert members is not None
    assert members.cells == (str(point[0]),)
    assert members.features == ("CD4",)
