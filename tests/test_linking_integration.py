from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from guanaco.linking.engine import state_for_target
from guanaco.linking.model import link, view
from guanaco.linking.runtime import LinkedView


def _adata(prefix: str = "cell") -> ad.AnnData:
    result = ad.AnnData(
        X=np.asarray(
            [[3.0, 0.0], [2.0, 0.1], [0.0, 4.0], [0.2, 3.0]],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {"cell_type": pd.Categorical(["T", "T", "B", "B"])},
            index=[f"{prefix}-{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["CD4", "MS4A1"]),
    )
    result.obsm["X_umap"] = np.asarray([[0.0, 0.0], [1.0, 0.0], [4.0, 4.0], [5.0, 4.0]])
    return result


def _follow(workspace, source_id, event, payload, target_id):
    source = workspace._view_by_id[source_id]
    members = workspace.registry.get(source.plot).decode_event(
        event, payload, source, workspace.store
    )
    workspace._record_members(source_id, members)
    target = workspace._view_by_id[target_id]
    state = state_for_target(
        workspace._state,
        target_id,
        workspace.links,
        workspace._table_views,
        workspace._row_key_maps,
    )
    rendered = workspace.registry.get(target.plot).render(
        target, workspace.store, state
    )
    return members, state, rendered


def _figure_row_ids(figure) -> set[str]:
    return {
        str(row[0]) for trace in figure.data for row in np.asarray(trace.customdata)
    }


def test_anndata_cell_selection_highlights_violin_as_selected_vs_others() -> None:
    workspace = LinkedView(
        _adata(),
        views=[
            view("umap", id="cells", color="cell_type"),
            view("violin", id="expression", keys=["CD4"], groupby="cell_type"),
        ],
        links=[link("cells", "expression")],
    )

    _, state, figure = _follow(
        workspace,
        "cells",
        "select",
        {"points": [{"customdata": ["cell-0"]}, {"customdata": ["cell-2"]}]},
        "expression",
    )

    assert state.is_highlighted("cell")
    assert {trace.name for trace in figure.data if trace.type == "violin"} == {
        "Selected",
        "Others",
    }


def test_anndata_feature_tile_changes_embedding_feature_context() -> None:
    workspace = LinkedView(
        _adata(),
        views=[
            view(
                "matrixplot",
                id="summary",
                var_names=["CD4", "MS4A1"],
                groupby="cell_type",
            ),
            view("umap", id="expression", color="CD4"),
        ],
        links=[link("summary", "expression", by="feature")],
    )

    members, state, figure = _follow(
        workspace,
        "summary",
        "click",
        {"points": [{"x": "MS4A1", "y": "B"}]},
        "expression",
    )

    assert members.features == ("MS4A1",)
    assert state.features == ("MS4A1",)
    assert figure.layout.title.text == "MS4A1"


def test_external_cell_index_highlights_matching_anndata_cells() -> None:
    adata = _adata()
    scores = pd.DataFrame(
        {"attribute_1": [0, 1, 2, 3], "attribute_2": [3, 2, 1, 0]},
        index=adata.obs_names,
    )
    workspace = LinkedView(
        {"scores": scores, "adata": adata},
        views=[
            view(
                "plotly.scatter",
                id="scores",
                data="scores",
                x="attribute_1",
                y="attribute_2",
            ),
            view("umap", id="cells", data="adata", color="cell_type"),
        ],
        links=[link("scores", "cells", by="cell")],
    )

    _, _, figure = _follow(
        workspace,
        "scores",
        "select",
        {"points": [{"customdata": ["cell-1"]}, {"customdata": ["cell-3"]}]},
        "cells",
    )
    selected = {
        str(adata.obs_names[int(trace.customdata[position])])
        for trace in figure.data
        for position in (trace.selectedpoints or ())
    }

    assert selected == {"cell-1", "cell-3"}


def test_table_key_maps_aggregated_marks_to_anndata_cells() -> None:
    adata = _adata()
    memberships = pd.DataFrame(
        {
            "cell_id": ["cell-0", "cell-1", "cell-2"],
            "source_group": ["T", "T", "B"],
            "target_group": ["B", "B", "Mono"],
            "enrichment": [2.0, 2.0, 1.0],
        },
        index=["pair1:cell-0", "pair1:cell-1", "pair2:cell-2"],
    )
    workspace = LinkedView(
        {"pairs": memberships, "spatial": adata},
        views=[
            view(
                "plotly.heatmap",
                id="neighborhoods",
                data="pairs",
                x="target_group",
                y="source_group",
                value="enrichment",
            ),
            view("umap", id="locations", data="spatial", color="cell_type"),
        ],
        links=[
            link(
                "neighborhoods",
                "locations",
                by="cell",
                key="cell_id",
                action="filter",
            )
        ],
    )

    _, state, figure = _follow(
        workspace,
        "neighborhoods",
        "click",
        {"points": [{"x": "B", "y": "T"}]},
        "locations",
    )

    assert state.cells == ("cell-0", "cell-1")
    assert figure.layout.title.text == "cell_type"
    assert sum(len(trace.x) if trace.x is not None else 0 for trace in figure.data) == 2


@pytest.mark.parametrize("feature", ["MS4A1", "quality_score"])
def test_external_feature_index_changes_anndata_feature_context(feature) -> None:
    adata = _adata()
    adata.obs["quality_score"] = [0.1, 0.7, 0.3, 0.9]
    feature_table = pd.DataFrame(
        {"importance": [2.5, 1.4], "stability": [0.95, 0.92]},
        index=["MS4A1", "quality_score"],
    )
    workspace = LinkedView(
        {"features": feature_table, "adata": adata},
        views=[
            view(
                "plotly.scatter",
                id="feature_summary",
                data="features",
                x="importance",
                y="stability",
            ),
            view("umap", id="embedding", data="adata", color="cell_type"),
        ],
        links=[link("feature_summary", "embedding", by="feature")],
    )

    _, state, figure = _follow(
        workspace,
        "feature_summary",
        "click",
        {"points": [{"customdata": [feature]}]},
        "embedding",
    )

    assert state.features == (feature,)
    assert figure.layout.title.text == feature


def _network_edges(component):
    cytoscape = component.children[1]
    return [element for element in cytoscape.elements if "source" in element["data"]]


def test_network_filters_pair_detail_through_long_table_rows() -> None:
    interactions = pd.DataFrame(
        {
            "sender": ["NK", "NK", "Mono"],
            "receiver": ["B", "B", "T"],
            "ligand": ["CD4", "IL7R", "LST1"],
            "receptor": ["HLA-DRA", "CD127", "S100A8"],
            "magnitude": [0.9, 0.8, 0.7],
        },
        index=["id1", "id2", "id3"],
    )
    workspace = LinkedView(
        interactions,
        views=[
            view("network", id="network", source="sender", target="receiver"),
            view(
                "plotly.scatter",
                id="pairs",
                x="ligand",
                y="receptor",
                color="magnitude",
            ),
        ],
        links=[link("network", "pairs")],
    )
    network = workspace.registry.get("network")
    component = network.render(
        workspace._view_by_id["network"], workspace.store, component_id="network"
    )
    edge = next(
        edge for edge in _network_edges(component) if edge["data"]["weight"] == 2
    )

    _, _, detail = _follow(workspace, "network", "click", edge["data"], "pairs")

    assert _figure_row_ids(detail) == {"id1", "id2"}


def test_aggregated_spatial_tile_filters_exact_location_rows() -> None:
    locations = pd.DataFrame(
        {
            "type_a": ["T", "T", "T", "B"],
            "type_b": ["B", "B", "B", "Mono"],
            "enrichment": [2.4, 2.4, 2.4, 1.2],
            "spatial_x": [10.0, 20.0, 30.0, 80.0],
            "spatial_y": [15.0, 25.0, 35.0, 70.0],
        },
        index=["pair1:spot1", "pair1:spot2", "pair1:spot3", "pair2:spot4"],
    )
    workspace = LinkedView(
        locations,
        views=[
            view(
                "plotly.heatmap",
                id="enrichment",
                x="type_a",
                y="type_b",
                value="enrichment",
            ),
            view("plotly.scatter", id="locations", x="spatial_x", y="spatial_y"),
        ],
        links=[link("enrichment", "locations")],
    )

    _, _, detail = _follow(
        workspace,
        "enrichment",
        "click",
        {"points": [{"x": "T", "y": "B"}]},
        "locations",
    )

    assert _figure_row_ids(detail) == {
        "pair1:spot1",
        "pair1:spot2",
        "pair1:spot3",
    }


def test_mudata_modalities_link_directly_by_shared_cell_ids() -> None:
    mdata = SimpleNamespace(mod={"rna": _adata("shared"), "protein": _adata("shared")})
    workspace = LinkedView(
        mdata,
        views=[
            view("umap", id="rna", data="rna", color="CD4"),
            view("umap", id="protein", data="protein", color="MS4A1"),
        ],
        links=[link("rna", "protein")],
    )

    assert workspace.links[0].selection_by == "cell"
    assert workspace.links[0].resolved_action == "highlight"
