import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

from guanaco.data.multiomics import MultiOmicsSource, try_build_multiomics_source
from guanaco.pages.matrix.callbacks.unpaired_multiomics_callbacks import (
    _panel_color_options,
    build_unpaired_embedding_figure,
)
from guanaco.pages.matrix.layouts.embedding_layout import (
    generate_embedding_plots,
    initialize_scatter_components,
)
from guanaco.pages.visualizations.layout import generate_visualization_sections


def _paired_mudata():
    cells = pd.Index(["c1", "c2", "c3", "c4"])
    rna = ad.AnnData(
        X=np.ones((4, 3)),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    adt = ad.AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=["p1", "p2"]),
    )
    rna.obsm["X_umap"] = np.arange(8, dtype=float).reshape(4, 2)
    adt.obsm["X_pca"] = np.arange(12, dtype=float).reshape(4, 3)
    adt.obsm["X_lsi"] = sparse.csc_matrix(
        np.arange(12, dtype=float).reshape(4, 3)
    )
    adt.obsm["X_isotypes"] = sparse.csc_matrix(np.ones((4, 4)))
    return mu.MuData({"rna": rna, "adt": adt})


def _unpaired_mudata():
    rna = ad.AnnData(
        X=np.asarray([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]]),
        obs=pd.DataFrame(
            {"cell_type": pd.Categorical(["T", "B", "T"])},
            index=["rna-1", "rna-2", "rna-3"],
        ),
        var=pd.DataFrame(index=["CD3D", "MS4A1"]),
    )
    adt = ad.AnnData(
        X=np.asarray([[1.0, 0.0], [3.0, 2.0]]),
        obs=pd.DataFrame(
            {"population": pd.Categorical(["T", "B"])},
            index=["adt-1", "adt-2"],
        ),
        var=pd.DataFrame(index=["CD3", "CD19"]),
    )
    rna.obsm["X_umap"] = np.asarray([[0, 0], [1, 1], [2, 0]], dtype=float)
    adt.obsm["X_umap"] = np.asarray([[10, 0], [11, 1]], dtype=float)
    return mu.MuData({"rna": rna, "adt": adt})


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


def test_multiomics_embedding_catalog_skips_sparse_feature_matrices():
    source = MultiOmicsSource(_paired_mudata())

    assert source.embedding_names == [
        "RNA · UMAP",
        "ADT · PCA",
        "ADT · LSI",
    ]
    assert "ADT · ISOTYPES" not in source.embedding_names
    assert source.base_adata.obsm["ADT · LSI"].shape == (4, 3)


def test_single_modality_embedding_controls_skip_isotype_matrix():
    adt = _paired_mudata().mod["adt"]

    embeddings, _columns, default, _default_columns = (
        initialize_scatter_components(adt)
    )

    assert embeddings == ["X_pca", "X_lsi"]
    assert default == "X_lsi"


def test_multiomics_feature_set_score_uses_scanpy_and_is_cached():
    mdata = _paired_mudata()
    mdata.mod["rna"].X = np.asarray(
        [
            [0.0, 1.0, 4.0],
            [1.0, 2.0, 2.0],
            [3.0, 1.0, 0.0],
            [5.0, 4.0, 1.0],
        ]
    )
    source = MultiOmicsSource(mdata)

    expected = mdata.mod["rna"].copy()
    sc.tl.score_genes(
        expected,
        ["g1", "g2"],
        score_name="expected",
        random_state=0,
        use_raw=False,
    )
    first = source.score_features(["RNA · g1", "RNA · g2"])
    second = source.score_features(["RNA · g2", "RNA · g1"])

    assert first == pytest.approx(expected.obs["expected"].to_numpy())
    assert second is first


def test_single_feature_uses_its_raw_values_and_is_cached():
    mdata = _paired_mudata()
    raw_values = np.asarray([0.0, 2.0, 5.0, 9.0])
    mdata.mod["rna"].X = np.column_stack(
        [raw_values, np.ones(4), np.zeros(4)]
    )
    source = MultiOmicsSource(mdata)

    first = source.score_features(["RNA · g1"])
    second = source.score_features(["RNA · g1"])

    assert first == pytest.approx(raw_values)
    assert second is first


def test_feature_set_score_rejects_features_from_different_matrices():
    source = MultiOmicsSource(_paired_mudata())

    with pytest.raises(ValueError, match="one matrix"):
        source.score_features(["RNA · g1", "ADT · p1"])


def test_unpaired_multiomics_source_keeps_modality_rows_separate():
    source, reason = try_build_multiomics_source(_unpaired_mudata())

    assert reason is None
    assert source is not None
    assert source.is_paired is False
    assert source.base_adata is None
    assert source.embedding_names == ["RNA · UMAP", "ADT · UMAP"]
    assert source.embedding_context("RNA · UMAP")[2].n_obs == 3
    assert source.embedding_context("ADT · UMAP")[2].n_obs == 2
    with pytest.raises(ValueError, match="cannot be materialized"):
        source.materialize(["RNA · CD3D"])


def test_unpaired_feature_score_stays_in_native_modality_rows():
    source = MultiOmicsSource(_unpaired_mudata())

    rna = source.score_modality_features(["RNA · CD3D"])
    adt = source.score_modality_features(["ADT · CD3"])

    assert rna.tolist() == [0.0, 2.0, 4.0]
    assert adt.tolist() == [1.0, 3.0]
    assert source.score_modality_features(["RNA · CD3D"]) is rna


def test_multiomics_feature_search_ranks_exact_raw_feature_first():
    source = MultiOmicsSource(_unpaired_mudata())

    assert source.search_all_features("CD3", limit=2) == [
        "ADT · CD3",
        "RNA · CD3D",
    ]


def test_paired_joint_view_exposes_virtual_features_to_marker_plots():
    mdata = _paired_mudata()
    cell_type = pd.Categorical(["T", "B", "T", "B"])
    for adata in mdata.mod.values():
        adata.obs["cell_type"] = cell_type
    source = MultiOmicsSource(mdata)

    sections = generate_visualization_sections(
        source.base_adata,
        source.default_features(),
        source.discrete_obs_names,
        "joint",
        multiomics_source=source,
        modality_name="multiomics",
    )

    marker_tabs = _by_id(sections, "joint-marker-tabs")
    assert [tab.value for tab in marker_tabs.children] == [
        "dotplot-tab",
        "heatmap-tab",
        "violin-tab",
    ]


def test_unpaired_scatter_layout_defaults_to_different_modalities():
    source = MultiOmicsSource(_unpaired_mudata())
    layout = generate_embedding_plots(
        None,
        "joint",
        scatter_defaults={
            "embedding_left": "RNA · UMAP",
            "embedding_right": "ADT · UMAP",
        },
        multiomics_source=source,
    )

    assert _by_id(layout, "joint-clustering-dropdown").value == "RNA · UMAP"
    assert _by_id(layout, "joint-right-clustering-dropdown").value == "ADT · UMAP"
    assert _by_id(layout, "joint-annotation-dropdown").value == "cell_type"
    assert _by_id(layout, "joint-scatter-gene-selection").value == "population"
    left_values = {
        option["value"]
        for option in _by_id(layout, "joint-annotation-dropdown").options
    }
    right_values = {
        option["value"]
        for option in _by_id(layout, "joint-scatter-gene-selection").options
    }
    assert "RNA · CD3D" in left_values
    assert "ADT · CD3" in right_values
    assert "ADT · CD3" not in left_values
    assert "RNA · CD3D" not in right_values
    assert _by_id(layout, "joint-controls-container") is not None
    assert _by_id(layout, "joint-annotation-scatter") is not None
    assert _by_id(layout, "joint-gene-scatter") is not None
    assert not any(
        str(getattr(item, "id", "")).startswith("joint-unpaired-")
        for item in _walk(layout)
    )


def test_unpaired_panel_switch_scopes_features_and_rows_to_embedding_modality():
    source = MultiOmicsSource(_unpaired_mudata())

    options, selected = _panel_color_options(
        source,
        "ADT · UMAP",
        None,
        "RNA · CD3D",
    )
    values = {option["value"] for option in options}
    assert selected == "ADT · CD3"
    assert "ADT · CD3" in values
    assert "RNA · CD3D" not in values

    figure = build_unpaired_embedding_figure(
        source,
        "ADT · UMAP",
        "ADT · CD3",
        continuous_colormap="viridis",
        discrete_colormap="cc/glasbey",
        marker_size=5,
        opacity=1,
        render_backend="scattergl",
        color_config=None,
    )
    assert len(figure.data[0].x) == 2
