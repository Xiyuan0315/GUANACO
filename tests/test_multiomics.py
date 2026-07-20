import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

from guanaco.data.multiomics import MultiOmicsSource
from guanaco.pages.matrix.layouts.embedding_layout import (
    initialize_scatter_components,
)


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
