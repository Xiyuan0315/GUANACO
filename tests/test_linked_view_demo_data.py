from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pytest

from guanaco.linking import link, linked_view, view
from guanaco.pages.matrix.plots.dotmatrix import plot_dot_matrix


_DEMO_DATA_PATH = (
    Path(__file__).resolve().parents[1] / "examples" / "linked_views" / "demo_data.py"
)
_SPEC = spec_from_file_location("guanaco_linked_demo_data", _DEMO_DATA_PATH)
assert _SPEC is not None and _SPEC.loader is not None
_DEMO_DATA = module_from_spec(_SPEC)
_SPEC.loader.exec_module(_DEMO_DATA)
make_single_cell = _DEMO_DATA.make_single_cell
load_pbmc_cd4_relationship = _DEMO_DATA.load_pbmc_cd4_relationship
load_spatial_relationship_demo = _DEMO_DATA.load_spatial_relationship_demo
make_mudata_demo = _DEMO_DATA.make_mudata_demo
make_pathway_demo = _DEMO_DATA.make_pathway_demo
significant_liana = _DEMO_DATA.significant_liana
_PBMC_H5MU = Path("/Users/xiyuanzhang/Documents/GUANACO_v2/data/mdata_pbmc.h5mu")
_SPATIAL_H5AD = Path(
    "/Users/xiyuanzhang/Documents/GUANACO_v2/data/visium_hne_spatial.h5ad"
)


def test_single_cell_demo_produces_informative_dotplot_fractions():
    adata = make_single_cell(seed=1006)
    genes = ["CD3D", "IL7R", "CCL5", "MS4A1", "NKG7", "LST1"]

    matrix = np.asarray(adata.X)
    assert 0.25 < np.mean(matrix == 0) < 0.8

    figure = plot_dot_matrix(
        adata,
        genes,
        "cell_type",
        selected_labels=list(adata.obs["cell_type"].cat.categories),
        plot_type="dotplot",
    )
    fractions = np.asarray(figure.data[0].customdata, dtype=float)[:, 1]

    assert np.unique(fractions).size > 5
    assert np.any(fractions < 0.5)
    assert not np.allclose(fractions, 1.0)


def test_mudata_fallback_uses_different_rna_and_protein_coordinates():
    mdata, *_ = make_mudata_demo()

    rna_umap = mdata.mod["rna"].obsm["X_umap"]
    protein_pca = mdata.mod["protein"].obsm["X_pca"]
    assert rna_umap.shape == protein_pca.shape
    assert not np.allclose(rna_umap, protein_pca)

    embedding, relationship, *_ = load_pbmc_cd4_relationship(mdata=mdata)
    assert relationship.index.equals(embedding.obs_names)


def test_liana_significance_filter_never_substitutes_non_significant_rows():
    interactions = _DEMO_DATA._fallback_liana()
    with pytest.raises(ValueError, match="No LIANA interactions pass"):
        significant_liana(interactions, pvalue=0.0, rank=0.0)


def test_pathway_demo_has_repeated_gene_keys_for_feature_linking():
    pathways, genes = make_pathway_demo()

    assert pathways.index.is_unique
    assert genes.index.is_unique
    assert genes["gene"].isin(_DEMO_DATA.MARKERS).all()
    assert genes.groupby("pathway")["gene"].size().min() > 1


@pytest.mark.skipif(not _SPATIAL_H5AD.is_file(), reason="Enriched spatial data unavailable")
def test_spatial_relationship_demo_keeps_full_matrix_and_child_keys():
    pairs, spatial, curves, selected = load_spatial_relationship_demo(
        _SPATIAL_H5AD
    )

    assert selected == _SPATIAL_H5AD.resolve()
    assert pairs["pair_id"].nunique() == 15 * 15
    assert pairs.index.is_unique
    assert curves["pair_id"].nunique() > 1
    assert curves.index.is_unique
    assert pairs["cell_id"].isin(spatial.obs_names).all()
    assert spatial.obsm["spatial"].shape[1] >= 2


@pytest.mark.skipif(not _PBMC_H5MU.is_file(), reason="PBMC MuData is unavailable")
def test_real_pbmc_cd4_demo_uses_one_shared_cell_index():
    embedding, relationship, summary, dataset_path, rna_feature, protein_label = (
        load_pbmc_cd4_relationship(_PBMC_H5MU)
    )

    assert dataset_path == _PBMC_H5MU.resolve()
    assert len(relationship) == embedding.n_obs == 8_000
    assert relationship.index.equals(embedding.obs_names)
    assert rna_feature == "CD4"
    assert protein_label.startswith("mean(adt_CD4-")
    assert np.isfinite(relationship["level1_spearman"]).all()
    assert {"T", "Mono"}.issubset(set(summary["level1"]))

    workspace = linked_view(
        {"embedding": embedding, "relationship": relationship},
        views=[
            view("umap", id="umap", data="embedding", color="level1"),
            view(
                "plotly.scatter",
                id="correlation",
                data="relationship",
                x="rna_cd4",
                y="protein_cd4",
                color="level1",
                group="correlation_label",
            ),
        ],
        links=[link("umap", "correlation", by="cell")],
    )

    assert workspace.describe()["links"] == [
        {
            "source": "umap",
            "target": "correlation",
            "by": "cell",
            "action": "highlight",
        }
    ]
