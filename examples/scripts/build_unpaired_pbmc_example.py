"""Build a compact, genuinely unpaired PBMC RNA+ATAC MuData example.

The two modalities come from separate public 10x Genomics experiments:

- RNA: PBMC 3k, Cell Ranger 1.1.0
- ATAC: PBMC 5k Next GEM, Cell Ranger ATAC 1.2.0

The script expects the two extracted matrix directories as arguments. It keeps
the modality matrices and embeddings independent and prefixes every barcode so
that no cell identity can accidentally match across modalities.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD


RNA_SOURCE = (
    "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
    "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)
ATAC_SOURCE = (
    "https://cf.10xgenomics.com/samples/cell-atac/1.2.0/"
    "atac_pbmc_5k_nextgem/"
    "atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.tar.gz"
)


def _cluster(adata: ad.AnnData, representation: str) -> None:
    """Add a stable exploratory cluster without requiring leidenalg."""
    values = np.asarray(adata.obsm[representation])
    n_clusters = max(2, min(12, int(round(np.sqrt(adata.n_obs / 50)))))
    labels = KMeans(
        n_clusters=n_clusters,
        random_state=0,
        n_init=10,
    ).fit_predict(values[:, : min(values.shape[1], 30)])
    adata.obs["cluster"] = pd.Categorical(labels.astype(str))


def prepare_rna(matrix_dir: Path) -> ad.AnnData:
    rna = sc.read_10x_mtx(
        matrix_dir,
        var_names="gene_symbols",
        make_unique=True,
        cache=False,
    )
    rna.var_names_make_unique()
    rna.obs_names = pd.Index(
        [f"rna:{barcode}" for barcode in rna.obs_names],
        name="cell_id",
    )
    rna.obs["dataset"] = pd.Categorical(["PBMC3k RNA"] * rna.n_obs)
    rna.obs["modality"] = pd.Categorical(["RNA"] * rna.n_obs)

    rna.var["mt"] = rna.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        rna,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    rna = rna[
        (rna.obs["n_genes_by_counts"] >= 200)
        & (rna.obs["n_genes_by_counts"] < 2500)
        & (rna.obs["pct_counts_mt"] < 5)
    ].copy()
    sc.pp.filter_genes(rna, min_cells=3)
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(
        rna,
        n_comps=50,
        mask_var="highly_variable",
        random_state=0,
    )
    sc.pp.neighbors(rna, use_rep="X_pca", metric="cosine", random_state=0)
    sc.tl.umap(rna, random_state=0)
    _cluster(rna, "X_pca")
    return rna


def prepare_atac(matrix_dir: Path) -> ad.AnnData:
    matrix = sparse.csr_matrix(mmread(matrix_dir / "matrix.mtx").T)
    barcodes = pd.read_csv(
        matrix_dir / "barcodes.tsv",
        header=None,
        sep="\t",
    )[0].astype(str)
    peaks = pd.read_csv(
        matrix_dir / "peaks.bed",
        header=None,
        sep="\t",
        names=["chrom", "start", "end"],
    )
    peak_names = (
        peaks["chrom"].astype(str)
        + ":"
        + peaks["start"].astype(str)
        + "-"
        + peaks["end"].astype(str)
    )
    atac = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            index=pd.Index(
                [f"atac:{barcode}" for barcode in barcodes],
                name="cell_id",
            )
        ),
        var=peaks.set_axis(pd.Index(peak_names, name="peak")).copy(),
    )
    atac.obs["dataset"] = pd.Categorical(["PBMC5k ATAC"] * atac.n_obs)
    atac.obs["modality"] = pd.Categorical(["ATAC"] * atac.n_obs)
    sc.pp.calculate_qc_metrics(
        atac,
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Retain broadly informative peaks while keeping the demo quick to load.
    sc.pp.filter_genes(atac, min_cells=max(10, int(atac.n_obs * 0.05)))
    counts = sparse.csr_matrix(atac.X, dtype=np.float32)
    row_sums = np.asarray(counts.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1
    tf = sparse.diags(1.0 / row_sums) @ counts
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    idf = np.log1p(atac.n_obs / np.maximum(detected, 1))
    tfidf = tf.multiply(idf).tocsr()
    lsi = TruncatedSVD(
        n_components=50,
        n_iter=10,
        random_state=0,
    ).fit_transform(tfidf)
    atac.obsm["X_lsi"] = np.asarray(lsi, dtype=np.float32)
    sc.pp.neighbors(
        atac,
        use_rep="X_lsi",
        n_pcs=30,
        metric="cosine",
        random_state=0,
    )
    sc.tl.umap(atac, random_state=0)
    _cluster(atac, "X_lsi")
    return atac


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna-dir", type=Path, required=True)
    parser.add_argument("--atac-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rna = prepare_rna(args.rna_dir)
    atac = prepare_atac(args.atac_dir)
    overlap = rna.obs_names.intersection(atac.obs_names)
    if len(overlap):
        raise RuntimeError("The unpaired example unexpectedly contains shared cells.")

    mdata = mu.MuData({"rna": rna, "atac": atac})
    mdata.uns["guanaco_example"] = {
        "paired": False,
        "description": "Independent 10x PBMC RNA and ATAC experiments.",
        "rna_source": RNA_SOURCE,
        "atac_source": ATAC_SOURCE,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(args.output, compression="gzip")
    print(
        {
            "output": str(args.output),
            "rna_shape": rna.shape,
            "atac_shape": atac.shape,
            "shared_cells": len(overlap),
        }
    )


if __name__ == "__main__":
    main()
