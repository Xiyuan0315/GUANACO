"""Prepare the Squidpy Visium H&E tutorial data for the Guanaco example."""

from __future__ import annotations

import argparse
from pathlib import Path

import squidpy as sq


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()

    adata = sq.datasets.visium_hne_adata()
    adata.obs["cluster"] = adata.obs["cluster"].astype("category")

    sq.gr.spatial_neighbors(adata)
    sq.gr.nhood_enrichment(
        adata,
        cluster_key="cluster",
        n_perms=1000,
        n_jobs=args.n_jobs,
        seed=0,
    )
    sq.gr.co_occurrence(adata, cluster_key="cluster")

    adata.uns["guanaco_spatial_example"] = {
        "source": "squidpy.datasets.visium_hne_adata",
        "tutorial": (
            "https://squidpy.readthedocs.io/en/stable/notebooks/"
            "tutorials/tutorial_visium_hne.html"
        ),
        "squidpy_version": sq.__version__,
        "cluster_key": "cluster",
        "nhood_permutations": 1000,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")

    print(f"Saved {adata.n_obs:,} spots × {adata.n_vars:,} genes to {args.output}")
    print(
        "Stored:",
        "spatial_connectivities, cluster_nhood_enrichment, cluster_co_occurrence",
    )


if __name__ == "__main__":
    main()
