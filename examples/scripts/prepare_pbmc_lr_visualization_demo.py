#!/usr/bin/env python3
"""Attach a small precomputed LR result table to the PBMC demo AnnData."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = Path(
    "/Users/xiyuanzhang/Documents/GUANACO_v2/data/pbmc_rna_demo.h5ad"
)
DEFAULT_OUTPUT = Path(
    "/Users/xiyuanzhang/Documents/GUANACO_v2/data/pbmc_lr_visualization_demo.h5ad"
)
DEFAULT_RESULTS = REPOSITORY / "examples/data/pbmc_ligand_receptor_demo.csv"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.source)
    interactions = pd.read_csv(args.results)
    adata.uns["ligand_receptor_demo"] = interactions
    adata.uns["guanaco_ligand_receptor_demo"] = {
        "description": (
            "Synthetic PBMC interactions for testing visualization only; "
            "these rows are not statistical inference results."
        ),
        "source_table": str(args.results),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(
        f"Wrote {args.output} with {len(interactions)} visualization-demo "
        "interactions."
    )


if __name__ == "__main__":
    main()
