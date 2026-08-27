#!/usr/bin/env python3
"""Run LIANA consensus inference on the GUANACO PBMC RNA example."""

from __future__ import annotations

import argparse
import platform
from importlib.metadata import version
from pathlib import Path

import liana as li
import scanpy as sc


DEFAULT_SOURCE = Path(
    "/Users/xiyuanzhang/Documents/GUANACO_v2/data/pbmc_rna_demo.h5ad"
)
DEFAULT_OUTPUT = Path(
    "/Users/xiyuanzhang/Documents/GUANACO_v2/data/pbmc_liana_demo.h5ad"
)
DEFAULT_RESULTS = Path(
    "/Users/xiyuanzhang/Documents/GUANACO_v2/data/pbmc_liana_results.csv"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--groupby", default="CellType")
    parser.add_argument("--resource", default="consensus")
    parser.add_argument("--expr-prop", type=float, default=0.1)
    parser.add_argument("--min-cells", type=int, default=5)
    parser.add_argument("--n-perms", type=int, default=1000)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.source)
    if args.groupby not in adata.obs:
        raise KeyError(f"{args.groupby!r} is not present in adata.obs.")

    li.mt.rank_aggregate(
        adata,
        groupby=args.groupby,
        resource_name=args.resource,
        expr_prop=args.expr_prop,
        min_cells=args.min_cells,
        n_perms=args.n_perms,
        n_jobs=args.n_jobs,
        seed=1337,
        use_raw=False,
        key_added="liana_res",
        inplace=True,
        verbose=True,
    )

    results = adata.uns["liana_res"]
    args.results.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.results, index=False)

    adata.uns["guanaco_liana"] = {
        "description": "LIANA consensus inference on normalized log1p RNA.",
        "groupby": args.groupby,
        "resource": args.resource,
        "expr_prop": args.expr_prop,
        "min_cells": args.min_cells,
        "n_perms": args.n_perms,
        "seed": 1337,
        "use_raw": False,
        "source": str(args.source),
        "results_csv": str(args.results),
        "python": platform.python_version(),
        "liana": version("liana"),
        "scanpy": version("scanpy"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Wrote {len(results):,} LIANA interactions to {args.results}")
    print(f"Wrote LIANA result AnnData to {args.output}")


if __name__ == "__main__":
    main()
