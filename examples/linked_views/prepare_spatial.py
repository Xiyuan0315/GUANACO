"""Export the original notebook's spatial views without its unused gene matrix."""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def export_spatial(source: Path, destination: Path):
    """Keep all spots, original images, coordinates, labels and computed results."""
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite {destination}")
    original = ad.read_h5ad(source, backed="r")
    try:
        result = ad.AnnData(
            X=sparse.csr_matrix((original.n_obs, 1), dtype=np.float32),
            obs=original.obs[["cluster"]].copy(),
            var=pd.DataFrame(index=["spatial_placeholder"]),
            obsm={"spatial": original.obsm["spatial"].copy()},
            obsp={"spatial_distances": original.obsp["spatial_distances"].copy()},
            uns={key: original.uns[key] for key in (
                "spatial", "cluster_colors", "cluster_nhood_enrichment",
                "cluster_co_occurrence", "guanaco_spatial_example",
            )},
        )
        result.uns["guanaco_spatial_example"]["gallery_export"] = (
            "All spots, images, scale factors and spatial results retained unchanged. "
            "Gene expression and unrelated metadata omitted; X is an empty placeholder."
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix="spatial-export-", dir=destination.parent) as stage:
            staged = Path(stage) / destination.name
            with ad.settings.override(allow_write_nullable_strings=True):
                result.write_h5ad(staged, compression="gzip")
            staged.replace(destination)
    finally:
        original.file.close()
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("data") / "visium_hne_spatial.h5ad")
    args = parser.parse_args()
    output = export_spatial(args.source, args.output)
    print(f"Saved spatial-only example: {output} ({output.stat().st_size / 1024**2:.2f} MiB)")
