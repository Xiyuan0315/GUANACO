"""Build a self-contained HTML demo using existing GUANACO plotting functions."""

import argparse
from pathlib import Path
import sys

import guanaco as gc
from guanaco.static_export import StaticPanel, freeze_violin, write_dashboard

# Reuse the deterministic simulated data already used by the linked-view gallery.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "linked_views"))
from demo_data import make_single_cell  # noqa: E402


def build_demo(path, *, overwrite=False):
    adata = make_single_cell(n_cells=540, seed=202)
    genes = ["CD3D", "IL7R", "CCL5", "MS4A1", "NKG7", "LST1"]
    groups = {"Cell type": "cell_type", "Condition": "condition"}
    common = {"show": False, "return_fig": True, "height": 440}
    dots, violins, heatmaps = {}, {}, {}
    for label, groupby in groups.items():
        dots[label] = gc.pl.dotplot(adata, genes, groupby, **common)
        violins[label] = freeze_violin(
            gc.pl.violin(adata, genes[:3], groupby, **common)
        )
        heatmaps[label] = gc.pl.heatmap(adata, genes, groupby, n_bins=540, **common)
    return write_dashboard(
        path,
        [
            StaticPanel(
                "dotplot",
                "Expression overview",
                dots,
                "Dot colour represents mean expression; dot size represents the fraction of cells expressing each gene.",
            ),
            StaticPanel(
                "violin",
                "Expression distributions",
                violins,
                "Three marker genes. Gaussian density shapes are precomputed from GUANACO's violin samples; dotted lines mark sample means.",
            ),
            StaticPanel(
                "heatmap",
                "Expression across cells",
                heatmaps,
                "Six markers across 540 simulated cells, ordered by the selected grouping. This is not a group-mean heatmap.",
                wide=True,
            ),
        ],
        title="A small atlas, ready to share.",
        description=(
            "GUANACO static dashboard demo · 540 simulated cells · 6 cell types · 6 marker genes. "
            "Change each chart's prepared grouping, zoom in, or save a figure. The panels are independent; no linking is configured."
        ),
        overwrite=overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("demo.html")
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    result = build_demo(args.output, overwrite=args.overwrite)
    print(
        f"Open in your browser: {result}\nSelf-contained HTML: {result.stat().st_size / 1024**2:.2f} MiB"
    )
