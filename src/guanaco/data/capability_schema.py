"""Shared visualization capability identifiers.

This module intentionally has no scientific-computing imports.  The config wizard
can import it without pulling the runtime data loader, Muon, or Scanpy into the
process.
"""

PLOT_KEYS = (
    "dotplot",
    "heatmap",
    "violin",
    "ridge",
    "pseudotime",
    "split-violin",
    "stacked-bar",
    "paga",
    "volcano",
    "ligand-receptor",
    "spatial-relationships",
    "peak-browser",
    "igv",
    "multiomics-composition",
    "cross-modal-concordance",
)
