"""Canonical metadata and capability checks for visualization plots."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from guanaco.pages.matrix.plots.atac_browser import has_genomic_peak_features
from guanaco.pages.matrix.plots.volcano import has_volcano_data


Workspace = Literal["feature-analysis", "dataset-exploration"]

FEATURE_WORKSPACE: Workspace = "feature-analysis"
EXPLORATION_WORKSPACE: Workspace = "dataset-exploration"


@dataclass(frozen=True)
class PlotSpec:
    """Stable user-facing identity and availability policy for one plot type."""

    key: str
    label: str
    workspace: Workspace
    default_enabled: bool = False


PLOT_SPECS = (
    PlotSpec("dotplot", "Dot plot", FEATURE_WORKSPACE, True),
    PlotSpec("heatmap", "Heatmap", FEATURE_WORKSPACE, True),
    PlotSpec("violin", "Violin Plot", FEATURE_WORKSPACE, True),
    PlotSpec("pseudotime", "Expression Trend", FEATURE_WORKSPACE),
    PlotSpec("split-violin", "Comparative Violin", EXPLORATION_WORKSPACE, True),
    PlotSpec("stacked-bar", "Composition", EXPLORATION_WORKSPACE, True),
    PlotSpec("paga", "PAGA", EXPLORATION_WORKSPACE),
    PlotSpec("volcano", "Volcano Plot", EXPLORATION_WORKSPACE),
    PlotSpec("grn", "GRN", EXPLORATION_WORKSPACE),
    PlotSpec("peak-browser", "Peak Browser", EXPLORATION_WORKSPACE, True),
    PlotSpec("igv", "IGV", EXPLORATION_WORKSPACE),
    PlotSpec(
        "cross-modal-concordance",
        "Cross-modal concordance",
        EXPLORATION_WORKSPACE,
    ),
)

PLOT_SPECS_BY_KEY = {spec.key: spec for spec in PLOT_SPECS}
MARKER_PLOTS = tuple(
    spec.key for spec in PLOT_SPECS if spec.workspace == FEATURE_WORKSPACE
)
EXPLORATORY_PLOTS = tuple(
    spec.key for spec in PLOT_SPECS if spec.workspace == EXPLORATION_WORKSPACE
)

_COMPONENT_ALIASES = {
    "vocano": "volcano",
    "grn-demo": "grn",
    "expression-trend": "pseudotime",
    "expression_trend": "pseudotime",
    "stacked-violin": "violin",
    "violin2": "split-violin",
    "split_violin": "split-violin",
    "splitviolin": "split-violin",
    "grouped-violin": "split-violin",
    "group-violin": "split-violin",
    "trackplot": "peak-browser",
    "track-plot": "peak-browser",
    "genome-browser": "peak-browser",
    "atac_browser": "peak-browser",
    "atac-browser": "peak-browser",
}


def _is_available(key, adata, *, has_igv):
    if key == "igv":
        return has_igv
    if adata is None:
        return False
    if key == "peak-browser":
        return has_genomic_peak_features(adata)
    if key == "paga":
        return "paga" in adata.uns and "connectivities" in adata.uns["paga"]
    if key == "volcano":
        return has_volcano_data(adata)
    return True


def resolve_plot_components(adata, optional_plot_components=None, *, has_igv=False):
    """Resolve configured plot aliases and remove unavailable plot types."""
    if optional_plot_components is None:
        selected = [
            spec.key
            for spec in PLOT_SPECS
            if spec.default_enabled and _is_available(spec.key, adata, has_igv=has_igv)
        ]
    else:
        selected = []
        for component in optional_plot_components:
            canonical = _COMPONENT_ALIASES.get(component, component)
            if (
                canonical != "igv"
                and canonical in PLOT_SPECS_BY_KEY
                and canonical not in selected
                and _is_available(canonical, adata, has_igv=has_igv)
            ):
                selected.append(canonical)

    # IGV is capability-driven rather than a wizard option.
    if has_igv and "igv" not in selected:
        selected.append("igv")
    return tuple(selected)
