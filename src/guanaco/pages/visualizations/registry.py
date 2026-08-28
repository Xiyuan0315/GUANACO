"""Canonical metadata and capability checks for visualization plots."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from guanaco.data.capabilities import PLOT_KEYS, has_plot_capability


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
    PlotSpec("network", "Network", EXPLORATION_WORKSPACE),
    PlotSpec("ligand-receptor", "Ligand–receptor", EXPLORATION_WORKSPACE),
    PlotSpec(
        "spatial-relationships",
        "Spatial relationships",
        EXPLORATION_WORKSPACE,
    ),
    PlotSpec("peak-browser", "Peak Browser", EXPLORATION_WORKSPACE, True),
    PlotSpec("igv", "IGV", EXPLORATION_WORKSPACE),
    PlotSpec(
        "multiomics-composition",
        "Multi-omics coverage",
        EXPLORATION_WORKSPACE,
        True,
    ),
    PlotSpec(
        "cross-modal-concordance",
        "Omics comparison",
        EXPLORATION_WORKSPACE,
    ),
)

PLOT_SPECS_BY_KEY = {spec.key: spec for spec in PLOT_SPECS}
if tuple(PLOT_SPECS_BY_KEY) != PLOT_KEYS:
    raise RuntimeError("Visualization specs and data capabilities are out of sync.")
MARKER_PLOTS = tuple(
    spec.key for spec in PLOT_SPECS if spec.workspace == FEATURE_WORKSPACE
)
EXPLORATORY_PLOTS = tuple(
    spec.key for spec in PLOT_SPECS if spec.workspace == EXPLORATION_WORKSPACE
)


def multiomics_plot_components(enabled) -> tuple[str, ...]:
    """Keep valid joint-view plots in one stable, shared order."""
    enabled = set(enabled)
    return (
        *(key for key in MARKER_PLOTS if key in enabled),
        *(key for key in ("multiomics-composition",) if key in enabled),
        "cross-modal-concordance",
    )


_COMPONENT_ALIASES = {
    "vocano": "volcano",
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
    "spatial-neighborhood": "spatial-relationships",
    "spatial-relations": "spatial-relationships",
    "omics-comparison": "cross-modal-concordance",
    "cross-modal-comparison": "cross-modal-concordance",
    "multiomics-coverage": "multiomics-composition",
    "omics-coverage": "multiomics-composition",
}

def is_plot_available(
    key,
    adata,
    *,
    has_igv=False,
    modality_name: str | None = None,
    feature_data_available: bool | None = None,
    discrete_data_available: bool | None = None,
) -> bool:
    """Return whether one plot can produce a useful view from this data."""
    key = _COMPONENT_ALIASES.get(key, key)
    return key in PLOT_SPECS_BY_KEY and has_plot_capability(
        key,
        adata,
        has_igv=has_igv,
        modality_name=modality_name,
        feature_data_available=feature_data_available,
        discrete_data_available=discrete_data_available,
    )


def available_plot_components(
    adata,
    *,
    has_igv=False,
    modality_name: str | None = None,
    feature_data_available: bool | None = None,
    discrete_data_available: bool | None = None,
) -> tuple[str, ...]:
    """Return every compatible plot in stable display order."""
    return tuple(
        spec.key
        for spec in PLOT_SPECS
        if is_plot_available(
            spec.key,
            adata,
            has_igv=has_igv,
            modality_name=modality_name,
            feature_data_available=feature_data_available,
            discrete_data_available=discrete_data_available,
        )
    )


def resolve_plot_components(
    adata,
    optional_plot_components=None,
    *,
    has_igv=False,
    modality_name: str | None = None,
    feature_data_available: bool | None = None,
    discrete_data_available: bool | None = None,
):
    """Resolve configured plot aliases and remove unavailable plot types."""
    if optional_plot_components is None:
        selected = [
            spec.key
            for spec in PLOT_SPECS
            if spec.default_enabled
            and is_plot_available(
                spec.key,
                adata,
                has_igv=has_igv,
                modality_name=modality_name,
                feature_data_available=feature_data_available,
                discrete_data_available=discrete_data_available,
            )
        ]
    else:
        selected = []
        for component in optional_plot_components:
            canonical = _COMPONENT_ALIASES.get(component, component)
            if (
                canonical != "igv"
                and canonical in PLOT_SPECS_BY_KEY
                and canonical not in selected
                and is_plot_available(
                    canonical,
                    adata,
                    has_igv=has_igv,
                    modality_name=modality_name,
                    feature_data_available=feature_data_available,
                    discrete_data_available=discrete_data_available,
                )
            ):
                selected.append(canonical)

    # IGV is capability-driven rather than a wizard option.
    if has_igv and "igv" not in selected:
        selected.append("igv")
    return tuple(selected)
