"""Composition and registry metadata for all visualization plots."""

from .layout import (
    generate_left_control,
    generate_visualization_sections,
)
from .registry import EXPLORATORY_PLOTS, MARKER_PLOTS, resolve_plot_components

__all__ = [
    "EXPLORATORY_PLOTS",
    "MARKER_PLOTS",
    "generate_left_control",
    "generate_visualization_sections",
    "resolve_plot_components",
]
