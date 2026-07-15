"""Callback registration for matrix-backed visualizations."""

from .register import (
    apply_relayout,
    filter_data,
    generate_embedding_plots,
    is_continuous_annotation,
    matrix_callbacks,
    register_exploratory_plot_callbacks,
    register_marker_plot_callbacks,
)

__all__ = [
    "apply_relayout",
    "filter_data",
    "generate_embedding_plots",
    "is_continuous_annotation",
    "matrix_callbacks",
    "register_exploratory_plot_callbacks",
    "register_marker_plot_callbacks",
]
