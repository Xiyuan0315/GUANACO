"""Compatibility re-exports for matrix plot layouts.

Layout implementations live in `guanaco.pages.matrix.layout_modules`.
"""

from guanaco.pages.matrix.layout_modules.dotplot_layout import generate_dotplot_layout
from guanaco.pages.matrix.layout_modules.heatmap_layout import generate_heatmap_layout
from guanaco.pages.matrix.layout_modules.pseudotime_layout import generate_pseudotime_layout
from guanaco.pages.matrix.layout_modules.stacked_bar_layout import generate_stacked_bar_layout
from guanaco.pages.matrix.layout_modules.violin_layout import generate_violin_layout

__all__ = [
    "generate_heatmap_layout",
    "generate_violin_layout",
    "generate_dotplot_layout",
    "generate_stacked_bar_layout",
    "generate_pseudotime_layout",
]
