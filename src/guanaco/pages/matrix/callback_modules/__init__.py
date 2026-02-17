"""Per-plot callback registration modules for matrix page."""

from .scatter_callbacks import register_scatter_callbacks
from .heatmap_callbacks import register_heatmap_callbacks
from .dotplot_callbacks import register_dotplot_callbacks
from .pseudotime_callbacks import register_pseudotime_callbacks
from .violin_callbacks import register_violin_callbacks
from .stacked_bar_callbacks import register_stacked_bar_callbacks

__all__ = [
    "register_scatter_callbacks",
    "register_heatmap_callbacks",
    "register_dotplot_callbacks",
    "register_pseudotime_callbacks",
    "register_violin_callbacks",
    "register_stacked_bar_callbacks",
]
