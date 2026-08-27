"""Optional, removable natural-language visualization demo.

The feature is intentionally kept outside the canonical visualization registry.
Add ``"ai-explorer"`` to ``optional_plot_components`` to expose it.
"""

from .callbacks import register_ai_explorer_callbacks
from .layout import generate_ai_explorer_layout


DEMO_COMPONENT_KEY = "ai-explorer"


def is_ai_explorer_requested(optional_plot_components) -> bool:
    """Return whether a dataset config explicitly opts into the demo."""
    return DEMO_COMPONENT_KEY in (optional_plot_components or ())


__all__ = [
    "DEMO_COMPONENT_KEY",
    "generate_ai_explorer_layout",
    "is_ai_explorer_requested",
    "register_ai_explorer_callbacks",
]
