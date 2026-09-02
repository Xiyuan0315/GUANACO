"""Shared visual defaults for GUANACO Plotly figures.

Plot builders remain free to choose domain-specific axes and encodings.  This
module owns the small set of presentation decisions that should not drift
between dashboard, notebook, and linked-view surfaces.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import plotly.graph_objects as go


GUANACO_QUALITATIVE = (
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
)

GUANACO_FONT_FAMILY = (
    "Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "
    "'Segoe UI', sans-serif"
)


def categorical_color_map(
    labels: Sequence[Any],
    palette: str | Sequence[str] | Mapping[Any, str] | None = None,
) -> dict[Any, str]:
    """Return a stable category-to-color mapping for arbitrary labels."""

    ordered = list(dict.fromkeys(labels))
    if isinstance(palette, Mapping):
        fallback = list(GUANACO_QUALITATIVE)
        return {
            label: str(palette.get(label, fallback[index % len(fallback)]))
            for index, label in enumerate(ordered)
        }
    if isinstance(palette, str):
        from guanaco.utils.colors import resolve_discrete_palette

        colors = list(resolve_discrete_palette(palette, len(ordered)) or ())
    elif palette is None:
        colors = list(GUANACO_QUALITATIVE)
    else:
        colors = [str(color) for color in palette]
    if not colors:
        colors = list(GUANACO_QUALITATIVE)
    return {label: colors[index % len(colors)] for index, label in enumerate(ordered)}


def apply_guanaco_figure_style(
    figure: go.Figure,
    *,
    title: str | None = None,
) -> go.Figure:
    """Apply GUANACO's shared Plotly surface without changing data encodings."""

    title_update: dict[str, Any] = {
        "x": 0.5,
        "xanchor": "center",
        "y": 0.98,
        "yanchor": "top",
        "automargin": True,
    }
    if title is not None:
        title_update["text"] = title
    figure.update_layout(
        title=title_update,
        paper_bgcolor="white",
        plot_bgcolor="white",
        font={"family": GUANACO_FONT_FAMILY, "color": "#334155"},
        hoverlabel={
            "bgcolor": "white",
            "font": {"family": GUANACO_FONT_FAMILY, "color": "#0F172A"},
        },
    )
    return figure


__all__ = [
    "GUANACO_FONT_FAMILY",
    "GUANACO_QUALITATIVE",
    "apply_guanaco_figure_style",
    "categorical_color_map",
]
