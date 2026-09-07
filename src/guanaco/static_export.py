"""Small offline dashboards of precomputed Plotly figures; no Dash or linking.

Notebook and app callers pass the same ordinary figures. Each panel owns its
presets independently; exporting never invents cross-panel selection behaviour.
"""

from dataclasses import dataclass
import html
import json
from pathlib import Path
import re
from typing import Mapping, Sequence

import numpy as np
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs
from plotly.utils import PlotlyJSONEncoder


@dataclass(frozen=True)
class StaticPanel:
    id: str
    title: str
    figures: Mapping[str, go.Figure]
    description: str = ""
    wide: bool = False


def freeze_violin(figure: go.Figure, *, resolution: int = 192) -> go.Figure:
    """Freeze GUANACO's vertical, width-scaled violins into sampled KDE polygons.

    Uses the existing figure's samples and shared bandwidth, not AnnData again.
    Gaussian densities are evaluated here; no raw violin samples, KDE traces or
    box calculations are sent to the reader. Density curves are approximated by
    ``resolution`` vertices per side. Means refer to the exported KDE sample,
    which GUANACO may already have downsampled for large groups.

    This minimal adapter intentionally rejects horizontal/split/box violins.
    Other kinds of precomputed plots can go directly to ``write_dashboard``.
    """
    if resolution < 16 or resolution > 4096:
        raise ValueError("resolution must be between 16 and 4096.")
    groups = {}
    for trace in figure.data:
        if trace.type != "violin":
            continue
        if (
            trace.orientation not in (None, "v")
            or trace.side not in (None, "both")
            or trace.scalemode != "width"
            or trace.box.visible
            or trace.points is not False
            or trace.spanmode not in ("hard", "soft")
        ):
            raise ValueError(
                "Only vertical, width-scaled GUANACO violins without points/boxes are supported."
            )
        labels = list(dict.fromkeys(str(value) for value in trace.x))
        if len(labels) != 1:
            raise ValueError("Each violin trace must represent one group.")
        axis = trace.xaxis or "x"
        groups.setdefault(axis, [])
        if labels[0] not in groups[axis]:
            groups[axis].append(labels[0])

    result = go.Figure(layout=figure.layout)
    for trace in figure.data:
        axis = trace.xaxis or "x"
        if trace.type != "violin":
            copied = trace.to_plotly_json()
            # GUANACO uses invisible scatter anchors for its secondary y axes.
            if axis in groups and copied.get("x") is not None:
                copied["x"] = [groups[axis].index(str(value)) for value in copied["x"]]
            result.add_trace(copied)
            continue
        values = np.asarray(trace.y, dtype=float)
        values = values[np.isfinite(values)]
        if not values.size:
            continue
        label = str(trace.x[0])
        center = groups[axis].index(label)
        lo, hi = float(values.min()), float(values.max())
        width = float(trace.width or 0.8) / 2
        flat = hi == lo and not trace.bandwidth
        if flat:
            grid, offsets = np.array([lo]), np.array([width])
        else:
            if trace.bandwidth is None:
                raise ValueError(
                    "Nonconstant violin traces require an explicit precomputed bandwidth."
                )
            bandwidth = max(float(trace.bandwidth), (hi - lo) / 10000)
            pad = 2 * bandwidth if trace.spanmode == "soft" else 0
            grid = np.linspace(lo - pad, hi + pad, resolution)
            density = np.zeros(resolution)
            # Bounded intermediate allocation even for caller-supplied large traces.
            for start in range(0, len(values), 1024):
                z = (grid[:, None] - values[None, start : start + 1024]) / bandwidth
                density += np.exp(-0.5 * z * z).sum(axis=1)
            offsets = width * density / density.max()
        mean = float(values.mean())
        metadata = {
            "group": label,
            "sample_count": len(values),
            "sample_mean": mean,
            "sample_median": float(np.median(values)),
            "precomputed": True,
        }
        x = np.concatenate([center - offsets, (center + offsets)[::-1]])
        y = np.concatenate([grid, grid[::-1]])
        result.add_trace(
            go.Scatter(
                x=x,
                y=y,
                xaxis=axis,
                yaxis=trace.yaxis,
                mode="lines",
                fill="toself" if not flat else None,
                fillcolor=trace.fillcolor,
                line={"color": trace.line.color, "width": 1},
                name=label,
                showlegend=False,
                meta=metadata,
                hoveron="fills" if not flat else "points",
                hovertemplate=(
                    "%{meta.group}<br>KDE sample: %{meta.sample_count} values"
                    "<br>Sample mean: %{meta.sample_mean:.3f}"
                    "<br>Sample median: %{meta.sample_median:.3f}<extra></extra>"
                ),
            )
        )
        if trace.meanline.visible and not flat:
            half = float(np.interp(mean, grid, offsets))
            result.add_trace(
                go.Scatter(
                    x=[center - half, center + half],
                    y=[mean, mean],
                    xaxis=axis,
                    yaxis=trace.yaxis,
                    mode="lines",
                    line={"color": "#334155", "width": 1, "dash": "dot"},
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
    for axis, labels in groups.items():
        result.update_layout(
            **{
                "xaxis" + axis[1:]: {
                    "type": "linear",
                    "tickmode": "array",
                    "tickvals": list(range(len(labels))),
                    "ticktext": labels,
                    "range": [-0.6, len(labels) - 0.4],
                }
            }
        )
    return result


def write_dashboard(
    path,
    panels: Sequence[StaticPanel],
    *,
    title="GUANACO · Static viewer",
    description="",
    overwrite=False,
) -> Path:
    """Write one self-contained HTML with independent, precomputed panel presets.

    Supports ordinary scatter, heatmap and bar figures in this first version.
    Statistical traces that calculate distributions in Plotly.js are rejected;
    pass ``freeze_violin(gc.pl.violin(..., return_fig=True))`` for violins.
    The export includes plotted data and is NOT an access-control mechanism.
    """
    if not panels:
        raise ValueError("At least one panel is required.")
    payload, seen = [], set()
    for panel in panels:
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", panel.id) or panel.id in seen:
            raise ValueError(
                "Panel IDs must be unique and contain only letters, digits, '_' or '-'."
            )
        seen.add(panel.id)
        if not panel.figures:
            raise ValueError(f"Panel {panel.id} needs at least one precomputed preset.")
        variants = []
        for label, figure in panel.figures.items():
            if figure.frames:
                raise ValueError(
                    "Animation frames are not supported in this minimal exporter."
                )
            unsupported = {t.type for t in figure.data} - {
                "scatter",
                "scattergl",
                "heatmap",
                "bar",
            }
            if unsupported:
                raise ValueError(
                    f"Precompute unsupported traces first: {sorted(unsupported)}"
                )
            images = list(figure.layout.images)
            if figure.layout.template:
                images += list(figure.layout.template.layout.images)
            if any(not str(image.source).startswith("data:image/") for image in images):
                raise ValueError(
                    "Offline figures require embedded data:image resources, not external images."
                )
            spec = figure.to_plotly_json()
            # Viewer controls are explicit presets, not hidden figure-side actions.
            if spec.get("layout", {}).get("updatemenus") or spec.get("layout", {}).get(
                "sliders"
            ):
                raise ValueError(
                    "Use panel presets instead of figure sliders/updatemenus."
                )
            variants.append({"label": str(label), "figure": spec})
        payload.append(
            {
                "id": panel.id,
                "title": panel.title,
                "description": panel.description,
                "wide": panel.wide,
                "presets": variants,
            }
        )
    # Keep resources OUTSIDE Dash's auto-loaded assets directory: this viewer
    # must not inject its scripts or global CSS into existing GUANACO dashboards.
    assets = Path(__file__).parent / "_static_viewer"
    # Escape script delimiters in untrusted labels/data; visible strings use textContent.
    data = json.dumps(payload, cls=PlotlyJSONEncoder).replace("<", "\\u003c")
    plotly_js = get_plotlyjs().replace("</script", "<\\/script")
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline' 'unsafe-eval'; style-src 'unsafe-inline'; img-src data: blob:; font-src data:; connect-src 'none'; worker-src blob:">
<title>{html.escape(title)}</title>
<style>{(assets / "viewer.css").read_text(encoding="utf-8")}</style>
</head><body><main>
<header><div class="eyebrow">GUANACO <span>DATA DISSEMINATION</span></div>
<h1>{html.escape(title)}</h1><p class="intro">{html.escape(description)}</p>
<div class="badges"><span>Offline HTML</span><span>Precomputed results</span><span>No linking</span></div></header>
<aside><strong>Explore prepared results.</strong> Each chart has its own grouping selector.
Switching a preset displays a saved result; it does not regroup cells or run an analysis.</aside>
<div id="dashboard" class="dashboard"></div>
<footer>This file contains the displayed data. Share only with permitted recipients.
No Python server, remote scripts or network requests are needed. Modern JavaScript-enabled browser required.</footer>
<noscript>This interactive viewer requires JavaScript to be enabled.</noscript>
</main><script>{plotly_js}</script>
<script type="application/json" id="guanaco-static-data">{data}</script>
<script>{(assets / "viewer.js").read_text(encoding="utf-8")}</script></body></html>"""
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w" if overwrite else "x", encoding="utf-8") as stream:
        stream.write(document)
    return path
