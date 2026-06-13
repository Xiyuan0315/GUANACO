"""Small matrix-backed ATAC genome browser utilities."""

from __future__ import annotations

import hashlib
import json
import re
import time
from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


_LOCUS_RE = re.compile(r"^\s*([^:\s]+)\s*:\s*([0-9,]+)\s*-\s*([0-9,]+)\s*$")
_PEAK_RE = re.compile(r"^\s*([^:\s]+)\s*:\s*([0-9,]+)\s*-\s*([0-9,]+)\s*$")


@dataclass
class PeakIndex:
    chroms: tuple[str, ...]
    starts: dict[str, np.ndarray]
    ends: dict[str, np.ndarray]
    var_indices: dict[str, np.ndarray]
    names: dict[str, np.ndarray]

    def bounds(self, chrom: str) -> tuple[int, int] | None:
        if chrom not in self.starts:
            return None
        return int(self.starts[chrom][0]), int(self.ends[chrom].max())


class _SimpleTTLCache:
    def __init__(self, max_items: int = 32, ttl_seconds: int = 180):
        self.max_items = max_items
        self.ttl_seconds = ttl_seconds
        self._store: OrderedDict[str, tuple[object, float]] = OrderedDict()

    def get(self, key: str):
        item = self._store.get(key)
        if item is None:
            return None
        value, ts = item
        if time.time() - ts > self.ttl_seconds:
            self._store.pop(key, None)
            return None
        self._store.move_to_end(key)
        return value

    def set(self, key: str, value) -> None:
        self._store[key] = (value, time.time())
        self._store.move_to_end(key)
        while len(self._store) > self.max_items:
            self._store.popitem(last=False)


_peak_index_cache: dict[tuple[int, int], PeakIndex] = {}
_signal_cache = _SimpleTTLCache(max_items=48, ttl_seconds=180)


def parse_locus(value: str | None) -> tuple[str, int, int] | None:
    if not value:
        return None
    match = _LOCUS_RE.match(str(value))
    if match is None:
        return None
    chrom = match.group(1)
    start = int(match.group(2).replace(",", ""))
    end = int(match.group(3).replace(",", ""))
    if end < start:
        start, end = end, start
    return chrom, max(0, start), max(0, end)


def format_locus(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{int(start):,}-{int(end):,}"


def _parse_peak_name(name: object) -> tuple[str, int, int] | None:
    match = _PEAK_RE.match(str(name))
    if match is None:
        return None
    chrom = match.group(1)
    start = int(match.group(2).replace(",", ""))
    end = int(match.group(3).replace(",", ""))
    if end < start:
        start, end = end, start
    return chrom, start, end


def has_genomic_peak_features(adata, *, sample_size: int = 2000, min_hits: int = 5) -> bool:
    if adata is None or getattr(adata, "n_vars", 0) == 0:
        return False
    hits = 0
    for name in list(adata.var_names[: min(sample_size, adata.n_vars)]):
        if _parse_peak_name(name) is not None:
            hits += 1
            if hits >= min_hits:
                return True
    var = getattr(adata, "var", None)
    if var is not None and {"chrom", "start", "end"}.issubset(set(var.columns)):
        return True
    return False


def build_peak_index(adata) -> PeakIndex:
    cache_key = (id(adata), getattr(adata, "n_vars", 0))
    cached = _peak_index_cache.get(cache_key)
    if cached is not None:
        return cached

    by_chrom: dict[str, list[tuple[int, int, int, str]]] = {}
    var = getattr(adata, "var", None)
    has_coord_cols = var is not None and {"chrom", "start", "end"}.issubset(set(var.columns))

    # Pull the coordinate columns into numpy arrays once instead of doing a
    # per-row ``var.iloc[idx]`` lookup -- the latter is an O(n_vars) sequence of
    # pandas row indexes that dominates index-building on 100k+ peak datasets.
    names_arr = np.asarray(adata.var_names, dtype=object)
    chrom_arr = start_arr = end_arr = None
    if has_coord_cols:
        chrom_arr = var["chrom"].astype(str).to_numpy()
        start_arr = pd.to_numeric(var["start"], errors="coerce").to_numpy()
        end_arr = pd.to_numeric(var["end"], errors="coerce").to_numpy()

    for idx in range(names_arr.size):
        name = names_arr[idx]
        parsed = None
        if has_coord_cols:
            start_val = start_arr[idx]
            end_val = end_arr[idx]
            if not (np.isnan(start_val) or np.isnan(end_val)):
                parsed = (chrom_arr[idx], int(start_val), int(end_val))
        if parsed is None:
            parsed = _parse_peak_name(name)
        if parsed is None:
            continue
        chrom, start, end = parsed
        by_chrom.setdefault(chrom, []).append((start, end, idx, str(name)))

    starts: dict[str, np.ndarray] = {}
    ends: dict[str, np.ndarray] = {}
    var_indices: dict[str, np.ndarray] = {}
    names: dict[str, np.ndarray] = {}

    def chrom_sort_key(chrom: str):
        normalized = chrom.removeprefix("chr")
        if normalized.isdigit():
            return (0, int(normalized))
        if normalized == "X":
            return (0, 23)
        if normalized == "Y":
            return (0, 24)
        if normalized in {"M", "MT"}:
            return (0, 25)
        return (1, normalized)

    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda item: (item[0], item[1]))
        starts[chrom] = np.asarray([r[0] for r in rows], dtype=np.int64)
        ends[chrom] = np.asarray([r[1] for r in rows], dtype=np.int64)
        var_indices[chrom] = np.asarray([r[2] for r in rows], dtype=np.int64)
        names[chrom] = np.asarray([r[3] for r in rows], dtype=object)

    index = PeakIndex(
        chroms=tuple(sorted(by_chrom.keys(), key=chrom_sort_key)),
        starts=starts,
        ends=ends,
        var_indices=var_indices,
        names=names,
    )
    _peak_index_cache[cache_key] = index
    return index


def default_region(adata, *, width: int = 1_000_000) -> dict[str, int | str]:
    index = build_peak_index(adata)
    preferred = next((c for c in ("chr1", "1") if c in index.starts), None)
    chrom = preferred or (index.chroms[0] if index.chroms else "chr1")
    starts = index.starts.get(chrom)
    if starts is None or starts.size == 0:
        return {"chrom": chrom, "start": 0, "end": width}
    center = int(starts[min(len(starts) // 2, len(starts) - 1)])
    start = max(0, center - width // 2)
    return {"chrom": chrom, "start": start, "end": start + width}


def coerce_region(index: PeakIndex, region: dict, *, min_span: int = 3_000, max_span: int = 5_000_000) -> dict[str, int | str]:
    chrom = str(region.get("chrom") or (index.chroms[0] if index.chroms else "chr1"))
    start = int(float(region.get("start", 0) or 0))
    end = int(float(region.get("end", start + 1_000_000) or start + 1_000_000))
    if end < start:
        start, end = end, start
    span = max(1, end - start)
    center = (start + end) // 2
    if span < min_span:
        span = min_span
    if span > max_span:
        span = max_span
    start = max(0, center - span // 2)
    end = start + span

    bounds = index.bounds(chrom)
    if bounds is not None:
        chrom_min, chrom_max = bounds
        if end < chrom_min:
            end = min(chrom_max, chrom_min + span)
            start = max(0, end - span)
        elif start > chrom_max:
            start = max(0, chrom_max - span)
            end = chrom_max
    return {"chrom": chrom, "start": int(start), "end": int(end)}


def peaks_in_region(index: PeakIndex, chrom: str, start: int, end: int, *, max_peaks: int = 400):
    if chrom not in index.starts:
        return {
            "starts": np.asarray([], dtype=np.int64),
            "ends": np.asarray([], dtype=np.int64),
            "var_indices": np.asarray([], dtype=np.int64),
            "names": np.asarray([], dtype=object),
            "total": 0,
            "downsampled": False,
        }

    starts = index.starts[chrom]
    stop = np.searchsorted(starts, end, side="right")
    candidate = np.arange(stop, dtype=np.int64)
    if candidate.size:
        candidate = candidate[index.ends[chrom][candidate] >= start]
    total = int(candidate.size)
    downsampled = False
    if total > max_peaks:
        take = np.linspace(0, total - 1, max_peaks).round().astype(np.int64)
        candidate = candidate[take]
        downsampled = True

    return {
        "starts": index.starts[chrom][candidate],
        "ends": index.ends[chrom][candidate],
        "var_indices": index.var_indices[chrom][candidate],
        "names": index.names[chrom][candidate],
        "total": total,
        "downsampled": downsampled,
    }


def _selected_signature(selected_cells) -> str:
    if not selected_cells:
        return "all"
    payload = json.dumps(list(selected_cells), separators=(",", ":"), default=str)
    digest = hashlib.md5(payload.encode("utf-8")).hexdigest()
    return f"{len(selected_cells)}:{digest}"


def _mean_or_detection(matrix, metric: str) -> np.ndarray:
    if matrix.shape[0] == 0:
        return np.zeros(matrix.shape[1], dtype=np.float32)
    if metric == "detection":
        return np.asarray((matrix > 0).mean(axis=0)).ravel().astype(np.float32, copy=False)
    return np.asarray(matrix.mean(axis=0)).ravel().astype(np.float32, copy=False)


def compute_atac_signal(
    adata,
    region: dict,
    *,
    selected_cells=None,
    groupby: str | None = None,
    labels=None,
    group_order=None,
    metric: str = "mean",
    max_peaks: int = 400,
):
    index = build_peak_index(adata)
    region = coerce_region(index, region)
    chrom = str(region["chrom"])
    start = int(region["start"])
    end = int(region["end"])
    peak_data = peaks_in_region(index, chrom, start, end, max_peaks=max_peaks)
    cols = peak_data["var_indices"]
    selected_sig = _selected_signature(selected_cells)
    label_sig = sorted(str(label) for label in labels) if labels else None
    cache_key = json.dumps(
        {
            "adata": id(adata),
            "chrom": chrom,
            "start": start,
            "end": end,
            "cols": cols.tolist(),
            "selected": selected_sig,
            "groupby": groupby,
            "labels": label_sig,
            "metric": metric,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    cached = _signal_cache.get(cache_key)
    if cached is not None:
        return cached

    if selected_cells:
        rows = adata.obs_names.get_indexer(selected_cells)
        rows = np.sort(rows[rows >= 0])
    else:
        rows = None
    n_cells = adata.n_obs if rows is None else int(rows.size)

    # Note: we deliberately do NOT bail out when the window has no peaks. The
    # grouping below still emits one (empty, flat) track per selected label, so
    # panning into a peak-free stretch keeps every cell-type track on screen
    # instead of collapsing the figure down to just the gene track.
    if rows is None:
        sub = adata.X[:, cols]
        row_labels = None
    else:
        sub = adata.X[rows, :][:, cols]
        row_labels = rows

    signals = []
    if groupby and groupby in adata.obs.columns:
        labels_source = adata.obs[groupby].astype(str).to_numpy()
        row_label_values = labels_source if row_labels is None else labels_source[row_labels]
        # Which groups get a track is driven by the left panel's selected labels;
        # selected_cells only narrows the cells inside each track. Order follows the
        # app-wide canonical order (group_order) so the tracks line up with the
        # heatmap/violin/dotplot rather than being sorted by cell count.
        if labels:
            wanted = [str(label) for label in labels]
        else:
            wanted = list(dict.fromkeys(row_label_values.tolist()))
        if group_order:
            rank = {str(group): i for i, group in enumerate(group_order)}
            wanted = sorted(dict.fromkeys(wanted), key=lambda group: rank.get(group, len(rank)))
        else:
            wanted = list(dict.fromkeys(wanted))
        for label in wanted:
            mask = row_label_values == label
            count = int(mask.sum())
            if count == 0:
                continue
            values = _mean_or_detection(sub[mask, :], metric)
            signals.append({"name": str(label), "values": values, "n_cells": count})
    else:
        values = _mean_or_detection(sub, metric)
        label = "Selected cells" if selected_cells else "All cells"
        signals.append({"name": label, "values": values, "n_cells": int(n_cells)})

    result = {
        "region": region,
        "peaks": peak_data,
        "signals": signals,
        "n_cells": int(n_cells),
        "metric": metric,
    }
    _signal_cache.set(cache_key, result)
    return result


def _pack_into_lane(lane_ends: list[int], start: int, end: int, padding: int) -> int:
    """Greedy interval packing: reuse the first lane that has cleared ``start``.

    Non-overlapping features share a row instead of each grabbing its own, which is
    what keeps the gene track compact (and readable) the way IGV's collapsed view is,
    rather than ballooning to one row per transcript.
    """
    for lane_idx, last_end in enumerate(lane_ends):
        if start > last_end + padding:
            lane_ends[lane_idx] = end
            return lane_idx
    lane_ends.append(end)
    return len(lane_ends) - 1


def _add_gene_model_track(
    fig,
    gene_models: list[dict],
    *,
    row: int,
    col: int,
    show_labels: bool = True,
    collapse: bool = False,
    view_span: int = 1,
) -> None:
    tx_x: list[int | None] = []
    tx_y: list[float | None] = []
    tx_text: list[str | None] = []
    exon_x: list[int | None] = []
    exon_y: list[float | None] = []
    exon_text: list[str | None] = []
    label_x: list[int] = []
    label_y: list[float] = []
    label_text: list[str] = []

    # Horizontal breathing room between packed features; wider when labels are drawn
    # so a gene's name doesn't sit on top of the next gene to its right.
    padding = max(1, int(view_span * (0.05 if show_labels else 0.006)))
    lane_ends: list[int] = []

    for gene in gene_models:
        gene_start = int(gene.get("start"))
        gene_end = int(gene.get("end"))
        direction = "<" if gene.get("strand") == "-" else ">"

        if collapse:
            # Zoomed-out / dense: one packed row per gene, gene body only (no exons).
            lane = _pack_into_lane(lane_ends, gene_start, gene_end, padding)
            y = float(lane)
            hover = (
                f"Gene: {gene.get('gene_name')}<br>"
                f"Gene ID: {gene.get('gene_id')}<br>"
                f"Type: {gene.get('gene_type', '')}<br>"
                f"Strand: {gene.get('strand')}<br>"
                f"{gene_start:,}-{gene_end:,}"
            )
            tx_x.extend([gene_start, gene_end, None])
            tx_y.extend([y, y, None])
            tx_text.extend([hover, hover, None])
            if show_labels:
                label_x.append(gene_start)
                label_y.append(y)
                label_text.append(f"{gene.get('gene_name')} {direction}")
            continue

        transcripts = gene.get("transcripts") or [
            {
                "transcript_id": gene.get("gene_id"),
                "transcript_name": gene.get("gene_name"),
                "start": gene.get("start"),
                "end": gene.get("end"),
                "strand": gene.get("strand"),
                "transcript_type": gene.get("gene_type", ""),
                "exons": [],
            }
        ]
        gene_label_added = False
        for tx in transcripts:
            tx_start = int(tx["start"])
            tx_end = int(tx["end"])
            lane = _pack_into_lane(lane_ends, tx_start, tx_end, padding)
            y = float(lane)
            hover = (
                f"Gene: {gene.get('gene_name')}<br>"
                f"Gene ID: {gene.get('gene_id')}<br>"
                f"Transcript: {tx.get('transcript_name') or tx.get('transcript_id')}<br>"
                f"Transcript ID: {tx.get('transcript_id')}<br>"
                f"Type: {tx.get('transcript_type') or gene.get('gene_type', '')}<br>"
                f"Strand: {tx.get('strand')}<br>"
                f"{tx_start:,}-{tx_end:,}"
            )
            tx_x.extend([tx_start, tx_end, None])
            tx_y.extend([y, y, None])
            tx_text.extend([hover, hover, None])
            for exon_start, exon_end in tx.get("exons", []):
                exon_x.extend([int(exon_start), int(exon_end), None])
                exon_y.extend([y, y, None])
                exon_text.extend([hover, hover, None])
            if show_labels and not gene_label_added:
                label_x.append(tx_start)
                label_y.append(y)
                label_text.append(f"{gene.get('gene_name')} {direction}")
                gene_label_added = True

    if tx_x:
        fig.add_trace(
            go.Scattergl(
                x=tx_x,
                y=tx_y,
                mode="lines",
                line={"color": "#555555", "width": 1.5},
                text=tx_text,
                hovertemplate="%{text}<extra>Transcript</extra>",
                showlegend=False,
                name="Transcripts",
            ),
            row=row,
            col=col,
        )
    if exon_x:
        fig.add_trace(
            go.Scattergl(
                x=exon_x,
                y=exon_y,
                mode="lines",
                line={"color": "#222222", "width": 7},
                text=exon_text,
                hovertemplate="%{text}<extra>Exon</extra>",
                showlegend=False,
                name="Exons",
            ),
            row=row,
            col=col,
        )
    if label_x:
        fig.add_trace(
            go.Scattergl(
                x=label_x,
                y=label_y,
                mode="text",
                text=label_text,
                textposition="top right",
                textfont={"size": 10, "color": "#222222"},
                hoverinfo="skip",
                showlegend=False,
                name="Gene labels",
            ),
            row=row,
            col=col,
        )
    lane_count = len(lane_ends)
    fig.update_yaxes(
        showticklabels=False,
        range=[max(lane_count, 1) - 0.35, -0.8],
        row=row,
        col=col,
    )


def _subplot_domain_ref(row: int) -> str:
    return "y domain" if row == 1 else f"y{row} domain"


def _add_track_badge(fig, *, row: int, text: str) -> None:
    fig.add_annotation(
        text=str(text),
        xref="paper",
        x=0.012,
        xanchor="left",
        yref=_subplot_domain_ref(row),
        y=0.94,
        yanchor="top",
        showarrow=False,
        font={"size": 11, "color": "#2f2f2f"},
        align="left",
        bgcolor="rgba(255,255,255,0.88)",
        bordercolor="rgba(90,90,90,0.35)",
        borderwidth=1,
        borderpad=4,
    )


def plot_atac_browser(
    signal_payload,
    gene_models: list[dict] | None = None,
    color_map: dict | None = None,
    gene_track_label: str = "Genes",
    y_mode: str = "auto",
    uirevision: str | None = None,
) -> go.Figure:
    region = signal_payload["region"]
    peaks = signal_payload["peaks"]
    signals = signal_payload["signals"]
    chrom = str(region["chrom"])
    start = int(region["start"])
    end = int(region["end"])
    has_gene_track = gene_models is not None

    if not signals and not has_gene_track:
        fig = go.Figure()
        fig.update_layout(
            template="plotly_white",
            height=260,
            xaxis={"range": [start, end], "title": "Genomic position", "showgrid": False, "zeroline": False},
            yaxis={"showgrid": False, "zeroline": False},
            showlegend=False,
        )
        return fig

    rows = len(signals) + (1 if has_gene_track else 0)
    row_heights = []
    if has_gene_track:
        row_heights.append(0.28)
    remaining = max(0.72, 1.0 - sum(row_heights))
    if signals:
        row_heights.extend([remaining / len(signals)] * len(signals))
    fig = make_subplots(
        rows=rows,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.025,
        row_heights=row_heights,
    )

    current_row = 1
    if has_gene_track:
        # IGV-like simplification when zoomed out: collapse isoforms to a single
        # gene body and drop the names once the view gets wide / gene-dense, so the
        # track stays legible instead of turning into a wall of overlapping labels.
        span = end - start
        n_genes = len(gene_models or [])
        collapse = span > 400_000 or n_genes > 10
        show_labels = n_genes <= 20
        _add_gene_model_track(
            fig,
            gene_models or [],
            row=current_row,
            col=1,
            show_labels=show_labels,
            collapse=collapse,
            view_span=max(1, span),
        )
        _add_track_badge(fig, row=current_row, text=gene_track_label)
        current_row += 1

    centers = (peaks["starts"].astype(np.float64) + peaks["ends"].astype(np.float64)) / 2.0
    widths = np.maximum(peaks["ends"] - peaks["starts"], 400)
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#5F4690", "#38A6A5"]
    metric_label = "Detection fraction" if signal_payload["metric"] == "detection" else "Mean accessibility"

    # In "shared" mode every track gets the same y-range (0 .. global max) so bar
    # heights are directly comparable between groups; "auto" lets each track scale
    # to its own peak.
    shared_max = None
    if y_mode == "shared":
        track_maxes = [float(np.nanmax(t["values"])) for t in signals if t["values"].size]
        shared_max = max(track_maxes) if track_maxes else None

    for track_idx, track in enumerate(signals):
        if color_map and track["name"] in color_map:
            track_color = color_map[track["name"]]
        else:
            track_color = palette[track_idx % len(palette)]
        fig.add_trace(
            go.Bar(
                x=centers,
                y=track["values"],
                width=widths,
                name=f"{track['name']} ({track['n_cells']})",
                marker={"color": track_color, "line": {"width": 0}},
                customdata=peaks["names"],
                hovertemplate=(
                    "Peak: %{customdata}<br>"
                    "Center: %{x:,.0f}<br>"
                    f"{metric_label}: " + "%{y:.3g}<br>"
                    f"Cells: {track['n_cells']}<extra>{track['name']}</extra>"
                ),
                showlegend=False,
            ),
            row=current_row,
            col=1,
        )
        # One y-tick per track: the peak height, at the top. Labelling 0 on every
        # track as well stacks a column of "0"s that collide between adjacent narrow
        # tracks (plotly then drops some, leaving zeros that read at a track's top on
        # one row and its bottom on the next). The baseline is obviously zero, so the
        # max alone is the informative number -- IGV labels tracks the same way.
        if shared_max is not None:
            track_top = shared_max
        else:
            track_top = float(np.nanmax(track["values"])) if track["values"].size else 0.0
        if track_top > 0:
            fig.update_yaxes(
                row=current_row,
                col=1,
                range=[0, track_top * 1.05],
                tickmode="array",
                tickvals=[track_top],
                ticktext=[f"{track_top:.2g}"],
            )
        else:
            fig.update_yaxes(row=current_row, col=1, range=[0, 1], tickmode="array", tickvals=[])
        _add_track_badge(fig, row=current_row, text=str(track["name"]))
        current_row += 1

    chrom_label = chrom[3:] if chrom.lower().startswith("chr") else chrom
    x_title = f"Chromosome {chrom_label}: {start:,}-{end:,}"
    fig.update_xaxes(range=[start, end], row=rows, col=1, title_text=x_title)
    fig.update_xaxes(showgrid=False, zeroline=False)
    # Lock the y-axes: navigation is horizontal-only (like IGV). dragmode="pan"
    # otherwise drags the y-axis too, and with a stable uirevision plotly would then
    # freeze each track's y-range across re-renders -- so panning into a new region
    # left the new bars outside the stuck range and they silently vanished.
    fig.update_yaxes(showgrid=False, zeroline=False, fixedrange=True)
    fig.update_layout(
        template="plotly_white",
        height=max(300, 130 + 105 * len(signals) + (120 if has_gene_track else 0)),
        margin={"l": 58, "r": 24, "t": 20, "b": 48},
        barmode="overlay",
        bargap=0,
        dragmode="pan",
        hovermode="x unified",
        uirevision=uirevision or "peak-browser",
        showlegend=False,
    )
    return fig


def signal_status(signal_payload) -> str:
    peaks = signal_payload["peaks"]
    selected = signal_payload["n_cells"]
    peak_text = f"{peaks['total']:,} peaks"
    if peaks["downsampled"]:
        peak_text += f" (showing {len(peaks['var_indices']):,})"
    return f"{peak_text} · {selected:,} cells"
