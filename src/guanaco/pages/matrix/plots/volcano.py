from __future__ import annotations

import csv
import io
import re
from typing import Any

import numpy as np
import plotly.graph_objects as go


DEFAULT_PADJ_THRESHOLD = 0.05
DEFAULT_X_THRESHOLD = 1.0
DEFAULT_TOP_N = 12
CORE_ENTRY_FIELDS = {"gene", "logfoldchange", "score", "padj", "neg_log10_padj"}


def decode_value(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def decode_extra_value(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def clipped_neglog10(padj: np.ndarray) -> tuple[np.ndarray, float]:
    padj = np.asarray(padj, dtype=float)
    finite_positive = padj[np.isfinite(padj) & (padj > 0)]
    floor = float(finite_positive.min()) if finite_positive.size else float(np.finfo(float).tiny)
    clipped = np.where(np.isfinite(padj) & (padj > 0), padj, floor)
    return -np.log10(clipped), floor


def normalize_lengths(entry_name: str, arrays: dict[str, Any]) -> None:
    lengths = {key: len(value) for key, value in arrays.items()}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"Entry '{entry_name}' has mismatched lengths: {lengths}")


def extract_extra_columns(raw_entry: dict[str, Any], gene_count: int) -> dict[str, list[Any]]:
    extra_columns: dict[str, list[Any]] = {}
    for key, value in raw_entry.items():
        decoded_key = decode_value(key)
        if decoded_key in CORE_ENTRY_FIELDS:
            continue
        try:
            if len(value) != gene_count:
                continue
        except TypeError:
            continue
        extra_columns[decoded_key] = [decode_extra_value(item) for item in value]
    return extra_columns


def has_volcano_data(adata: Any) -> bool:
    return bool(adata is not None and ("volcano" in adata.uns or "rank_genes_groups" in adata.uns))


def load_volcano_payload(adata: Any) -> dict[str, Any]:
    if "volcano" in adata.uns:
        return extract_generic_volcano(adata, key="volcano")
    if "rank_genes_groups" in adata.uns:
        return extract_scanpy_rank_genes_groups(adata, key="rank_genes_groups")
    raise ValueError("No supported DE result found. Expected adata.uns['volcano'] or adata.uns['rank_genes_groups'].")


def extract_scanpy_rank_genes_groups(adata: Any, key: str = "rank_genes_groups") -> dict[str, Any]:
    rank_genes_groups = adata.uns[key]
    names = rank_genes_groups["names"]
    groups = list(names.dtype.names or [])
    pval_key = "pvals_adj" if "pvals_adj" in rank_genes_groups else "pvals"

    params = {}
    if "params" in rank_genes_groups:
        for param_key, param_value in rank_genes_groups["params"].items():
            params[decode_value(param_key)] = decode_value(param_value)

    entries: dict[str, Any] = {}
    for group in groups:
        gene = [decode_value(item) for item in names[group]]
        logfoldchange = np.asarray(rank_genes_groups["logfoldchanges"][group], dtype=float)
        score = (
            np.asarray(rank_genes_groups["scores"][group], dtype=float)
            if "scores" in rank_genes_groups
            else np.full(len(gene), np.nan)
        )
        padj = np.asarray(rank_genes_groups[pval_key][group], dtype=float)
        neg_log10_padj, clipped_floor = clipped_neglog10(padj)

        entry = {
            "gene": gene,
            "logfoldchange": logfoldchange.tolist(),
            "score": score.tolist(),
            "padj": padj.tolist(),
            "neg_log10_padj": neg_log10_padj.tolist(),
            "extra_columns": {},
            "clipped_floor": clipped_floor,
            "meta": {
                "source": "scanpy",
                "group": decode_value(group),
                "key": key,
                "pval_key": pval_key,
                "params": params,
            },
        }
        normalize_lengths(
            decode_value(group),
            {
                "gene": entry["gene"],
                "logfoldchange": entry["logfoldchange"],
                "score": entry["score"],
                "padj": entry["padj"],
                "neg_log10_padj": entry["neg_log10_padj"],
            },
        )
        entries[decode_value(group)] = entry

    return {
        "source_type": "scanpy",
        "entries": entries,
        "default_entry": next(iter(entries), None),
    }


def extract_generic_volcano(adata: Any, key: str = "volcano") -> dict[str, Any]:
    volcano = adata.uns[key]
    if "entries" not in volcano:
        raise ValueError("adata.uns['volcano'] must contain an 'entries' mapping.")

    entries: dict[str, Any] = {}
    for entry_name, raw_entry in volcano["entries"].items():
        decoded_entry_name = decode_value(entry_name)
        gene = [decode_value(item) for item in raw_entry["gene"]]
        logfoldchange = np.asarray(raw_entry["logfoldchange"], dtype=float)
        score = (
            np.asarray(raw_entry["score"], dtype=float)
            if "score" in raw_entry
            else np.full(len(gene), np.nan)
        )

        if "padj" in raw_entry:
            pval_key = "padj"
            padj = np.asarray(raw_entry["padj"], dtype=float)
            neg_log10_padj, clipped_floor = clipped_neglog10(padj)
        elif "neg_log10_padj" in raw_entry:
            pval_key = "neg_log10_padj"
            neg_log10_padj = np.asarray(raw_entry["neg_log10_padj"], dtype=float)
            padj = np.power(10.0, -neg_log10_padj)
            finite_padj = padj[np.isfinite(padj)]
            clipped_floor = float(finite_padj.min()) if finite_padj.size else float("nan")
        else:
            raise ValueError(f"Entry '{decoded_entry_name}' must contain either 'padj' or 'neg_log10_padj'.")

        entry = {
            "gene": gene,
            "logfoldchange": logfoldchange.tolist(),
            "score": score.tolist(),
            "padj": padj.tolist(),
            "neg_log10_padj": neg_log10_padj.tolist(),
            "extra_columns": extract_extra_columns(raw_entry, len(gene)),
            "clipped_floor": clipped_floor,
            "meta": {
                "source": "volcano",
                "entry": decoded_entry_name,
                "pval_key": pval_key,
            },
        }
        normalize_lengths(
            decoded_entry_name,
            {
                "gene": entry["gene"],
                "logfoldchange": entry["logfoldchange"],
                "score": entry["score"],
                "padj": entry["padj"],
                "neg_log10_padj": entry["neg_log10_padj"],
            },
        )
        entries[decoded_entry_name] = entry

    default_entry = decode_value(volcano.get("default_entry", next(iter(entries), "")))
    if default_entry not in entries:
        default_entry = next(iter(entries), None)
    return {
        "source_type": "volcano",
        "entries": entries,
        "default_entry": default_entry,
    }


def volcano_entry_options(adata: Any) -> tuple[list[dict[str, str]], str | None]:
    if not has_volcano_data(adata):
        return [], None
    try:
        payload = load_volcano_payload(adata)
    except Exception:
        return [], None
    options = [{"label": entry_name, "value": entry_name} for entry_name in payload["entries"]]
    return options, payload["default_entry"]


def pvalue_key(entry: dict[str, Any]) -> str:
    meta = entry.get("meta", {}) or {}
    return decode_value(meta.get("pval_key", "padj"))


def pvalue_axis_label(entry: dict[str, Any]) -> str:
    key = pvalue_key(entry)
    if key == "neg_log10_padj":
        return key
    return f"-log10({key})"


def x_axis_options(entry: dict[str, Any] | None = None) -> list[dict[str, str]]:
    options = [{"label": "log fold change", "value": "logfoldchange"}]
    if entry is None or np.isfinite(np.asarray(entry.get("score", []), dtype=float)).any():
        options.append({"label": "score", "value": "score"})
    return options


def classify_points(
    x_values: np.ndarray,
    padj: np.ndarray,
    padj_threshold: float,
    x_threshold: float,
) -> np.ndarray:
    categories = np.full(len(x_values), "Not significant", dtype=object)
    significant = np.isfinite(padj) & (padj <= padj_threshold)
    categories[significant & np.isfinite(x_values) & (x_values >= x_threshold)] = "Upregulated"
    categories[significant & np.isfinite(x_values) & (x_values <= -x_threshold)] = "Downregulated"
    return categories


def pick_label_indices(
    entry: dict[str, Any],
    x_field: str,
    padj_threshold: float,
    x_threshold: float,
    top_n: int,
) -> np.ndarray:
    if top_n <= 0:
        return np.array([], dtype=int)

    x_values = np.asarray(entry[x_field], dtype=float)
    padj = np.asarray(entry["padj"], dtype=float)
    candidate = np.where(
        np.isfinite(padj)
        & (padj <= padj_threshold)
        & np.isfinite(x_values)
        & (np.abs(x_values) >= x_threshold)
    )[0]
    if candidate.size == 0:
        return np.array([], dtype=int)

    order = np.lexsort((-np.abs(x_values[candidate]), padj[candidate]))
    return candidate[order][:top_n]


def selected_deg_indices(
    entry: dict[str, Any],
    x_field: str,
    padj_threshold: float,
    x_threshold: float,
) -> np.ndarray:
    x_values = np.asarray(entry[x_field], dtype=float)
    padj = np.asarray(entry["padj"], dtype=float)
    return np.where(
        np.isfinite(padj)
        & (padj <= padj_threshold)
        & np.isfinite(x_values)
        & (np.abs(x_values) >= x_threshold)
    )[0]


def clean_csv_value(value: Any) -> Any:
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return ""
    return value


def deg_rows(
    entry: dict[str, Any],
    x_field: str,
    padj_threshold: float,
    x_threshold: float,
) -> list[dict[str, Any]]:
    indices = selected_deg_indices(entry, x_field, padj_threshold, x_threshold)
    category = classify_points(
        np.asarray(entry[x_field], dtype=float),
        np.asarray(entry["padj"], dtype=float),
        padj_threshold,
        x_threshold,
    )
    extra_columns = entry.get("extra_columns", {}) or {}

    rows: list[dict[str, Any]] = []
    for idx in indices:
        row = {
            "gene": entry["gene"][idx],
            "selected_x_axis": x_field,
            "selected_x_value": clean_csv_value(entry[x_field][idx]),
            "category": category[idx],
            "logfoldchange": clean_csv_value(entry["logfoldchange"][idx]),
            "score": clean_csv_value(entry["score"][idx]),
            "padj": clean_csv_value(entry["padj"][idx]),
            "neg_log10_padj": clean_csv_value(entry["neg_log10_padj"][idx]),
        }
        for key, values in extra_columns.items():
            row[key] = clean_csv_value(values[idx])
        rows.append(row)

    rows.sort(
        key=lambda item: (
            item["padj"] if item["padj"] != "" else float("inf"),
            -abs(float(item["selected_x_value"] or 0)),
        )
    )
    return rows


def deg_summary(entry: dict[str, Any], x_field: str, padj_threshold: float, x_threshold: float) -> dict[str, Any]:
    rows = deg_rows(entry, x_field, padj_threshold, x_threshold)
    up_count = sum(1 for row in rows if row["category"] == "Upregulated")
    down_count = sum(1 for row in rows if row["category"] == "Downregulated")
    return {
        "total": len(rows),
        "up": up_count,
        "down": down_count,
        "criteria": f"{pvalue_key(entry)} <= {padj_threshold:g}, |{x_field}| >= {x_threshold:g}",
    }


def deg_csv(entry: dict[str, Any], x_field: str, padj_threshold: float, x_threshold: float) -> str:
    rows = deg_rows(entry, x_field, padj_threshold, x_threshold)
    default_columns = [
        "gene",
        "selected_x_axis",
        "selected_x_value",
        "category",
        "logfoldchange",
        "score",
        "padj",
        "neg_log10_padj",
    ]
    extra_columns = list((entry.get("extra_columns", {}) or {}).keys())
    columns = default_columns + [column for column in extra_columns if column not in default_columns]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return cleaned or "degs"


def _compact_threshold(value: float) -> str:
    return f"{float(value):g}".replace("-", "neg").replace(".", "")


def _filename_pvalue_label(entry: dict[str, Any]) -> str:
    key = pvalue_key(entry).lower()
    if key in {"pvals_adj", "padj", "p_adj", "adjusted_pvalue", "adjusted_p_value"}:
        return "padj"
    if key in {"pvals", "pval", "pvalue", "p_value", "p"}:
        return "p"
    return safe_filename(key)


def _filename_x_axis_label(x_field: str) -> str:
    if x_field == "logfoldchange":
        return "log2fc"
    return safe_filename(x_field)


def volcano_degs_filename(
    entry_name: str,
    entry: dict[str, Any],
    x_field: str,
    padj_threshold: float,
    x_threshold: float,
) -> str:
    return (
        f"{safe_filename(entry_name)}_"
        f"{_filename_pvalue_label(entry)}{_compact_threshold(padj_threshold)}_"
        f"{_filename_x_axis_label(x_field)}{_compact_threshold(x_threshold)}_deg.csv"
    )


def finite_extent(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -1.0, 1.0

    vmin = float(finite.min())
    vmax = float(finite.max())
    if vmin == vmax:
        return vmin - 1.0, vmax + 1.0

    pad = (vmax - vmin) * 0.08
    return vmin - pad, vmax + pad


def plot_volcano(
    entry_name: str,
    entry: dict[str, Any],
    x_field: str = "logfoldchange",
    padj_threshold: float = DEFAULT_PADJ_THRESHOLD,
    x_threshold: float = DEFAULT_X_THRESHOLD,
    top_n: int = DEFAULT_TOP_N,
) -> go.Figure:
    colors = {
        "Upregulated": "#d1495b",
        "Downregulated": "#2b59c3",
        "Not significant": "#9aa5b1",
    }
    axis_labels = {
        "logfoldchange": "log fold change",
        "score": "score",
    }
    if x_field not in axis_labels:
        x_field = "logfoldchange"

    x = np.asarray(entry[x_field], dtype=float)
    y = np.asarray(entry["neg_log10_padj"], dtype=float)
    y_label = pvalue_axis_label(entry)
    logfoldchange = np.asarray(entry["logfoldchange"], dtype=float)
    score = np.asarray(entry["score"], dtype=float)
    padj = np.asarray(entry["padj"], dtype=float)
    gene = np.asarray(entry["gene"], dtype=object)

    category = classify_points(x, padj, padj_threshold, x_threshold)
    fig = go.Figure()

    for category_name in ["Upregulated", "Downregulated", "Not significant"]:
        mask = category == category_name
        fig.add_trace(
            go.Scattergl(
                x=x[mask],
                y=y[mask],
                mode="markers",
                name=category_name,
                marker={"size": 7, "opacity": 0.72, "color": colors[category_name]},
                customdata=np.column_stack([gene[mask], logfoldchange[mask], score[mask], padj[mask]]),
                hovertemplate=(
                    "Gene: %{customdata[0]}<br>"
                    f"{axis_labels[x_field]}: %{{x:.3f}}<br>"
                    f"{y_label}: %{{y:.3f}}<br>"
                    "logFC: %{customdata[1]:.3f}<br>"
                    "score: %{customdata[2]:.3f}<br>"
                    f"{pvalue_key(entry)}: %{{customdata[3]:.3e}}<extra></extra>"
                ),
            )
        )

    label_idx = pick_label_indices(entry, x_field, padj_threshold, x_threshold, top_n)
    if label_idx.size:
        fig.add_trace(
            go.Scattergl(
                x=x[label_idx],
                y=y[label_idx],
                mode="markers+text",
                text=gene[label_idx],
                textposition="top center",
                textfont={"size": 10, "color": "#1f2933"},
                marker={"size": 9, "color": "#111111"},
                name=f"Top {label_idx.size} labels",
                hoverinfo="skip",
            )
        )

    y_threshold = -np.log10(max(padj_threshold, np.nextafter(0.0, 1.0)))
    xmin, xmax = finite_extent(x)
    ymin, ymax = finite_extent(y)
    shapes = [
        {
            "type": "line",
            "x0": xmin,
            "x1": xmax,
            "y0": y_threshold,
            "y1": y_threshold,
            "line": {"color": "#555555", "width": 1.5, "dash": "dash"},
        },
        {
            "type": "line",
            "x0": x_threshold,
            "x1": x_threshold,
            "y0": ymin,
            "y1": ymax,
            "line": {"color": "#555555", "width": 1.5, "dash": "dot"},
        },
        {
            "type": "line",
            "x0": -x_threshold,
            "x1": -x_threshold,
            "y0": ymin,
            "y1": ymax,
            "line": {"color": "#555555", "width": 1.5, "dash": "dot"},
        },
    ]

    fig.update_layout(
        template="plotly_white",
        height=620,
        title=f"Volcano Plot: {entry_name}",
        xaxis_title=axis_labels[x_field],
        yaxis_title=y_label,
        legend_title_text="Category",
        margin={"l": 50, "r": 20, "t": 60, "b": 45},
        shapes=shapes,
    )
    return fig


def empty_volcano_figure(message: str) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(
        template="plotly_white",
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[
            {
                "text": message,
                "xref": "paper",
                "yref": "paper",
                "x": 0.5,
                "y": 0.5,
                "showarrow": False,
                "font": {"size": 14},
            }
        ],
    )
    return fig
