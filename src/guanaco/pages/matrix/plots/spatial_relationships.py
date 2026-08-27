"""Linked views of precomputed Squidpy spatial relationship statistics."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist

from guanaco.data.loader import obs_col


NHOOD_SUFFIX = "_nhood_enrichment"
CO_OCCURRENCE_SUFFIX = "_co_occurrence"


def _categories(adata, cluster_key: str) -> list[str]:
    values = obs_col(adata.obs, cluster_key)
    if not isinstance(values.dtype, pd.CategoricalDtype):
        return []
    return [str(value) for value in values.cat.categories]


def _result_mapping(adata, key: str) -> Mapping[str, Any] | None:
    result = adata.uns.get(key)
    return result if isinstance(result, Mapping) else None


def spatial_relationship_keys(adata) -> list[str]:
    """Return categorical obs keys with both compatible Squidpy results."""
    if adata is None:
        return []
    keys = []
    for uns_key in adata.uns:
        if not str(uns_key).endswith(NHOOD_SUFFIX):
            continue
        cluster_key = str(uns_key)[: -len(NHOOD_SUFFIX)]
        if cluster_key not in adata.obs.columns:
            continue
        categories = _categories(adata, cluster_key)
        n_clusters = len(categories)
        if n_clusters == 0:
            continue
        nhood = _result_mapping(adata, str(uns_key))
        occurrence = _result_mapping(
            adata, f"{cluster_key}{CO_OCCURRENCE_SUFFIX}"
        )
        if nhood is None or occurrence is None:
            continue
        zscore = np.asarray(nhood.get("zscore"))
        count = np.asarray(nhood.get("count"))
        occ = np.asarray(occurrence.get("occ"))
        interval = np.asarray(occurrence.get("interval"))
        if (
            zscore.shape == (n_clusters, n_clusters)
            and count.shape == zscore.shape
            and occ.ndim == 3
            and occ.shape[:2] == zscore.shape
            and interval.ndim == 1
            and occ.shape[2] == max(0, interval.size - 1)
        ):
            keys.append(cluster_key)
    return sorted(keys, key=str.casefold)


def load_spatial_relationship_payload(adata, cluster_key: str) -> dict[str, Any]:
    """Load and validate one pair of precomputed Squidpy results."""
    if cluster_key not in spatial_relationship_keys(adata):
        raise ValueError(
            f"No compatible neighborhood enrichment and co-occurrence results "
            f"were found for {cluster_key!r}."
        )
    nhood = adata.uns[f"{cluster_key}{NHOOD_SUFFIX}"]
    occurrence = adata.uns[f"{cluster_key}{CO_OCCURRENCE_SUFFIX}"]
    return {
        "cluster_key": cluster_key,
        "categories": _categories(adata, cluster_key),
        "zscore": np.asarray(nhood["zscore"], dtype=float),
        "count": np.asarray(nhood["count"], dtype=float),
        "occ": np.asarray(occurrence["occ"], dtype=float),
        "distance": np.asarray(occurrence["interval"], dtype=float)[1:],
    }


def default_spatial_pair(payload: Mapping[str, Any]) -> tuple[str, str]:
    """Choose the strongest positive off-diagonal neighborhood relationship."""
    zscore = np.asarray(payload["zscore"], dtype=float)
    categories = list(payload["categories"])
    candidates = zscore.copy()
    np.fill_diagonal(candidates, np.nan)
    if np.isfinite(candidates).any():
        index = int(np.nanargmax(candidates))
    else:
        index = 0
    row, column = np.unravel_index(index, zscore.shape)
    return categories[row], categories[column]


def spatial_pair_from_click(
    click_data: Mapping[str, Any] | None,
    payload: Mapping[str, Any],
) -> tuple[str, str]:
    """Resolve heatmap click data as (conditional group, observed group)."""
    categories = list(payload["categories"])
    if click_data:
        points = click_data.get("points") or []
        if points:
            conditional = str(points[0].get("y", ""))
            observed = str(points[0].get("x", ""))
            if conditional in categories and observed in categories:
                return conditional, observed
    return default_spatial_pair(payload)


def empty_spatial_figure(message: str) -> go.Figure:
    figure = go.Figure()
    figure.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font={"color": "#6c757d", "size": 14},
    )
    figure.update_xaxes(visible=False)
    figure.update_yaxes(visible=False)
    figure.update_layout(
        margin={"l": 30, "r": 20, "t": 50, "b": 30},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    return figure


def _cluster_columns(values: np.ndarray):
    """Return a stable hierarchical column order and its linkage matrix."""
    if values.shape[1] < 2:
        return np.arange(values.shape[1]), None
    profiles = np.nan_to_num(values.T, nan=0.0, posinf=0.0, neginf=0.0)
    distances = pdist(profiles, metric="correlation")
    if not np.isfinite(distances).all() or np.allclose(distances, 0):
        distances = pdist(profiles, metric="euclidean")
    if distances.size == 0 or np.allclose(distances, 0):
        return np.arange(values.shape[1]), None
    column_linkage = linkage(distances, method="average")
    return leaves_list(column_linkage), column_linkage


def _add_column_dendrogram(
    figure: go.Figure,
    column_linkage: np.ndarray | None,
) -> None:
    if column_linkage is None:
        return
    tree = dendrogram(column_linkage, no_plot=True)
    max_height = max(max(values) for values in tree["dcoord"]) or 1.0
    for x_values, y_values in zip(
        tree["icoord"],
        tree["dcoord"],
        strict=True,
    ):
        figure.add_trace(
            go.Scatter(
                x=[(value - 5.0) / 10.0 for value in x_values],
                y=[value / max_height for value in y_values],
                xaxis="x2",
                yaxis="y2",
                mode="lines",
                line={"color": "#555555", "width": 1},
                hoverinfo="skip",
                showlegend=False,
            )
        )


def plot_neighborhood_enrichment(
    payload: Mapping[str, Any],
    *,
    cluster_columns: bool = True,
) -> go.Figure:
    categories = list(payload["categories"])
    zscore = np.asarray(payload["zscore"], dtype=float)
    count = np.asarray(payload["count"], dtype=float)
    if cluster_columns:
        column_order, column_linkage = _cluster_columns(zscore)
    else:
        column_order = np.arange(len(categories))
        column_linkage = None
    ordered_columns = [categories[index] for index in column_order]
    zscore = zscore[:, column_order]
    count = count[:, column_order]
    finite = np.abs(zscore[np.isfinite(zscore)])
    color_limit = float(finite.max()) if finite.size else 1.0
    color_limit = max(color_limit, 1.0)

    figure = go.Figure(
        go.Heatmap(
            x=ordered_columns,
            y=categories,
            z=zscore,
            customdata=count,
            colorscale="RdBu_r",
            zmin=-color_limit,
            zmax=color_limit,
            zmid=0,
            colorbar={
                "title": "Z-score",
                "thickness": 14,
                "x": 0.91,
                "len": 0.76,
                "y": 0.4,
            },
            hovertemplate=(
                "Conditional: %{y}<br>"
                "Neighbor: %{x}<br>"
                "Z-score: %{z:.3f}<br>"
                "Observed edges: %{customdata:,.0f}<extra></extra>"
            ),
        )
    )
    show_dendrogram = column_linkage is not None
    if show_dendrogram:
        _add_column_dendrogram(figure, column_linkage)

    layout = {
        "title": {
            "text": "Neighborhood enrichment",
            "x": 0.02,
            "xanchor": "left",
        },
        "clickmode": "event+select",
        "margin": {
            "l": 135,
            "r": 35,
            "t": 72 if show_dendrogram else 62,
            "b": 115,
        },
        "paper_bgcolor": "white",
        "plot_bgcolor": "white",
        "xaxis": {
            "domain": [0, 0.88],
            "tickangle": -45,
            "title": "Neighbor group",
            "categoryorder": "array",
            "categoryarray": ordered_columns,
        },
        "yaxis": {
            "domain": [0, 0.8] if show_dendrogram else [0, 1],
            "autorange": "reversed",
            "title": "Conditional group",
        },
        "uirevision": (
            f"nhood-{payload['cluster_key']}-"
            f"{'clustered' if cluster_columns else 'original'}"
        ),
    }
    if show_dendrogram:
        layout.update(
            {
                "xaxis2": {
                    "domain": [0, 0.88],
                    "range": [-0.5, len(ordered_columns) - 0.5],
                    "showgrid": False,
                    "zeroline": False,
                    "showticklabels": False,
                    "fixedrange": True,
                },
                "yaxis2": {
                    "domain": [0.84, 0.97],
                    "range": [0, 1.05],
                    "showgrid": False,
                    "zeroline": False,
                    "showticklabels": False,
                    "fixedrange": True,
                },
            }
        )
    figure.update_layout(**layout)
    return figure


def plot_co_occurrence_group(
    payload: Mapping[str, Any],
    conditional: str,
    emphasized: str | None = None,
    colors: list[str] | None = None,
) -> go.Figure:
    """Plot every observed group around one conditional group.

    This mirrors ``squidpy.pl.co_occurrence(..., clusters=conditional)`` while
    allowing the heatmap-selected observed group to remain visually emphasized.
    """
    categories = list(payload["categories"])
    conditional_index = categories.index(conditional)
    distance = np.asarray(payload["distance"], dtype=float)
    occurrence = np.asarray(payload["occ"], dtype=float)
    palette = colors or [
        "#636EFA",
        "#EF553B",
        "#00CC96",
        "#AB63FA",
        "#FFA15A",
        "#19D3F3",
        "#FF6692",
        "#B6E880",
        "#FF97FF",
        "#FECB52",
    ]

    figure = go.Figure()
    for observed_index, observed in enumerate(categories):
        is_emphasized = observed == emphasized
        figure.add_trace(
            go.Scatter(
                x=distance,
                y=occurrence[conditional_index, observed_index, :],
                mode="lines",
                name=observed,
                line={
                    "color": palette[observed_index % len(palette)],
                    "width": 3.5 if is_emphasized else 1.8,
                },
                opacity=1.0 if is_emphasized else 0.72,
                hoverinfo="skip",
            )
        )
    figure.add_hline(
        y=1,
        line={"color": "#9a9a9a", "dash": "dash", "width": 1},
        annotation_text="Expected",
        annotation_position="bottom right",
    )
    figure.update_layout(
        title={
            "text": (
                f"<i>p</i>(group | {conditional}) / "
                "<i>p</i>(group)"
            ),
            "x": 0.02,
            "xanchor": "left",
        },
        margin={"l": 70, "r": 145, "t": 62, "b": 65},
        paper_bgcolor="white",
        plot_bgcolor="white",
        hovermode=False,
        legend={
            "title": payload["cluster_key"],
            "orientation": "v",
            "x": 1.01,
            "xanchor": "left",
            "y": 0.5,
            "yanchor": "middle",
        },
        uirevision=f"co-occurrence-{payload['cluster_key']}-{conditional}",
    )
    figure.update_xaxes(
        title="Distance (spatial coordinate units)",
        showgrid=False,
        zeroline=False,
    )
    figure.update_yaxes(
        title="Co-occurrence ratio",
        gridcolor="#ececec",
        zeroline=False,
    )
    return figure


def plot_co_occurrence_pair(
    payload: Mapping[str, Any],
    conditional: str,
    observed: str,
) -> go.Figure:
    """Backward-compatible wrapper for callers using the former pair API."""
    return plot_co_occurrence_group(
        payload,
        conditional,
        emphasized=observed,
    )
