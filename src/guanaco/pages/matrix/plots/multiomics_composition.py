"""Observation-level modality coverage for multi-omics datasets."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from guanaco.utils.colors import (
    DEFAULT_DISCRETE_COLORMAP,
    resolve_discrete_palette,
)


MISSING_LABEL = "Missing"
MISSING_COLOR = "#E1E4E8"
SEPARATOR_COLOR = "#111111"
MAX_DISPLAY_COLUMNS = 5_000


def _metadata_labels(source, group_by: str | None) -> pd.Series:
    if not group_by:
        return pd.Series("All observations", index=source.coverage_ids, dtype="string")
    values = source.coverage_metadata(group_by)
    return values.astype("string").fillna(MISSING_LABEL)


def _ordered_categories(values: pd.Series) -> list[str]:
    if isinstance(values.dtype, pd.CategoricalDtype):
        present = set(values.dropna().astype(str))
        categories = [str(value) for value in values.cat.categories if str(value) in present]
        if MISSING_LABEL in set(values.astype(str)) and MISSING_LABEL not in categories:
            categories.append(MISSING_LABEL)
        return categories
    return sorted(values.dropna().astype(str).unique().tolist())


def normalize_required_modalities(source, required_modalities) -> list[str]:
    requested = [
        str(value)
        for value in (required_modalities or [])
        if str(value) in source.modalities
    ]
    return list(dict.fromkeys(requested)) or list(source.modalities)


def multiomics_coverage_summary(
    source,
    *,
    group_by: str | None = None,
    required_modalities=None,
) -> pd.DataFrame:
    """Count modality coverage and complete cases for every metadata group."""
    coverage = source.coverage_matrix
    groups = _metadata_labels(source, group_by)
    required = normalize_required_modalities(source, required_modalities)
    ready = coverage[required].all(axis=1)
    rows = []
    for group in _ordered_categories(groups):
        mask = groups.astype(str) == group
        n_group = int(mask.sum())
        row = {
            "Group": group,
            "N": n_group,
        }
        for modality in required:
            row[str(modality).upper()] = int(coverage.loc[mask, modality].sum())
        complete = int(ready.loc[mask].sum())
        row["Complete"] = complete
        row["Complete %"] = round(100 * complete / n_group, 1) if n_group else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def _app_palette(n_colors: int, color_config=None) -> list[str]:
    """Use the dataset/app palette, with GUANACO's default only as fallback."""
    palette = list(color_config or [])
    if not palette:
        palette = list(
            resolve_discrete_palette(DEFAULT_DISCRETE_COLORMAP, n_colors) or []
        )
    if not palette:
        palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44"]
    repeats = (max(1, n_colors) // len(palette)) + 1
    return (palette * repeats)[: max(1, n_colors)]


def _discrete_colorscale(colors: list[str]) -> list[list[float | str]]:
    if len(colors) == 1:
        return [[0.0, colors[0]], [1.0, colors[0]]]
    scale: list[list[float | str]] = []
    n_colors = len(colors)
    for index, color in enumerate(colors):
        start = index / n_colors
        end = (index + 1) / n_colors
        scale.extend([[start, color], [end, color]])
    return scale


def _display_order(
    coverage: pd.DataFrame,
    groups: pd.Series,
    modalities,
) -> pd.Index:
    if coverage.empty:
        return coverage.index
    category_order = {
        value: index for index, value in enumerate(_ordered_categories(groups))
    }
    pattern_weights = 2 ** np.arange(len(modalities) - 1, -1, -1)
    pattern_score = coverage[list(modalities)].to_numpy(dtype=np.int8) @ pattern_weights
    ordering = pd.DataFrame(
        {
            "group": groups.astype(str).map(category_order).fillna(len(category_order)),
            "pattern": pattern_score,
            "observation": coverage.index.astype(str),
        },
        index=coverage.index,
    )
    return ordering.sort_values(
        ["group", "pattern", "observation"],
        ascending=[True, False, True],
        kind="stable",
    ).index


def _downsample_order(order: pd.Index, limit: int) -> tuple[pd.Index, bool]:
    if len(order) <= limit:
        return order, False
    positions = np.linspace(0, len(order) - 1, num=limit, dtype=np.int64)
    return order.take(np.unique(positions)), True


def empty_multiomics_composition_figure(message: str) -> go.Figure:
    figure = go.Figure()
    figure.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
    )
    figure.update_layout(
        template="plotly_white",
        height=430,
        margin={"l": 40, "r": 20, "t": 40, "b": 30},
        xaxis={"visible": False},
        yaxis={"visible": False},
    )
    return figure


def _block_counts(
    coverage: pd.DataFrame,
    groups: pd.Series,
    modalities,
) -> dict[tuple[str, str, bool], int]:
    group_values = groups.astype(str)
    counts = {}
    for modality in modalities:
        measured = coverage[modality].astype(bool)
        for group in _ordered_categories(groups):
            group_mask = group_values.eq(group)
            counts[(modality, group, True)] = int((group_mask & measured).sum())
            counts[(modality, group, False)] = int((group_mask & ~measured).sum())
    return counts


def build_multiomics_composition_figure(
    source,
    *,
    group_by: str | None = None,
    required_modalities=None,
    color_config=None,
    max_display_columns: int = MAX_DISPLAY_COLUMNS,
) -> go.Figure:
    """Build the observation-level modality coverage heatmap."""
    full_coverage = source.coverage_matrix
    full_groups = _metadata_labels(source, group_by)
    if full_coverage.empty:
        return empty_multiomics_composition_figure("No observations are available.")

    modalities = normalize_required_modalities(source, required_modalities)
    order = _display_order(full_coverage, full_groups, modalities)
    display_order, sampled = _downsample_order(order, max_display_columns)
    coverage = full_coverage.loc[display_order]
    groups = full_groups.loc[display_order]

    categories = _ordered_categories(full_groups) if group_by else []
    palette = _app_palette(len(modalities) + len(categories), color_config)
    modality_colors = {
        modality: palette[index] for index, modality in enumerate(modalities)
    }
    category_colors = {
        category: palette[len(modalities) + index]
        for index, category in enumerate(categories)
    }

    colors = [MISSING_COLOR]
    code_by_modality = {}
    for modality in modalities:
        code_by_modality[modality] = len(colors)
        colors.append(modality_colors[modality])
    code_by_category = {}
    for category in categories:
        code_by_category[category] = len(colors)
        colors.append(category_colors[category])

    row_labels = []
    z_rows = []
    hover_rows = []
    dimensions = source.modality_dimensions
    group_values = groups.astype(str)
    metadata_counts = full_groups.astype(str).value_counts().to_dict()
    block_counts = _block_counts(full_coverage, full_groups, modalities)

    if group_by:
        row_labels.append(str(group_by))
        z_rows.append([code_by_category[value] for value in group_values])
        hover_rows.append(
            [
                "<br>".join(
                    [
                        f"Metadata: {group_by}",
                        f"Group: {value}",
                        f"Number of samples: {metadata_counts[value]:,}",
                    ]
                )
                for value in group_values
            ]
        )

    for modality in modalities:
        n_obs, n_vars = dimensions[modality]
        modality_label = str(modality).upper()
        row_labels.append(f"{modality_label}  ·  N={n_obs:,}  ·  D={n_vars:,}")
        present = coverage[modality].to_numpy(dtype=bool)
        z_rows.append(
            np.where(present, code_by_modality[modality], 0).astype(int).tolist()
        )
        hover_rows.append(
            [
                "<br>".join(
                    [
                        f"Modality: {modality_label}",
                        f"Measured: {'Yes' if is_present else 'No'}",
                        *([f"{group_by}: {group_value}"] if group_by else []),
                        (
                            "Number of samples: "
                            f"{block_counts[(modality, group_value, bool(is_present))]:,}"
                        ),
                        f"Features: {n_vars:,}",
                    ]
                )
                for is_present, group_value in zip(
                    present,
                    group_values,
                    strict=True,
                )
            ]
        )

    x_positions = np.arange(len(display_order))
    figure = go.Figure(
        go.Heatmap(
            x=x_positions,
            y=row_labels,
            z=z_rows,
            customdata=hover_rows,
            colorscale=_discrete_colorscale(colors),
            zmin=-0.5,
            zmax=len(colors) - 0.5,
            showscale=False,
            hovertemplate="%{customdata}<extra></extra>",
            xgap=0,
            ygap=2,
        )
    )

    if group_by and len(categories) <= 20:
        for category in categories:
            figure.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker={
                        "size": 10,
                        "color": category_colors[category],
                        "symbol": "square",
                    },
                    name=category,
                    legendgroup="metadata",
                    hoverinfo="skip",
                )
            )

    if group_by:
        display_groups = group_values.to_numpy()
        boundaries = np.flatnonzero(display_groups[1:] != display_groups[:-1]) + 0.5
        for boundary in boundaries:
            figure.add_vline(
                x=float(boundary),
                line_width=2,
                line_color=SEPARATOR_COLOR,
            )

    subtitle = (
        f"Showing {len(display_order):,} of {len(order):,} observations"
        if sampled
        else f"N={len(order):,} observations"
    )
    figure.update_layout(
        template="plotly_white",
        title={
            "text": f"Multi-omics coverage<br><sup>{subtitle}</sup>",
            "x": 0.5,
            "xanchor": "center",
        },
        height=max(410, 150 + 66 * len(row_labels)),
        margin={"l": 210, "r": 40, "t": 75, "b": 25},
        hovermode="closest",
        dragmode="zoom",
        plot_bgcolor="#FFFFFF",
        legend={
            "title": {"text": str(group_by) if group_by else ""},
            "orientation": "h",
            "y": -0.12,
            "x": 0,
        },
        xaxis={
            "title": None,
            "showticklabels": False,
            "showgrid": False,
            "zeroline": False,
            "range": [-0.5, len(display_order) - 0.5],
        },
        yaxis={
            "autorange": "reversed",
            "fixedrange": True,
            "showgrid": False,
            "zeroline": False,
            "automargin": True,
        },
        uirevision=(
            f"multiomics-composition::{group_by or 'none'}::"
            f"{','.join(modalities)}"
        ),
    )
    return figure
