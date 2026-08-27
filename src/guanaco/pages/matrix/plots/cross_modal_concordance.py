"""Analysis and Plotly figures for paired cross-modal feature concordance."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.graph_objects as go


MAX_DISPLAY_POINTS = 8_000
MIN_GROUP_CELLS = 30
GROUP_CORRELATION = "group_correlation"
RELATIVE_SKEW = "relative_skew"
SKEW_COLOR_LIMIT = 2.0
DISAGREEMENT_COLORSCALE = [
    [0.0, "#2166AC"],
    [0.5, "#D4D2CC"],
    [1.0, "#D55E00"],
]
CONCORDANCE_COLORSCALE = [
    [0.0, "#B2182B"],
    [0.5, "#D4D2CC"],
    [1.0, "#1B7837"],
]


@dataclass(frozen=True)
class ConcordanceAnalysis:
    valid_mask: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z_x: np.ndarray
    z_y: np.ndarray
    disagreement: np.ndarray
    spearman: float


@dataclass(frozen=True)
class GroupCorrelationAnalysis:
    groups: np.ndarray
    correlations: np.ndarray
    counts: np.ndarray


@dataclass(frozen=True)
class UnpairedGroupComparison:
    groups: np.ndarray
    mean_a: np.ndarray
    mean_b: np.ndarray
    z_a: np.ndarray
    z_b: np.ndarray
    counts_a: np.ndarray
    counts_b: np.ndarray
    spearman: float


def _correlation(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def analyze_concordance(x, y) -> ConcordanceAnalysis:
    """Measure global concordance and standardized relative modality skew."""
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    if len(x) != len(y):
        raise ValueError("Paired features must contain the same number of cells.")

    valid_mask = np.isfinite(x) & np.isfinite(y)
    paired_x = x[valid_mask]
    paired_y = y[valid_mask]
    if len(paired_x) < 3:
        raise ValueError("At least three finite paired cells are required.")

    x_std = float(np.std(paired_x))
    y_std = float(np.std(paired_y))
    if x_std <= np.finfo(float).eps or y_std <= np.finfo(float).eps:
        raise ValueError("Both selected features must vary across the selected cells.")

    z_x = (paired_x - float(np.mean(paired_x))) / x_std
    z_y = (paired_y - float(np.mean(paired_y))) / y_std
    relative_skew = z_x - z_y

    ranked_x = pd.Series(paired_x).rank(method="average").to_numpy()
    ranked_y = pd.Series(paired_y).rank(method="average").to_numpy()
    spearman = _correlation(ranked_x, ranked_y)

    return ConcordanceAnalysis(
        valid_mask=valid_mask,
        x=paired_x,
        y=paired_y,
        z_x=z_x,
        z_y=z_y,
        disagreement=relative_skew,
        spearman=spearman,
    )


def calculate_group_correlations(
    analysis: ConcordanceAnalysis,
    groups,
    *,
    min_cells: int = MIN_GROUP_CELLS,
) -> GroupCorrelationAnalysis:
    """Map each metadata group to its within-group Spearman correlation."""
    raw_groups = pd.Series(groups)
    if len(raw_groups) != len(analysis.valid_mask):
        raise ValueError("Metadata groups do not match the paired cells.")

    valid_groups = raw_groups[analysis.valid_mask].astype("string").reset_index(drop=True)
    labels = valid_groups.fillna("Missing").to_numpy(dtype=str)
    correlations = np.full(len(valid_groups), np.nan, dtype=np.float64)
    counts = np.zeros(len(valid_groups), dtype=np.int64)

    for group in valid_groups.dropna().unique():
        mask = valid_groups.eq(group).fillna(False).to_numpy(dtype=bool)
        count = int(mask.sum())
        counts[mask] = count
        if count < min_cells:
            continue
        ranked_x = pd.Series(analysis.x[mask]).rank(method="average").to_numpy()
        ranked_y = pd.Series(analysis.y[mask]).rank(method="average").to_numpy()
        correlations[mask] = _correlation(ranked_x, ranked_y)

    if not np.any(np.isfinite(correlations)):
        raise ValueError(
            f"No groups contain at least {min_cells} cells with variable paired "
            "features."
        )
    return GroupCorrelationAnalysis(labels, correlations, counts)


def analyze_unpaired_group_comparison(
    values_a,
    groups_a,
    values_b,
    groups_b,
) -> UnpairedGroupComparison:
    """Compare independent modalities after aggregating shared metadata groups."""

    def summarize(values, groups):
        frame = pd.DataFrame(
            {
                "value": np.asarray(values, dtype=np.float64).reshape(-1),
                "group": pd.Series(
                    np.asarray(groups).reshape(-1),
                    dtype="string",
                ),
            }
        )
        if len(frame["value"]) != len(frame["group"]):
            raise ValueError("Feature values and metadata groups do not align.")
        frame = frame[np.isfinite(frame["value"]) & frame["group"].notna()]
        frame["group"] = frame["group"].astype(str)
        return frame.groupby("group", sort=False)["value"].agg(["mean", "count"])

    summary_a = summarize(values_a, groups_a)
    summary_b = summarize(values_b, groups_b)
    shared = [
        group for group in summary_a.index if group in set(summary_b.index)
    ]
    if len(shared) < 2:
        raise ValueError(
            "The selected metadata columns need at least two shared group labels."
        )

    mean_a = summary_a.loc[shared, "mean"].to_numpy(dtype=np.float64)
    mean_b = summary_b.loc[shared, "mean"].to_numpy(dtype=np.float64)
    std_a = float(np.std(mean_a))
    std_b = float(np.std(mean_b))
    if std_a <= np.finfo(float).eps or std_b <= np.finfo(float).eps:
        raise ValueError(
            "Both feature profiles must vary across the shared metadata groups."
        )
    z_a = (mean_a - float(np.mean(mean_a))) / std_a
    z_b = (mean_b - float(np.mean(mean_b))) / std_b
    ranked_a = pd.Series(mean_a).rank(method="average").to_numpy()
    ranked_b = pd.Series(mean_b).rank(method="average").to_numpy()
    return UnpairedGroupComparison(
        groups=np.asarray(shared, dtype=str),
        mean_a=mean_a,
        mean_b=mean_b,
        z_a=z_a,
        z_b=z_b,
        counts_a=summary_a.loc[shared, "count"].to_numpy(dtype=np.int64),
        counts_b=summary_b.loc[shared, "count"].to_numpy(dtype=np.int64),
        spearman=_correlation(ranked_a, ranked_b),
    )


def empty_concordance_figure(message: str, *, height: int = 520) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font={"color": "#6c757d", "size": 14},
    )
    fig.update_layout(
        template="plotly_white",
        height=height,
        margin={"l": 30, "r": 30, "t": 30, "b": 30},
        xaxis={"visible": False},
        yaxis={"visible": False},
    )
    return fig


def _display_indices(size: int) -> np.ndarray:
    if size <= MAX_DISPLAY_POINTS:
        return np.arange(size)
    rng = np.random.default_rng(0)
    return np.sort(rng.choice(size, MAX_DISPLAY_POINTS, replace=False))


def _view_values(analysis, view_mode, group_correlation):
    if view_mode == RELATIVE_SKEW:
        return (
            analysis.disagreement,
            DISAGREEMENT_COLORSCALE,
            -SKEW_COLOR_LIMIT,
            SKEW_COLOR_LIMIT,
        )
    if view_mode == GROUP_CORRELATION and group_correlation is not None:
        values = np.asarray(group_correlation.correlations, dtype=np.float64)
        if len(values) != len(analysis.x):
            raise ValueError("Group correlations do not match the paired cells.")
        return values, CONCORDANCE_COLORSCALE, -1.0, 1.0
    raise ValueError(f"Unsupported concordance view: {view_mode}")


def _point_selection_style():
    return {
        "selected": {"marker": {"opacity": 1.0, "size": 6}},
        "unselected": {"marker": {"opacity": 0.08}},
    }


def build_concordance_embedding(
    coordinates,
    analysis: ConcordanceAnalysis,
    cell_ids,
    feature_a: str,
    feature_b: str,
    embedding: str,
    *,
    view_mode: str = RELATIVE_SKEW,
    group_correlation: GroupCorrelationAnalysis | None = None,
    group_by: str | None = None,
) -> go.Figure:
    coords = np.asarray(coordinates)[analysis.valid_mask, :2]
    cell_ids = np.asarray(cell_ids, dtype=str)[analysis.valid_mask]
    score_values, colorscale, cmin, cmax = _view_values(
        analysis, view_mode, group_correlation
    )
    marker_values = np.nan_to_num(score_values, nan=0.0)
    selected = _display_indices(len(score_values))
    customdata = cell_ids[selected].astype(object).reshape(-1, 1)

    if view_mode == RELATIVE_SKEW:
        colorbar = {
            "title": "Relative skew",
            "tickmode": "array",
            "tickvals": [cmin, 0, cmax],
            "ticktext": ["B-skewed", "0", "A-skewed"],
        }
    else:
        if group_correlation is None or not group_by:
            raise ValueError("Select categorical metadata for Correlation view.")
        colorbar = {
            "title": "Within-group<br>Spearman ρ",
            "tickmode": "array",
            "tickvals": [-1, 0, 1],
            "ticktext": ["−1", "0", "+1"],
        }

    fig = go.Figure(
        go.Scattergl(
            x=coords[selected, 0],
            y=coords[selected, 1],
            mode="markers",
            customdata=customdata,
            marker={
                "size": 4,
                "opacity": 0.75,
                "color": marker_values[selected],
                "colorscale": colorscale,
                "cmin": cmin,
                "cmax": cmax,
                "cmid": 0,
                "colorbar": colorbar,
            },
            hoverinfo="skip",
            selectedpoints=None,
            **_point_selection_style(),
        )
    )
    fig.update_layout(
        template="plotly_white",
        meta={"spearman": analysis.spearman},
        dragmode="lasso",
        height=520,
        margin={"l": 55, "r": 80, "t": 25, "b": 50},
        xaxis={"title": f"{embedding} 1", "showgrid": False, "zeroline": False},
        yaxis={
            "title": f"{embedding} 2",
            "showgrid": False,
            "zeroline": False,
            "scaleanchor": "x",
            "scaleratio": 1,
        },
    )
    return fig


def build_disagreement_embedding(
    coordinates,
    analysis: ConcordanceAnalysis,
    cell_ids,
    feature_a: str,
    feature_b: str,
    embedding: str,
) -> go.Figure:
    """Backward-compatible relative-skew embedding builder."""
    return build_concordance_embedding(
        coordinates,
        analysis,
        cell_ids,
        feature_a,
        feature_b,
        embedding,
        view_mode=RELATIVE_SKEW,
    )


def _linear_fit_line(x, y):
    """Return the endpoints of the least-squares line through paired values."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    centered_x = x - float(np.mean(x))
    slope = float(np.sum(centered_x * (y - float(np.mean(y))))) / float(
        np.sum(centered_x * centered_x)
    )
    line_x = np.asarray([np.min(x), np.max(x)], dtype=np.float64)
    intercept = float(np.mean(y)) - slope * float(np.mean(x))
    return line_x, intercept + slope * line_x


def build_feature_scatter(
    analysis: ConcordanceAnalysis,
    cell_ids,
    feature_a: str,
    feature_b: str,
    *,
    view_mode: str = RELATIVE_SKEW,
    group_correlation: GroupCorrelationAnalysis | None = None,
    group_by: str | None = None,
) -> go.Figure:
    """Build the explanatory feature-pair scatter for the active score."""
    cell_ids = np.asarray(cell_ids, dtype=str)[analysis.valid_mask]
    score_values, colorscale, cmin, cmax = _view_values(
        analysis, view_mode, group_correlation
    )
    marker_values = np.nan_to_num(score_values, nan=0.0)
    selected = _display_indices(len(score_values))
    customdata = cell_ids[selected].astype(object).reshape(-1, 1)

    if view_mode == RELATIVE_SKEW:
        x_values = analysis.z_x
        y_values = analysis.z_y
        x_title = f"{feature_a} (z-score)"
        y_title = f"{feature_b} (z-score)"
    else:
        if group_correlation is None or not group_by:
            raise ValueError("Select categorical metadata for Correlation view.")
        x_values = analysis.x
        y_values = analysis.y
        x_title = f"{feature_a} expression"
        y_title = f"{feature_b} expression"

    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=x_values[selected],
            y=y_values[selected],
            mode="markers",
            customdata=customdata,
            marker={
                "size": 4,
                "opacity": 0.6,
                "color": marker_values[selected],
                "colorscale": colorscale,
                "cmin": cmin,
                "cmax": cmax,
                "cmid": 0,
                "showscale": False,
            },
            hoverinfo="skip",
            selectedpoints=None,
            name="Cells",
            **_point_selection_style(),
        )
    )

    if view_mode == RELATIVE_SKEW:
        line_min = float(min(np.min(analysis.z_x), np.min(analysis.z_y)))
        line_max = float(max(np.max(analysis.z_x), np.max(analysis.z_y)))
        fig.add_trace(
            go.Scatter(
                x=[line_min, line_max],
                y=[line_min, line_max],
                mode="lines",
                line={"color": "#6b6a64", "dash": "dash", "width": 1.5},
                hoverinfo="skip",
                name="Equal standardized expression",
            )
        )
    else:
        line_x, line_y = _linear_fit_line(analysis.x, analysis.y)
        fig.add_trace(
            go.Scatter(
                x=line_x,
                y=line_y,
                mode="lines",
                line={"color": "#4e5d6c", "width": 2},
                hoverinfo="skip",
                name="Linear fit",
            )
        )
        if np.isfinite(analysis.spearman):
            fig.add_annotation(
                text=f"Overall Spearman ρ = {analysis.spearman:.3f}",
                x=0.98,
                y=0.98,
                xref="paper",
                yref="paper",
                xanchor="right",
                yanchor="top",
                showarrow=False,
                font={"color": "#343a40", "size": 15},
            )

    fig.update_layout(
        template="plotly_white",
        meta={"spearman": analysis.spearman},
        dragmode="lasso",
        height=520,
        margin={"l": 70, "r": 25, "t": 25, "b": 65},
        xaxis={"title": x_title, "showgrid": False, "zeroline": True},
        yaxis={"title": y_title, "showgrid": False, "zeroline": True},
        showlegend=False,
    )
    return fig


def build_group_summary_heatmap(
    analysis: ConcordanceAnalysis,
    group_by: str,
    *,
    view_mode: str,
    groups,
    group_correlation: GroupCorrelationAnalysis | None = None,
) -> go.Figure:
    """Summarize the active concordance score for each metadata group."""
    raw_groups = pd.Series(groups)
    if len(raw_groups) != len(analysis.valid_mask):
        raise ValueError("Metadata groups do not match the paired cells.")
    valid_groups = (
        raw_groups[analysis.valid_mask]
        .astype("string")
        .fillna("Missing")
        .reset_index(drop=True)
        .to_numpy(dtype=str)
    )

    if view_mode == RELATIVE_SKEW:
        point_scores = analysis.disagreement
        column_label = "Mean relative skew"
        colorscale = DISAGREEMENT_COLORSCALE
        score_min, score_max = -SKEW_COLOR_LIMIT, SKEW_COLOR_LIMIT
        value_format = "+.2f"
    elif view_mode == GROUP_CORRELATION and group_correlation is not None:
        if len(group_correlation.groups) != len(valid_groups):
            raise ValueError("Group correlations do not match the paired cells.")
        point_scores = group_correlation.correlations
        column_label = "Spearman ρ"
        colorscale = CONCORDANCE_COLORSCALE
        score_min, score_max = -1.0, 1.0
        value_format = ".2f"
    else:
        raise ValueError("Group scores are unavailable for the selected view.")

    labels = []
    scores = []
    for group in pd.unique(valid_groups):
        mask = valid_groups == group
        finite_scores = point_scores[mask]
        finite_scores = finite_scores[np.isfinite(finite_scores)]
        labels.append(str(group))
        scores.append(float(np.mean(finite_scores)) if len(finite_scores) else np.nan)

    score_values = np.asarray(scores, dtype=np.float64)
    score_text = np.asarray(
        [
            [format(value, value_format) if np.isfinite(value) else "NA"]
            for value in score_values
        ]
    )
    fig = go.Figure(
        go.Heatmap(
            z=score_values[:, None],
            x=[column_label],
            y=labels,
            zmin=score_min,
            zmax=score_max,
            zmid=0,
            colorscale=colorscale,
            showscale=False,
            text=score_text,
            texttemplate="%{text}",
            hovertemplate=(
                f"{group_by}: %{{y}}<br>"
                f"{column_label}: %{{z:.2f}}<extra></extra>"
            ),
            hoverongaps=False,
        )
    )
    fig.update_layout(
        template="plotly_white",
        height=520,
        margin={"l": 80, "r": 10, "t": 50, "b": 25},
        showlegend=False,
    )
    fig.update_xaxes(side="top", showgrid=False)
    fig.update_yaxes(
        title=group_by,
        autorange="reversed",
        automargin=True,
        showgrid=False,
    )
    return fig


def build_unpaired_group_scatter(
    comparison: UnpairedGroupComparison,
    label_a: str,
    label_b: str,
) -> go.Figure:
    """Show concordance of independent modalities across matched groups."""
    customdata = np.column_stack(
        [
            comparison.mean_a,
            comparison.mean_b,
            comparison.counts_a,
            comparison.counts_b,
        ]
    )
    limit = float(
        max(
            1.0,
            np.max(np.abs(comparison.z_a)),
            np.max(np.abs(comparison.z_b)),
        )
    )
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=comparison.z_a,
            y=comparison.z_b,
            text=comparison.groups,
            customdata=customdata,
            mode="markers+text",
            textposition="top center",
            marker={"size": 10, "color": "#4e5d6c"},
            hovertemplate=(
                "%{text}<br>"
                f"{label_a}: %{{customdata[0]:.3g}} (n=%{{customdata[2]:.0f}})<br>"
                f"{label_b}: %{{customdata[1]:.3g}} (n=%{{customdata[3]:.0f}})"
                "<extra></extra>"
            ),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[-limit, limit],
            y=[-limit, limit],
            mode="lines",
            line={"color": "#a0a0a0", "dash": "dash", "width": 1.5},
            hoverinfo="skip",
            showlegend=False,
        )
    )
    correlation = (
        f"{comparison.spearman:.3f}"
        if np.isfinite(comparison.spearman)
        else "NA"
    )
    fig.add_annotation(
        text=f"Across-group Spearman ρ = {correlation}",
        x=0.98,
        y=0.98,
        xref="paper",
        yref="paper",
        xanchor="right",
        yanchor="top",
        showarrow=False,
        font={"color": "#343a40", "size": 14},
    )
    fig.update_layout(
        template="plotly_white",
        height=520,
        margin={"l": 70, "r": 25, "t": 35, "b": 70},
        xaxis={
            "title": f"{label_a} group mean (z-score)",
            "range": [-limit * 1.1, limit * 1.1],
            "showgrid": False,
            "zeroline": True,
        },
        yaxis={
            "title": f"{label_b} group mean (z-score)",
            "range": [-limit * 1.1, limit * 1.1],
            "showgrid": False,
            "zeroline": True,
            "scaleanchor": "x",
            "scaleratio": 1,
        },
    )
    return fig


def build_unpaired_group_heatmap(
    comparison: UnpairedGroupComparison,
    label_a: str,
    label_b: str,
) -> go.Figure:
    """Show standardized modality profiles and their group-level difference."""
    values = np.column_stack(
        [
            comparison.z_a,
            comparison.z_b,
            comparison.z_a - comparison.z_b,
        ]
    )
    hover = np.empty(values.shape, dtype=object)
    for index, group in enumerate(comparison.groups):
        hover[index, 0] = (
            f"{group}<br>{label_a}<br>Mean: {comparison.mean_a[index]:.3g}"
            f"<br>Cells: {comparison.counts_a[index]}"
        )
        hover[index, 1] = (
            f"{group}<br>{label_b}<br>Mean: {comparison.mean_b[index]:.3g}"
            f"<br>Cells: {comparison.counts_b[index]}"
        )
        hover[index, 2] = (
            f"{group}<br>Relative profile: "
            f"{comparison.z_a[index] - comparison.z_b[index]:+.2f}"
        )
    fig = go.Figure(
        go.Heatmap(
            z=values,
            x=[label_a, label_b, "A − B"],
            y=comparison.groups,
            zmin=-2,
            zmax=2,
            zmid=0,
            colorscale=DISAGREEMENT_COLORSCALE,
            text=hover,
            hovertemplate="%{text}<extra></extra>",
            colorbar={"title": "Group profile<br>z-score"},
        )
    )
    fig.update_layout(
        template="plotly_white",
        height=520,
        margin={"l": 90, "r": 70, "t": 35, "b": 80},
    )
    fig.update_xaxes(side="top", automargin=True)
    fig.update_yaxes(autorange="reversed", automargin=True)
    return fig
