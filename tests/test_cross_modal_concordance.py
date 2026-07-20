import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from dash import Dash, no_update

from guanaco.pages.matrix.callbacks.cross_modal_concordance_callbacks import (
    _pair_data,
    register_cross_modal_concordance_callbacks,
)
from guanaco.pages.matrix.plots.cross_modal_concordance import (
    GROUP_CORRELATION,
    MAX_DISPLAY_POINTS,
    RELATIVE_SKEW,
    SKEW_COLOR_LIMIT,
    analyze_concordance,
    build_concordance_embedding,
    build_disagreement_embedding,
    build_feature_scatter,
    build_group_summary_heatmap,
    calculate_group_correlations,
)
from guanaco.pages.visualizations.registry import EXPLORATION_WORKSPACE


def test_concordance_recovers_perfect_linear_relation():
    x = np.arange(50, dtype=float)
    result = analyze_concordance(x, 2 * x + 3)

    assert result.spearman == pytest.approx(1.0)
    assert result.disagreement == pytest.approx(np.zeros(50), abs=1e-12)


def test_disagreement_is_symmetric_when_features_are_swapped():
    x = np.asarray([0.0, 0.3, 0.8, 1.4, 2.5, 3.1, 4.7])
    y = np.asarray([0.2, 0.9, 0.5, 2.0, 1.8, 4.2, 3.9])

    forward = analyze_concordance(x, y)
    reverse = analyze_concordance(y, x)

    assert np.abs(forward.disagreement) == pytest.approx(
        np.abs(reverse.disagreement)
    )
    assert forward.disagreement == pytest.approx(-reverse.disagreement)


def test_disagreement_is_feature_a_z_score_minus_feature_b_z_score():
    x = np.asarray([0.0, 0.5, 2.0, 3.0, 8.0])
    y = np.asarray([4.0, 3.0, 2.0, 2.0, 1.0])

    result = analyze_concordance(x, y)

    assert result.disagreement == pytest.approx(result.z_x - result.z_y)


def test_concordance_rejects_a_constant_feature():
    with pytest.raises(ValueError, match="must vary"):
        analyze_concordance(np.ones(20), np.arange(20))


def test_pair_data_supports_selecting_the_same_feature_twice():
    cell_ids = pd.Index(["a", "b", "c", "d"])
    values = np.arange(4, dtype=float)
    source = SimpleNamespace(
        cell_ids=cell_ids,
        base_adata=SimpleNamespace(),
        score_features=lambda _features: values,
    )

    _frame, positions, result = _pair_data(
        source,
        "RNA · gene",
        "RNA · gene",
        None,
        None,
    )

    assert positions.tolist() == [0, 1, 2, 3]
    assert result.spearman == pytest.approx(1.0)
    assert result.disagreement == pytest.approx(np.zeros(4))


def test_embedding_uses_a_fixed_skew_color_range():
    result = analyze_concordance(
        np.asarray([0.0, 1.0, 2.0, 4.0]),
        np.asarray([0.0, 0.5, 3.0, 2.0]),
    )

    figure = build_disagreement_embedding(
        np.arange(8, dtype=float).reshape(4, 2),
        result,
        np.asarray(["a", "b", "c", "d"]),
        "RNA · gene",
        "PROTEIN · marker",
        "RNA · UMAP",
    )

    assert figure.data[0].marker.cmin == -SKEW_COLOR_LIMIT
    assert figure.data[0].marker.cmax == SKEW_COLOR_LIMIT
    assert figure.data[0].type == "scattergl"
    assert figure.layout.title.text is None
    assert figure.layout.height == 520


def test_cross_modal_point_plots_share_an_eight_thousand_cell_sample():
    size = MAX_DISPLAY_POINTS + 500
    x = np.arange(size, dtype=float)
    result = analyze_concordance(x, np.sqrt(x + 1))
    cell_ids = np.asarray([f"cell-{index}" for index in range(size)])

    embedding = build_concordance_embedding(
        np.column_stack([x, x]),
        result,
        cell_ids,
        "RNA · gene",
        "PROTEIN · marker",
        "RNA · UMAP",
    )
    relationship = build_feature_scatter(
        result,
        cell_ids,
        "RNA · gene",
        "PROTEIN · marker",
    )

    assert len(embedding.data[0].x) == 8_000
    assert len(relationship.data[0].x) == 8_000
    assert embedding.data[0].customdata[:, 0].tolist() == (
        relationship.data[0].customdata[:, 0].tolist()
    )


def test_group_correlation_maps_each_group_spearman_back_to_its_cells():
    increasing = np.arange(40, dtype=float)
    decreasing = np.arange(40, dtype=float)
    small = np.arange(10, dtype=float)
    x = np.concatenate([increasing, decreasing, small])
    y = np.concatenate([increasing, decreasing[::-1], small])
    groups = np.asarray(["A"] * 40 + ["B"] * 40 + ["small"] * 10)
    result = analyze_concordance(x, y)

    scores = calculate_group_correlations(result, groups, min_cells=30)

    assert scores.correlations[:40] == pytest.approx(np.ones(40))
    assert scores.correlations[40:80] == pytest.approx(-np.ones(40))
    assert np.all(np.isnan(scores.correlations[80:]))
    assert scores.groups[[0, 40, 80]].tolist() == ["A", "B", "small"]


def test_group_correlation_view_uses_raw_axes_and_group_values():
    x = np.tile(np.arange(40, dtype=float), 2)
    y = np.concatenate(
        [np.arange(40, dtype=float), np.arange(40, dtype=float)[::-1]]
    )
    groups = np.asarray(["A"] * 40 + ["B"] * 40)
    result = analyze_concordance(x, y)
    group_scores = calculate_group_correlations(result, groups)
    cell_ids = np.asarray([f"cell-{index}" for index in range(len(x))])

    embedding = build_concordance_embedding(
        np.arange(160, dtype=float).reshape(80, 2),
        result,
        cell_ids,
        "RNA · gene",
        "PROTEIN · marker",
        "RNA · UMAP",
        view_mode=GROUP_CORRELATION,
        group_correlation=group_scores,
        group_by="cell_type",
    )
    relationship = build_feature_scatter(
        result,
        cell_ids,
        "RNA · gene",
        "PROTEIN · marker",
        view_mode=GROUP_CORRELATION,
        group_correlation=group_scores,
        group_by="cell_type",
    )

    assert np.asarray(embedding.data[0].marker.color[:40]) == pytest.approx(
        np.ones(40)
    )
    assert np.asarray(embedding.data[0].marker.color[40:]) == pytest.approx(
        -np.ones(40)
    )
    assert embedding.data[0].hoverinfo == "skip"
    assert embedding.data[0].customdata.shape == (80, 1)
    assert embedding.data[0].customdata[0][0] == "cell-0"
    assert np.asarray(relationship.data[0].x) == pytest.approx(x)
    assert np.asarray(relationship.data[0].y) == pytest.approx(y)
    assert relationship.data[1].name == "Linear fit"
    assert len(relationship.data[1].x) == 2
    assert relationship.layout.xaxis.title.text.endswith("expression")
    assert relationship.data[0].hoverinfo == "skip"
    assert relationship.data[1].hoverinfo == "skip"
    assert [annotation.text for annotation in relationship.layout.annotations] == [
        f"Overall Spearman ρ = {result.spearman:.3f}"
    ]


def test_group_summary_heatmap_shows_the_active_score_on_hover():
    x = np.tile(np.arange(40, dtype=float), 2)
    y = np.concatenate(
        [np.arange(40, dtype=float), np.arange(40, dtype=float)[::-1]]
    )
    groups = np.asarray(["A"] * 40 + ["B"] * 40)
    result = analyze_concordance(x, y)
    group_scores = calculate_group_correlations(result, groups)

    correlation_figure = build_group_summary_heatmap(
        result,
        "cell_type",
        view_mode=GROUP_CORRELATION,
        groups=groups,
        group_correlation=group_scores,
    )
    skew_figure = build_group_summary_heatmap(
        result,
        "cell_type",
        view_mode=RELATIVE_SKEW,
        groups=groups,
    )

    assert len(correlation_figure.data) == 1
    assert list(correlation_figure.data[0].x) == ["Spearman ρ"]
    assert list(correlation_figure.data[0].y) == ["A", "B"]
    assert correlation_figure.data[0].text[:, 0].tolist() == ["1.00", "-1.00"]
    assert "cell_type: %{y}" in correlation_figure.data[0].hovertemplate

    assert len(skew_figure.data) == 1
    assert list(skew_figure.data[0].x) == ["Mean relative skew"]
    assert list(skew_figure.data[0].y) == ["A", "B"]
    expected_skew = [
        f"{np.mean(result.disagreement[:40]):+.2f}",
        f"{np.mean(result.disagreement[40:]):+.2f}",
    ]
    assert skew_figure.data[0].text[:, 0].tolist() == expected_skew
    assert "Mean relative skew: %{z:.2f}" in skew_figure.data[0].hovertemplate


def test_relative_skew_scatter_uses_z_scores_and_a_diagonal():
    result = analyze_concordance(
        np.asarray([0.0, 1.0, 2.0, 4.0]),
        np.asarray([0.0, 0.5, 3.0, 2.0]),
    )

    figure = build_feature_scatter(
        result,
        np.asarray(["a", "b", "c", "d"]),
        "RNA · gene",
        "PROTEIN · marker",
        view_mode=RELATIVE_SKEW,
    )

    assert np.asarray(figure.data[0].x) == pytest.approx(result.z_x)
    assert np.asarray(figure.data[0].y) == pytest.approx(result.z_y)
    assert np.asarray(figure.data[1].x) == pytest.approx(figure.data[1].y)
    assert figure.data[0].hoverinfo == "skip"
    assert not figure.layout.annotations


def test_embedding_only_updates_the_embedding_plot():
    app = Dash(__name__)
    register_cross_modal_concordance_callbacks(app, object(), "test")

    left_key = next(
        key
        for key in app.callback_map
        if "test-concordance-embedding-plot.figure" in key
    )
    right_key = next(
        key for key in app.callback_map if "test-concordance-feature-plot.figure" in key
    )
    summary_key = next(
        key
        for key in app.callback_map
        if "test-concordance-group-summary-plot.figure" in key
    )
    left_inputs = app.callback_map[left_key]["inputs"]
    right_inputs = app.callback_map[right_key]["inputs"]
    summary_inputs = app.callback_map[summary_key]["inputs"]
    left_states = app.callback_map[left_key]["state"]
    right_states = app.callback_map[right_key]["state"]
    embedding_id = "test-concordance-embedding"
    view_mode_id = "test-concordance-view-mode"
    group_by_id = "test-concordance-group-by"
    main_view_trigger_id = "test-concordance-main-view-trigger"
    committed_features_id = "test-concordance-committed-features"

    assert embedding_id in {item["id"] for item in left_inputs}
    assert embedding_id not in {item["id"] for item in right_inputs}
    assert group_by_id not in {item["id"] for item in left_inputs}
    assert group_by_id not in {item["id"] for item in right_inputs}
    assert view_mode_id not in {item["id"] for item in left_inputs}
    assert view_mode_id not in {item["id"] for item in right_inputs}
    assert main_view_trigger_id in {item["id"] for item in left_inputs}
    assert main_view_trigger_id in {item["id"] for item in right_inputs}
    assert committed_features_id in {item["id"] for item in left_inputs}
    assert committed_features_id in {item["id"] for item in right_inputs}
    assert committed_features_id in {item["id"] for item in summary_inputs}
    assert "test-concordance-feature-a" not in {
        item["id"] for item in (*left_inputs, *right_inputs, *summary_inputs)
    }
    assert view_mode_id in {item["id"] for item in left_states}
    assert view_mode_id in {item["id"] for item in right_states}
    assert group_by_id in {item["id"] for item in left_states}
    assert group_by_id in {item["id"] for item in right_states}
    assert group_by_id in {item["id"] for item in summary_inputs}
    assert main_view_trigger_id not in {item["id"] for item in summary_inputs}
    assert embedding_id not in {item["id"] for item in summary_inputs}
    assert not any("test-concordance-group-plot" in key for key in app.callback_map)


def test_feature_sets_are_only_committed_after_successful_scoring():
    valid_features = {"RNA · a1", "RNA · a2", "ADT · b1"}

    def score_features(features):
        if features[0].startswith("RNA"):
            return np.asarray([0.0, 1.0, 2.0, 3.0])
        return np.asarray([0.0, 1.0, 4.0, 9.0])

    source = SimpleNamespace(
        is_feature=valid_features.__contains__,
        score_features=score_features,
    )
    app = Dash(__name__)
    register_cross_modal_concordance_callbacks(app, source, "test")
    commit_key = next(
        key
        for key in app.callback_map
        if "test-concordance-committed-features.data" in key
        and "test-concordance-update-status.children" in key
    )
    commit_callback = app.callback_map[commit_key]["callback"].__wrapped__

    comparison, message = commit_callback(
        1,
        ["RNA · a1", "RNA · a2"],
        ["ADT · b1"],
    )
    assert comparison == {
        "feature_a": ["RNA · a1", "RNA · a2"],
        "feature_b": ["ADT · b1"],
    }
    assert message == ""

    comparison, message = commit_callback(2, ["missing"], ["ADT · b1"])
    assert comparison is no_update
    assert message.startswith("Could not update:")


def test_cross_modal_figures_are_reused_when_returning_to_the_tab():
    class Source:
        def __init__(self):
            self.cell_ids = pd.Index([f"cell-{index}" for index in range(60)])
            obs = pd.DataFrame(
                {"cell_type": ["A"] * 30 + ["B"] * 30}, index=self.cell_ids
            )
            self.base_adata = SimpleNamespace(
                obsm={"RNA · UMAP": np.arange(120, dtype=float).reshape(60, 2)},
                obs=obs,
                obs_names=self.cell_ids,
            )
            self.embedding_names = ["RNA · UMAP"]
            self.obs_names = ["cell_type"]
            self.score_calls = 0
            self._scores = {}
            self.values = {
                "RNA · gene": np.tile(np.arange(30, dtype=float), 2),
                "PROTEIN · marker": np.concatenate(
                    [np.arange(30, dtype=float), np.arange(30, dtype=float)[::-1]]
                ),
            }

        def score_features(self, features):
            key = tuple(sorted(features))
            if key not in self._scores:
                self.score_calls += 1
                self._scores[key] = np.mean(
                    np.column_stack([self.values[feature] for feature in key]), axis=1
                )
            return self._scores[key]

    source = Source()
    figures = {}

    def make_cache_key(kind, _adata, **kwargs):
        return f"{kind}:{json.dumps(kwargs, sort_keys=True)}"

    app = Dash(__name__)
    register_cross_modal_concordance_callbacks(
        app,
        source,
        "test",
        make_cache_key=make_cache_key,
        hash_list_signature=lambda values: values,
        cached_figure_get=figures.get,
        cached_figure_set=figures.__setitem__,
    )
    left_key = next(
        key
        for key in app.callback_map
        if "test-concordance-embedding-plot.figure" in key
    )
    right_key = next(
        key for key in app.callback_map if "test-concordance-feature-plot.figure" in key
    )
    summary_key = next(
        key
        for key in app.callback_map
        if "test-concordance-group-summary-plot.figure" in key
    )
    left_callback = app.callback_map[left_key]["callback"].__wrapped__
    right_callback = app.callback_map[right_key]["callback"].__wrapped__
    summary_callback = app.callback_map[summary_key]["callback"].__wrapped__
    comparison = {
        "feature_a": ["RNA · gene"],
        "feature_b": ["PROTEIN · marker"],
    }
    filtered_data = None
    cells_hash = {"len": 0, "hash": None}
    active_workspace = EXPLORATION_WORKSPACE
    active_tab = "cross-modal-concordance-tab"
    relative_trigger = {"view_mode": RELATIVE_SKEW, "group_by": None}
    correlation_trigger = {
        "view_mode": GROUP_CORRELATION,
        "group_by": "cell_type",
    }

    first_left = left_callback(
        comparison,
        relative_trigger,
        "RNA · UMAP",
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        RELATIVE_SKEW,
        None,
        None,
    )
    second_left = left_callback(
        comparison,
        relative_trigger,
        "RNA · UMAP",
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        RELATIVE_SKEW,
        None,
        None,
    )
    assert second_left is first_left
    assert source.score_calls == 2

    first_right = right_callback(
        comparison,
        relative_trigger,
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        RELATIVE_SKEW,
        None,
        None,
    )
    second_right = right_callback(
        comparison,
        relative_trigger,
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        RELATIVE_SKEW,
        None,
        None,
    )
    assert second_right is first_right
    assert source.score_calls == 2

    first_correlation = left_callback(
        comparison,
        correlation_trigger,
        "RNA · UMAP",
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        GROUP_CORRELATION,
        "cell_type",
        None,
    )
    second_correlation = left_callback(
        comparison,
        correlation_trigger,
        "RNA · UMAP",
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        GROUP_CORRELATION,
        "cell_type",
        None,
    )
    assert second_correlation is first_correlation

    right_correlation = right_callback(
        comparison,
        correlation_trigger,
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        GROUP_CORRELATION,
        "cell_type",
        None,
    )
    assert right_correlation.data[0].customdata[0][0] == "cell-0"

    summary = summary_callback(
        comparison,
        GROUP_CORRELATION,
        "cell_type",
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        None,
    )
    assert list(summary.data[0].y) == ["A", "B"]
    relative_summary = summary_callback(
        comparison,
        RELATIVE_SKEW,
        "cell_type",
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        None,
    )
    assert list(relative_summary.data[0].x) == ["Mean relative skew"]
    assert list(relative_summary.data[0].y) == ["A", "B"]
    assert source.score_calls == 2
