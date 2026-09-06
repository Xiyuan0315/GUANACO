import numpy as np
import pandas as pd
import pytest
import plotly.graph_objects as go
from anndata import AnnData

from guanaco.pages.matrix.plots import violin2 as violin2_module
from guanaco.pages.matrix.plots.violin2 import plot_violin2_new
from guanaco.utils.gene_extraction_utils import clear_gene_cache
from guanaco.utils.obs_utils import SELECTION_GROUP, selection_group_values
from guanaco.widget import violin_grouped


def _violin2_adata(n_per_group=4):
    groups = []
    conditions = []
    expression = []
    for group, offset in [("A", 0.0), ("B", 10.0)]:
        for condition, condition_offset in [("ctrl", 0.0), ("stim", 2.0)]:
            for i in range(n_per_group):
                groups.append(group)
                conditions.append(condition)
                expression.append(offset + condition_offset + i)

    obs = pd.DataFrame(
        {
            "cell_type": groups,
            "condition": conditions,
        },
        index=[f"c{i}" for i in range(len(groups))],
    )
    x = np.asarray(expression, dtype=np.float32).reshape(-1, 1)
    var = pd.DataFrame(index=["GeneA"])
    return AnnData(X=x, obs=obs, var=var)


@pytest.fixture(autouse=True)
def _clear_gene_expression_cache():
    clear_gene_cache()
    violin2_module._statistics_cache.clear()
    yield
    clear_gene_cache()
    violin2_module._statistics_cache.clear()


def test_violin2_empty_label_filter_returns_message_figure():
    fig = plot_violin2_new(
        _violin2_adata(),
        "GeneA",
        "cell_type",
        None,
        "mode1",
        labels=["Missing"],
    )

    assert len(fig.data) == 0
    assert fig.layout.annotations[0].text == "No data available for selected filters"


def test_violin2_accepts_show_points_from_widget_wrapper():
    adata = _violin2_adata()

    direct_fig = plot_violin2_new(
        adata,
        "GeneA",
        "cell_type",
        None,
        "mode1",
        show_points=True,
        test_method="none",
    )
    widget_fig = violin_grouped(
        adata,
        "GeneA",
        "cell_type",
        show_points=True,
        show=False,
        return_fig=True,
    )

    assert all(trace.points == "all" for trace in direct_fig.data)
    assert all(trace.points == "all" for trace in widget_fig.data)


def test_violin2_downsampling_only_affects_drawn_traces(monkeypatch):
    adata = _violin2_adata(n_per_group=6)
    captured = {}

    def fake_downsample(df, group_cols):
        return df.groupby(group_cols, observed=True, sort=False).head(2)

    def fake_p_values(df, *args, **kwargs):
        captured["rows_used_for_stats"] = len(df)
        return {"overall": 0.01}

    monkeypatch.setattr(violin2_module, "_downsample_for_violins", fake_downsample)
    monkeypatch.setattr(violin2_module, "calculate_p_values_by_mode", fake_p_values)

    fig = plot_violin2_new(
        adata,
        "GeneA",
        "cell_type",
        None,
        "mode1",
        test_method="mwu-test",
    )

    assert captured["rows_used_for_stats"] == adata.n_obs
    assert sorted(len(trace.y) for trace in fig.data) == [2, 2]


@pytest.mark.parametrize("mode", ["mode1", "mode2", "mode3", "mode4"])
def test_style_changes_reuse_full_data_statistics(monkeypatch, mode):
    adata = _violin2_adata()
    calls = []

    def calculate(df, *args):
        calls.append((len(df), args))
        return {"overall": 0.01}

    monkeypatch.setattr(violin2_module, "calculate_p_values_by_mode", calculate)
    options = dict(
        key="GeneA", meta1="cell_type",
        meta2=None if mode == "mode1" else "condition", mode=mode,
    )
    first = plot_violin2_new(adata, **options, palette=["red", "blue"])
    restyled = plot_violin2_new(
        adata, **options, palette=["green", "orange"], show_box=True,
    )
    assert len(calls) == 1
    assert calls[0][0] == adata.n_obs
    assert first.data[0].fillcolor != restyled.data[0].fillcolor
    assert all(trace.box.visible for trace in restyled.data)
    assert [a.text for a in first.layout.annotations] == [
        a.text for a in restyled.layout.annotations
    ]


@pytest.mark.parametrize("mode", ["mode1", "mode2", "mode3"])
def test_cached_statistics_match_direct_calculation(mode):
    adata = _violin2_adata()
    meta2 = None if mode == "mode1" else "condition"
    df = violin2_module._build_expression_frame(
        adata, "GeneA", "cell_type", meta2, mode,
    )
    method, _ = violin2_module.determine_test_method(2, 2, mode)
    args = (df, "cell_type", meta2, mode, method)
    expected = violin2_module.calculate_p_values_by_mode(*args)
    assert violin2_module._cached_p_values_by_mode(*args) == expected
    assert violin2_module._cached_p_values_by_mode(*args) == expected


def test_statistics_cache_tracks_data_groups_categories_and_test(monkeypatch):
    df = pd.DataFrame({
        "Expression": [1., 2., 3., 4.],
        "group": pd.Categorical(["A", "A", "B", "B"]),
    })
    calls = []

    def calculate(frame, *args):
        calls.append(len(frame))
        return {"overall": 0.01, "model_summary": {"meta1_p": 0.01}}

    monkeypatch.setattr(violin2_module, "calculate_p_values_by_mode", calculate)

    def cached(frame, **kwargs):
        return violin2_module._cached_p_values_by_mode(
            frame, "group", None, "mode1", kwargs.pop("method", "mwu-test"),
            **kwargs,
        )

    result = cached(df)
    result["model_summary"]["meta1_p"] = 1.0
    assert cached(df.copy())["model_summary"]["meta1_p"] == 0.01
    assert len(calls) == 1
    cached(df.iloc[:3])  # changed cell selection
    changed = df.copy()
    changed.loc[0, "Expression"] = 10.0  # different layer/transformed values
    cached(changed)
    changed = df.copy()
    changed.loc[0, "group"] = "B"  # changed lasso/metadata groups
    cached(changed)
    changed = df.copy()
    changed["group"] = changed["group"].cat.reorder_categories(["B", "A"])
    cached(changed)  # same values, different model reference level
    cached(df, method="ttest")
    cached(df, labels=["A"])
    assert len(calls) == 7


def test_statistics_cache_is_bounded_and_does_not_retain_failed_fits(monkeypatch):
    monkeypatch.setattr(violin2_module, "_STATISTICS_CACHE_SIZE", 2)
    calls = []

    def calculate(df, *args):
        calls.append(1)
        return {"overall": 0.1}

    monkeypatch.setattr(violin2_module, "calculate_p_values_by_mode", calculate)
    frames = [pd.DataFrame({"Expression": [i], "group": ["A"]}) for i in range(3)]

    def cached(df):
        return violin2_module._cached_p_values_by_mode(
            df, "group", None, "mode1", "ttest",
        )

    for frame in frames:
        cached(frame)
    assert len(violin2_module._statistics_cache) == 2
    cached(frames[0])
    assert len(calls) == 4
    violin2_module._statistics_cache.clear()

    def failed(df, *args):
        calls.append(1)
        return {"error": "fit failed"}

    monkeypatch.setattr(violin2_module, "calculate_p_values_by_mode", failed)
    cached(frames[0])
    cached(frames[0])
    assert len(calls) == 6
    assert not violin2_module._statistics_cache


def test_violin2_zscore_layout_keeps_negative_values_visible():
    fig = plot_violin2_new(
        _violin2_adata(),
        "GeneA",
        "cell_type",
        None,
        "mode1",
        transformation="zscore",
        test_method="none",
    )

    y_range = tuple(fig.layout.yaxis.range)
    assert y_range[0] < 0
    assert y_range[1] > 0


def test_violin2_accepts_session_local_lasso_groups():
    adata = _violin2_adata()
    fig = plot_violin2_new(
        adata,
        "GeneA",
        SELECTION_GROUP,
        None,
        "mode1",
        test_method="none",
        group_values={
            SELECTION_GROUP: selection_group_values(
                adata,
                ["c1", "c2"],
            )
        },
    )

    assert {trace.name for trace in fig.data} == {"Selected", "Others"}


def test_model_summary_annotation_uses_public_helper_signature():
    fig = go.Figure()
    df = pd.DataFrame({"Expression": [0.0, 1.0]})

    violin2_module.add_p_value_annotations_new(
        fig,
        {"model_summary": {"meta1_p": 0.02, "meta2_p": 0.5}},
        df,
        "mode3",
        meta1="cell_type",
        meta2="condition",
    )

    assert "cell_type" in fig.layout.annotations[0].text
    assert "condition" in fig.layout.annotations[0].text


def test_downsample_for_violins_caps_each_drawn_group():
    df = pd.DataFrame(
        {
            "Expression": np.arange(20),
            "cell_type": ["A"] * 10 + ["B"] * 10,
            "condition": ["ctrl", "stim"] * 10,
        }
    )

    sampled = violin2_module._downsample_for_violins(
        df,
        ["cell_type", "condition"],
        cap=3,
        seed=1,
    )

    counts = sampled.groupby(["cell_type", "condition"], observed=True).size()
    assert set(counts.index) == set(df.groupby(["cell_type", "condition"], observed=True).size().index)
    assert counts.max() == 3
