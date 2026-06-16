import numpy as np
import pandas as pd
import pytest
import plotly.graph_objects as go
from anndata import AnnData

from guanaco.pages.matrix.plots import violin2 as violin2_module
from guanaco.pages.matrix.plots.violin2 import plot_violin2_new
from guanaco.utils.gene_extraction_utils import clear_gene_cache
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
    yield
    clear_gene_cache()


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
