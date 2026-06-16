import numpy as np
import pandas as pd
import pytest
import plotly.graph_objects as go
from anndata import AnnData
from dash.exceptions import PreventUpdate

from guanaco.pages.matrix.plots import violin1 as violin1_module
from guanaco.pages.matrix.plots.violin1 import plot_violin1
from guanaco.utils.gene_extraction_utils import clear_gene_cache


def _violin_adata():
    x = np.array(
        [
            [0.0, 1.0, 5.0],
            [2.0, 0.0, 5.0],
            [3.0, 2.0, 5.0],
            [4.0, 3.0, 5.0],
        ],
        dtype=np.float32,
    )
    obs = pd.DataFrame(
        {"cell_type": ["A", "B", "A", "B"]},
        index=["c1", "c2", "c3", "c4"],
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneFlat"])
    return AnnData(X=x, obs=obs, var=var)


@pytest.fixture(autouse=True)
def _clear_violin_cache():
    violin1_module._violin_data_cache.clear()
    clear_gene_cache()
    yield
    violin1_module._violin_data_cache.clear()
    clear_gene_cache()


def test_violin_empty_inputs_prevent_update():
    adata = _violin_adata()

    with pytest.raises(PreventUpdate):
        plot_violin1(adata, [], "cell_type", labels=["A"])

    with pytest.raises(PreventUpdate):
        plot_violin1(adata, ["GeneA"], "cell_type", labels=[])


def test_violin_missing_gene_or_empty_label_filter_returns_empty_figure():
    adata = _violin_adata()

    missing_gene_fig = plot_violin1(adata, ["Missing"], "cell_type", labels=["A"])
    empty_label_fig = plot_violin1(adata, ["GeneA"], "cell_type", labels=["Missing"])

    assert isinstance(missing_gene_fig, go.Figure)
    assert isinstance(empty_label_fig, go.Figure)
    assert len(missing_gene_fig.data) == 0
    assert len(empty_label_fig.data) == 0


def test_violin_filters_labels_and_keeps_gene_labels():
    fig = plot_violin1(
        _violin_adata(),
        ["GeneA", "GeneB"],
        "cell_type",
        labels=["B"],
        show_box=True,
        groupby_label_color_map={"A": "#E69F00", "B": "#56B4E9"},
    )

    # One violin per gene for label B; the other traces are invisible anchors that
    # force each row's right-hand numeric axis to render.
    violins = [trace for trace in fig.data if trace.type == "violin"]
    assert [trace.name for trace in violins] == ["B", "B"]
    assert all(trace.box.visible is True for trace in violins)
    # Gene names are the left y-axis tick labels (heatmap-style, automargin-sized),
    # not annotations -- so long names can't be clipped.
    gene_labels = [ax.ticktext[0] for ax in fig.select_yaxes() if ax.ticktext]
    assert gene_labels == ["GeneA", "GeneB"]
    np.testing.assert_allclose(np.asarray(violins[0].y), [2.0, 4.0])


def test_violin_log_transform_is_applied_before_plotting():
    fig = plot_violin1(
        _violin_adata(),
        ["GeneA"],
        "cell_type",
        labels=["A"],
        transformation="log1p",
    )

    np.testing.assert_allclose(np.asarray(fig.data[0].y), np.log1p([0.0, 3.0]))


def test_violin_constant_gene_gets_padded_y_range():
    fig = plot_violin1(
        _violin_adata(),
        ["GeneFlat"],
        "cell_type",
        labels=["A", "B"],
    )

    assert tuple(fig.layout.yaxis.range) == pytest.approx((4.75, 5.25))


def test_violin_cache_limit_keeps_most_recent_fifty_entries():
    for i in range(51):
        violin1_module._cache_violin_data(f"k{i}", {"i": i})

    assert len(violin1_module._violin_data_cache) == 50
    assert "k0" not in violin1_module._violin_data_cache
    assert violin1_module._violin_data_cache["k50"] == {"i": 50}
