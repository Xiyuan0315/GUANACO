import numpy as np
import pandas as pd
import pytest
import plotly.graph_objects as go
from anndata import AnnData
from dash.exceptions import PreventUpdate

from guanaco.pages.matrix.plots import violin1 as violin1_module
from guanaco.pages.matrix.plots.violin1 import plot_ridge, plot_violin1
from guanaco.pages.matrix.callbacks.violin_callbacks import (
    _ridge_gene_options,
    _violin1_graph_style,
)
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
    assert fig.layout.autosize is True
    assert fig.layout.width is None
    assert fig.layout.height == 400


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


def test_ridge_plot_uses_horizontal_group_densities_and_violin_colors():
    fig = plot_ridge(
        _violin_adata(),
        "GeneA",
        "cell_type",
        labels=["A", "B"],
        show_box=True,
        groupby_label_color_map={"A": "#E69F00", "B": "#56B4E9"},
    )

    ridges = [trace for trace in fig.data if trace.type == "violin"]
    assert [trace.name for trace in ridges] == ["A", "B"]
    assert all(trace.orientation == "h" for trace in ridges)
    assert all(trace.side == "positive" for trace in ridges)
    assert all(trace.box.visible is True for trace in ridges)
    np.testing.assert_allclose(np.asarray(ridges[0].x), [0.0, 3.0])
    np.testing.assert_allclose(np.asarray(ridges[1].x), [2.0, 4.0])
    assert set(ridges[0].y) == {"A"}
    assert set(ridges[1].y) == {"B"}
    assert ridges[0].fillcolor == "rgba(230,159,0,0.6)"
    assert ridges[1].fillcolor == "rgba(86,180,233,0.6)"
    # Plotly lays y-axis categories out bottom-to-top, so reversing the axis
    # array displays the included-group selection order from top-to-bottom.
    assert tuple(fig.layout.yaxis.categoryarray) == ("B", "A")
    assert fig.layout.title.text == "<b>GeneA</b>"
    assert fig.layout.showlegend is True
    assert fig.layout.legend.title.text == "cell_type"
    assert fig.layout.legend.itemclick == "toggle"
    assert fig.layout.legend.itemdoubleclick == "toggleothers"


def test_ridge_gene_options_default_to_first_selected_gene_and_preserve_switch():
    options, value = _ridge_gene_options(["GeneB", "GeneA", "GeneB"], None)
    assert options == [
        {"label": "GeneB", "value": "GeneB"},
        {"label": "GeneA", "value": "GeneA"},
    ]
    assert value == "GeneB"

    _, value = _ridge_gene_options(["GeneB", "GeneA"], "GeneA")
    assert value == "GeneA"


def test_violin_cache_limit_keeps_most_recent_fifty_entries():
    for i in range(51):
        violin1_module._cache_violin_data(f"k{i}", {"i": i})

    assert len(violin1_module._violin_data_cache) == 50
    assert "k0" not in violin1_module._violin_data_cache
    assert violin1_module._violin_data_cache["k50"] == {"i": 50}


def test_violin_graph_container_preserves_dynamic_height():
    assert _violin1_graph_style({"layout": {"height": 840}}) == {
        "width": "100%",
        "height": "840px",
        "minHeight": "400px",
    }
