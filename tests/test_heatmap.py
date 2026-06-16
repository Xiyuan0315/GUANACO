import numpy as np
import pandas as pd
from anndata import AnnData

from guanaco.pages.matrix.plots import heatmap as heatmap_module
from guanaco.pages.matrix.plots.heatmap import plot_heatmap2_continuous, plot_unified_heatmap


PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def _heatmap_adata(n_obs=6):
    x = np.arange(n_obs * 3, dtype=np.float32).reshape(n_obs, 3)
    groups = np.array(["A", "B"] * ((n_obs + 1) // 2))[:n_obs]
    batches = np.array(["x", "y", "x"] * ((n_obs + 2) // 3))[:n_obs]
    obs = pd.DataFrame(
        {
            "cell_type": groups,
            "batch": batches,
            "score": np.linspace(0, 1, n_obs, dtype=np.float32),
        },
        index=[f"c{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])
    return AnnData(X=x, obs=obs, var=var)


def test_heatmap_no_valid_genes_returns_message_figure():
    fig = plot_unified_heatmap(
        _heatmap_adata(),
        ["Missing"],
        "cell_type",
        color_config=PALETTE,
    )

    assert fig.layout.height == 400
    assert fig.layout.annotations[0].text == "No valid genes found in the dataset"


def test_categorical_heatmap_filters_labels_and_counts_original_cells():
    fig = plot_unified_heatmap(
        _heatmap_adata(),
        ["GeneA", "GeneB"],
        "cell_type",
        labels=["B"],
        color_config=PALETTE,
    )

    heatmap = fig.data[0]
    primary_bar = fig.data[1]
    assert np.asarray(heatmap.z).shape == (2, 3)
    assert primary_bar.name == "B"
    assert primary_bar.x == (3,)
    assert fig.layout.xaxis.range == (0, 3)


def test_binned_categorical_heatmap_uses_weighted_column_centers():
    fig = plot_unified_heatmap(
        _heatmap_adata(n_obs=24),
        ["GeneA"],
        "cell_type",
        max_cells=4,
        n_bins=4,
        color_config=PALETTE,
    )

    heatmap = fig.data[0]
    assert np.asarray(heatmap.z).shape == (1, 4)
    assert np.asarray(heatmap.x).tolist() == [3.0, 9.0, 15.0, 21.0]
    assert fig.layout.xaxis.range == (0, 24)


def test_secondary_categorical_binning_restores_group_columns():
    adata = _heatmap_adata(n_obs=22)
    fig = plot_unified_heatmap(
        adata,
        ["GeneA"],
        "cell_type",
        groupby2="batch",
        max_cells=4,
        n_bins=4,
        color_config=PALETTE,
    )

    secondary_bar_names = [trace.name for trace in fig.data[3:]]
    assert set(secondary_bar_names).issubset({"x", "y"})
    assert len(fig.layout.annotations) >= 4


def test_continuous_secondary_heatmap_uses_shared_primary_legend():
    fig = plot_heatmap2_continuous(
        _heatmap_adata(),
        ["GeneA"],
        "cell_type",
        "score",
        color_config=PALETTE,
    )

    annotation_texts = [annotation.text for annotation in fig.layout.annotations]
    assert "<b>cell_type</b>" in annotation_texts
    assert any("A" in text for text in annotation_texts)
    assert any("B" in text for text in annotation_texts)


def test_category_colorscale_maps_each_category_to_a_step():
    colorscale = heatmap_module._category_colorscale(["A", "B"], {"A": "red", "B": "blue"})

    assert colorscale == [[0.0, "red"], [0.5, "red"], [0.5, "blue"], [1.0, "blue"]]
    assert heatmap_module._category_colorscale([], {}) == []
