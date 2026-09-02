import numpy as np
import pandas as pd
from anndata import AnnData

from guanaco.pages.matrix.plots import pseudotime as pseudotime_module
from guanaco.pages.matrix.plots.pseudotime import plot_genes_in_pseudotime


def _pseudotime_adata(n_obs=60):
    rng = np.random.default_rng(1)
    x = rng.random((n_obs, 3), dtype=np.float64).astype(np.float32)
    obs = pd.DataFrame(
        {
            "pseudotime": np.linspace(0, 1, n_obs, dtype=np.float32),
            "lineage": np.array(["L1", "L2"] * (n_obs // 2))[:n_obs],
        },
        index=[f"c{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])
    return AnnData(X=x, obs=obs, var=var)


def test_no_valid_genes_returns_message_figure():
    fig = plot_genes_in_pseudotime(_pseudotime_adata(), ["Missing"])

    assert fig.layout.height == 400
    assert fig.layout.annotations[0].text == "No valid genes found in the dataset"


def test_missing_pseudotime_key_returns_message_figure():
    fig = plot_genes_in_pseudotime(_pseudotime_adata(), ["GeneA"], pseudotime_key="absent_key")

    assert fig.layout.height == 400
    assert fig.layout.annotations[0].text == "'absent_key' not found in adata.obs"


def test_over_strict_min_expr_returns_no_cells_message():
    fig = plot_genes_in_pseudotime(_pseudotime_adata(), ["GeneA"], min_expr=1e9)

    assert fig.layout.annotations[0].text == "No cells pass the filtering criteria"


def test_basic_plot_builds_one_subplot_row_per_gene():
    fig = plot_genes_in_pseudotime(_pseudotime_adata(), ["GeneA", "GeneB"], min_expr=0.0)

    # subplot_titles carries one annotation per gene.
    titles = [a.text for a in fig.layout.annotations]
    assert "GeneA" in titles and "GeneB" in titles
    assert len(fig.data) > 0


def test_each_gene_trend_uses_its_own_expression_filter():
    adata = _pseudotime_adata()
    pseudotime = adata.obs["pseudotime"].to_numpy()
    adata.X[:, 0] = np.where(pseudotime < 0.5, 0.2, 1.0)
    adata.X[:, 1] = np.where(pseudotime < 0.5, 1.0, 0.0)

    single_gene = plot_genes_in_pseudotime(adata, ["GeneA"], min_expr=0.5)
    multiple_genes = plot_genes_in_pseudotime(adata, ["GeneA", "GeneB"], min_expr=0.5)
    single_trend = next(trace for trace in single_gene.data if trace.name == "Smoothed")
    multiple_trend = next(trace for trace in multiple_genes.data if trace.name == "Smoothed")

    np.testing.assert_allclose(single_trend.x, multiple_trend.x)
    np.testing.assert_allclose(single_trend.y, multiple_trend.y)


def test_message_figure_helper_centers_text():
    fig = pseudotime_module._message_figure("hello")

    assert fig.layout.height == 400
    assert fig.layout.annotations[0].text == "hello"
    assert fig.layout.annotations[0].x == 0.5
