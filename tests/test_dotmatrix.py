import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from dash.exceptions import PreventUpdate

from guanaco.pages.matrix.plots import dotmatrix as dotmatrix_module
from guanaco.pages.matrix.plots.dotmatrix import plot_dot_matrix


def _dot_adata(n_obs=12):
    rng = np.random.default_rng(0)
    x = rng.random((n_obs, 3), dtype=np.float64).astype(np.float32)
    groups = np.array(["A", "B", "C"] * ((n_obs + 2) // 3))[:n_obs]
    obs = pd.DataFrame(
        {"cell_type": groups},
        index=[f"c{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])
    return AnnData(X=x, obs=obs, var=var)


def test_no_valid_genes_raises_prevent_update():
    with pytest.raises(PreventUpdate):
        plot_dot_matrix(_dot_adata(), ["Missing"], "cell_type", selected_labels=None)


def test_dotplot_builds_scatter_with_one_marker_per_group_gene_pair():
    adata = _dot_adata()
    fig = plot_dot_matrix(adata, ["GeneA", "GeneB"], "cell_type", selected_labels=None,
                          plot_type="dotplot")

    main_scatter = fig.data[0]
    # 3 groups x 2 genes = 6 dots in the primary trace.
    assert len(main_scatter.x) == 6
    assert main_scatter.mode == "markers"


def test_matrixplot_builds_heatmap_with_group_by_gene_shape():
    adata = _dot_adata()
    fig = plot_dot_matrix(adata, ["GeneA", "GeneB", "GeneC"], "cell_type", selected_labels=None,
                          plot_type="matrixplot")

    heatmap = fig.data[0]
    assert heatmap.type == "heatmap"
    # 3 groups (rows) x 3 genes (cols)
    assert np.asarray(heatmap.z).shape == (3, 3)


def test_colorbar_helper_is_shared_styling():
    bar = dotmatrix_module._colorbar("Z-score")
    assert bar["title"] == "Z-score"
    assert bar["len"] == 0.6
    assert bar["x"] == 0.98


def test_dendro_axis_layout_only_emits_requested_axes():
    none = dotmatrix_module._dendro_axis_layout(0.7, 1.0, False, False)
    assert none == {}

    both = dotmatrix_module._dendro_axis_layout(0.7, 0.86, True, True)
    assert set(both) == {"xaxis3", "yaxis3", "xaxis4", "yaxis4"}


def test_dendro_axis_layout_clamps_right_domain_for_matrixplot():
    # main_x_right close to the paper edge: matrixplot clamps the upper bound to 0.99.
    clamped = dotmatrix_module._dendro_axis_layout(0.96, 1.0, True, False, clamp_right_domain=True)
    assert clamped["xaxis3"]["domain"][1] == 0.99

    unclamped = dotmatrix_module._dendro_axis_layout(0.96, 1.0, True, False, clamp_right_domain=False)
    assert unclamped["xaxis3"]["domain"][1] == pytest.approx(1.0)
