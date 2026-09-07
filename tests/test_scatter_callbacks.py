import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad

from guanaco.pages.matrix.layouts.embedding_layout import (
    selectable_scatter_annotations,
)
from guanaco.utils.obs_utils import (
    OTHERS_LABEL,
    SELECTED_LABEL,
    selected_cell_view,
    selection_group_context,
)


def test_scatter_annotations_use_global_filter_rule_for_categorical_metadata():
    n_obs = 60
    obs = pd.DataFrame(
        {
            "condition": ["control", "treated"] * (n_obs // 2),
            "cell_type": pd.Categorical(["B", "T", "Mono"] * (n_obs // 3)),
            "cell_id": pd.Categorical([f"cell-{index}" for index in range(n_obs)]),
            "score": np.linspace(0, 1, n_obs),
        }
    )
    adata = AnnData(
        X=np.ones((n_obs, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1"]),
    )

    assert selectable_scatter_annotations(adata) == [
        "condition",
        "cell_type",
        "score",
    ]


def test_highlighting_groups_lasso_cells_within_current_filter_universe():
    adata = AnnData(
        X=np.ones((4, 1), dtype=np.float32),
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["G1"]),
    )

    plot_adata, groups = selection_group_context(
        adata,
        {
            "selected_cells": ["c2"],
            "universe_cells": ["c1", "c2", "c3"],
        },
    )

    assert plot_adata.obs_names.tolist() == ["c1", "c2", "c3"]
    assert groups.tolist() == [OTHERS_LABEL, SELECTED_LABEL, OTHERS_LABEL]
    assert "Selected / Others" not in adata.obs.columns


def test_selected_cells_compose_with_backed_filter_view(tmp_path):
    path = tmp_path / "backed-selection.h5ad"
    AnnData(
        X=np.arange(8, dtype=np.float32).reshape(4, 2),
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    ).write_h5ad(path)

    backed = read_h5ad(path, backed="r")
    try:
        globally_filtered = backed[[0, 2, 3]]
        selected = selected_cell_view(globally_filtered, ["c3", "c4"])

        assert selected.obs_names.tolist() == ["c3", "c4"]
        assert selected.isbacked
        np.testing.assert_array_equal(
            np.asarray(selected.X),
            np.array([[4, 5], [6, 7]], dtype=np.float32),
        )
    finally:
        backed.file.close()
