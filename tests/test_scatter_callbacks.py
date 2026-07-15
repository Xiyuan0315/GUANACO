import numpy as np
import pandas as pd
from anndata import AnnData

from guanaco.pages.matrix.callbacks.scatter_callbacks import _is_reset_relayout
from guanaco.pages.matrix.layouts.embedding_layout import (
    selectable_scatter_annotations,
)


def test_is_reset_relayout_true_for_xaxis_autorange():
    # Double-click "reset axes" emits an autorange relayout event.
    assert _is_reset_relayout({"xaxis.autorange": True}) is True


def test_is_reset_relayout_true_for_yaxis_autorange():
    assert _is_reset_relayout({"yaxis.autorange": True}) is True


def test_is_reset_relayout_false_for_zoom_range_event():
    # Zoom/pan emits explicit axis ranges, not autorange -> not a reset.
    zoom = {"xaxis.range[0]": 0.0, "xaxis.range[1]": 5.0}
    assert _is_reset_relayout(zoom) is False


def test_is_reset_relayout_false_for_empty_or_none():
    assert _is_reset_relayout(None) is False
    assert _is_reset_relayout({}) is False


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
