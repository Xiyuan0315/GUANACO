import numpy as np
import pandas as pd
from anndata import AnnData

from guanaco.pages.matrix.callbacks.stacked_bar_callbacks import (
    _composition_control_state,
)
from guanaco.pages.matrix.plots.stacked_bar import (
    plot_composition_hierarchy,
    plot_stacked_bar,
)
from guanaco.utils.obs_utils import SELECTION_GROUP, selection_group_values


def _adata():
    obs = pd.DataFrame(
        {
            "lineage": pd.Categorical(["Lymphoid", "Lymphoid", "Myeloid", "Myeloid"]),
            "cell_type": pd.Categorical(["B", "T", "Mono", "Mono"]),
        },
        index=[f"cell-{index}" for index in range(4)],
    )
    return AnnData(
        X=np.ones((4, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1"]),
    )


def test_hierarchy_builds_parent_and_child_nodes_with_counts():
    fig = plot_composition_hierarchy(
        "lineage",
        "cell_type",
        _adata(),
        parent_color_map={"Lymphoid": "red", "Myeloid": "blue"},
        child_color_map={"B": "gold", "T": "green", "Mono": "purple"},
        parent_order=["Myeloid", "Lymphoid"],
    )

    trace = fig.data[0]
    values_by_id = dict(zip(trace.ids, trace.values, strict=True))
    colors_by_id = dict(zip(trace.ids, trace.marker.colors, strict=True))

    assert trace.type == "icicle"
    assert list(trace.ids[:3]) == [
        "all-cells",
        "lineage::Myeloid",
        "lineage::Myeloid/cell_type::Mono",
    ]
    assert values_by_id["all-cells"] == 4
    assert values_by_id["lineage::Myeloid"] == 2
    assert values_by_id["lineage::Myeloid/cell_type::Mono"] == 2
    assert colors_by_id["lineage::Myeloid"] == "blue"
    assert colors_by_id["lineage::Myeloid/cell_type::Mono"] == "purple"


def test_same_metadata_uses_a_single_hierarchy_level():
    fig = plot_composition_hierarchy(
        "lineage",
        "lineage",
        _adata(),
    )

    assert list(fig.data[0].ids) == [
        "all-cells",
        "lineage::Lymphoid",
        "lineage::Myeloid",
    ]


def test_hierarchy_controls_use_hierarchy_terms_and_hide_bar_value():
    (
        parent_label,
        child_label,
        value_style,
        tooltip,
        da_controls_style,
        da_panel_style,
    ) = _composition_control_state("hierarchy")

    assert parent_label == "Parent level:"
    assert child_label == "Child level:"
    assert value_style == {"display": "none"}
    assert "nested beneath each parent" in tooltip
    assert da_controls_style == {"display": "none"}
    assert da_panel_style == {"display": "none"}


def test_bar_controls_keep_the_original_terms_and_value_control():
    (
        parent_label,
        child_label,
        value_style,
        _tooltip,
        da_controls_style,
        da_panel_style,
    ) = _composition_control_state("bars")

    assert parent_label == "Group bars by:"
    assert child_label == "Stack bars by:"
    assert value_style == {}
    assert da_controls_style == {"marginBottom": "10px"}
    assert da_panel_style == {"marginTop": "24px", "width": "100%"}


def test_stacked_bar_can_swap_to_horizontal_axes():
    fig = plot_stacked_bar(
        "lineage",
        "cell_type",
        "prop",
        _adata(),
        swap_axes=True,
    )

    assert all(trace.orientation == "h" for trace in fig.data)
    assert fig.layout.yaxis.autorange == "reversed"


def test_stacked_bar_accepts_session_local_lasso_groups():
    adata = _adata()
    fig = plot_stacked_bar(
        "lineage",
        SELECTION_GROUP,
        "prop",
        adata,
        group_values={
            SELECTION_GROUP: selection_group_values(
                adata,
                ["cell-1"],
            )
        },
    )

    assert {trace.name for trace in fig.data} == {"Selected", "Others"}
