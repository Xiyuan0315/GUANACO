import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from guanaco.pages.ai_explorer.planner import (
    build_data_profile,
    local_demo_plan,
    validate_plan,
)
from guanaco.pages.visualizations.layout import generate_visualization_sections


def _adata():
    obs = pd.DataFrame(
        {
            "IGHV_status": pd.Categorical(["Mutated", "Unmutated", "Mutated"]),
            "SF3B1_status": pd.Categorical(
                ["Not altered", "Altered", "Not altered"]
            ),
        },
        index=["patient-1", "patient-2", "patient-3"],
    )
    return AnnData(
        np.array([[0.9, 0.1], [0.7, 0.3], [1.0, 0.2]]),
        obs=obs,
        var=pd.DataFrame(index=["ibrutinib_3", "drug_b"]),
    )


def _walk(component):
    if component is None:
        return
    if isinstance(component, (list, tuple)):
        for child in component:
            yield from _walk(child)
        return
    yield component
    children = getattr(component, "children", None)
    if children is not None:
        yield from _walk(children)


def test_profile_does_not_read_or_include_matrix_values():
    profile = build_data_profile(
        _adata(),
        "Compare ibrutinib_3 by IGHV_status",
        ["ibrutinib_3"],
    )

    prompt = profile.for_prompt()
    assert prompt["n_observations"] == 3
    assert {item["key"] for item in prompt["available_fields"]} == {
        "obs:IGHV_status",
        "obs:SF3B1_status",
        "feature:ibrutinib_3",
    }
    assert set(prompt) == {
        "n_observations",
        "n_features",
        "available_fields",
        "privacy_note",
    }


def test_local_demo_builds_a_valid_comparative_violin():
    question = "How does ibrutinib_3 differ by IGHV_status and SF3B1_status?"
    profile = build_data_profile(_adata(), question, ["ibrutinib_3"])

    plan = validate_plan(local_demo_plan(question, profile), profile)

    assert plan["answerable"] is True
    assert plan["plots"][0]["type"] == "violin"
    assert plan["plots"][0]["x"] == "obs:IGHV_status"
    assert plan["plots"][0]["color"] == "obs:SF3B1_status"


def test_unknown_planner_field_is_rejected_instead_of_rendered():
    profile = build_data_profile(_adata(), "Compare survival", ["ibrutinib_3"])
    hallucinated = {
        "answerable": True,
        "answer_summary": "",
        "missing_data": [],
        "plots": [
            {
                "type": "violin",
                "title": "Invented result",
                "x": "obs:invented_group",
                "y": "feature:ibrutinib_3",
                "color": None,
                "features": [],
                "group_by": None,
                "reason": "",
            }
        ],
    }

    with pytest.raises(ValueError, match="unavailable fields"):
        validate_plan(hallucinated, profile)


def test_short_feature_is_not_matched_inside_an_unavailable_measurement_name():
    rna = AnnData(
        np.ones((3, 2)),
        obs=_adata().obs.copy(),
        var=pd.DataFrame(index=["T", "LPL"]),
    )
    question = "How does ibrutinib_3 differ by IGHV_status and SF3B1_status?"
    profile = build_data_profile(rna, question, ["LPL"])

    assert "feature:T" not in profile.by_key
    plan = validate_plan(local_demo_plan(question, profile), profile)
    assert plan["answerable"] is False
    assert plan["plots"] == []


def test_single_categorical_count_question_builds_a_count_bar():
    question = "How many people are mutated in IGHV_status group?"
    profile = build_data_profile(_adata(), question, ["ibrutinib_3"])

    plan = validate_plan(local_demo_plan(question, profile), profile)

    assert plan["answerable"] is True
    assert plan["plots"][0]["type"] == "count_bar"
    assert plan["plots"][0]["x"] == "obs:IGHV_status"


def test_demo_tab_is_opt_in_and_outside_the_plot_registry():
    sections = generate_visualization_sections(
        _adata(),
        ["ibrutinib_3"],
        ["IGHV_status", "SF3B1_status"],
        "cll-drug",
        optional_plot_components=["split-violin", "ai-explorer"],
    )
    tabs = next(
        item
        for item in _walk(sections)
        if getattr(item, "id", None) == "cll-drug-exploratory-tabs"
    )

    assert [tab.value for tab in tabs.children] == [
        "split-violin-tab",
        "ai-explorer-tab",
    ]
    assert tabs.children[-1].label == "Ask your data (demo)"
