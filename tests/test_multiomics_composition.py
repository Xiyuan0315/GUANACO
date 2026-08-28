import warnings

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd

from guanaco.data.multiomics import MultiOmicsSource
from guanaco.pages.matrix.plots.multiomics_composition import (
    build_multiomics_composition_figure,
    multiomics_coverage_summary,
)
from guanaco.pages.visualizations.layout import generate_visualization_sections


def _coverage_mudata():
    specs = {
        "rna": (["p1", "p2", "p3", "p4"], 5),
        "drug": (["p1", "p3", "p4"], 3),
        "methylation": (["p2", "p3", "p4"], 4),
    }
    modalities = {}
    for modality, (observations, n_vars) in specs.items():
        modalities[modality] = ad.AnnData(
            X=np.ones((len(observations), n_vars), dtype=np.float32),
            obs=pd.DataFrame(index=observations),
            var=pd.DataFrame(
                index=[f"{modality}-{index}" for index in range(n_vars)]
            ),
        )
    mdata = mu.MuData(modalities)
    mdata.obs["response"] = pd.Categorical(
        ["Relapsed", "Relapsed", "Stable", "Stable"]
    )
    mdata.obs["binary_flag"] = [1, 0, 1, 0]
    return mdata


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


def _by_id(component, component_id):
    return next(
        item for item in _walk(component) if getattr(item, "id", None) == component_id
    )


def test_coverage_source_uses_observation_membership_without_embeddings():
    source = MultiOmicsSource(_coverage_mudata())

    assert source.supports_embedding_view is False
    assert source.coverage_ids.tolist() == ["p1", "p2", "p3", "p4"]
    assert source.coverage_matrix.astype(int).to_dict("list") == {
        "rna": [1, 1, 1, 1],
        "drug": [1, 0, 1, 1],
        "methylation": [0, 1, 1, 1],
    }
    assert {"response", "binary_flag"} <= set(source.coverage_metadata_names)
    assert source.modality_dimensions == {
        "rna": (4, 5),
        "drug": (3, 3),
        "methylation": (3, 4),
    }


def test_coverage_summary_counts_model_ready_observations_by_clinical_group():
    source = MultiOmicsSource(_coverage_mudata())

    summary = multiomics_coverage_summary(
        source,
        group_by="response",
        required_modalities=["rna", "drug"],
    ).set_index("Group")

    assert summary.loc["Relapsed"].to_dict() == {
        "N": 2.0,
        "RNA": 2.0,
        "DRUG": 1.0,
        "Complete": 1.0,
        "Complete %": 50.0,
    }
    assert summary.loc["Stable", "Complete"] == 2
    assert summary.loc["Stable", "Complete %"] == 100.0


def test_coverage_figure_uses_grouped_heatmap_and_aggregated_block_hover():
    source = MultiOmicsSource(_coverage_mudata())

    figure = build_multiomics_composition_figure(
        source,
        group_by="response",
        required_modalities=["rna", "drug"],
    )

    heatmap = figure.data[0]
    assert heatmap.type == "heatmap"
    assert heatmap.y[0] == "response"
    assert [label.split()[0] for label in heatmap.y[1:]] == ["RNA", "DRUG"]
    assert "Model-ready" not in heatmap.y
    assert len(heatmap.x) == 4
    assert figure.layout.xaxis.title.text is None
    assert figure.layout.plot_bgcolor == "#FFFFFF"
    assert figure.layout.shapes[0].line.color == "#111111"

    drug_hover = list(heatmap.customdata[2])
    assert any(
        "Modality: DRUG" in text
        and "Measured: No" in text
        and "response: Relapsed" in text
        and "Number of samples: 1" in text
        for text in drug_hover
    )


def test_joint_exploratory_workspace_contains_multiomics_coverage_tab():
    source = MultiOmicsSource(_coverage_mudata())
    sections = generate_visualization_sections(
        None,
        [],
        [],
        "joint",
        optional_plot_components=["multiomics-composition"],
        multiomics_source=source,
        modality_name="multiomics",
    )

    tabs = _by_id(sections, "joint-exploratory-tabs")
    assert [tab.value for tab in tabs.children] == [
        "multiomics-composition-tab",
        "cross-modal-concordance-tab",
    ]
    assert tabs.children[0].label == "Multi-omics coverage"
    assert _by_id(sections, "joint-multiomics-composition-group-by") is not None
    assert _by_id(sections, "joint-multiomics-composition-required").value == [
        "rna",
        "drug",
        "methylation",
    ]
    details = _by_id(sections, "joint-multiomics-composition-summary-details")
    assert getattr(details, "open", None) in (None, False)
    assert all(
        getattr(component, "id", None) not in {
            "joint-multiomics-composition-palette",
            "joint-multiomics-composition-groups",
        }
        for component in _walk(sections)
    )
    assert all(
        getattr(component, "className", None) != "multiomics-composition-intro"
        for component in _walk(sections)
    )


def test_coverage_metadata_merge_avoids_empty_concat_futurewarning():
    mdata = _coverage_mudata()
    mdata.obs["shared"] = pd.Categorical([None, None, None, None])
    mdata.mod["rna"].obs["shared"] = pd.Categorical(["A", "A", "B", "B"])

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "error",
            message="The behavior of array concatenation with empty entries",
            category=FutureWarning,
        )
        source = MultiOmicsSource(mdata)

    assert source.coverage_metadata("shared").tolist() == ["A", "A", "B", "B"]
