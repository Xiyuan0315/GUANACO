import numpy as np
import pandas as pd
from anndata import AnnData

from guanaco.pages.visualizations.layout import (
    generate_visualization_sections,
    resolve_plot_components,
)
from guanaco.pages.matrix.layouts.heatmap_layout import generate_heatmap_layout
from guanaco.pages.matrix.layouts.paga_layout import generate_paga_layout


def _adata():
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["A", "B", "A", "B"]),
            "condition": pd.Categorical(["control", "control", "case", "case"]),
            "sample": ["S1", "S1", "S2", "S2"],
        },
        index=[f"cell-{index}" for index in range(4)],
    )
    return AnnData(
        X=np.ones((4, 3), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
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


def _by_id(component, component_id):
    return next(
        item for item in _walk(component) if getattr(item, "id", None) == component_id
    )


def test_plot_components_keep_aliases_and_add_igv_by_capability():
    enabled = resolve_plot_components(
        _adata(),
        ["dotplot", "split_violin", "genome-browser"],
        has_igv=True,
    )

    # The legacy genome-browser key continues to mean the matrix Peak Browser;
    # it is removed here because this fixture has no genomic peak variables.
    assert enabled == ("dotplot", "split-violin", "igv")


def test_result_driven_plots_are_removed_when_results_are_missing():
    adata = _adata()

    assert resolve_plot_components(
        adata,
        [
            "paga",
            "volcano",
            "ligand-receptor",
            "spatial-relationships",
            "peak-browser",
        ],
    ) == ()


def test_visualizations_share_one_workspace_with_feature_analysis_as_default():
    prefix = "dataset-atac"
    sections = generate_visualization_sections(
        _adata(),
        ["G1", "G2"],
        ["cell_type", "condition"],
        prefix,
        optional_plot_components=[
            "dotplot",
            "violin",
            "ridge",
            "split-violin",
            "stacked-bar",
        ],
        genome_tracks={"sample": [{"name": "track"}]},
        ref_track={"label": "hg38"},
    )

    assert len(sections) == 1
    assert sections[0].className == "visualization-section"
    assert sections[0].children.className == "card-elevated"
    assert sections[0].children.children.className == "plot-section"
    workspace = _by_id(sections, f"{prefix}-visualization-workspace-tabs")
    assert workspace.value == "feature-analysis"
    assert [tab.value for tab in workspace.children] == [
        "feature-analysis",
        "dataset-exploration",
    ]
    assert [tab.label for tab in workspace.children] == [
        "Markers visualization",
        "Exploratory visualization",
    ]
    assert all(
        tab.className == "visualization-workspace-tab"
        and tab.selected_className == "visualization-workspace-tab--selected"
        for tab in workspace.children
    )
    section_text = {item for item in _walk(sections) if isinstance(item, str)}
    assert "Feature Visualization" not in section_text
    assert "Exploratory Visualization" not in section_text
    marker_tabs = _by_id(sections, f"{prefix}-marker-tabs")
    exploratory_tabs = _by_id(sections, f"{prefix}-exploratory-tabs")

    assert marker_tabs.value == "dotplot-tab"
    assert marker_tabs.children[0].label == "Dot plot"
    assert [tab.value for tab in marker_tabs.children] == [
        "dotplot-tab",
        "violin-tab",
        "ridge-tab",
    ]
    assert _by_id(sections, f"{prefix}-ridge-gene-selection").value == "G1"
    assert _by_id(sections, f"{prefix}-ridge-show-box").value == []
    ridge_grid = _by_id(sections, f"{prefix}-ridge-grid")
    assert ridge_grid.isDraggable is True
    assert ridge_grid.isResizable is True
    assert _by_id(sections, f"{prefix}-ridge-grid-item") is not None
    assert all(
        "visualization-plot-tab" in tab.className
        and "visualization-plot-tab--selected" in tab.selected_className
        for tab in [*marker_tabs.children, *exploratory_tabs.children]
    )
    assert _by_id(sections, f"{prefix}-violin-plot1").responsive is True
    assert [tab.value for tab in exploratory_tabs.children] == [
        "split-violin-tab",
        "stacked-bar-tab",
        "igv-tab",
    ]
    assert [tab.label for tab in exploratory_tabs.children] == [
        "Comparative Violin",
        "Composition",
        "IGV",
    ]
    composition_view = _by_id(sections, f"{prefix}-composition-view")
    assert composition_view.value == []
    assert [option["value"] for option in composition_view.options] == ["hierarchy"]
    assert _by_id(sections, f"{prefix}-composition-swap-axes").value == []
    sample_unit = _by_id(sections, f"{prefix}-composition-sample-unit")
    assert sample_unit.value is None
    assert "sample" in {option["value"] for option in sample_unit.options}
    assert _by_id(
        sections, f"{prefix}-composition-alr-reference"
    ).value is None
    assert _by_id(sections, f"{prefix}-stacked-bar-x-group").value == "cell_type"
    assert _by_id(sections, f"{prefix}-stacked-bar-stack-by").value == "condition"
    composition_grid = _by_id(sections, f"{prefix}-stacked-bar-grid")
    assert all(
        breakpoint_layout[0]["minW"] == 3
        for breakpoint_layout in composition_grid.layouts.values()
    )
    assert _by_id(sections, f"{prefix}-composition-alr-reference") is not None
    differential_abundance = _by_id(
        sections, f"{prefix}-composition-da-panel"
    )
    assert _by_id(
        differential_abundance, f"{prefix}-composition-da-toggle-label"
    ).children == "▸ Differential abundance"
    da_tooltip = _by_id(
        differential_abundance, f"{prefix}-composition-da-info-tooltip"
    )
    assert "sample-level cell counts" in da_tooltip.children
    assert _by_id(
        differential_abundance, f"{prefix}-composition-da-collapse"
    ).style == {"display": "none", "width": "100%"}
    assert (
        _by_id(differential_abundance, f"{prefix}-composition-da-controls")
        is not None
    )
    assert (
        _by_id(differential_abundance, f"{prefix}-composition-da-plot")
        is not None
    )
    assert _by_id(
        differential_abundance, f"{prefix}-composition-da-results"
    ).style == {"display": "none"}
    assert (
        _by_id(workspace.children[0], f"{prefix}-single-cell-genes-selection")
        is not None
    )
    controls = _by_id(sections[0], f"{prefix}-feature-controls")
    assert controls.className == "feature-controls-panel"
    control_text = {item for item in _walk(controls) if isinstance(item, str)}
    assert {"Features", "Group cells by", "Included groups"} <= control_text
    assert any(
        getattr(item, "className", None) == "feature-controls-sidebar"
        for item in _walk(workspace.children[0])
    )
    assert _by_id(workspace.children[1], f"{prefix}-igv-genome-select") is not None
    assert not any(
        getattr(item, "id", None) == f"{prefix}-single-cell-genes-selection"
        for item in _walk(workspace.children[1])
    )


def test_tracks_only_modality_keeps_feature_analysis_as_default():
    prefix = "dataset-genome"
    sections = generate_visualization_sections(
        None,
        [],
        [],
        prefix,
        genome_tracks={"sample": [{"name": "track"}]},
        ref_track={"label": "hg38"},
    )

    assert len(sections) == 1
    workspace = _by_id(sections, f"{prefix}-visualization-workspace-tabs")
    assert workspace.value == "feature-analysis"
    assert "No feature-analysis plots are available for this modality." in {
        item for item in _walk(workspace.children[0]) if isinstance(item, str)
    }
    exploratory_tabs = _by_id(sections, f"{prefix}-exploratory-tabs")
    assert [tab.value for tab in exploratory_tabs.children] == ["igv-tab"]
    assert _by_id(sections, f"{prefix}-igv-genome-select") is not None


def test_paga_separates_selections_from_download_action():
    prefix = "dataset-rna"
    layout = generate_paga_layout(_adata(), prefix)
    controls = layout.children[1]
    selections_row, download_row = controls.children

    selection_ids = {
        getattr(item, "id", None) for item in _walk(selections_row)
    }
    assert {
        f"{prefix}-paga-color-mode",
        f"{prefix}-paga-obs-dropdown",
        f"{prefix}-paga-gene-dropdown",
        f"{prefix}-paga-threshold",
    } <= selection_ids
    assert _by_id(download_row, f"{prefix}-paga-download-svg") is not None
    assert f"{prefix}-paga-download-svg" not in selection_ids


def test_categorical_layout_selectors_exclude_high_cardinality_ids():
    n_obs = 60
    obs = pd.DataFrame(
        {
            "condition": ["control", "treated"] * (n_obs // 2),
            "cell_type": pd.Categorical(["B", "T", "Mono"] * (n_obs // 3)),
            "cell_id": pd.Categorical([f"cell-{index}" for index in range(n_obs)]),
        }
    )
    adata = AnnData(
        X=np.ones((n_obs, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1"]),
    )

    heatmap = generate_heatmap_layout(adata, "p")
    heatmap_values = {
        option["value"]
        for option in _by_id(heatmap, "p-heatmap-label-dropdown").options
    }
    paga = generate_paga_layout(adata, "p")
    paga_values = {
        option["value"]
        for option in _by_id(paga, "p-paga-obs-dropdown").options
    }

    assert heatmap_values == {"None", "condition", "cell_type"}
    assert paga_values == {"condition", "cell_type"}
