"""Modality-scoped composition for feature analysis and dataset exploration.

This module owns where registered plots appear and gives each one an explicit
control model: shared feature controls or controls contained in an exploratory tab.
"""

from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.pages.matrix.layouts.atac_browser_layout import (
    generate_atac_browser_layout,
)
from guanaco.pages.matrix.layouts.cross_modal_concordance_layout import (
    generate_cross_modal_concordance_layout,
)
from guanaco.pages.matrix.layouts.dotplot_layout import generate_dotplot_layout
from guanaco.pages.matrix.layouts.grn_demo_layout import generate_grn_demo_layout
from guanaco.pages.matrix.layouts.heatmap_layout import generate_heatmap_layout
from guanaco.pages.matrix.layouts.paga_layout import generate_paga_layout
from guanaco.pages.matrix.layouts.pseudotime_layout import generate_pseudotime_layout
from guanaco.pages.matrix.layouts.stacked_bar_layout import generate_stacked_bar_layout
from guanaco.pages.matrix.layouts.violin_layout import (
    generate_split_violin_layout,
    generate_violin_layout,
)
from guanaco.pages.matrix.layouts.volcano_layout import generate_volcano_layout
from guanaco.pages.visualizations.plots.igv import build_igv_layout
from guanaco.pages.visualizations.registry import (
    EXPLORATION_WORKSPACE,
    EXPLORATORY_PLOTS,
    FEATURE_WORKSPACE,
    MARKER_PLOTS,
    PLOT_SPECS_BY_KEY,
    resolve_plot_components,
)


def generate_left_control(default_gene_markers, label_list, prefix):
    """Shared feature/annotation controls used only by feature plots."""
    default_gene_markers = list(default_gene_markers or [])
    label_list = list(label_list or [])
    default_annotation = label_list[len(label_list) // 2] if label_list else None

    genes_selection = html.Div(
        [
            dbc.RadioItems(
                id=f"{prefix}-gene-input-mode",
                options=[
                    {"label": "Dropdown", "value": "dropdown"},
                    {"label": "Text", "value": "text"},
                ],
                value="dropdown",
                inline=True,
                style={"fontSize": "14px", "marginBottom": "10px"},
            ),
            dcc.Dropdown(
                id=f"{prefix}-single-cell-genes-selection",
                options=[
                    {"label": gene, "value": gene} for gene in default_gene_markers
                ],
                value=default_gene_markers,
                multi=True,
                style={"marginBottom": "15px", "fontSize": "12px"},
                className="custom-dropdown",
            ),
            dcc.Textarea(
                id=f"{prefix}-single-cell-genes-textarea",
                placeholder="Enter genes separated by commas (e.g., Gene1, Gene2, Gene3)",
                value=", ".join(default_gene_markers),
                style={
                    "width": "100%",
                    "height": "80px",
                    "marginBottom": "10px",
                    "display": "none",
                },
                className="custom-textarea",
            ),
            html.Div(
                id=f"{prefix}-gene-input-error",
                style={"color": "red", "fontSize": "12px", "marginBottom": "10px"},
            ),
        ]
    )

    return html.Div(
        [
            html.Label(
                "Features", style={"fontWeight": "bold", "marginBottom": "5px"}
            ),
            genes_selection,
            html.Label(
                "Group cells by",
                style={"fontWeight": "bold", "marginBottom": "5px"},
            ),
            dcc.Dropdown(
                id=f"{prefix}-single-cell-annotation-dropdown",
                options=[{"label": label, "value": label} for label in label_list],
                value=default_annotation,
                style={"marginBottom": "15px"},
                clearable=False,
                className="custom-dropdown",
            ),
            html.Label(
                "Included groups",
                style={"fontWeight": "bold", "marginBottom": "5px"},
            ),
            dcc.Dropdown(
                id=f"{prefix}-single-cell-label-selection",
                multi=True,
                style={"marginBottom": "15px", "fontSize": "12px"},
                className="custom-dropdown",
            ),
        ],
        id=f"{prefix}-feature-controls",
        className="feature-controls-panel",
    )


def _plot_tab_classes(key):
    """Keep plot-type styling and optional icons isolated from tab behavior."""
    base = f"visualization-plot-tab visualization-plot-tab--{key}"
    return base, f"{base} visualization-plot-tab--selected"


def _marker_tabs(adata, markers, labels, prefix, enabled):
    factories = {
        "dotplot": lambda: generate_dotplot_layout(prefix),
        "heatmap": lambda: generate_heatmap_layout(adata, prefix),
        "violin": lambda: generate_violin_layout(markers, labels, prefix),
        "pseudotime": lambda: generate_pseudotime_layout(prefix),
    }
    children = []
    for key in MARKER_PLOTS:
        if key not in enabled:
            continue
        class_name, selected_class_name = _plot_tab_classes(key)
        children.append(
            dcc.Tab(
                label=PLOT_SPECS_BY_KEY[key].label,
                value=f"{key}-tab",
                children=[factories[key]()],
                className=class_name,
                selected_className=selected_class_name,
            )
        )
    if not children:
        return None
    return dcc.Tabs(
        children,
        id=f"{prefix}-marker-tabs",
        value=children[0].value,
        className="custom-tabs",
    )


def _exploratory_tabs(
    adata,
    markers,
    labels,
    prefix,
    enabled,
    gene_annotation_path,
    genome_tracks,
    multiomics_source=None,
):
    factories = {
        "split-violin": lambda: generate_split_violin_layout(
            adata,
            markers,
            labels,
            prefix,
        ),
        "stacked-bar": lambda: generate_stacked_bar_layout(adata, labels, prefix),
        "paga": lambda: generate_paga_layout(adata, prefix),
        "volcano": lambda: generate_volcano_layout(adata, prefix),
        "grn": lambda: generate_grn_demo_layout(adata, prefix),
        "peak-browser": lambda: generate_atac_browser_layout(
            adata,
            prefix,
            gene_annotation_path=gene_annotation_path,
            discrete_label_list=labels,
        ),
        "igv": lambda: build_igv_layout(prefix, genome_tracks or {}),
        "cross-modal-concordance": lambda: generate_cross_modal_concordance_layout(
            multiomics_source,
            prefix,
        ),
    }
    children = []
    for key in EXPLORATORY_PLOTS:
        if key not in enabled:
            continue
        class_name, selected_class_name = _plot_tab_classes(key)
        children.append(
            dcc.Tab(
                label=PLOT_SPECS_BY_KEY[key].label,
                value=f"{key}-tab",
                children=[factories[key]()],
                className=class_name,
                selected_className=selected_class_name,
            )
        )

    if not children:
        return None
    return dcc.Tabs(
        children,
        id=f"{prefix}-exploratory-tabs",
        value=children[0].value,
        className="custom-tabs",
    )


def _section(body):
    return html.Section(
        dbc.Card(
            html.Div(body, className="plot-section"),
            className="card-elevated",
        ),
        className="visualization-section",
    )


def _empty_workspace(message):
    return html.Div(message, className="visualization-workspace-empty")


def _workspace_tab(label, value, content):
    return dcc.Tab(
        label=label,
        value=value,
        children=[html.Div(content, className="visualization-workspace-content")],
        className="visualization-workspace-tab",
        selected_className="visualization-workspace-tab--selected",
    )


def generate_visualization_sections(
    adata,
    default_gene_markers,
    discrete_label_list,
    prefix,
    optional_plot_components=None,
    gene_annotation_path=None,
    genome_tracks=None,
    ref_track=None,
    multiomics_source=None,
):
    """Build one modality-scoped workspace with two control models."""
    has_igv = bool(genome_tracks) and bool(ref_track)
    enabled = resolve_plot_components(
        adata,
        optional_plot_components,
        has_igv=has_igv,
    )
    if multiomics_source is not None:
        enabled = tuple(key for key in enabled if key in MARKER_PLOTS) + (
            "cross-modal-concordance",
        )
    marker_tabs = _marker_tabs(
        adata, default_gene_markers, discrete_label_list, prefix, enabled
    )
    exploratory_tabs = _exploratory_tabs(
        adata,
        default_gene_markers,
        discrete_label_list,
        prefix,
        enabled,
        gene_annotation_path,
        genome_tracks,
        multiomics_source,
    )

    if marker_tabs is None and exploratory_tabs is None:
        return []

    feature_content = (
        dbc.Row(
            [
                dbc.Col(
                    generate_left_control(
                        default_gene_markers, discrete_label_list, prefix
                    ),
                    xs=12,
                    md=4,
                    xl=2,
                    className="feature-controls-sidebar",
                ),
                dbc.Col(marker_tabs, xs=12, md=8, xl=10),
            ],
            className="g-0",
        )
        if marker_tabs is not None
        else _empty_workspace(
            "No feature-analysis plots are available for this modality."
        )
    )
    exploration_content = (
        exploratory_tabs
        if exploratory_tabs is not None
        else _empty_workspace(
            "No dataset-exploration plots are available for this modality."
        )
    )

    workspace = dcc.Tabs(
        [
            _workspace_tab(
                "Markers visualization",
                FEATURE_WORKSPACE,
                feature_content,
            ),
            _workspace_tab(
                "Exploratory visualization",
                EXPLORATION_WORKSPACE,
                exploration_content,
            ),
        ],
        id=f"{prefix}-visualization-workspace-tabs",
        value=FEATURE_WORKSPACE,
        className="visualization-workspace-tabs",
    )
    return [_section(workspace)]
