"""Layout for the interactive multi-omics coverage/composition view."""

from dash import dcc, html
import dash_bootstrap_components as dbc

from guanaco.utils.plot_config import common_config
from guanaco.utils.ui_helpers import (
    graph_flex_container,
    labeled_dropdown,
    responsive_graph_grid,
)


def generate_multiomics_composition_layout(source, prefix):
    """Build controls, coverage heatmap, and complete-case summary."""
    modalities = list(source.modalities)
    metadata_names = source.coverage_metadata_names

    controls = dbc.Row(
        [
            dbc.Col(
                labeled_dropdown(
                    "Group and color by:",
                    f"{prefix}-multiomics-composition-group-by",
                    [{"label": value, "value": value} for value in metadata_names],
                    value=None,
                    clearable=True,
                    placeholder="Select clinical or sample metadata",
                ),
                xs=12,
                lg=6,
            ),
            dbc.Col(
                labeled_dropdown(
                    "Required modalities:",
                    f"{prefix}-multiomics-composition-required",
                    [
                        {"label": str(modality).upper(), "value": modality}
                        for modality in modalities
                    ],
                    value=modalities,
                    clearable=False,
                    multi=True,
                    placeholder="Select model inputs",
                ),
                xs=12,
                lg=6,
            ),
        ],
        className="g-3 multiomics-composition-controls",
    )

    grid = responsive_graph_grid(
        f"{prefix}-multiomics-composition-grid",
        f"{prefix}-multiomics-composition-grid-item",
        graph_flex_container(
            f"{prefix}-multiomics-composition-plot",
            container_id=f"{prefix}-multiomics-composition-grid-item",
            config=common_config,
        ),
        w=12,
        h=16,
        min_w=6,
        min_h=10,
        max_w=12,
        max_h=24,
    )

    return html.Div(
        [
            controls,
            dcc.Loading(grid, type="circle"),
            html.Details(
                [
                    html.Summary(
                        "Coverage by group",
                        className="multiomics-composition-summary-title",
                    ),
                    html.Div(id=f"{prefix}-multiomics-composition-summary"),
                ],
                id=f"{prefix}-multiomics-composition-summary-details",
                className="multiomics-composition-summary",
            ),
        ],
        className="multiomics-composition",
    )
