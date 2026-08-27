"""Self-contained layout for linked precomputed spatial statistics."""

import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.pages.matrix.plots.spatial_relationships import (
    empty_spatial_figure,
    spatial_relationship_keys,
)
from guanaco.utils.plot_config import common_config
from guanaco.utils.ui_helpers import labeled_dropdown


def generate_spatial_relationships_layout(adata, prefix: str):
    cluster_keys = spatial_relationship_keys(adata)
    options = [{"label": key, "value": key} for key in cluster_keys]
    default_key = cluster_keys[0] if cluster_keys else None

    return html.Div(
        [
            html.Div(
                [
                    labeled_dropdown(
                        "Group cells by",
                        f"{prefix}-spatial-relationships-groupby",
                        options,
                        value=default_key,
                        clearable=False,
                        wrapper_style={"width": "240px"},
                        dropdown_style={"width": "240px"},
                    ),
                    html.Div(
                        [
                            html.Label(
                                "Appearance:",
                                style={
                                    "fontWeight": "bold",
                                    "marginBottom": "5px",
                                },
                            ),
                            dbc.Checklist(
                                id=f"{prefix}-spatial-cluster-columns",
                                options=[
                                    {
                                        "label": "Cluster columns",
                                        "value": "cluster",
                                    }
                                ],
                                value=["cluster"],
                                inline=True,
                                switch=True,
                            ),
                        ]
                    ),
                    html.Div(
                        "Select a heatmap cell to inspect that pair across distance.",
                        className="spatial-relationships-hint",
                    ),
                ],
                className="spatial-relationships-controls",
            ),
            html.Div(
                [
                    html.Div(
                        dcc.Graph(
                            id=f"{prefix}-spatial-nhood-heatmap",
                            figure=empty_spatial_figure(
                                "Open this tab to view neighborhood enrichment."
                            ),
                            config=common_config,
                            responsive=True,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="spatial-relationships-panel",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=f"{prefix}-spatial-co-occurrence-curve",
                            figure=empty_spatial_figure(
                                "Select a neighborhood pair to view co-occurrence."
                            ),
                            config=common_config,
                            responsive=True,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="spatial-relationships-panel",
                    ),
                ],
                className="spatial-relationships-figures",
            ),
        ],
        className="spatial-relationships-layout",
    )
