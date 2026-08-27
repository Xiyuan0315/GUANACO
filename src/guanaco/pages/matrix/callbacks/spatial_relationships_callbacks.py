"""Callbacks linking precomputed Squidpy neighborhood and co-occurrence views."""

from dash import Input, Output, no_update

from guanaco.pages.matrix.plots.spatial_relationships import (
    empty_spatial_figure,
    load_spatial_relationship_payload,
    plot_co_occurrence_group,
    plot_neighborhood_enrichment,
    spatial_pair_from_click,
)
from guanaco.utils.colors import resolve_discrete_palette


_WORKSPACE = "dataset-exploration"
_TAB = "spatial-relationships-tab"


def register_spatial_relationships_callbacks(
    app,
    adata,
    prefix: str,
    *,
    color_config=None,
):
    @app.callback(
        Output(f"{prefix}-spatial-nhood-heatmap", "figure"),
        Output(f"{prefix}-spatial-nhood-heatmap", "clickData"),
        Input(f"{prefix}-spatial-relationships-groupby", "value"),
        Input(f"{prefix}-spatial-cluster-columns", "value"),
        Input(f"{prefix}-visualization-workspace-tabs", "value"),
        Input(f"{prefix}-exploratory-tabs", "value"),
    )
    def update_neighborhood_heatmap(
        cluster_key,
        cluster_columns,
        active_workspace,
        active_tab,
    ):
        if active_workspace != _WORKSPACE or active_tab != _TAB:
            return no_update, no_update
        if not cluster_key:
            return (
                empty_spatial_figure(
                    "No compatible precomputed Squidpy results are available."
                ),
                None,
            )
        try:
            payload = load_spatial_relationship_payload(adata, cluster_key)
            return (
                plot_neighborhood_enrichment(
                    payload,
                    cluster_columns="cluster" in (cluster_columns or []),
                ),
                None,
            )
        except ValueError as exc:
            return empty_spatial_figure(str(exc)), None

    @app.callback(
        Output(f"{prefix}-spatial-co-occurrence-curve", "figure"),
        Input(f"{prefix}-spatial-nhood-heatmap", "clickData"),
        Input(f"{prefix}-spatial-relationships-groupby", "value"),
        Input(f"{prefix}-discrete-color-map-dropdown", "value"),
        Input(f"{prefix}-visualization-workspace-tabs", "value"),
        Input(f"{prefix}-exploratory-tabs", "value"),
    )
    def update_co_occurrence_curve(
        click_data,
        cluster_key,
        discrete_color_map,
        active_workspace,
        active_tab,
    ):
        if active_workspace != _WORKSPACE or active_tab != _TAB:
            return no_update
        if not cluster_key:
            return empty_spatial_figure(
                "No compatible precomputed Squidpy results are available."
            )
        try:
            payload = load_spatial_relationship_payload(adata, cluster_key)
            conditional, observed = spatial_pair_from_click(click_data, payload)
            colors = resolve_discrete_palette(
                discrete_color_map,
                len(payload["categories"]),
                default=color_config,
            )
            return plot_co_occurrence_group(
                payload,
                conditional,
                emphasized=observed,
                colors=colors,
            )
        except ValueError as exc:
            return empty_spatial_figure(str(exc))
