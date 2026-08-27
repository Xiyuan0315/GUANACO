"""Callbacks for independent scatter panels in unpaired multi-omics data."""

from __future__ import annotations

import pandas as pd
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from guanaco.pages.matrix.layouts.embedding_layout import (
    EMBEDDING_PREFIXES,
    selectable_scatter_annotations,
)
from guanaco.pages.matrix.plots.embedding import (
    plot_coexpression_embedding,
    plot_embedding,
)
from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.search import ranked_substring_matches


def _is_continuous_obs(adata, column: str) -> bool:
    if column not in adata.obs.columns:
        return False
    dtype = adata.obs.dtypes[column]
    return bool(
        pd.api.types.is_numeric_dtype(dtype)
        and adata.obs[column].nunique(dropna=True) >= 50
    )


def _panel_color_options(source, embedding, query, current):
    modality = source.embedding_modality(embedding)
    if modality is None:
        return [], None
    adata = source.modality_adata(modality)
    annotations = selectable_scatter_annotations(adata)
    needle = str(query or "").strip().lower()
    if needle:
        annotations = ranked_substring_matches(annotations, needle, limit=10)
        features = source.search_features(modality, needle, limit=20)
    else:
        features = source.search_features(modality, None, limit=20)

    current_is_valid = (
        current in adata.obs.columns
        or (
            source.is_feature(current)
            and source.feature_modality(current) == modality
        )
    )
    selected = current if current_is_valid else source.first_feature(modality)
    values = list(dict.fromkeys([*annotations, *features]))
    if selected and selected not in values:
        values.append(selected)
    return [{"label": value, "value": value} for value in values], selected


def build_unpaired_embedding_figure(
    source,
    embedding,
    color,
    *,
    continuous_colormap,
    discrete_colormap,
    marker_size,
    opacity,
    render_backend,
    color_config,
    x_axis=None,
    y_axis=None,
    order="max",
    legend_show="right",
    axis_show=True,
):
    """Render one panel from only the modality that owns its embedding."""
    modality, raw_embedding, adata = source.embedding_context(embedding)
    if source.is_feature(color):
        if source.feature_modality(color) != modality:
            raise ValueError("The selected feature belongs to another modality.")
        raw_color = source.feature_key(color)
        mode = "continuous"
    elif color in adata.obs.columns:
        raw_color = color
        mode = "continuous" if _is_continuous_obs(adata, color) else "categorical"
    else:
        raise ValueError("The selected color is unavailable in this modality.")

    palette = None
    if mode == "categorical":
        palette = resolve_discrete_palette(
            discrete_colormap,
            adata.obs[raw_color].nunique(dropna=True),
            default=color_config,
        )
    return plot_embedding(
        adata=adata,
        adata_full=adata,
        embedding_key=raw_embedding,
        color=raw_color,
        x_axis=x_axis,
        y_axis=y_axis,
        mode=mode,
        order=order if mode == "continuous" else None,
        continuous_color_map=continuous_colormap or "Viridis",
        discrete_color_map=palette,
        marker_size=marker_size,
        opacity=opacity,
        render_backend=render_backend,
        legend_show=legend_show,
        axis_show=axis_show,
        source_adata=adata,
    )


def register_unpaired_multiomics_callbacks(
    app,
    source,
    prefix: str,
    *,
    embedding_render_backend="scattergl",
    color_config=None,
):
    """Register two independent modality-scoped scatter callback sets."""

    @app.callback(
        Output(f"{prefix}-controls-container", "style"),
        Output(f"{prefix}-toggle-button", "children"),
        Input(f"{prefix}-toggle-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def toggle_controls(n_clicks):
        if n_clicks % 2:
            return {"display": "block"}, "Hide controls"
        return {"display": "none"}, "More controls"

    def coordinate_controls(embedding, id_prefix):
        _modality, raw_key, adata = source.embedding_context(embedding)
        label = EMBEDDING_PREFIXES.get(
            raw_key, raw_key.removeprefix("X_")
        )
        columns = [
            f"{label}{index + 1}"
            for index in range(adata.obsm[raw_key].shape[1])
        ]
        options = [{"label": column, "value": column} for column in columns]
        x_value = columns[0]
        y_value = columns[1]
        style = (
            {"display": "flex", "marginBottom": "15px"}
            if raw_key == "X_pca"
            else {"display": "none"}
        )
        return (
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("X-axis:"),
                            dcc.Dropdown(
                                id=f"{id_prefix}-x-axis",
                                options=options,
                                value=x_value,
                                clearable=False,
                                style={"fontSize": "14px"},
                            ),
                        ],
                        style={"flex": "1", "paddingRight": "10px"},
                    ),
                    html.Div(
                        [
                            html.Label("Y-axis:"),
                            dcc.Dropdown(
                                id=f"{id_prefix}-y-axis",
                                options=options,
                                value=y_value,
                                clearable=False,
                                style={"fontSize": "14px"},
                            ),
                        ],
                        style={"flex": "1", "paddingLeft": "10px"},
                    ),
                ],
                style=style,
            ),
            x_value,
            y_value,
        )

    @app.callback(
        Output(f"{prefix}-coordinates-dropdowns", "children"),
        Output(f"{prefix}-x-axis", "value"),
        Output(f"{prefix}-y-axis", "value"),
        Input(f"{prefix}-clustering-dropdown", "value"),
    )
    def update_left_coordinates(embedding):
        return coordinate_controls(embedding, prefix)

    @app.callback(
        Output(f"{prefix}-right-coordinates-dropdowns", "children"),
        Output(f"{prefix}-right-x-axis", "value"),
        Output(f"{prefix}-right-y-axis", "value"),
        Input(f"{prefix}-right-clustering-dropdown", "value"),
    )
    def update_right_coordinates(embedding):
        return coordinate_controls(embedding, f"{prefix}-right")

    def register_color_search(role: str, component_id: str, embedding_id: str):
        @app.callback(
            Output(component_id, "options"),
            Output(component_id, "value"),
            Input(embedding_id, "value"),
            Input(component_id, "search_value"),
            State(component_id, "value"),
        )
        def update_color_options(embedding, query, current):
            return _panel_color_options(source, embedding, query, current)

    register_color_search(
        "left",
        f"{prefix}-annotation-dropdown",
        f"{prefix}-clustering-dropdown",
    )
    register_color_search(
        "right",
        f"{prefix}-scatter-gene-selection",
        f"{prefix}-right-clustering-dropdown",
    )

    @app.callback(
        Output(f"{prefix}-scatter-gene2-selection", "options"),
        Input(f"{prefix}-scatter-gene2-selection", "search_value"),
        State(f"{prefix}-right-clustering-dropdown", "value"),
        State(f"{prefix}-scatter-gene2-selection", "value"),
    )
    def update_second_feature(query, embedding, current):
        modality = source.embedding_modality(embedding)
        if modality is None:
            raise PreventUpdate
        values = source.search_features(modality, query, limit=20)
        if current and current not in values:
            values.append(current)
        return [{"label": value, "value": value} for value in values]

    @app.callback(
        Output(f"{prefix}-gene2-container", "style"),
        Output(f"{prefix}-threshold-container", "style"),
        Input(f"{prefix}-coexpression-toggle", "value"),
    )
    def toggle_coexpression_controls(mode):
        if mode == "coexpression":
            return {"display": "block"}, {"display": "block"}
        return {"display": "none"}, {"display": "none"}

    @app.callback(
        Output(f"{prefix}-annotation-scatter", "figure"),
        Input(f"{prefix}-clustering-dropdown", "value"),
        Input(f"{prefix}-x-axis", "value"),
        Input(f"{prefix}-y-axis", "value"),
        Input(f"{prefix}-annotation-dropdown", "value"),
        Input(f"{prefix}-marker-size-slider", "value"),
        Input(f"{prefix}-opacity-slider", "value"),
        Input(f"{prefix}-scatter-legend-toggle", "value"),
        Input(f"{prefix}-axis-toggle", "value"),
        Input(f"{prefix}-discrete-color-map-dropdown", "value"),
        Input(f"{prefix}-plot-order", "value"),
        Input(f"{prefix}-scatter-color-map-dropdown", "value"),
    )
    def update_left_scatter(
        embedding,
        x_axis,
        y_axis,
        color,
        marker_size,
        opacity,
        legend_show,
        axis_show,
        discrete_colormap,
        order,
        continuous_colormap,
    ):
        return build_unpaired_embedding_figure(
            source,
            embedding,
            color,
            x_axis=x_axis,
            y_axis=y_axis,
            continuous_colormap=continuous_colormap,
            discrete_colormap=discrete_colormap,
            marker_size=marker_size,
            opacity=opacity,
            order=order,
            legend_show=legend_show,
            axis_show=axis_show,
            render_backend=embedding_render_backend,
            color_config=color_config,
        )

    @app.callback(
        Output(f"{prefix}-gene-scatter", "figure"),
        Input(f"{prefix}-scatter-gene-selection", "value"),
        Input(f"{prefix}-right-clustering-dropdown", "value"),
        Input(f"{prefix}-right-x-axis", "value"),
        Input(f"{prefix}-right-y-axis", "value"),
        Input(f"{prefix}-plot-order", "value"),
        Input(f"{prefix}-scatter-color-map-dropdown", "value"),
        Input(f"{prefix}-marker-size-slider", "value"),
        Input(f"{prefix}-opacity-slider", "value"),
        Input(f"{prefix}-axis-toggle", "value"),
        Input(f"{prefix}-coexpression-toggle", "value"),
        Input(f"{prefix}-scatter-gene2-selection", "value"),
        Input(f"{prefix}-gene1-threshold-slider", "value"),
        Input(f"{prefix}-gene2-threshold-slider", "value"),
        Input(f"{prefix}-scatter-legend-toggle", "value"),
        Input(f"{prefix}-discrete-color-map-dropdown", "value"),
    )
    def update_right_scatter(
        color,
        embedding,
        x_axis,
        y_axis,
        order,
        continuous_colormap,
        marker_size,
        opacity,
        axis_show,
        coexpression_mode,
        second_color,
        threshold1,
        threshold2,
        legend_show,
        discrete_colormap,
    ):
        modality, raw_embedding, adata = source.embedding_context(embedding)
        if (
            coexpression_mode == "coexpression"
            and source.is_feature(color)
            and source.is_feature(second_color)
            and source.feature_modality(color) == modality
            and source.feature_modality(second_color) == modality
        ):
            return plot_coexpression_embedding(
                adata=adata,
                embedding_key=raw_embedding,
                gene1=source.feature_key(color),
                gene2=source.feature_key(second_color),
                x_axis=x_axis,
                y_axis=y_axis,
                threshold1=threshold1,
                threshold2=threshold2,
                marker_size=marker_size,
                opacity=opacity,
                legend_show=legend_show,
                axis_show=axis_show,
                source_adata=adata,
            )
        return build_unpaired_embedding_figure(
            source,
            embedding,
            color,
            x_axis=x_axis,
            y_axis=y_axis,
            continuous_colormap=continuous_colormap,
            discrete_colormap=discrete_colormap,
            marker_size=marker_size,
            opacity=opacity,
            order=order,
            legend_show=legend_show,
            axis_show=axis_show,
            render_backend=embedding_render_backend,
            color_config=color_config,
        )
