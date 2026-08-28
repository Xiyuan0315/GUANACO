"""Callbacks for multi-omics observation coverage and cohort feasibility."""

import dash_bootstrap_components as dbc
from dash import Input, Output, html, no_update

from guanaco.pages.matrix.plots.multiomics_composition import (
    build_multiomics_composition_figure,
    multiomics_coverage_summary,
    normalize_required_modalities,
)
from guanaco.pages.visualizations.registry import EXPLORATION_WORKSPACE


_ACTIVE_TAB = "multiomics-composition-tab"


def _summary_table(frame, required_modalities):
    if frame.empty:
        return html.Div(
            "No observations remain for the selected groups.",
            className="visualization-workspace-empty",
        )
    required = " + ".join(str(value).upper() for value in required_modalities)
    return html.Div(
        [
            html.Div(
                f"Complete = observations measured for {required}.",
                className="multiomics-composition-summary-note",
            ),
            dbc.Table.from_dataframe(
                frame,
                striped=True,
                bordered=False,
                hover=True,
                responsive=True,
                index=False,
                class_name="multiomics-composition-table",
            ),
        ]
    )


def register_multiomics_composition_callbacks(
    app,
    source,
    prefix,
    *,
    color_config=None,
):
    """Register metadata grouping and lazy coverage rendering."""

    @app.callback(
        Output(f"{prefix}-multiomics-composition-plot", "figure"),
        Output(f"{prefix}-multiomics-composition-summary", "children"),
        Input(f"{prefix}-multiomics-composition-group-by", "value"),
        Input(f"{prefix}-multiomics-composition-required", "value"),
        Input(f"{prefix}-visualization-workspace-tabs", "value"),
        Input(f"{prefix}-exploratory-tabs", "value"),
    )
    def render_coverage(
        group_by,
        required_modalities,
        active_workspace,
        active_tab,
    ):
        if (
            active_workspace != EXPLORATION_WORKSPACE
            or active_tab != _ACTIVE_TAB
        ):
            return no_update, no_update

        required = normalize_required_modalities(source, required_modalities)
        figure = build_multiomics_composition_figure(
            source,
            group_by=group_by,
            required_modalities=required,
            color_config=color_config,
        )
        summary = multiomics_coverage_summary(
            source,
            group_by=group_by,
            required_modalities=required,
        )
        return figure, _summary_table(summary, required)
