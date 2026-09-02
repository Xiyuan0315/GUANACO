"""Self-contained layout for precomputed ligand–receptor results."""

import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.data.ligand_receptor import load_default_ligand_receptor_result
from guanaco.pages.matrix.plots.ligand_receptor import (
    _empty,
    empty_ligand_receptor_figure,
    metric_slider_settings,
)
from guanaco.utils.plot_config import common_config
from guanaco.utils.ui_helpers import labeled_dropdown


def _control_values(payload):
    if not payload:
        return (
            [],
            None,
            [],
            None,
            [],
            [],
        )
    records = payload.get("records", [])
    metric_options = [
        {"label": metric["label"], "value": metric["value"]}
        for metric in payload.get("metrics", [])
    ]
    sources = sorted(
        {str(row["source"]) for row in records if row.get("source") is not None},
        key=str.casefold,
    )
    targets = sorted(
        {str(row["target"]) for row in records if row.get("target") is not None},
        key=str.casefold,
    )
    return (
        metric_options,
        payload.get("default_magnitude"),
        metric_options,
        payload.get("default_specificity"),
        [{"label": value, "value": value} for value in sources],
        [{"label": value, "value": value} for value in targets],
    )


def generate_ligand_receptor_layout(adata, prefix: str):
    payload = load_default_ligand_receptor_result(adata)
    (
        magnitude_options,
        default_magnitude,
        specificity_options,
        default_specificity,
        sender_options,
        receiver_options,
    ) = _control_values(payload)
    magnitude_slider = metric_slider_settings(payload, default_magnitude)
    specificity_slider = metric_slider_settings(payload, default_specificity)
    controls = html.Div(
        [
            labeled_dropdown(
                "Sender groups",
                f"{prefix}-lr-senders",
                sender_options,
                value=[],
                multi=True,
                placeholder="All senders",
                wrapper_style={"width": "230px"},
                dropdown_style={"width": "230px"},
            ),
            labeled_dropdown(
                "Receiver groups",
                f"{prefix}-lr-receivers",
                receiver_options,
                value=[],
                multi=True,
                placeholder="All receivers",
                wrapper_style={"width": "230px"},
                dropdown_style={"width": "230px"},
            ),
        ],
        className="lr-controls",
    )
    color_controls = html.Div(
        [
            labeled_dropdown(
                "Color by",
                f"{prefix}-lr-magnitude",
                magnitude_options,
                value=default_magnitude,
                clearable=False,
                placeholder="Magnitude metric",
                wrapper_style={"width": "100%"},
                dropdown_style={"width": "100%"},
            ),
            html.Div(
                [
                    html.Label(
                        "Color range",
                        style={"fontWeight": "bold", "marginBottom": "5px"},
                    ),
                    dcc.RangeSlider(
                        id=f"{prefix}-lr-magnitude-range",
                        **magnitude_slider,
                        tooltip={"placement": "bottom"},
                        updatemode="mouseup",
                        className="dbc-slider",
                    ),
                ],
                className="lr-range-control",
            ),
        ],
        className="lr-metric-control",
    )
    size_controls = html.Div(
        [
            labeled_dropdown(
                "Size by",
                f"{prefix}-lr-specificity",
                specificity_options,
                value=default_specificity,
                clearable=True,
                placeholder="Optional specificity metric",
                wrapper_style={"width": "100%"},
                dropdown_style={"width": "100%"},
            ),
            html.Div(
                [
                    html.Label(
                        "Size range",
                        style={"fontWeight": "bold", "marginBottom": "5px"},
                    ),
                    dcc.RangeSlider(
                        id=f"{prefix}-lr-specificity-range",
                        **specificity_slider,
                        tooltip={"placement": "bottom"},
                        updatemode="mouseup",
                        className="dbc-slider",
                    ),
                ],
                id=f"{prefix}-lr-specificity-range-wrap",
                className="lr-range-control",
                style={
                    "display": "block" if default_specificity else "none",
                },
            ),
        ],
        className="lr-metric-control",
    )
    advanced_controls = html.Div(
        [color_controls, size_controls],
        className="lr-controls lr-metric-controls",
    )
    advanced_toggle = dbc.Button(
        "▸ More options",
        id=f"{prefix}-lr-options-toggle",
        color="link",
        size="sm",
        style={
            "padding": "2px 0",
            "textDecoration": "none",
            "fontWeight": "normal",
            "marginBottom": "10px",
        },
    )
    advanced_panel = dbc.Collapse(
        advanced_controls,
        id=f"{prefix}-lr-options-collapse",
        is_open=False,
    )
    download_button = html.Button(
        "⬇ Download SVG",
        id=f"{prefix}-lr-download-svg",
        n_clicks=0,
        className="lr-download-button lr-download-button--overlay",
    )
    return html.Div(
        [
            dcc.Store(id=f"{prefix}-lr-view-store"),
            dcc.Store(id=f"{prefix}-lr-selection-store"),
            controls,
            advanced_toggle,
            advanced_panel,
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                _empty("Choose an embedded result."),
                                id=f"{prefix}-lr-network",
                                className="lr-network-container",
                            ),
                            download_button,
                        ],
                        className="lr-graph-shell lr-linked-panel",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=f"{prefix}-lr-dotplot",
                            figure=empty_ligand_receptor_figure(
                                "Load a precomputed result to view interactions."
                            ),
                            config=common_config,
                            responsive=True,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="lr-linked-panel",
                    ),
                ],
                className="lr-linked-views",
            ),
            html.Div(
                _empty("Select a circle edge or dot to inspect its interactions."),
                id=f"{prefix}-lr-detail",
                className="lr-detail-container",
            ),
        ],
        className="lr-layout",
    )
