"""Self-contained layout for database-derived RNA networks."""

import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.pages.matrix.plots.network import _empty_component
from guanaco.utils.ui_helpers import labeled_dropdown


NETWORK_TYPE_OPTIONS = [
    {"label": "PPI", "value": "ppi"},
    {"label": "TF–gene", "value": "tf-gene"},
    {"label": "Metabolite–gene", "value": "metabolite"},
    {"label": "miRNA–mRNA", "value": "mirna"},
]

NETWORK_LAYOUT_OPTIONS = [
    {"label": "Force-directed", "value": "cose"},
    {"label": "Hierarchical", "value": "breadthfirst"},
    {"label": "Circle", "value": "circle"},
    {"label": "Concentric", "value": "concentric"},
]

NETWORK_VIEW_OPTIONS = [
    {"label": "First-order", "value": "first-order"},
    {"label": "Input only", "value": "input-only"},
    {"label": "Minimum network", "value": "minimum"},
]

TF_DIRECTION_OPTIONS = [
    {"label": "Targets", "value": "targets"},
    {"label": "Regulators", "value": "regulators"},
    {"label": "Both", "value": "both"},
]


def generate_network_layout(prefix: str):
    gene_control = html.Div(
        [
            html.Label("Genes", style={"fontWeight": "bold", "marginBottom": "5px"}),
            dcc.Textarea(
                id=f"{prefix}-network-genes",
                placeholder="TP53, EGFR, CDK2",
                className="network-gene-input",
            ),
        ],
        className="network-control network-control--genes",
    )
    network_type_control = labeled_dropdown(
        "Network type",
        f"{prefix}-network-type",
        NETWORK_TYPE_OPTIONS,
        value="ppi",
        clearable=False,
        dropdown_style={"width": "220px"},
        wrapper_style={"width": "220px"},
    )
    ppi_network_type_control = html.Div(
        [
            html.Div(
                [
                    html.Span("Functional"),
                    dbc.Checklist(
                        id=f"{prefix}-network-string-physical",
                        options=[{"label": "Physical", "value": "physical"}],
                        value=[],
                        inline=True,
                        switch=True,
                        className="network-string-mode-switch",
                    ),
                ],
                className="network-string-mode-options",
            ),
        ],
        id=f"{prefix}-network-string-type-wrap",
        className="network-string-mode-control",
    )
    tf_direction_control = html.Div(
        labeled_dropdown(
            "Direction",
            f"{prefix}-network-tf-direction",
            TF_DIRECTION_OPTIONS,
            value="targets",
            clearable=False,
            dropdown_style={"width": "180px"},
            wrapper_style={"width": "180px"},
        ),
        id=f"{prefix}-network-tf-direction-wrap",
        style={"display": "none"},
    )
    build_button = html.Button(
        "Build network",
        id=f"{prefix}-network-build",
        n_clicks=0,
        className="update-other-plots-button network-build-button",
    )
    primary_controls = html.Div(
        [
            gene_control,
            network_type_control,
            ppi_network_type_control,
            tf_direction_control,
            build_button,
        ],
        id=f"{prefix}-network-primary-controls",
        className="network-controls network-controls--primary",
    )

    view_control = labeled_dropdown(
        "Network view",
        f"{prefix}-network-view",
        NETWORK_VIEW_OPTIONS,
        value="first-order",
        clearable=False,
        dropdown_style={"width": "180px"},
        wrapper_style={"width": "180px"},
    )
    layout_control = labeled_dropdown(
        "Layout",
        f"{prefix}-network-layout",
        NETWORK_LAYOUT_OPTIONS,
        value="cose",
        clearable=False,
        dropdown_style={"width": "180px"},
        wrapper_style={"width": "180px"},
    )
    ppi_confidence_control = html.Div(
        [
            html.Label(
                "Confidence",
                style={"fontWeight": "bold", "marginBottom": "5px"},
            ),
            dcc.Slider(
                id=f"{prefix}-network-string-confidence",
                min=0.4,
                max=0.9,
                step=0.05,
                value=0.7,
                marks={0.4: "0.4", 0.7: "0.7", 0.9: "0.9"},
                tooltip={"placement": "bottom"},
                updatemode="mouseup",
                className="dbc-slider",
            ),
        ],
        id=f"{prefix}-network-string-confidence-wrap",
        className="network-confidence-control",
    )
    download_button = html.Button(
        "⬇ Download SVG",
        id=f"{prefix}-network-download-svg",
        n_clicks=0,
        className="network-download-button network-download-button--overlay",
    )
    secondary_controls = html.Div(
        [
            view_control,
            layout_control,
            ppi_confidence_control,
        ],
        id=f"{prefix}-network-secondary-controls",
        className="network-controls network-controls--secondary",
    )
    options_toggle = dbc.Button(
        "▸ More options",
        id=f"{prefix}-network-options-toggle",
        color="link",
        size="sm",
        style={
            "padding": "2px 0",
            "textDecoration": "none",
            "fontWeight": "normal",
        },
    )
    options_panel = dbc.Collapse(
        secondary_controls,
        id=f"{prefix}-network-options-collapse",
        is_open=False,
    )
    controls = html.Div(
        [primary_controls, options_toggle, options_panel],
        className="network-control-rows",
    )
    return html.Div(
        [
            dcc.Store(id=f"{prefix}-network-graph-store"),
            dcc.Store(id=f"{prefix}-network-source-store"),
            dcc.Store(id=f"{prefix}-network-rendered-key"),
            dcc.Store(id=f"{prefix}-network-highlight-store"),
            controls,
            html.Div(id=f"{prefix}-network-status", className="network-status"),
            html.Div(
                [
                    dcc.Loading(
                        html.Div(
                            _empty_component(
                                "Enter a gene list, choose a network type, and build the network."
                            ),
                            id=f"{prefix}-network-graph",
                            className="network-graph-container",
                        ),
                        id=f"{prefix}-network-loading",
                        type="circle",
                        # Node highlighting only updates the nested Cytoscape stylesheet.
                        # Restrict the spinner to graph replacement so a tap does not flash
                        # or visually resemble a full network rebuild.
                        target_components={f"{prefix}-network-graph": "children"},
                    ),
                    download_button,
                ],
                id=f"{prefix}-network-graph-shell",
                className="network-graph-shell",
            ),
            html.Details(
                [
                    html.Summary("Key regulator enrichment"),
                    html.Div(id=f"{prefix}-network-regulator-results"),
                ],
                id=f"{prefix}-network-regulator-section",
                open=True,
                style={"display": "none"},
                className="network-enrichment-section",
            ),
        ],
        className="network-layout",
    )
