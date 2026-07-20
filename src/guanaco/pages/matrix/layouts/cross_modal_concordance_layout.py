"""Self-contained controls and plots for paired cross-modal concordance."""

import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.pages.matrix.plots.cross_modal_concordance import (
    GROUP_CORRELATION,
    RELATIVE_SKEW,
)
from guanaco.utils.plot_config import common_config, concordance_scatter_config
from guanaco.utils.ui_helpers import labeled_dropdown


def _feature_control(prefix, role, feature):
    return labeled_dropdown(
        f"Feature {role.upper()}:",
        f"{prefix}-concordance-feature-{role}",
        ([{"label": feature, "value": feature}] if feature else []),
        value=[feature] if feature else [],
        clearable=False,
        multi=True,
        placeholder="Search features from any modality",
    )


def generate_cross_modal_concordance_layout(source, prefix):
    """Build the paired-feature concordance exploratory tab."""
    modalities = list(source.modalities)
    modality_a = modalities[0]
    modality_b = modalities[1]
    feature_a = source.first_feature(modality_a)
    feature_b = source.first_feature(modality_b)

    feature_controls = html.Div(
        [
            html.Div(
                _feature_control(prefix, "a", feature_a),
                className="concordance-control concordance-control--feature",
            ),
            html.Div(
                _feature_control(prefix, "b", feature_b),
                className="concordance-control concordance-control--feature",
            ),
            html.Div(
                labeled_dropdown(
                    "View:",
                    f"{prefix}-concordance-view-mode",
                    [
                        {"label": "Relative skew", "value": RELATIVE_SKEW},
                        {
                            "label": "Correlation",
                            "value": GROUP_CORRELATION,
                        },
                    ],
                    value=RELATIVE_SKEW,
                    clearable=False,
                ),
                className="concordance-control concordance-control--view",
            ),
            html.Div(
                dbc.Button(
                    "Update comparison",
                    id=f"{prefix}-concordance-update-comparison",
                    n_clicks=0,
                    className="update-other-plots-button concordance-update-button",
                ),
                className="concordance-control concordance-control--action",
            ),
        ],
        className="concordance-control-row concordance-controls",
    )

    analysis_controls = html.Div(
        [
            html.Div(
                labeled_dropdown(
                    "Embedding:",
                    f"{prefix}-concordance-embedding",
                    [
                        {"label": embedding, "value": embedding}
                        for embedding in source.embedding_names
                    ],
                    value=source.preferred_embedding(modality_a),
                    clearable=False,
                ),
                className="concordance-control concordance-control--embedding",
            ),
            html.Div(
                labeled_dropdown(
                    "Group cells by:",
                    f"{prefix}-concordance-group-by",
                    [
                        {"label": obs, "value": obs}
                        for obs in source.discrete_obs_names
                    ],
                    value=None,
                    clearable=True,
                    placeholder="Select categorical metadata",
                ),
                id=f"{prefix}-concordance-group-control",
                className="concordance-control concordance-control--group",
            ),
        ],
        className="concordance-control-row concordance-analysis-controls",
    )

    embedding_panel = html.Div(
        [
            html.Div(
                dcc.Loading(
                    dcc.Graph(
                        id=f"{prefix}-concordance-embedding-plot",
                        config=concordance_scatter_config,
                        responsive=True,
                        style={"height": "520px", "width": "100%"},
                    ),
                    type="circle",
                ),
                className="concordance-main-plot",
            ),
        ],
        className="concordance-panel",
    )

    expression_panel = html.Div(
        [
            html.Div(
                dcc.Loading(
                    dcc.Graph(
                        id=f"{prefix}-concordance-feature-plot",
                        config=concordance_scatter_config,
                        responsive=True,
                        style={"height": "520px", "width": "100%"},
                    ),
                    type="circle",
                ),
                className="concordance-main-plot",
            ),
        ],
        className="concordance-panel",
    )

    summary_panel = html.Div(
        dcc.Loading(
            dcc.Graph(
                id=f"{prefix}-concordance-group-summary-plot",
                config=common_config,
                responsive=True,
                style={"height": "520px", "width": "100%"},
            ),
            type="circle",
        ),
        className="concordance-panel concordance-summary-panel",
    )

    return html.Div(
        [
            dcc.Store(
                id=f"{prefix}-concordance-committed-features",
                data={
                    "feature_a": [feature_a] if feature_a else [],
                    "feature_b": [feature_b] if feature_b else [],
                },
            ),
            dcc.Store(id=f"{prefix}-concordance-linked-cells"),
            dcc.Store(id=f"{prefix}-concordance-highlight-link"),
            dcc.Store(id=f"{prefix}-concordance-main-view-trigger"),
            feature_controls,
            html.Div(
                id=f"{prefix}-concordance-update-status",
                className="concordance-update-status",
            ),
            analysis_controls,
            html.Div(
                [
                    embedding_panel,
                    expression_panel,
                    summary_panel,
                ],
                id=f"{prefix}-concordance-main-row",
                className="concordance-main-row concordance-main-row--with-summary",
            ),
        ],
        className="cross-modal-concordance",
    )
