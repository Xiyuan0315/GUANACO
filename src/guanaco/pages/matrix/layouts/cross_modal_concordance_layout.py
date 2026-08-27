"""Self-contained controls and plots for paired cross-modal concordance."""

import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.data.loader import obs_col
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
    """Build the paired or metadata-level unpaired comparison tab."""
    if not source.is_paired:
        return _generate_unpaired_comparison_layout(source, prefix)

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


def _generate_unpaired_comparison_layout(source, prefix):
    modality_a, modality_b = list(source.modalities[:2])
    feature_a = source.first_feature(modality_a)
    feature_b = source.first_feature(modality_b)
    groups_a = source.modality_discrete_obs_names(modality_a)
    groups_b = source.modality_discrete_obs_names(modality_b)
    shared_columns = [column for column in groups_a if column in set(groups_b)]
    adata_a = source.modality_adata(modality_a)
    adata_b = source.modality_adata(modality_b)
    comparable_columns = []
    for column in shared_columns:
        labels_a = set(
            obs_col(adata_a.obs, column).dropna().astype(str).unique()
        )
        labels_b = set(
            obs_col(adata_b.obs, column).dropna().astype(str).unique()
        )
        if len(labels_a & labels_b) >= 2:
            comparable_columns.append(column)
    default_a = (
        comparable_columns[0]
        if comparable_columns
        else (
            shared_columns[0]
            if shared_columns
            else (groups_a[0] if groups_a else None)
        )
    )
    default_b = default_a if default_a in groups_b else (groups_b[0] if groups_b else None)

    return html.Div(
        [
            dcc.Store(
                id=f"{prefix}-concordance-committed-features",
                data={
                    "feature_a": [feature_a] if feature_a else [],
                    "feature_b": [feature_b] if feature_b else [],
                },
            ),
            html.Div(
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
                        dbc.Button(
                            "Update comparison",
                            id=f"{prefix}-concordance-update-comparison",
                            n_clicks=0,
                            className=(
                                "update-other-plots-button "
                                "concordance-update-button"
                            ),
                        ),
                        className="concordance-control concordance-control--action",
                    ),
                ],
                className="concordance-control-row concordance-controls",
            ),
            html.Div(
                id=f"{prefix}-concordance-update-status",
                className="concordance-update-status",
            ),
            html.Div(
                [
                    html.Div(
                        labeled_dropdown(
                            f"{modality_a.upper()} groups:",
                            f"{prefix}-concordance-group-a",
                            [
                                {"label": column, "value": column}
                                for column in groups_a
                            ],
                            value=default_a,
                            clearable=False,
                        ),
                        className="concordance-control concordance-control--group",
                    ),
                    html.Div(
                        labeled_dropdown(
                            f"{modality_b.upper()} groups:",
                            f"{prefix}-concordance-group-b",
                            [
                                {"label": column, "value": column}
                                for column in groups_b
                            ],
                            value=default_b,
                            clearable=False,
                        ),
                        className="concordance-control concordance-control--group",
                    ),
                ],
                className=(
                    "concordance-control-row concordance-analysis-controls"
                ),
            ),
            html.Div(
                [
                    html.Div(
                        dcc.Loading(
                            dcc.Graph(
                                id=f"{prefix}-concordance-unpaired-scatter",
                                config=common_config,
                                responsive=True,
                                style={"height": "520px", "width": "100%"},
                            ),
                            type="circle",
                        ),
                        className="concordance-panel",
                    ),
                    html.Div(
                        dcc.Loading(
                            dcc.Graph(
                                id=f"{prefix}-concordance-unpaired-heatmap",
                                config=common_config,
                                responsive=True,
                                style={"height": "520px", "width": "100%"},
                            ),
                            type="circle",
                        ),
                        className="concordance-panel",
                    ),
                ],
                className="concordance-main-row",
            ),
        ],
        className="cross-modal-concordance",
    )
