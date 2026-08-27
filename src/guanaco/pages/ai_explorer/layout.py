"""Layout for the optional Ask your data demo."""

import dash_bootstrap_components as dbc
from dash import dcc, html


def generate_ai_explorer_layout(prefix, default_features=None):
    default_features = list(default_features or [])
    example = (
        f"How does {default_features[0]} differ between the main clinical groups?"
        if default_features
        else "How do the measured variables differ between the main groups?"
    )
    return html.Div(
        [
            html.Div(
                [
                    html.H5("Ask your data", className="mb-1"),
                    html.P(
                        "Describe a data question. GUANACO first checks whether "
                        "the required variables exist, then builds only validated plots.",
                        className="text-muted mb-3",
                    ),
                    dcc.Textarea(
                        id=f"{prefix}-ai-explorer-question",
                        placeholder=example,
                        style={"width": "100%", "minHeight": "92px"},
                        className="custom-textarea",
                    ),
                    html.Div(
                        [
                            dbc.Button(
                                "Generate visualization",
                                id=f"{prefix}-ai-explorer-submit",
                                color="primary",
                                n_clicks=0,
                            ),
                            html.Span(
                                "The matrix stays local; only field names and the "
                                "question are sent when an LLM is configured.",
                                className="text-muted ms-3 small",
                            ),
                        ],
                        className="d-flex align-items-center flex-wrap gap-2",
                    ),
                ],
                className="p-3",
            ),
            dcc.Loading(
                html.Div(
                    id=f"{prefix}-ai-explorer-results",
                    className="px-3 pb-3",
                ),
                type="circle",
            ),
        ],
        id=f"{prefix}-ai-explorer",
        className="ai-explorer-demo",
    )
