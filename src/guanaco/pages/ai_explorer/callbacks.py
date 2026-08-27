"""Dash callbacks and local rendering for the optional AI explorer demo."""

from __future__ import annotations

import pandas as pd
import plotly.express as px
from dash import Input, Output, State, dcc, html

from .planner import build_data_profile, field_values, plan_question


def _label(field_key):
    return field_key.split(":", 1)[1]


def _plot_figure(adata, spec):
    plot_type = spec["type"]
    if plot_type == "violin":
        x = field_values(adata, spec["x"])
        y = field_values(adata, spec["y"])
        frame = pd.DataFrame({_label(spec["x"]): x, _label(spec["y"]): y})
        color_name = None
        if spec.get("color"):
            color_name = _label(spec["color"])
            frame[color_name] = field_values(adata, spec["color"])
        figure = px.violin(
            frame.dropna(),
            x=_label(spec["x"]),
            y=_label(spec["y"]),
            color=color_name,
            box=True,
            points="outliers",
            title=spec["title"],
        )
    elif plot_type == "count_bar":
        x_name = _label(spec["x"])
        values = field_values(adata, spec["x"])
        counts = (
            values.astype("string")
            .fillna("Missing")
            .value_counts(dropna=False)
            .rename_axis(x_name)
            .rename("count")
            .reset_index()
        )
        figure = px.bar(
            counts,
            x=x_name,
            y="count",
            text="count",
            title=f"{spec['title']} (n={len(values)})",
            labels={"count": "Number of observations"},
        )
        figure.update_traces(textposition="outside", cliponaxis=False)
    elif plot_type == "stacked_bar":
        x_name, color_name = _label(spec["x"]), _label(spec["color"])
        frame = pd.DataFrame(
            {
                x_name: field_values(adata, spec["x"]),
                color_name: field_values(adata, spec["color"]),
            }
        ).dropna()
        counts = (
            frame.groupby([x_name, color_name], observed=True)
            .size()
            .rename("count")
            .reset_index()
        )
        totals = counts.groupby(x_name)["count"].transform("sum")
        counts["percent"] = 100 * counts["count"] / totals
        figure = px.bar(
            counts,
            x=x_name,
            y="percent",
            color=color_name,
            custom_data=["count"],
            title=spec["title"],
            labels={"percent": "Percent of observations"},
        )
        figure.update_traces(hovertemplate="%{y:.1f}% (n=%{customdata[0]})<extra></extra>")
    elif plot_type == "scatter":
        x_name, y_name = _label(spec["x"]), _label(spec["y"])
        frame = pd.DataFrame(
            {
                x_name: field_values(adata, spec["x"]),
                y_name: field_values(adata, spec["y"]),
            }
        )
        color_name = None
        if spec.get("color"):
            color_name = _label(spec["color"])
            frame[color_name] = field_values(adata, spec["color"])
        figure = px.scatter(
            frame.dropna(),
            x=x_name,
            y=y_name,
            color=color_name,
            trendline=None,
            title=spec["title"],
        )
    else:
        group_key = spec["group_by"]
        group_name = _label(group_key)
        frame = pd.DataFrame({group_name: field_values(adata, group_key)})
        feature_names = []
        for key in spec["features"]:
            name = _label(key)
            feature_names.append(name)
            frame[name] = field_values(adata, key)
        means = frame.dropna(subset=[group_name]).groupby(group_name, observed=True)[feature_names].mean()
        figure = px.imshow(
            means.T,
            aspect="auto",
            color_continuous_scale="Viridis",
            labels={"x": group_name, "y": "Feature", "color": "Mean"},
            title=spec["title"],
        )
    figure.update_layout(template="plotly_white", margin={"l": 60, "r": 30, "t": 70, "b": 55})
    return figure


def register_ai_explorer_callbacks(app, adata, prefix, default_features=None):
    """Register one isolated callback; deleting this package removes the demo."""

    @app.callback(
        Output(f"{prefix}-ai-explorer-results", "children"),
        Input(f"{prefix}-ai-explorer-submit", "n_clicks"),
        State(f"{prefix}-ai-explorer-question", "value"),
        prevent_initial_call=True,
    )
    def answer_data_question(_n_clicks, question):
        question = (question or "").strip()
        if not question:
            return html.Div("Enter a question first.", className="alert alert-warning")
        try:
            profile = build_data_profile(adata, question, default_features)
            plan, mode = plan_question(question, profile)
        except Exception as exc:
            return html.Div(
                [
                    html.Strong("The question could not be planned."),
                    html.Div(str(exc), className="small mt-1"),
                ],
                className="alert alert-danger",
            )

        mode_note = html.Div(f"Planner: {mode}", className="text-muted small mb-2")
        if not plan["answerable"]:
            missing = plan.get("missing_data") or []
            return html.Div(
                [
                    mode_note,
                    html.Div(
                        [
                            html.Strong("This data cannot answer that question."),
                            html.Div(plan["answer_summary"], className="mt-1"),
                            html.Ul([html.Li(item) for item in missing], className="mb-0 mt-2")
                            if missing
                            else None,
                        ],
                        className="alert alert-warning",
                    ),
                ]
            )

        plot_cards = []
        for index, spec in enumerate(plan["plots"]):
            try:
                figure = _plot_figure(adata, spec)
            except Exception as exc:
                return html.Div(
                    [
                        mode_note,
                        html.Div(
                            f"The validated plot could not be rendered: {exc}",
                            className="alert alert-danger",
                        ),
                    ]
                )
            plot_cards.append(
                html.Div(
                    [
                        html.P(spec["reason"], className="text-muted small mb-1"),
                        dcc.Graph(
                            id=f"{prefix}-ai-explorer-plot-{index}",
                            figure=figure,
                            config={"displaylogo": False, "responsive": True},
                        ),
                    ],
                    className="border rounded p-2 mb-3",
                )
            )
        return [
            mode_note,
            html.Div(plan["answer_summary"], className="alert alert-success"),
            *plot_cards,
        ]
