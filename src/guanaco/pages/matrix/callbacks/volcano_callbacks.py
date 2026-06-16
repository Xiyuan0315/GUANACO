from dash import Input, Output, State, dcc, html, no_update
from dash.exceptions import PreventUpdate

from guanaco.utils.render_guard import signature

from guanaco.pages.matrix.plots.volcano import (
    DEFAULT_PADJ_THRESHOLD,
    DEFAULT_TOP_N,
    DEFAULT_X_THRESHOLD,
    deg_csv,
    deg_summary,
    empty_volcano_figure,
    load_volcano_payload,
    plot_volcano,
    volcano_degs_filename,
    x_axis_options,
)


_DEFAULT_X_FIELD = "logfoldchange"


def _summary_component(summary):
    return html.Div(
        [
            html.H5("DEGs after threshold", style={"marginTop": 0}),
            html.Div(f"{summary['total']} genes left"),
            html.Div(f"Up: {summary['up']}  Down: {summary['down']}"),
            html.Div(
                f"Criteria: {summary['criteria']}",
                style={"color": "#52606d", "marginTop": "6px", "fontSize": "12px"},
            ),
        ]
    )


def _resolve_x_field(entry, x_field):
    valid_x_fields = {option["value"] for option in x_axis_options(entry)}
    return x_field if x_field in valid_x_fields else _DEFAULT_X_FIELD


def _resolve_deg_thresholds(padj_threshold, x_threshold):
    return (
        float(padj_threshold if padj_threshold is not None else DEFAULT_PADJ_THRESHOLD),
        float(x_threshold if x_threshold is not None else DEFAULT_X_THRESHOLD),
    )


def _resolve_top_n(top_n):
    return int(top_n if top_n is not None else DEFAULT_TOP_N)


def register_volcano_callbacks(app, adata, prefix):
    @app.callback(
        Output(f"{prefix}-volcano-x-axis-dropdown", "options"),
        Output(f"{prefix}-volcano-x-axis-dropdown", "value"),
        Input(f"{prefix}-volcano-entry-dropdown", "value"),
        State(f"{prefix}-volcano-x-axis-dropdown", "value"),
    )
    def update_volcano_x_axis_options(entry_name, current_x_field):
        if not entry_name:
            return x_axis_options(), _DEFAULT_X_FIELD

        try:
            payload = load_volcano_payload(adata)
            entry = payload["entries"][entry_name]
        except Exception:
            return x_axis_options(), _DEFAULT_X_FIELD

        options = x_axis_options(entry)
        valid_values = {option["value"] for option in options}
        value = current_x_field if current_x_field in valid_values else _DEFAULT_X_FIELD
        return options, value

    @app.callback(
        Output(f"{prefix}-volcano-plot", "figure"),
        Output(f"{prefix}-volcano-deg-summary", "children"),
        Output(f"{prefix}-volcano-rendered-key", "data"),
        Input(f"{prefix}-volcano-entry-dropdown", "value"),
        Input(f"{prefix}-volcano-x-axis-dropdown", "value"),
        Input(f"{prefix}-volcano-padj-threshold", "value"),
        Input(f"{prefix}-volcano-x-threshold", "value"),
        Input(f"{prefix}-volcano-top-n", "value"),
        Input(f"{prefix}-single-cell-tabs", "value"),
        State(f"{prefix}-volcano-plot", "figure"),
        State(f"{prefix}-volcano-rendered-key", "data"),
    )
    def update_volcano_plot(
        entry_name,
        x_field,
        padj_threshold,
        x_threshold,
        top_n,
        active_tab,
        current_figure,
        rendered_key,
    ):
        if active_tab != "volcano-tab":
            return no_update, no_update, no_update
        if not entry_name:
            return (
                empty_volcano_figure(
                    "No volcano or rank_genes_groups result is available."
                ),
                "No DE result selected.",
                None,
            )

        cache_key = signature(
            "volcano", entry_name, x_field, padj_threshold, x_threshold, top_n
        )
        if cache_key == rendered_key and current_figure:
            return no_update, no_update, no_update

        try:
            payload = load_volcano_payload(adata)
            entry = payload["entries"][entry_name]
            x_field = _resolve_x_field(entry, x_field)
            resolved_padj_threshold, resolved_x_threshold = _resolve_deg_thresholds(
                padj_threshold, x_threshold
            )
            resolved_top_n = _resolve_top_n(top_n)
            fig = plot_volcano(
                entry_name=entry_name,
                entry=entry,
                x_field=x_field,
                padj_threshold=resolved_padj_threshold,
                x_threshold=resolved_x_threshold,
                top_n=resolved_top_n,
            )
            summary = _summary_component(
                deg_summary(
                    entry, x_field, resolved_padj_threshold, resolved_x_threshold
                )
            )
            return fig, summary, cache_key
        except Exception as exc:
            return (
                empty_volcano_figure(str(exc)),
                html.Div(str(exc), style={"color": "#b42318"}),
                None,
            )

    @app.callback(
        Output(f"{prefix}-volcano-degs-download", "data"),
        Input(f"{prefix}-volcano-download-button", "n_clicks"),
        State(f"{prefix}-volcano-entry-dropdown", "value"),
        State(f"{prefix}-volcano-x-axis-dropdown", "value"),
        State(f"{prefix}-volcano-padj-threshold", "value"),
        State(f"{prefix}-volcano-x-threshold", "value"),
        prevent_initial_call=True,
    )
    def download_volcano_degs(
        n_clicks, entry_name, x_field, padj_threshold, x_threshold
    ):
        if not n_clicks or not entry_name:
            raise PreventUpdate

        payload = load_volcano_payload(adata)
        entry = payload["entries"][entry_name]
        x_field = _resolve_x_field(entry, x_field)
        resolved_padj_threshold, resolved_x_threshold = _resolve_deg_thresholds(
            padj_threshold, x_threshold
        )
        filename = volcano_degs_filename(
            entry_name,
            entry,
            x_field,
            resolved_padj_threshold,
            resolved_x_threshold,
        )
        return dcc.send_string(
            deg_csv(entry, x_field, resolved_padj_threshold, resolved_x_threshold),
            filename,
        )
