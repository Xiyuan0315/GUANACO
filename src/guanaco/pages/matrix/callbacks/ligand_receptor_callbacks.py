"""Callbacks for precomputed ligand–receptor result visualization."""

from __future__ import annotations

from dash import Input, Output, State, ctx
from dash.exceptions import PreventUpdate

from guanaco.data.ligand_receptor import (
    LigandReceptorResultError,
    load_default_ligand_receptor_result,
)
from guanaco.pages.matrix.plots.ligand_receptor import (
    _empty,
    build_ligand_receptor_network,
    build_ligand_receptor_view,
    ligand_receptor_detail,
    ligand_receptor_highlight_stylesheet,
    ligand_receptor_stylesheet,
    metric_slider_settings,
    plot_ligand_receptor_dotplot,
)
from guanaco.utils.colors import (
    resolve_continuous_colorscale,
    resolve_discrete_palette,
)


def _triggered_property() -> str:
    try:
        return next(iter(ctx.triggered_prop_ids), "")
    except Exception:
        return ""


def _selection_candidate(triggered, node, edge, click_data):
    if triggered.endswith(".tapNodeData") and node:
        return {"kind": "node", "node_id": str(node.get("id"))}
    if triggered.endswith(".tapEdgeData") and edge:
        return {"kind": "edge", "edge_id": str(edge.get("id"))}
    if triggered.endswith(".clickData") and click_data:
        points = click_data.get("points") or []
        custom = points[0].get("customdata") if points else None
        try:
            has_interaction = custom is not None and len(custom) > 1
        except TypeError:
            has_interaction = False
        if has_interaction:
            return {
                "kind": "interaction",
                "interaction_id": str(custom[0]),
                "edge_id": str(custom[1]),
            }
    return None


def register_ligand_receptor_callbacks(
    app,
    adata,
    prefix: str,
    *,
    color_config,
):
    payload = load_default_ligand_receptor_result(adata)

    @app.callback(
        Output(f"{prefix}-lr-options-collapse", "is_open"),
        Output(f"{prefix}-lr-options-toggle", "children"),
        Input(f"{prefix}-lr-options-toggle", "n_clicks"),
        State(f"{prefix}-lr-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_options(_n_clicks, is_open):
        now_open = not is_open
        label = "▾ More options" if now_open else "▸ More options"
        return now_open, label

    @app.callback(
        Output(f"{prefix}-lr-magnitude-range", "min"),
        Output(f"{prefix}-lr-magnitude-range", "max"),
        Output(f"{prefix}-lr-magnitude-range", "step"),
        Output(f"{prefix}-lr-magnitude-range", "value"),
        Output(f"{prefix}-lr-magnitude-range", "marks"),
        Output(f"{prefix}-lr-specificity-range-wrap", "style"),
        Output(f"{prefix}-lr-specificity-range", "min"),
        Output(f"{prefix}-lr-specificity-range", "max"),
        Output(f"{prefix}-lr-specificity-range", "step"),
        Output(f"{prefix}-lr-specificity-range", "value"),
        Output(f"{prefix}-lr-specificity-range", "marks"),
        Input(f"{prefix}-lr-magnitude", "value"),
        Input(f"{prefix}-lr-specificity", "value"),
    )
    def update_metric_ranges(magnitude, specificity):
        magnitude_settings = metric_slider_settings(payload, magnitude)
        specificity_settings = metric_slider_settings(payload, specificity)
        return (
            magnitude_settings["min"],
            magnitude_settings["max"],
            magnitude_settings["step"],
            magnitude_settings["value"],
            magnitude_settings["marks"],
            {"display": "block"} if specificity else {"display": "none"},
            specificity_settings["min"],
            specificity_settings["max"],
            specificity_settings["step"],
            specificity_settings["value"],
            specificity_settings["marks"],
        )

    @app.callback(
        Output(f"{prefix}-lr-network", "children"),
        Output(f"{prefix}-lr-view-store", "data"),
        Input(f"{prefix}-lr-magnitude", "value"),
        Input(f"{prefix}-lr-specificity", "value"),
        Input(f"{prefix}-lr-magnitude-range", "value"),
        Input(f"{prefix}-lr-specificity-range", "value"),
        Input(f"{prefix}-lr-senders", "value"),
        Input(f"{prefix}-lr-receivers", "value"),
        Input(f"{prefix}-discrete-color-map-dropdown", "value"),
    )
    def render_circle(
        magnitude,
        specificity,
        magnitude_range,
        specificity_range,
        senders,
        receivers,
        discrete_color_map,
    ):
        if not payload or not magnitude:
            return (
                _empty("Load a precomputed interaction result."),
                None,
            )
        try:
            view = build_ligand_receptor_view(
                payload,
                magnitude=magnitude,
                specificity=specificity,
                magnitude_range=magnitude_range,
                specificity_range=specificity_range,
                senders=senders,
                receivers=receivers,
            )
        except (ValueError, LigandReceptorResultError) as exc:
            return _empty(str(exc), error=True), None

        palette = resolve_discrete_palette(
            discrete_color_map,
            max(1, len(view["groups"])),
            default=color_config,
        ) or ["#808080"]
        view["node_colors"] = {
            group: palette[index % len(palette)]
            for index, group in enumerate(view["groups"])
        }
        component = build_ligand_receptor_network(
            view,
            component_id=f"{prefix}-lr-cytoscape-view",
        )
        return component, view

    @app.callback(
        Output(f"{prefix}-lr-selection-store", "data"),
        Input(f"{prefix}-lr-view-store", "data"),
        Input(f"{prefix}-lr-cytoscape-view", "tapNodeData"),
        Input(f"{prefix}-lr-cytoscape-view", "tapEdgeData"),
        Input(f"{prefix}-lr-dotplot", "clickData"),
        State(f"{prefix}-lr-selection-store", "data"),
        prevent_initial_call=True,
    )
    def update_selection(_view, node, edge, click_data, current):
        triggered = _triggered_property()
        if triggered == f"{prefix}-lr-view-store.data":
            return None
        selection = _selection_candidate(triggered, node, edge, click_data)
        if not selection:
            raise PreventUpdate
        return None if selection == current else selection

    @app.callback(
        Output(f"{prefix}-lr-dotplot", "figure"),
        Input(f"{prefix}-lr-view-store", "data"),
        Input(f"{prefix}-lr-selection-store", "data"),
        Input(f"{prefix}-scatter-color-map-dropdown", "value"),
    )
    def update_dotplot(view, selection, color_map):
        selection = selection or {}
        if (
            _triggered_property() == f"{prefix}-lr-selection-store.data"
            and selection.get("kind") == "interaction"
        ):
            raise PreventUpdate
        selected_node = (
            selection.get("node_id") if selection.get("kind") == "node" else None
        )
        selected_edge = (
            selection.get("edge_id")
            if selection.get("kind") == "edge"
            else None
        )
        return plot_ligand_receptor_dotplot(
            view,
            selected_node=selected_node,
            selected_edge=selected_edge,
            colorscale=resolve_continuous_colorscale(color_map),
        )

    @app.callback(
        Output(f"{prefix}-lr-detail", "children"),
        Input(f"{prefix}-lr-selection-store", "data"),
        Input(f"{prefix}-lr-view-store", "data"),
    )
    def inspect_interaction(selection, view):
        if not view:
            return _empty("Load a result to inspect ligand–receptor pairs.")
        selection = selection or {}
        if selection.get("kind") == "interaction":
            return ligand_receptor_detail(
                view,
                interaction_id=selection.get("interaction_id"),
            )
        if selection.get("kind") == "edge":
            return ligand_receptor_detail(
                view,
                edge_id=selection.get("edge_id"),
            )
        if selection.get("kind") == "node":
            return ligand_receptor_detail(
                view,
                node_id=selection.get("node_id"),
            )
        return ligand_receptor_detail(view)

    @app.callback(
        Output(f"{prefix}-lr-cytoscape-view", "stylesheet"),
        Input(f"{prefix}-lr-selection-store", "data"),
        Input(f"{prefix}-lr-view-store", "data"),
        prevent_initial_call=True,
    )
    def highlight_selection(selection, view):
        if not view:
            raise PreventUpdate
        selection = selection or {}
        if selection.get("kind") == "node":
            return ligand_receptor_highlight_stylesheet(
                view,
                node_id=selection.get("node_id"),
            )
        if selection.get("kind") in {"edge", "interaction"}:
            return ligand_receptor_highlight_stylesheet(
                view,
                edge_id=selection.get("edge_id"),
            )
        return ligand_receptor_stylesheet()

    @app.callback(
        Output(f"{prefix}-lr-cytoscape-view", "generateImage"),
        Input(f"{prefix}-lr-download-svg", "n_clicks"),
        prevent_initial_call=True,
    )
    def download_svg(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        return {"type": "svg", "action": "download", "filename": "ligand_receptor"}
