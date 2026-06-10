from dash import Input, Output, State, no_update

from guanaco.pages.matrix.plots.grn_demo import (
    build_grn_cytoscape,
    default_grn_edge_threshold,
    grn_dataframe,
    grn_has_weight,
)
from guanaco.utils.render_guard import signature


def register_grn_demo_callbacks(app, adata, prefix):
    @app.callback(
        Output(f"{prefix}-grn-demo", "children"),
        Output(f"{prefix}-grn-demo-rendered-key", "data"),
        Input(f"{prefix}-grn-demo-context", "value"),
        Input(f"{prefix}-grn-demo-layout", "value"),
        Input(f"{prefix}-grn-demo-threshold", "value"),
        Input(f"{prefix}-single-cell-tabs", "value"),
        State(f"{prefix}-grn-demo-rendered-key", "data"),
    )
    def update_grn_demo(selected_context, layout_name, edge_threshold, active_tab, rendered_key):
        if active_tab != "grn-tab":
            return no_update, no_update

        cache_key = signature("grn", selected_context, layout_name, edge_threshold)
        if cache_key == rendered_key:
            return no_update, no_update

        try:
            grn_df = grn_dataframe(adata)
            has_weight = grn_has_weight(grn_df)
            resolved_threshold = edge_threshold if edge_threshold is not None else default_grn_edge_threshold(grn_df)
        except Exception:
            has_weight = False
            resolved_threshold = None

        component = build_grn_cytoscape(
            adata,
            component_id=f"{prefix}-grn-demo-cytoscape-view",
            selected_context=selected_context,
            edge_threshold=float(resolved_threshold) if has_weight and resolved_threshold is not None else None,
            layout_name=layout_name or "cose",
        )
        return component, cache_key

    # Download the GRN network as SVG (the Cytoscape equivalent of the Plotly
    # camera icon). Setting generateImage triggers a client-side download.
    @app.callback(
        Output(f"{prefix}-grn-demo-cytoscape-view", "generateImage"),
        Input(f"{prefix}-grn-demo-download-svg", "n_clicks"),
        prevent_initial_call=True,
    )
    def download_grn_svg(n_clicks):
        if not n_clicks:
            return no_update
        return {"type": "svg", "action": "download", "filename": "grn_network"}
