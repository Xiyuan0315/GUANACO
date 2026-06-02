from dash import Input, Output, no_update

from guanaco.pages.matrix.plots.grn_demo import (
    build_grn_cytoscape,
    default_grn_edge_threshold,
    grn_dataframe,
    grn_has_weight,
)


def register_grn_demo_callbacks(app, adata, prefix):
    @app.callback(
        Output(f"{prefix}-grn-demo", "children"),
        Input(f"{prefix}-grn-demo-context", "value"),
        Input(f"{prefix}-grn-demo-layout", "value"),
        Input(f"{prefix}-grn-demo-threshold", "value"),
        Input(f"{prefix}-single-cell-tabs", "value"),
    )
    def update_grn_demo(selected_context, layout_name, edge_threshold, active_tab):
        if active_tab != "grn-tab":
            return no_update

        try:
            grn_df = grn_dataframe(adata)
            has_weight = grn_has_weight(grn_df)
            resolved_threshold = edge_threshold if edge_threshold is not None else default_grn_edge_threshold(grn_df)
        except Exception:
            has_weight = False
            resolved_threshold = None

        return build_grn_cytoscape(
            adata,
            component_id=f"{prefix}-grn-demo-cytoscape-view",
            selected_context=selected_context,
            edge_threshold=float(resolved_threshold) if has_weight and resolved_threshold is not None else None,
            layout_name=layout_name or "cose",
        )
