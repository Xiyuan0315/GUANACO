"""Callbacks for the self-contained exploratory Network tab."""

from __future__ import annotations

from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from guanaco.data.network_sources import (
    NetworkEdge,
    NetworkSourceError,
    fetch_network_edges,
    fetch_trrust_table,
)
from guanaco.pages.matrix.analysis.network import (
    MAX_NETWORK_EDGES,
    MAX_NETWORK_NODES,
    build_network_payload,
    trrust_key_regulator_enrichment,
    validate_gene_list,
)
from guanaco.pages.matrix.plots.network import (
    build_key_regulator_component,
    build_network_cytoscape,
    network_highlight_stylesheet,
    network_stylesheet,
)
from guanaco.utils.obs_utils import selected_cell_view
from guanaco.utils.render_guard import signature


def _network_status(graph):
    notices = []
    if graph.get("network_type") == "mirna":
        retained = graph.get("retained_regulator_count", 0)
        eligible = graph.get("eligible_regulator_count", 0)
        if retained == 0:
            notices.append("No shared miRNA regulators were found for these genes.")
        else:
            notices.append(
                f"Showing the top {retained} of {eligible} shared miRNA regulators."
            )
    if graph["truncated"]:
        notices.append(
            "Showing "
            f"{len(graph['nodes'])} of {graph['available_node_count']} nodes and "
            f"{len(graph['edges'])} of {graph['available_edge_count']} interactions. "
            "The view was limited to keep the network responsive."
        )
    missing_count = len(graph["missing_expression"])
    if missing_count:
        notices.append(
            f"{missing_count} database gene node(s) are not measured in this RNA matrix."
        )
    return " ".join(notices)


def _build_view(
    adata,
    source,
    view_mode,
    ppi_confidence=0.7,
    tf_direction="targets",
):
    edges = [NetworkEdge(**edge) for edge in source["edges"]]
    threshold = min(0.9, max(0.4, float(ppi_confidence or 0.7)))
    if source["network_type"] == "ppi":
        edges = [
            edge
            for edge in edges
            if edge.score is not None and float(edge.score) >= threshold
        ]
    resolved_tf_direction = (
        tf_direction if tf_direction in {"targets", "regulators", "both"} else "targets"
    )
    if source["network_type"] == "tf-gene" and resolved_tf_direction != "both":
        query = {gene.casefold() for gene in source["query_genes"]}
        if resolved_tf_direction == "targets":
            edges = [edge for edge in edges if edge.source.casefold() in query]
        else:
            edges = [edge for edge in edges if edge.target.casefold() in query]
    graph = build_network_payload(
        adata,
        source["query_genes"],
        source["network_type"],
        edges,
        view_mode=view_mode or "first-order",
        max_nodes=MAX_NETWORK_NODES,
        max_edges=MAX_NETWORK_EDGES,
        mirna_regulator_filter="shared-20",
    )
    graph["ppi_network_type"] = source.get("ppi_network_type", "functional")
    return graph


def register_network_callbacks(
    app,
    adata,
    prefix: str,
    *,
    organism: str = "human",
    resolve_plot_adata_from_filter,
    hash_list_signature,
    edge_fetcher=fetch_network_edges,
    regulator_fetcher=fetch_trrust_table,
):
    def render_graph(graph, layout_name, continuous_colormap):
        return build_network_cytoscape(
            graph,
            component_id=f"{prefix}-network-cytoscape-view",
            layout_name=layout_name or "cose",
            continuous_colormap=continuous_colormap or "viridis",
        )

    @app.callback(
        Output(f"{prefix}-network-options-collapse", "is_open"),
        Output(f"{prefix}-network-options-toggle", "children"),
        Input(f"{prefix}-network-options-toggle", "n_clicks"),
        State(f"{prefix}-network-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_options(_n_clicks, is_open):
        now_open = not is_open
        return now_open, "▾ More options" if now_open else "▸ More options"

    @app.callback(
        Output(f"{prefix}-network-graph", "children"),
        Output(f"{prefix}-network-source-store", "data"),
        Output(f"{prefix}-network-graph-store", "data"),
        Output(f"{prefix}-network-rendered-key", "data"),
        Output(f"{prefix}-network-highlight-store", "data"),
        Output(f"{prefix}-network-status", "children"),
        Output(f"{prefix}-network-regulator-section", "style"),
        Output(f"{prefix}-network-regulator-results", "children"),
        Input(f"{prefix}-network-build", "n_clicks"),
        State(f"{prefix}-network-genes", "value"),
        State(f"{prefix}-network-type", "value"),
        State(f"{prefix}-network-string-physical", "value"),
        State(f"{prefix}-network-tf-direction", "value"),
        State(f"{prefix}-network-view", "value"),
        State(f"{prefix}-network-string-confidence", "value"),
        State(f"{prefix}-network-layout", "value"),
        State(f"{prefix}-scatter-color-map-dropdown", "value"),
        State(f"{prefix}-global-filtered-data", "data"),
        State(f"{prefix}-selected-cells-store", "data"),
        State(f"{prefix}-selected-cells-hash", "data"),
        State(f"{prefix}-visualization-workspace-tabs", "value"),
        State(f"{prefix}-exploratory-tabs", "value"),
        State(f"{prefix}-network-rendered-key", "data"),
        prevent_initial_call=True,
    )
    def build_network(
        n_clicks,
        gene_text,
        network_type,
        string_physical,
        tf_direction,
        view_mode,
        ppi_confidence,
        layout_name,
        continuous_colormap,
        filtered_data,
        selected_cells,
        selected_cells_hash,
        active_workspace,
        active_tab,
        rendered_key,
    ):
        if not n_clicks:
            raise PreventUpdate
        if active_workspace != "dataset-exploration" or active_tab != "network-tab":
            raise PreventUpdate
        try:
            genes = validate_gene_list(gene_text)
            plot_adata = resolve_plot_adata_from_filter(filtered_data)
            plot_adata = selected_cell_view(plot_adata, selected_cells)
            if plot_adata.n_obs == 0:
                raise ValueError(
                    "No cells remain after applying Global Data Filter and lasso selection."
                )
            ppi_network_type = (
                "physical" if "physical" in (string_physical or []) else "functional"
            )
            selection_key = selected_cells_hash or (
                hash_list_signature(selected_cells) if selected_cells else None
            )
            cache_key = signature(
                "network",
                organism,
                network_type,
                ppi_network_type if network_type == "ppi" else None,
                sorted(g.casefold() for g in genes),
                hash_list_signature((filtered_data or {}).get("cell_indices")),
                selection_key,
            )
            if cache_key == rendered_key:
                return (no_update,) * 8
            if network_type == "ppi":
                edges = edge_fetcher(
                    network_type,
                    genes,
                    organism,
                    ppi_network_type=ppi_network_type,
                )
            else:
                edges = edge_fetcher(network_type, genes, organism)
            source = {
                "network_type": network_type,
                "ppi_network_type": ppi_network_type,
                "query_genes": genes,
                "edges": [edge.to_dict() for edge in edges],
            }
            graph = _build_view(
                plot_adata,
                source,
                view_mode,
                ppi_confidence,
                tf_direction,
            )
            component = render_graph(graph, layout_name, continuous_colormap)
            regulator_style = {"display": "none"}
            regulator_component = None
            if network_type == "tf-gene":
                enrichment = trrust_key_regulator_enrichment(
                    adata,
                    genes,
                    regulator_fetcher(organism),
                )
                regulator_style = {"display": "block"}
                regulator_component = build_key_regulator_component(enrichment)
        except (ValueError, NetworkSourceError) as exc:
            # Keep the last successful graph visible when a new remote query fails.
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                str(exc),
                no_update,
                no_update,
            )

        return (
            component,
            source,
            graph,
            cache_key,
            None,
            _network_status(graph),
            regulator_style,
            regulator_component,
        )

    @app.callback(
        Output(f"{prefix}-network-graph", "children", allow_duplicate=True),
        Output(f"{prefix}-network-graph-store", "data", allow_duplicate=True),
        Output(f"{prefix}-network-highlight-store", "data", allow_duplicate=True),
        Output(f"{prefix}-network-status", "children", allow_duplicate=True),
        Input(f"{prefix}-network-view", "value"),
        Input(f"{prefix}-network-string-confidence", "value"),
        Input(f"{prefix}-network-tf-direction", "value"),
        State(f"{prefix}-network-source-store", "data"),
        State(f"{prefix}-network-layout", "value"),
        State(f"{prefix}-scatter-color-map-dropdown", "value"),
        State(f"{prefix}-global-filtered-data", "data"),
        State(f"{prefix}-selected-cells-store", "data"),
        prevent_initial_call=True,
    )
    def update_network_view(
        view_mode,
        ppi_confidence,
        tf_direction,
        source,
        layout_name,
        continuous_colormap,
        filtered_data,
        selected_cells,
    ):
        if not source:
            raise PreventUpdate
        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        plot_adata = selected_cell_view(plot_adata, selected_cells)
        if plot_adata.n_obs == 0:
            raise PreventUpdate
        graph = _build_view(
            plot_adata,
            source,
            view_mode,
            ppi_confidence,
            tf_direction,
        )
        return (
            render_graph(graph, layout_name, continuous_colormap),
            graph,
            None,
            _network_status(graph),
        )

    @app.callback(
        Output(f"{prefix}-network-string-confidence-wrap", "style"),
        Output(f"{prefix}-network-string-type-wrap", "style"),
        Output(f"{prefix}-network-tf-direction-wrap", "style"),
        Input(f"{prefix}-network-type", "value"),
    )
    def toggle_network_specific_controls(network_type):
        return (
            {"display": "block"} if network_type == "ppi" else {"display": "none"},
            {"display": "block"} if network_type == "ppi" else {"display": "none"},
            {"display": "block"} if network_type == "tf-gene" else {"display": "none"},
        )

    @app.callback(
        Output(f"{prefix}-network-graph", "children", allow_duplicate=True),
        Output(f"{prefix}-network-highlight-store", "data", allow_duplicate=True),
        Input(f"{prefix}-network-layout", "value"),
        Input(f"{prefix}-scatter-color-map-dropdown", "value"),
        State(f"{prefix}-network-graph-store", "data"),
        prevent_initial_call=True,
    )
    def update_network_appearance(layout_name, continuous_colormap, graph):
        if not graph:
            raise PreventUpdate
        return render_graph(graph, layout_name, continuous_colormap), None

    @app.callback(
        Output(f"{prefix}-network-cytoscape-view", "stylesheet"),
        Output(f"{prefix}-network-highlight-store", "data", allow_duplicate=True),
        Input(f"{prefix}-network-cytoscape-view", "tapNodeData"),
        State(f"{prefix}-network-graph-store", "data"),
        State(f"{prefix}-network-highlight-store", "data"),
        prevent_initial_call=True,
    )
    def highlight_network_node(node_data, graph, current_highlight):
        if not node_data or not graph:
            raise PreventUpdate
        node_id = str(node_data.get("id"))
        if current_highlight and current_highlight.get("node") == node_id:
            return (
                network_stylesheet(),
                None,
            )
        return network_highlight_stylesheet(graph, node_id), {"node": node_id}

    @app.callback(
        Output(f"{prefix}-network-cytoscape-view", "generateImage"),
        Input(f"{prefix}-network-download-svg", "n_clicks"),
        prevent_initial_call=True,
    )
    def download_network_svg(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        return {"type": "svg", "action": "download", "filename": "network"}
