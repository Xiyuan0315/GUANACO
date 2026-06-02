from dash import Input, Output, State, html, no_update

from guanaco.utils.colors import resolve_discrete_palette


def _empty_paga_component(message):
    return html.Div(
        message,
        style={
            "padding": "16px",
            "margin": "8px",
            "color": "#4A4A4A",
            "backgroundColor": "#F7F8FA",
            "border": "1px solid #D8DDE6",
            "borderRadius": "6px",
        },
    )


def register_paga_callbacks(
    app,
    adata,
    prefix,
    *,
    build_paga_cytoscape,
    color_config,
):
    @app.callback(
        Output(f"{prefix}-paga-gene-dropdown", "options"),
        Input(f"{prefix}-paga-gene-dropdown", "search_value"),
        State(f"{prefix}-paga-gene-dropdown", "value"),
    )
    def update_paga_gene_selection(search_value, current_value):
        if not search_value:
            values = [current_value] if current_value else []
            return [{"label": value, "value": value} for value in values]

        matches = [gene for gene in adata.var_names if search_value.lower() in gene.lower()][:25]
        if current_value and current_value not in matches:
            matches = [current_value] + matches
        return [{"label": gene, "value": gene} for gene in matches]

    @app.callback(
        Output(f"{prefix}-paga-hover-detail", "children"),
        Input(f"{prefix}-paga-cytoscape-view", "mouseoverNodeData"),
    )
    def update_paga_hover_detail(node_data):
        if not node_data:
            return ""

        hover_text = node_data.get("hover_text")
        if not hover_text:
            return ""

        lines = []
        for line in str(hover_text).split("<br>"):
            lines.append(html.Div(line))
        return lines

    @app.callback(
        Output(f"{prefix}-paga-obs-wrapper", "style"),
        Output(f"{prefix}-paga-gene-wrapper", "style"),
        Input(f"{prefix}-paga-color-mode", "value"),
        Input(f"{prefix}-paga-obs-dropdown", "value"),
    )
    def toggle_paga_controls(color_mode, obs_key):
        obs_style = {"marginBottom": "10px"} if color_mode == "obs" else {"display": "none", "marginBottom": "10px"}
        gene_style = {"display": "none", "marginBottom": "10px"} if color_mode == "obs" else {"marginBottom": "10px"}
        return obs_style, gene_style

    @app.callback(
        Output(f"{prefix}-paga", "children"),
        [
            Input(f"{prefix}-paga-color-mode", "value"),
            Input(f"{prefix}-paga-obs-dropdown", "value"),
            Input(f"{prefix}-paga-gene-dropdown", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-paga-threshold", "value"),
            Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
            Input(f"{prefix}-single-cell-label-selection", "value"),
            Input(f"{prefix}-selected-cells-store", "data"),
            Input(f"{prefix}-single-cell-tabs", "value"),
        ],
    )
    def update_paga(
        color_mode,
        obs_key,
        gene,
        continuous_colormap,
        discrete_colormap,
        edge_threshold,
        selected_annotation,
        selected_labels,
        selected_cells,
        active_tab,
    ):
        if active_tab != "paga-tab":
            return no_update

        color_mode = color_mode or "obs"

        if color_mode == "gene" and not gene:
            return _empty_paga_component("Select a gene to color the PAGA graph.")
        if color_mode == "obs" and not obs_key:
            return _empty_paga_component("Select an obs column to color the PAGA graph.")

        discrete_palette = None
        if color_mode == "obs":
            n_colors = adata.obs[obs_key].nunique() if obs_key in adata.obs.columns else 0
            discrete_palette = resolve_discrete_palette(
                discrete_colormap,
                n_colors,
                default=color_config,
            )

        edge_threshold = float(edge_threshold if edge_threshold is not None else 0.03)

        return build_paga_cytoscape(
            adata,
            component_id=f"{prefix}-paga-cytoscape-view",
            color_mode=color_mode,
            obs_key=obs_key,
            gene=gene,
            continuous_color_map=continuous_colormap or "Viridis",
            discrete_palette=discrete_palette,
            edge_threshold=edge_threshold,
            selected_annotation=selected_annotation,
            selected_labels=selected_labels,
            selected_cells=selected_cells,
        )
