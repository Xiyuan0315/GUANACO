import json
from pathlib import Path

import dash_bootstrap_components as dbc
import plotly.express as px
from dash import dcc, html

from guanaco.config import scatter_config, gene_scatter_config


def _palette_names():
    cvd_color_path = Path(__file__).resolve().parents[1] / "cvd_color.json"
    with open(cvd_color_path, "r") as f:
        palette_json = json.load(f)

    plotly_qualitative_palettes = {}
    for name in dir(px.colors.qualitative):
        if name.startswith("_"):
            continue
        palette = getattr(px.colors.qualitative, name, None)
        if isinstance(palette, (list, tuple)) and palette and all(isinstance(c, str) for c in palette):
            plotly_qualitative_palettes[name] = list(palette)
    for k, v in plotly_qualitative_palettes.items():
        if k not in palette_json.get("color_palettes", {}):
            palette_json["color_palettes"][k] = v
    try:
        import colorcet as cc

        for name, palette in cc.palette.items():
            if "glasbey" not in name.lower():
                continue
            key = f"colorcet/{name}"
            if key not in palette_json.get("color_palettes", {}):
                palette_json["color_palettes"][key] = list(palette)
    except Exception:
        # colorcet is optional; keep existing palettes if unavailable.
        pass
    return list(palette_json["color_palettes"].keys())


def _continuous_colormap_options():
    """Plotly + colorcet continuous colormap options."""
    options = [{"label": scale, "value": scale} for scale in px.colors.named_colorscales()]
    try:
        import colorcet as cc

        cc_names = sorted(
            name for name in cc.palette.keys()
            if "glasbey" not in name.lower()
        )
        options.extend(
            {"label": f"colorcet/{name}", "value": f"cc:{name}"}
            for name in cc_names
        )
    except Exception:
        # colorcet is optional; keep Plotly maps if unavailable.
        pass
    return options


def initialize_scatter_components(adata):
    embedding_prefixes = {
        "X_umap": "UMAP",
        "X_pca": "PCA",
        "X_tsne": "t-SNE",
        "X_diffmap": "DiffMap",
        "X_phate": "PHATE",
        "X_draw_graph_fa": "FA",
    }

    obsm_list = list(adata.obsm.keys())
    embedding_columns = {
        key: [f"{embedding_prefixes.get(key, key.upper())}{i+1}" for i in range(adata.obsm[key].shape[1])]
        for key in obsm_list
    }
    default_embedding = obsm_list[-1]
    default_columns = embedding_columns[default_embedding]
    return obsm_list, embedding_columns, default_embedding, default_columns


def create_control_components(adata, prefix):
    obsm_list, _, default_embedding, default_columns = initialize_scatter_components(adata)
    clustering_dropdown = dcc.Dropdown(
        id=f"{prefix}-clustering-dropdown",
        options=[{"label": key.replace("X_", "").upper(), "value": key} for key in obsm_list],
        value="X_umap" if "X_umap" in obsm_list else default_embedding,
        placeholder="Select Clustering Method",
        style={"marginBottom": "15px", "fontSize": "14px"},
        clearable=False,
    )

    coordinates_dropdowns = html.Div(
        [
            html.Label("X-axis:", style={"fontWeight": "bold", "marginBottom": "5px"}),
            dcc.Dropdown(
                id=f"{prefix}-x-axis",
                options=[{"label": col, "value": col} for col in default_columns],
                value=default_columns[0],
            ),
            html.Label("Y-axis:", style={"fontWeight": "bold", "marginTop": "15px", "marginBottom": "5px"}),
            dcc.Dropdown(
                id=f"{prefix}-y-axis",
                options=[{"label": col, "value": col} for col in default_columns],
                value=default_columns[1],
            ),
        ],
        style={"marginBottom": "15px", "fontSize": "14px"},
    )

    return clustering_dropdown, coordinates_dropdowns


def generate_annotation_dropdown(anno_list, prefix, default_value=None):
    return dcc.Dropdown(
        id=f"{prefix}-annotation-dropdown",
        options=[{"label": label, "value": label} for label in anno_list],
        placeholder="Search annotations or genes...",
        value=default_value if default_value else (anno_list[0] if anno_list else None),
        style={"marginBottom": "15px"},
    )


def generate_scatter_gene_selection(combined_list, prefix, default_value=None):
    return dcc.Dropdown(
        id=f"{prefix}-scatter-gene-selection",
        options=[{"label": label, "value": label} for label in combined_list],
        value=default_value if default_value else (combined_list[0] if combined_list else None),
        placeholder="Search annotations or genes...",
        style={"marginBottom": "15px"},
    )


def create_global_metadata_filter(adata, prefix):
    categorical_columns = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype == "category" or adata.obs[col].dtype == "object":
            unique_vals = adata.obs[col].unique()
            if len(unique_vals) < 100:
                categorical_columns.append(col)

    filter_components = []
    for col in categorical_columns:
        unique_values = sorted([str(val) for val in adata.obs[col].unique()])
        filter_components.append(
            html.Div(
                [
                    html.Label(f"{col}:", style={"fontWeight": "bold", "marginRight": "10px", "minWidth": "120px"}),
                    dcc.Dropdown(
                        id={"type": f"{prefix}-global-metadata-filter", "column": col},
                        options=[{"label": val, "value": val} for val in unique_values],
                        value=unique_values,
                        multi=True,
                        placeholder=f"Select {col}...",
                        style={"flex": "1", "minWidth": "200px"},
                    ),
                ],
                style={"display": "flex", "alignItems": "center", "marginBottom": "8px", "padding": "5px"},
            )
        )

    return html.Div(
        [
            html.Div(
                [
                    html.H5("🔍 Global Data Filter", style={"margin": "0", "color": "#2c3e50", "display": "inline-block"}),
                    html.Div(
                        [
                            html.Span("Active cells: ", style={"fontWeight": "bold"}),
                            html.Span(
                                id=f"{prefix}-global-cell-count",
                                children=f"{adata.n_obs:,}",
                                style={"color": "#27ae60", "fontWeight": "bold", "fontSize": "16px"},
                            ),
                            html.Span(f" / {adata.n_obs:,} total", style={"color": "#7f8c8d"}),
                            html.Span(
                                id=f"{prefix}-filter-preview",
                                children="",
                                style={
                                    "marginLeft": "15px",
                                    "color": "#e67e22",
                                    "fontSize": "14px",
                                    "fontStyle": "italic",
                                },
                            ),
                        ],
                        style={"display": "inline-block", "marginLeft": "20px"},
                    ),
                ],
                style={"marginBottom": "15px"},
            ),
            dbc.Collapse(
                [
                    html.Div(filter_components, style={"maxHeight": "300px", "overflowY": "auto"}),
                    html.Div(
                        [
                            dbc.Button("Select All", id=f"{prefix}-select-all-filters", color="success", size="sm", style={"marginRight": "10px"}),
                            dbc.Button("Clear All", id=f"{prefix}-clear-all-filters", color="warning", size="sm", style={"marginRight": "10px"}),
                            dbc.Button(
                                "Apply Filter",
                                id=f"{prefix}-apply-global-filter",
                                color="primary",
                                size="sm",
                                style={"fontWeight": "bold", "minWidth": "100px"},
                            ),
                        ],
                        style={"textAlign": "center", "marginTop": "15px"},
                    ),
                ],
                id=f"{prefix}-global-filter-collapse",
                is_open=False,
            ),
            html.Div(
                [
                    dbc.Button(
                        "▼ Show Filters",
                        id=f"{prefix}-toggle-global-filter",
                        color="link",
                        size="sm",
                        style={"padding": "5px 10px", "textDecoration": "none"},
                    )
                ],
                style={"textAlign": "center", "marginTop": "10px"},
            ),
            dcc.Store(id=f"{prefix}-global-filtered-data", data={"cell_indices": None, "n_cells": adata.n_obs}),
        ],
        style={
            "backgroundColor": "#f8f9fa",
            "border": "2px solid #dee2e6",
            "borderRadius": "10px",
            "padding": "15px",
            "marginBottom": "20px",
            "boxShadow": "0 2px 4px rgba(0,0,0,0.1)",
        },
    )


def generate_embedding_plots(adata, prefix):
    palette_names = _palette_names()
    scatter_transformation_selection = html.Div(
        [
            dbc.RadioItems(
                id=f"{prefix}-scatter-log-or-zscore",
                options=[{"label": "None", "value": None}, {"label": "Log", "value": "log"}],
                value=None,
                inline=True,
                style={"fontSize": "14px"},
            )
        ]
    )

    scatter_order_selection = html.Div(
        [
            dbc.RadioItems(
                id=f"{prefix}-plot-order",
                options=[
                    {"label": "Max-1st", "value": "max"},
                    {"label": "Min-1st", "value": "min"},
                    {"label": "Original", "value": "original"},
                    {"label": "Random", "value": "random"},
                ],
                value="max",
                inline=True,
                style={"fontSize": "14px"},
            )
        ]
    )

    colorscales = _continuous_colormap_options()
    color_map_dropdown = dcc.Dropdown(
        id=f"{prefix}-scatter-color-map-dropdown",
        options=colorscales,
        value="viridis",
        style={"marginBottom": "10px"},
        clearable=False,
    )

    color_map_discrete_dropdown = dcc.Dropdown(
        id=f"{prefix}-discrete-color-map-dropdown",
        options=[{"label": name, "value": name} for name in palette_names],
        value=None,
        placeholder="Default color",
        style={"marginBottom": "10px"},
        clearable=True,
        className="custom-dropdown",
    )

    marker_size_slider = dcc.Slider(
        id=f"{prefix}-marker-size-slider",
        min=1,
        max=10,
        value=5,
        marks=None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider",
    )
    opacity_slider = dcc.Slider(
        id=f"{prefix}-opacity-slider",
        min=0.1,
        max=1.0,
        value=1,
        marks=None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider",
    )
    scatter_legend_toggle = dbc.RadioItems(
        id=f"{prefix}-scatter-legend-toggle",
        options=[{"label": "on data", "value": "on data"}, {"label": "right", "value": "right"}],
        value="right",
        inline=True,
        style={"fontSize": "14px"},
    )
    axis_toggle = dbc.RadioItems(
        id=f"{prefix}-axis-toggle",
        options=[{"label": "Show", "value": True}, {"label": "Hide", "value": False}],
        value=True,
        inline=True,
        style={"fontSize": "14px"},
    )
    graphic_control = html.Div(
        id=f"{prefix}-controls-container",
        children=[
            html.Div([html.Label("Plot order: ", className="control-label"), scatter_order_selection], style={"marginBottom": "15px"}),
            html.Div([html.Label("Continous ColorMap:", className="control-label"), color_map_dropdown], className="dbc", style={"marginBottom": "15px"}),
            html.Div([html.Label("Discrete ColorMap:", className="control-label"), color_map_discrete_dropdown], className="dbc", style={"marginBottom": "15px"}),
            html.Div([html.Label("Marker Size:", className="control-label"), marker_size_slider], style={"marginBottom": "15px"}),
            html.Div([html.Label("Opacity:", className="control-label"), opacity_slider], style={"marginBottom": "15px"}),
            html.Div([html.Label("Legend Location:", className="control-label"), scatter_legend_toggle], style={"marginBottom": "15px"}),
            html.Div([html.Label("Axis:", className="control-label"), axis_toggle], style={"marginBottom": "15px"}),
        ],
        style={"display": "none"},
    )

    annotations = adata.obs.columns.tolist()
    sample_genes = adata.var_names[:20].tolist()
    combined_list = annotations + sample_genes
    default_annotation = annotations[0] if annotations else None
    # Use an obs annotation as default on the right plot to avoid heavy
    # gene-expression/datashader work at initial page load.
    default_gene = annotations[0] if annotations else (sample_genes[0] if sample_genes else None)

    clustering_dropdown, coordinates_dropdowns = create_control_components(adata, prefix)
    spatial_imgkey_control = html.Div(
        [
            html.Label("Spatial Image:", className="control-label"),
            dcc.Dropdown(
                id=f"{prefix}-spatial-imgkey-dropdown",
                options=[],
                value=None,
                clearable=False,
                placeholder="Select spatial image key",
                style={"fontSize": "14px"},
            ),
        ],
        id=f"{prefix}-spatial-imgkey-container",
        style={"display": "none", "marginBottom": "15px"},
    )

    return html.Div(
        [
            dcc.Store(id=f"{prefix}-selected-cells-store"),
            create_global_metadata_filter(adata, prefix),
            dbc.Row(
                [
                    dbc.Col(
                        html.Div(
                            [
                                html.Div([html.Label("Dimension Reduction:", className="control-label"), clustering_dropdown], className="dbc", style={"marginBottom": "15px"}),
                                spatial_imgkey_control,
                                html.Div([html.Div(coordinates_dropdowns, id=f"{prefix}-coordinates-dropdowns")], className="dbc", style={"marginBottom": "15px"}),
                                html.Div([html.Label("Transformation:", className="control-label"), scatter_transformation_selection], style={"marginBottom": "15px"}),
                                html.Button("More controls", id=f"{prefix}-toggle-button", n_clicks=0, style={"marginBottom": "10px", "border": "1px solid", "borderRadius": "5px"}),
                                graphic_control,
                            ]
                        ),
                        xs=12,
                        sm=12,
                        md=4,
                        lg=4,
                        xl=2,
                        style={"borderRight": "1px solid #ddd", "padding": "10px"},
                    ),
                    dbc.Col(
                        html.Div(
                            [
                                html.Label("Select Annotation/Variable", style={"fontWeight": "bold", "marginBottom": "5px"}),
                                generate_annotation_dropdown(anno_list=combined_list, prefix=prefix, default_value=default_annotation),
                                html.Div(style={"height": "25px", "marginBottom": "5px"}),
                                dcc.Loading(
                                    id=f"{prefix}-loading-annotaion-scatter",
                                    type="circle",
                                    children=dcc.Graph(id=f"{prefix}-annotation-scatter", config=scatter_config, style={"height": "60vh", "width": "100%"}),
                                    style={"height": "60vh", "width": "100%"},
                                ),
                            ],
                            className="dbc",
                            style={"marginBottom": "20px"},
                        ),
                        xs=12,
                        sm=12,
                        md=4,
                        lg=4,
                        xl=5,
                    ),
                    dbc.Col(
                        html.Div(
                            [
                                html.Label("Select Annotation/Variable:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                                generate_scatter_gene_selection(combined_list=combined_list, prefix=prefix, default_value=default_gene),
                                dbc.RadioItems(
                                    id=f"{prefix}-coexpression-toggle",
                                    options=[{"label": "Single Gene", "value": "single"}, {"label": "Co-expression", "value": "coexpression"}],
                                    value="single",
                                    inline=True,
                                    style={"marginBottom": "10px", "fontSize": "14px"},
                                ),
                                html.Div(
                                    [
                                        html.Label("Second Gene:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                                        dcc.Dropdown(
                                            id=f"{prefix}-scatter-gene2-selection",
                                            options=[{"label": label, "value": label} for label in adata.var_names.to_list()[:10]],
                                            value=adata.var_names.to_list()[1],
                                            placeholder="Search and select second gene...",
                                            style={"marginBottom": "10px"},
                                        ),
                                    ],
                                    id=f"{prefix}-gene2-container",
                                    style={"display": "none"},
                                ),
                                html.Div(
                                    [
                                        html.Label("Expression Thresholds:", style={"fontWeight": "bold", "marginBottom": "5px"}),
                                        html.Div(
                                            [
                                                html.Label("Gene 1 Threshold:", style={"fontSize": "12px"}),
                                                dcc.Slider(
                                                    id=f"{prefix}-gene1-threshold-slider",
                                                    min=0,
                                                    max=1,
                                                    value=0.5,
                                                    marks=None,
                                                    tooltip={"placement": "bottom", "always_visible": True},
                                                    className="dbc-slider",
                                                ),
                                            ],
                                            style={"marginBottom": "10px"},
                                        ),
                                        html.Div(
                                            [
                                                html.Label("Gene 2 Threshold:", style={"fontSize": "12px"}),
                                                dcc.Slider(
                                                    id=f"{prefix}-gene2-threshold-slider",
                                                    min=0,
                                                    max=1,
                                                    value=0.5,
                                                    marks=None,
                                                    tooltip={"placement": "bottom", "always_visible": True},
                                                    className="dbc-slider",
                                                ),
                                            ],
                                            style={"marginBottom": "10px"},
                                        ),
                                    ],
                                    id=f"{prefix}-threshold-container",
                                    style={"display": "none"},
                                ),
                                dcc.Loading(
                                    id=f"{prefix}-loading-gene-scatter",
                                    type="circle",
                                    children=dcc.Graph(id=f"{prefix}-gene-scatter", config=gene_scatter_config, style={"height": "60vh", "width": "100%"}),
                                    style={"height": "60vh", "width": "100%"},
                                ),
                            ],
                            className="dbc",
                            style={"marginBottom": "20px"},
                        ),
                        xs=12,
                        sm=12,
                        md=4,
                        lg=4,
                        xl=5,
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(width={"size": 0, "offset": 0, "md": 4, "lg": 4, "xl": 2}),
                    dbc.Col(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            dbc.Button(
                                                "Update other Plots",
                                                id=f"{prefix}-update-plots-button",
                                                color="primary",
                                                n_clicks=0,
                                                style={"marginRight": "10px", "display": "inline-block"},
                                            ),
                                            dbc.DropdownMenu(
                                                [dbc.DropdownMenuItem("Cell IDs (.txt)", id=f"{prefix}-download-cellids")],
                                                label="Download",
                                                color="secondary",
                                                id=f"{prefix}-download-menu",
                                                disabled=True,
                                                style={"display": "inline-block"},
                                            ),
                                            dcc.Download(id=f"{prefix}-download-cells-data"),
                                        ],
                                        style={"textAlign": "left"},
                                    ),
                                    html.Div(id=f"{prefix}-selection-status", style={"textAlign": "left", "marginTop": "5px"}),
                                ],
                                style={"marginTop": "10px"},
                            )
                        ],
                        xs=12,
                        sm=12,
                        md=4,
                        lg=4,
                        xl=5,
                    ),
                    dbc.Col(width={"size": 0, "md": 4, "lg": 4, "xl": 5}),
                ]
            ),
        ]
    )
