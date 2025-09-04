import json
from pathlib import Path
import warnings
import pandas as pd
import dash
from dash import dcc, html, Input, Output, exceptions, State, callback_context, ALL, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from guanaco.pages.matrix.cellplotly.embedding import plot_continuous_embedding, plot_coexpression_embedding,plot_categorical_embedding_with_fixed_colors,plot_continuous_annotation
from guanaco.pages.matrix.cellplotly.heatmap import plot_unified_heatmap
from guanaco.pages.matrix.cellplotly.violin1 import plot_violin1
from guanaco.pages.matrix.cellplotly.violin2 import plot_violin2_new
from guanaco.pages.matrix.cellplotly.stacked_bar import plot_stacked_bar
from guanaco.pages.matrix.cellplotly.dotmatrix import plot_dot_matrix
from guanaco.pages.matrix.cellplotly.pseudotime import plot_genes_in_pseudotime
from guanaco.pages.matrix.layout import generate_heatmap_layout,generate_violin_layout,generate_dotplot_layout,generate_stacked_bar_layout,generate_pseudotime_layout
# Import configs
from guanaco.config import scatter_config, gene_scatter_config
from guanaco.data_loader import color_config
warnings.filterwarnings('ignore', message='.*observed=False.*')

# helper functions
cvd_color_path = Path(__file__).parent / "cvd_color.json"
with open(cvd_color_path, "r") as f:
    palette_json = json.load(f)
palette_names = list(palette_json["color_palettes"].keys())

def apply_relayout(fig, relayout):
    if not relayout:
        return fig

    if all(k in relayout for k in ["xaxis.range[0]", "xaxis.range[1]", "yaxis.range[0]", "yaxis.range[1]"]):
        fig.update_layout(
            xaxis=dict(range=[relayout["xaxis.range[0]"], relayout["xaxis.range[1]"]]),
            yaxis=dict(range=[relayout["yaxis.range[0]"], relayout["yaxis.range[1]"]]),
        )
        return fig

    if "xaxis.range" in relayout and "yaxis.range" in relayout:
        fig.update_layout(
            xaxis=dict(range=relayout["xaxis.range"]),
            yaxis=dict(range=relayout["yaxis.range"]),
        )
        return fig

    # 3) Reset axes (autorange)
    if "xaxis.autorange" in relayout or "autosize" in relayout:
        fig.update_layout(
            xaxis=dict(autorange=True),
            yaxis=dict(autorange=True),
        )
        return fig

    return fig


def initialize_scatter_components(adata):
    embedding_prefixes = {
        "X_umap": "UMAP", "X_pca": "PCA", "X_tsne": "t-SNE", "X_diffmap": "DiffMap",
        "X_phate": "PHATE", "X_draw_graph_fa": "FA"
    }
    
    obsm_list = list(adata.obsm.keys())
    embedding_columns = {
        key: [f'{embedding_prefixes.get(key, key.upper())}{i+1}' for i in range(adata.obsm[key].shape[1])]
        for key in obsm_list
    }
    default_embedding = obsm_list[-1]
    default_columns = embedding_columns[default_embedding]
    
    return obsm_list, embedding_columns, default_embedding, default_columns


def create_control_components(adata, prefix):
    obsm_list, _, default_embedding, default_columns = initialize_scatter_components(adata)
    clustering_dropdown = dcc.Dropdown(
        id=f'{prefix}-clustering-dropdown',
        options=[{'label': key.replace("X_", "").upper(), 'value': key} for key in obsm_list],
        value='X_umap' if 'X_umap' in obsm_list else default_embedding,
        placeholder="Select Clustering Method",
        style={'marginBottom': '15px','fontSize': '14px'},
        clearable=False
    )
    
    coordinates_dropdowns = html.Div([
        html.Label("X-axis:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dcc.Dropdown(id=f'{prefix}-x-axis', options=[{'label': col, 'value': col} for col in default_columns], value=default_columns[0]),
        html.Label("Y-axis:", style={'fontWeight': 'bold', 'marginTop': '15px', 'marginBottom': '5px'}),
        dcc.Dropdown(id=f'{prefix}-y-axis', options=[{'label': col, 'value': col} for col in default_columns], value=default_columns[1])
    ], style={'marginBottom': '15px','fontSize': '14px'})
    
    return clustering_dropdown, coordinates_dropdowns


def generate_annotation_dropdown(anno_list, prefix, default_value=None):
    # Generate dropdown with both annotations and genes
    return dcc.Dropdown(
        id=f'{prefix}-annotation-dropdown', 
        options=[{'label': label, 'value': label} for label in anno_list],
        placeholder="Search annotations or genes...", 
        value=default_value if default_value else (anno_list[0] if anno_list else None),
        style={'marginBottom': '15px'}
    )


def generate_scatter_gene_selection(combined_list, prefix, default_value=None):
    return dcc.Dropdown(
        id=f'{prefix}-scatter-gene-selection', 
        options=[{'label': label, 'value': label} for label in combined_list], 
        value=default_value if default_value else (combined_list[0] if combined_list else None), 
        placeholder="Search annotations or genes...", 
        style={'marginBottom': '15px'}
    )

def create_global_metadata_filter(adata, prefix):
    """Create global metadata filter component"""
    
    # Get all categorical metadata columns
    categorical_columns = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'category' or adata.obs[col].dtype == 'object':
            unique_vals = adata.obs[col].unique()
            if len(unique_vals) < 100:  # Reasonable limit for categorical
                categorical_columns.append(col)
    
    # Create filter components
    filter_components = []
    
    for col in categorical_columns:
        unique_values = sorted([str(val) for val in adata.obs[col].unique()])
        
        filter_component = html.Div([
            html.Label(f"{col}:", style={'fontWeight': 'bold', 'marginRight': '10px', 'minWidth': '120px'}),
            dcc.Dropdown(
                id={'type': f'{prefix}-global-metadata-filter', 'column': col},
                options=[{'label': val, 'value': val} for val in unique_values],
                value=unique_values,  # All selected by default
                multi=True,
                placeholder=f"Select {col}...",
                style={'flex': '1', 'minWidth': '200px'}
            )
        ], style={
            'display': 'flex',
            'alignItems': 'center',
            'marginBottom': '8px',
            'padding': '5px'
        })
        
        filter_components.append(filter_component)
    
    # Global filter panel
    global_filter_panel = html.Div([
        # Header with cell count
        html.Div([
            html.H5("🔍 Global Data Filter", 
                   style={'margin': '0', 'color': '#2c3e50', 'display': 'inline-block'}),
            html.Div([
                html.Span("Active cells: ", style={'fontWeight': 'bold'}),
                html.Span(
                    id=f'{prefix}-global-cell-count',
                    children=f"{adata.n_obs:,}",
                    style={'color': '#27ae60', 'fontWeight': 'bold', 'fontSize': '16px'}
                ),
                html.Span(f" / {adata.n_obs:,} total", style={'color': '#7f8c8d'}),
                html.Span(
                    id=f'{prefix}-filter-preview',
                    children="",
                    style={'marginLeft': '15px', 'color': '#e67e22', 'fontSize': '14px', 'fontStyle': 'italic'}
                )
            ], style={'display': 'inline-block', 'marginLeft': '20px'})
        ], style={'marginBottom': '15px'}),
        
        # Collapsible filter section
        dbc.Collapse([
            html.Div(filter_components, style={'maxHeight': '300px', 'overflowY': 'auto'}),
            
            # Quick action buttons
            html.Div([
                dbc.Button("Select All", id=f'{prefix}-select-all-filters', 
                          color="success", size="sm", style={'marginRight': '10px'}),
                dbc.Button("Clear All", id=f'{prefix}-clear-all-filters', 
                          color="warning", size="sm", style={'marginRight': '10px'}),
                dbc.Button("Apply Filter", id=f'{prefix}-apply-global-filter', 
                          color="primary", size="sm", 
                          style={'fontWeight': 'bold', 'minWidth': '100px'})
            ], style={'textAlign': 'center', 'marginTop': '15px'})
            
        ], id=f'{prefix}-global-filter-collapse', is_open=False),
        
        # Toggle button
        html.Div([
            dbc.Button(
                "▼ Show Filters",
                id=f'{prefix}-toggle-global-filter',
                color="link",
                size="sm",
                style={'padding': '5px 10px', 'textDecoration': 'none'}
            )
        ], style={'textAlign': 'center', 'marginTop': '10px'}),
        
        # Hidden store for filtered data (start empty for better performance)
        dcc.Store(id=f'{prefix}-global-filtered-data', data={'cell_indices': None, 'n_cells': adata.n_obs})
        
    ], style={
        'backgroundColor': '#f8f9fa',
        'border': '2px solid #dee2e6',
        'borderRadius': '10px',
        'padding': '15px',
        'marginBottom': '20px',
        'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
    })
    
    return global_filter_panel



def generate_embedding_plots(adata,prefix):

    scatter_transformation_selection = html.Div([
        dbc.RadioItems(
            id=f'{prefix}-scatter-log-or-zscore',
            options=[
                {'label': 'None', 'value':None},
                {'label': 'Log', 'value': 'log'}
            ],
            value=None,
            inline=True,
            style={'fontSize': '14px'} 
        )
    ])

    scatter_order_selection = html.Div([
        dbc.RadioItems(
            id=f'{prefix}-plot-order',
            options=[
                {'label': 'Max-1st', 'value': 'max'},  # None option with 'False' value
                {'label': 'Min-1st', 'value': 'min'},
                {'label': 'Original', 'value': 'original'},
                {'label': 'Random', 'value': 'random'},
            ],
            value='max',
            inline=True,
            style={'fontSize': '14px'}

        )
    ])

    colorscales = px.colors.named_colorscales()
    # Dropdown for color map selection
    color_map_dropdown = dcc.Dropdown(
        id=f"{prefix}-scatter-color-map-dropdown",
        options=[{"label": scale, "value": scale} for scale in colorscales],
        value="viridis",  # Default colorscale
        style={"marginBottom": "10px"},
        clearable = False
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


    # Slider for marker size
    marker_size_slider = dcc.Slider(
        id=f'{prefix}-marker-size-slider',
        min=1,
        max=10,
        value=5,  # Default marker size
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )

    # Slider for opacity
    opacity_slider = dcc.Slider(
        id=f"{prefix}-opacity-slider",
        min=0.1,
        max=1.0,
        value=1,  # Default opacity
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )
    # toggle for scatter plot legend
    scatter_legend_toggle = dbc.RadioItems(
        id=f'{prefix}-scatter-legend-toggle',
        options=[
            {'label': 'on data', 'value': 'on data'},
            {'label': 'right', 'value': 'right'}
        ],
        value='right',
        inline=True,
        style={'fontSize': '14px'} 
    )
    # Toggle for axis show
    axis_toggle = dbc.RadioItems(
        id=f'{prefix}-axis-toggle',
        options=[
            {'label': 'Show', 'value': True},
            {'label': 'Hide', 'value': False}
        ],
        value=True,
        inline=True,
        style={'fontSize': '14px'}
    )

    graphic_control = html.Div(
        id=f"{prefix}-controls-container",
        children=[
            html.Div(
                [
                    html.Label("Plot order: ", className="control-label"),
                    scatter_order_selection,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Continous ColorMap:",  className="control-label"),
                    color_map_dropdown,
                ],
                className="dbc",
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Discrete ColorMap:",  className="control-label"),
                    color_map_discrete_dropdown,
                ],
                className="dbc",
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Marker Size:",  className="control-label"),
                    marker_size_slider,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Opacity:",  className="control-label"),
                    opacity_slider,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Legend Location:",  className="control-label"),
                    scatter_legend_toggle,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Axis:",  className="control-label"),
                    axis_toggle,
                ],
                style={'marginBottom': '15px'}
            ),
        ],
        style={'display': 'none'}  # hidden by default
    )
    # Get annotations and genes
    annotations = adata.obs.columns.tolist()
    sample_genes = adata.var_names[:20].tolist()
    
    # Create combined list for both dropdowns
    combined_list = annotations + sample_genes
    
    # Set default values: first annotation for left, first gene for right
    default_annotation = annotations[0] if annotations else None
    default_gene = sample_genes[0] if sample_genes else None

    clustering_dropdown, coordinates_dropdowns = create_control_components(adata, prefix)
    layout = html.Div([
        # Add Store component to hold selected cells data
        dcc.Store(id=f'{prefix}-selected-cells-store'),
        
        create_global_metadata_filter(adata, prefix),
        
        dbc.Row([
            # Left column: Controls
            dbc.Col(
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Dimension Reduction:", className = "control-label"),
                            clustering_dropdown,
                        ],
                        className="dbc",
                        style={'marginBottom': '15px'}
                    ),
                    html.Div(
                        [
                            html.Div(coordinates_dropdowns, id=f'{prefix}-coordinates-dropdowns')
                        ],
                        className="dbc",
                        style={'marginBottom': '15px'}
                    ),
                    html.Div(
                        [
                            html.Label("Transformation:",  className = "control-label"),
                            scatter_transformation_selection,
                        ],
                        style={'marginBottom': '15px'}
                    ),
                    html.Button("More controls", id=f"{prefix}-toggle-button", n_clicks=0, style={'marginBottom': '10px','border':'1px solid','borderRadius': '5px'}),

                    graphic_control,
                ],
            ),
            xs=12, sm=12, md=4, lg=4, xl=2, 
            style={"borderRight": "1px solid #ddd", "padding": "10px"}
        ),

    dbc.Col(
        html.Div(
            [
                html.Label("Select Annotation/check:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                generate_annotation_dropdown(anno_list=combined_list, prefix=prefix, default_value=default_annotation),
                # Add invisible spacer to align with gene scatter's radio buttons
                html.Div(style={'height': '25px', 'marginBottom': '5px'}),  # Same height as RadioItems
                dcc.Loading(
                    id=f"{prefix}-loading-annotaion-scatter",
                    type="circle",
                    children=dcc.Graph(
                        id=f'{prefix}-annotation-scatter', 
                        config=scatter_config,
                        style={'height': '60vh', 'width': '100%'}
                    ),
                    style={'height': '60vh', 'width': '100%'}
                ),
            ],
            className="dbc",
            style={'marginBottom': '20px'}
        ),
        xs=12, sm=12, md=4, lg=4, xl=5 # Full width on small screens, equal split of remaining space on larger screens
    ),
    dbc.Col(
        html.Div(
            [
                html.Label("Select Annotation/Gene:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                generate_scatter_gene_selection(combined_list=combined_list, prefix=prefix, default_value=default_gene),
                # Add toggle for co-expression mode
                dbc.RadioItems(
                    id=f'{prefix}-coexpression-toggle',
                    options=[
                        {'label': 'Single Gene', 'value': 'single'},
                        {'label': 'Co-expression', 'value': 'coexpression'}
                    ],
                    value='single',
                    inline=True,
                    style={'marginBottom': '10px', 'fontSize': '14px'}
                ),
                # Second gene selection (hidden by default)
                html.Div(
                    [
                        html.Label("Second Gene:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id=f'{prefix}-scatter-gene2-selection',
                            options=[{'label': label, 'value': label} for label in adata.var_names.to_list()[:10]],
                            value=adata.var_names.to_list()[1],
                            placeholder="Search and select second gene...",
                            style={'marginBottom': '10px'}
                        ),
                    ],
                    id=f'{prefix}-gene2-container',
                    style={'display': 'none'}
                ),
                # Threshold sliders 
                html.Div(
                    [
                        html.Label("Expression Thresholds:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                        html.Div([
                            html.Label("Gene 1 Threshold:", style={'fontSize': '12px'}),
                            dcc.Slider(
                                id=f'{prefix}-gene1-threshold-slider',
                                min=0,
                                max=1,
                                value=0.5,
                                marks=None,
                                tooltip={"placement": "bottom", "always_visible": True},
                                className="dbc-slider"
                            ),
                        ], style={'marginBottom': '10px'}),
                        html.Div([
                            html.Label("Gene 2 Threshold:", style={'fontSize': '12px'}),
                            dcc.Slider(
                                id=f'{prefix}-gene2-threshold-slider',
                                min=0,
                                max=1,
                                value=0.5,
                                marks=None,
                                tooltip={"placement": "bottom", "always_visible": True},
                                className="dbc-slider"
                            ),
                        ], style={'marginBottom': '10px'}),
                    ],
                    id=f'{prefix}-threshold-container',
                    style={'display': 'none'}
                ),
                dcc.Loading(
                    id=f"{prefix}-loading-gene-scatter",
                    type="circle",
                    children=dcc.Graph(
                        id=f'{prefix}-gene-scatter', 
                        config=gene_scatter_config,
                        style={'height': '60vh', 'width': '100%'}
                    ),
                    style={'height': '60vh', 'width': '100%'}
                ),
            ],
            className="dbc",
            style={'marginBottom': '20px'}
        ),
        xs=12, sm=12, md=4, lg=4, xl=5  # Full width on small screens, equal split of remaining space on larger screens
    ),
        ]),
        
        # Add control buttons as a separate row below the scatter plots
        dbc.Row([
            dbc.Col(width={"size": 0, "offset": 0, "md": 4, "lg": 4, "xl": 2}),  # Empty column to align with control column
            dbc.Col([
                html.Div([
                    html.Div([
                        dbc.Button(
                            "Update other Plots",
                            id=f"{prefix}-update-plots-button",
                            color="primary",
                            n_clicks=0,
                            style={'marginRight': '10px', 'display': 'inline-block'}
                        ),
                        dbc.DropdownMenu(
                            [
                                dbc.DropdownMenuItem("Cell IDs (.txt)", id=f"{prefix}-download-cellids"),
                                dbc.DropdownMenuItem("Subset AnnData (.h5ad)", id=f"{prefix}-download-adata"),
                            ],
                            label="Download",
                            color="secondary",
                            id=f"{prefix}-download-menu",
                            disabled=True,
                            style={'display': 'inline-block'}
                        ),
                        dcc.Download(id=f"{prefix}-download-cells-data")
                    ], style={'textAlign': 'left'}),
                    html.Div(id=f"{prefix}-selection-status", style={'textAlign': 'left', 'marginTop': '5px'})
                ], style={'marginTop': '10px'})
            ], xs=12, sm=12, md=4, lg=4, xl=5),  # Align with annotation scatter column
            dbc.Col(width={"size": 0, "md": 4, "lg": 4, "xl": 5})  # Empty column to balance layout
        ])
    ])

    return layout


def filter_data(adata, annotation, selected_labels, selected_cells=None):
    if selected_cells is not None and len(selected_cells) > 0:
        adata_filtered = adata[selected_cells]
    elif selected_labels and annotation:
        # Apply label filtering
        cell_indices = adata.obs[annotation].isin(selected_labels)
        adata_filtered = adata[cell_indices]
    else:
        # No filtering
        adata_filtered = adata
    
    return adata_filtered


def generate_left_control(default_gene_markers, label_list, prefix):
    # Radio buttons for input mode selection
    input_mode_radio = dbc.RadioItems(
        id=f'{prefix}-gene-input-mode',
        options=[
            {'label': 'Dropdown', 'value': 'dropdown'},
            {'label': 'Text', 'value': 'text'}
        ],
        value='dropdown',
        inline=True,
        style={'fontSize': '14px', 'marginBottom': '10px'}
    )
    
    # Dropdown for gene selection
    genes_dropdown = dcc.Dropdown(
        id=f'{prefix}-single-cell-genes-selection',
        options=[{'label': gene, 'value': gene} for gene in default_gene_markers],
        value=default_gene_markers,
        multi=True,
        style={'marginBottom': '15px', 'font-size': '12px'},
        className='custom-dropdown'
    )
    
    # Text area for gene input (initially hidden)
    genes_textarea = dcc.Textarea(
        id=f'{prefix}-single-cell-genes-textarea',
        placeholder='Enter genes separated by commas (e.g., Gene1, Gene2, Gene3)',
        value=', '.join(default_gene_markers),
        style={'width': '100%', 'height': '80px', 'marginBottom': '10px', 'display': 'none'},
        className='custom-textarea'
    )
    
    # Error message div
    error_message = html.Div(
        id=f'{prefix}-gene-input-error',
        style={'color': 'red', 'fontSize': '12px', 'marginBottom': '10px'}
    )
    
    # Container for gene selection
    genes_selection = html.Div([
        input_mode_radio,
        genes_dropdown,
        genes_textarea,
        error_message
    ])
    
    annotation_filter = dcc.Dropdown(
        id=f'{prefix}-single-cell-annotation-dropdown',
        options=[{'label': label, 'value': label} for label in label_list],
        value=label_list[len(label_list) // 2],
        style={'marginBottom': '15px'},
        clearable=False,
        className='custom-dropdown'
    )
    
    label_list_selection = dcc.Dropdown(
        id=f'{prefix}-single-cell-label-selection',
        multi=True,
        style={'marginBottom': '15px', 'font-size': '12px'},
        className='custom-dropdown'
    )

    
    return html.Div([
        html.Label('Select Variables:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        genes_selection,
        html.Label('Select Annotation:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        annotation_filter,
        html.Label('Select Labels:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        label_list_selection,
    ])


def generate_other_plots(adata, default_gene_markers, discrete_label_list, prefix):

    tabs = dcc.Tabs([
        dcc.Tab(label='Heatmap', value='heatmap-tab', children=[
            html.Div(generate_heatmap_layout(adata, prefix))
        ]),
        dcc.Tab(label='Violin Plot', value='violin-tab', children=[
            html.Div(generate_violin_layout(default_gene_markers, discrete_label_list, prefix))
        ]),
        dcc.Tab(label='Dotplot', value='dotplot-tab', children=[
            html.Div(generate_dotplot_layout(prefix))
        ]),
        dcc.Tab(label='Stacked Bar', value='stacked-bar-tab', children=[
            html.Div(generate_stacked_bar_layout(discrete_label_list,prefix))
        ]),
        dcc.Tab(label='Pseudotime Plot', value='pseudotime-tab', children=[
            html.Div(generate_pseudotime_layout(prefix))
        ]),
    ], id=f'{prefix}-single-cell-tabs', value='dotplot-tab', className='custom-tabs')

    return dbc.Row([
        dbc.Col(
            generate_left_control(default_gene_markers, discrete_label_list, prefix),
            xs=12, sm=12, md=4, lg=4, xl=2
        ),
        dbc.Col(tabs, xs=12, sm=12, md=8, lg=8, xl=10)
    ], style={'marginBottom': '50px'})

# ============= Helper Functions =============

def is_continuous_annotation(adata, annotation, threshold=50):
    """Check if an annotation is continuous based on unique value count and data type."""
    if annotation not in adata.obs.columns:
        return False
    
    # Check data type
    dtype = adata.obs[annotation].dtype
    if dtype in ['float32', 'float64', 'int32', 'int64']:
        # Numeric type - check unique values
        n_unique = adata.obs[annotation].nunique()
        return n_unique >= threshold
    return False

# ============= Main Callback Functions =============

def matrix_callbacks(app, adata, prefix):
    """Combined callback registration for both scatter and other plots"""
    
    # ===== Global Filter Callbacks =====
    
    @app.callback(
        [Output(f'{prefix}-global-filter-collapse', 'is_open'),
         Output(f'{prefix}-toggle-global-filter', 'children')],
        Input(f'{prefix}-toggle-global-filter', 'n_clicks'),
        State(f'{prefix}-global-filter-collapse', 'is_open'),
        prevent_initial_call=True
    )
    def toggle_global_filter(n_clicks, is_open):
        if is_open:
            return False, "▼ Show Filters"
        else:
            return True, "▲ Hide Filters"
    
    @app.callback(
        [Output({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'value'),
         Output(f'{prefix}-filter-preview', 'children')],
        [Input(f'{prefix}-select-all-filters', 'n_clicks'),
         Input(f'{prefix}-clear-all-filters', 'n_clicks'),
         Input({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'value')],
        [State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'options'),
         State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'id')],
        prevent_initial_call=True
    )
    def update_all_filters_and_preview(select_clicks, clear_clicks, current_values, all_options, all_ids):
        ctx = callback_context
        if not ctx.triggered:
            raise PreventUpdate
        
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        values = []
        
        if f'{prefix}-select-all-filters' in trigger_id:
            # Select all values for each filter
            for options in all_options:
                if options:
                    values.append([opt['value'] for opt in options])
                else:
                    values.append([])
        elif f'{prefix}-clear-all-filters' in trigger_id:
            # Clear all filters
            values = [[] for _ in all_options]
        else:
            # Use current values for real-time preview
            values = current_values or []
        
        # Calculate preview cell count in real-time
        if values and all_ids:
            mask = pd.Series(True, index=adata.obs.index)
            
            # Apply each filter
            for i, (filter_values, filter_id) in enumerate(zip(values, all_ids)):
                if filter_values:  # Only apply if values are selected
                    column = filter_id['column']
                    col_mask = adata.obs[column].astype(str).isin(filter_values)
                    mask = mask & col_mask
            
            # Get preview count
            preview_count = mask.sum()
            percentage = (preview_count / adata.n_obs) * 100
            preview_text = f"Preview: {preview_count:,} cells will be selected"
        else:
            preview_text = ""
        
        return values, preview_text
    
    @app.callback(
        [Output(f'{prefix}-global-filtered-data', 'data'),
         Output(f'{prefix}-global-cell-count', 'children'),
         Output(f'{prefix}-filter-preview', 'children', allow_duplicate=True)],
        Input(f'{prefix}-apply-global-filter', 'n_clicks'),
        [State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'value'),
         State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'id')],
        prevent_initial_call=True
    )
    def apply_global_filter(n_clicks, filter_values, filter_ids):
        if not n_clicks:
            raise PreventUpdate
        
        # Start with all cells
        mask = pd.Series(True, index=adata.obs.index)
        
        # Apply each filter
        for i, (values, filter_id) in enumerate(zip(filter_values, filter_ids)):
            if values:  # Only apply if values are selected
                column = filter_id['column']
                col_mask = adata.obs[column].astype(str).isin(values)
                mask = mask & col_mask
        
        # Get filtered cell indices
        filtered_indices = adata.obs.index[mask].tolist()
        n_filtered = len(filtered_indices)
        
        # Simple status message without "Filtered by:" details
        preview_text = f"Applied: {n_filtered:,} cells selected"
        
        # Update cell count display
        cell_count_text = f"{n_filtered:,}"
        
        return {
            'cell_indices': filtered_indices,
            'n_cells': n_filtered
        }, cell_count_text, preview_text
    
    # ===== Gene Selection Callbacks =====
    
    @app.callback(
        [Output(f'{prefix}-single-cell-genes-selection', 'style'),
         Output(f'{prefix}-single-cell-genes-textarea', 'style')],
        Input(f'{prefix}-gene-input-mode', 'value')
    )
    def toggle_gene_input_display(input_mode):
        if input_mode == 'dropdown':
            return {'marginBottom': '15px', 'font-size': '12px'}, {'display': 'none'}
        else:
            return {'display': 'none'}, {'width': '100%', 'height': '80px', 'marginBottom': '10px'}
    
    @app.callback(
        [Output(f'{prefix}-single-cell-genes-selection', 'value', allow_duplicate=True),
         Output(f'{prefix}-gene-input-error', 'children')],
        [Input(f'{prefix}-single-cell-genes-textarea', 'value'),
         Input(f'{prefix}-gene-input-mode', 'value')],
        prevent_initial_call=True
    )
    def validate_gene_input(textarea_value, input_mode):
        if input_mode != 'text' or not textarea_value:
            return no_update, ''
        
        # Parse the textarea input - handle various formats
        # First, check if it looks like a Python list (has brackets)
        text = textarea_value.strip()
        if text.startswith('[') and text.endswith(']'):
            text = text[1:-1]  # Remove brackets
        elif text.startswith('(') and text.endswith(')'):
            text = text[1:-1]  # Remove parentheses
        
        # Split by comma and clean each gene name
        input_genes = []
        for item in text.split(','):
            # Remove various quote types and whitespace
            gene = item.strip()
            # Remove quotes (single, double, backticks, smart quotes)
            gene = gene.strip('"\'`''""')
            if gene:  # Only add non-empty strings
                input_genes.append(gene)
        
        # Get all available genes (case-insensitive mapping)
        available_genes = list(adata.var_names)
        gene_map = {gene.upper(): gene for gene in available_genes}
        
        valid_genes = []
        invalid_genes = []
        seen_genes = set()
        duplicate_genes = []
        
        for gene in input_genes:
            gene_upper = gene.upper()
            if gene_upper in gene_map:
                actual_gene = gene_map[gene_upper]
                if actual_gene not in seen_genes:
                    valid_genes.append(actual_gene)
                    seen_genes.add(actual_gene)
                else:
                    duplicate_genes.append(gene)
            else:
                invalid_genes.append(gene)
        
        error_messages = []
        if invalid_genes:
            error_messages.append(f"Invalid genes not found: {', '.join(invalid_genes)}")
        if duplicate_genes:
            error_messages.append(f"Duplicate genes removed: {', '.join(duplicate_genes)}")
        
        error_message = ' | '.join(error_messages)
        
        return valid_genes, error_message
    
    @app.callback(
        Output(f'{prefix}-single-cell-genes-textarea', 'value'),
        Input(f'{prefix}-single-cell-genes-selection', 'value'),
        State(f'{prefix}-gene-input-mode', 'value'),
        prevent_initial_call=True
    )
    def sync_dropdown_to_textarea(dropdown_value, input_mode):
        if input_mode == 'dropdown' and dropdown_value:
            return ', '.join(dropdown_value)
        return no_update
    
    # ===== Scatter Plot Callbacks =====
    
    @app.callback(
        Output(f"{prefix}-controls-container", "style"),
        Output(f"{prefix}-toggle-button", "children"),
        Input(f"{prefix}-toggle-button", "n_clicks"),
        prevent_initial_call=True
    )
    def toggle_controls(n_clicks):
        if n_clicks % 2 == 1:
            return {'display': 'block'}, "Hide controls"
        else:
            return {'display': 'none'}, "More controls"
    
    @app.callback(
    [Output(f'{prefix}-coordinates-dropdowns', 'children'),
    Output(f'{prefix}-x-axis', 'value'),
    Output(f'{prefix}-y-axis', 'value')],
    Input(f'{prefix}-clustering-dropdown', 'value')
    )
    def update_coordinates_dropdowns(selected_clustering):
        _, embedding_columns, _, _ = initialize_scatter_components(adata)
        selected_columns = embedding_columns[selected_clustering]
        options = [{'label': col, 'value': col} for col in selected_columns]
        x_value = selected_columns[0]
        y_value = selected_columns[1] if len(selected_columns) > 1 else selected_columns[0]
        
        return (
            html.Div([
                html.Div([
                    html.Label("X-axis:"),
                    dcc.Dropdown(id=f'{prefix}-x-axis', options=options, value=x_value, clearable=False,style={'fontSize': '14px'})
                ], style={'flex': '1', 'paddingRight': '10px'}),
                
                html.Div([
                    html.Label("Y-axis:"),
                    dcc.Dropdown(id=f'{prefix}-y-axis', options=options, value=y_value,clearable=False,style={'fontSize': '14px'})
                ], style={'flex': '1', 'paddingLeft': '10px'})
            ], style={'display': 'flex', 'marginBottom': '15px'}),
            x_value,
            y_value
        )
    
    @app.callback(
        Output(f'{prefix}-annotation-dropdown', 'options'),
        Input(f'{prefix}-annotation-dropdown', 'search_value')
    )
    def update_annotation_dropdown(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        
        # Include both metadata and genes
        label_list = adata.obs.keys().to_list()
        gene_list = adata.var_names.to_list()
        
        # Search in metadata first, then genes
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        matching_genes = [gene for gene in gene_list if search_value.lower() in gene.lower()]
        
        # Combine results with metadata first, then genes (limited to 10 total)
        all_matches = matching_labels + matching_genes
        return [{'label': item, 'value': item} for item in all_matches[:10]]
    
    @app.callback(
        Output(f'{prefix}-scatter-gene-selection', 'options'),
        Input(f'{prefix}-scatter-gene-selection', 'search_value')
    )
    def update_scatter_gene_selection(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
                # Include both metadata and genes
        label_list = adata.obs.keys().to_list()
        gene_list = adata.var_names.to_list()
        
        # Search in metadata first, then genes
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        matching_genes = [gene for gene in gene_list if search_value.lower() in gene.lower()]
        
        # Combine results with metadata first, then genes (limited to 10 total)
        all_matches = matching_labels + matching_genes
        return [{'label': item, 'value': item} for item in all_matches[:10]]
        # gene_list = adata.var_names.to_list()
        # matching_genes = [gene for gene in gene_list if search_value.lower() in gene.lower()]
        # return [{'label': gene, 'value': gene} for gene in matching_genes[:20]]
    
    @app.callback(
        Output(f'{prefix}-scatter-gene2-selection', 'options'),
        Input(f'{prefix}-scatter-gene2-selection', 'search_value')
    )
    def update_scatter_gene2_selection(search_value):
        if not search_value:
            raise exceptions.PreventUpdate
        gene_list = adata.var_names.to_list()
        matching_genes = [gene for gene in gene_list if search_value.lower() in gene.lower()]
        return [{'label': gene, 'value': gene} for gene in matching_genes[:20]]
    
    @app.callback(
        [Output(f'{prefix}-gene2-container', 'style'),
         Output(f'{prefix}-threshold-container', 'style')],
        Input(f'{prefix}-coexpression-toggle', 'value')
    )
    def toggle_coexpression_controls(mode):
        if mode == 'coexpression':
            return {'display': 'block'}, {'display': 'block'}
        else:
            return {'display': 'none'}, {'display': 'none'}
    
    @app.callback(
        Output(f'{prefix}-annotation-scatter', 'figure'),
        [Input(f'{prefix}-clustering-dropdown', 'value'),
         Input(f'{prefix}-x-axis', 'value'),
         Input(f'{prefix}-y-axis', 'value'),
         Input(f'{prefix}-annotation-dropdown', 'value'),
         Input(f'{prefix}-marker-size-slider', 'value'),
         Input(f'{prefix}-opacity-slider', 'value'),
         Input(f'{prefix}-scatter-legend-toggle', 'value'),
         Input(f'{prefix}-axis-toggle', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-scatter-log-or-zscore', 'value'),  # Add for continuous transformations
         Input(f'{prefix}-plot-order', 'value'),  # Add for continuous ordering
         Input(f'{prefix}-scatter-color-map-dropdown', 'value'),  # Add for continuous color maps
         Input(f'{prefix}-global-filtered-data', 'data'),  # Add global filtered data
         Input(f'{prefix}-gene-scatter', 'relayoutData'),  # Add gene scatter relayout for bidirectional sync
         ]
    )
    def update_annotation_scatter(clustering_method, x_axis, y_axis, annotation, 
                                marker_size, opacity, legend_show, axis_show, 
                                discrete_color_map, transformation, order, continuous_color_map, 
                                filtered_data, gene_relayout):
        if not annotation:
            raise exceptions.PreventUpdate
        
        # Only use filtered data if filter has actually been applied (not default state)
        if (filtered_data and 
            filtered_data.get('cell_indices') is not None and 
            filtered_data.get('n_cells', adata.n_obs) < adata.n_obs):
            plot_adata = adata[filtered_data['cell_indices']]
        else:
            plot_adata = adata
        
        # Check if annotation is a gene or metadata
        if annotation in adata.var_names:
            # This is a gene - use continuous gene plotting
            color_map = continuous_color_map or 'Viridis'
            
            fig = plot_continuous_embedding(
                adata=plot_adata,
                embedding_key=clustering_method,
                color=annotation,
                x_axis=x_axis,
                y_axis=y_axis,
                transformation=transformation,
                order=order,
                color_map=color_map,
                marker_size=marker_size,
                opacity=opacity,
                axis_show=axis_show,
            )
        elif is_continuous_annotation(adata, annotation):
            # Use continuous annotation plotting
            color_map = continuous_color_map or 'Viridis'
            
            fig = plot_continuous_annotation(
                adata=plot_adata,
                embedding_key=clustering_method,
                annotation=annotation,
                x_axis=x_axis,
                y_axis=y_axis,
                transformation=None,  # Disable transformation for annotation data
                order=order,
                color_map=color_map,
                marker_size=marker_size,
                opacity=opacity,
                axis_show=axis_show,
            )
        else:
            # Use categorical plotting
            if discrete_color_map is None:
                color_map = color_config
            else:
                color_map = palette_json["color_palettes"][discrete_color_map]
            
            fig = plot_categorical_embedding_with_fixed_colors(
                adata=plot_adata,
                adata_full=adata,  
                gene=None,  
                embedding_key=clustering_method,
                color=annotation,
                x_axis=x_axis,
                y_axis=y_axis,
                color_map=color_map,
                marker_size=marker_size,
                opacity=opacity,
                legend_show=legend_show,
                axis_show=axis_show,
            )
        
        fig = apply_relayout(fig, gene_relayout)
        return fig
    
    @app.callback(
        Output(f'{prefix}-gene-scatter', 'figure'),
        [Input(f'{prefix}-scatter-gene-selection', 'value'),
         Input(f'{prefix}-clustering-dropdown', 'value'),
         Input(f'{prefix}-x-axis', 'value'),
         Input(f'{prefix}-y-axis', 'value'),
         Input(f'{prefix}-scatter-log-or-zscore', 'value'),
         Input(f'{prefix}-plot-order', 'value'),
         Input(f'{prefix}-scatter-color-map-dropdown', 'value'),
         Input(f'{prefix}-marker-size-slider', 'value'),
         Input(f'{prefix}-opacity-slider', 'value'),
         Input(f'{prefix}-annotation-scatter', 'relayoutData'),
         Input(f'{prefix}-axis-toggle', 'value'),
         Input(f'{prefix}-coexpression-toggle', 'value'),
         Input(f'{prefix}-scatter-gene2-selection', 'value'),
         Input(f'{prefix}-gene1-threshold-slider', 'value'),
         Input(f'{prefix}-gene2-threshold-slider', 'value'),
         Input(f'{prefix}-scatter-legend-toggle', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),  # Add discrete color map
         Input(f'{prefix}-global-filtered-data', 'data'),  # Add global filtered data
         ]
    )
    def update_gene_scatter(gene_name, clustering, x_axis, y_axis, transformation, order, 
                           color_map, marker_size, opacity, annotation_relayout, axis_show,
                           coexpression_mode, gene2_name, threshold1, threshold2, legend_show, 
                           discrete_color_map, filtered_data):
        if not gene_name:
            raise exceptions.PreventUpdate
        
        # Only use filtered data if filter has actually been applied (not default state)
        if (filtered_data and 
            filtered_data.get('cell_indices') is not None and 
            filtered_data.get('n_cells', adata.n_obs) < adata.n_obs):
            plot_adata = adata[filtered_data['cell_indices']]
        else:
            plot_adata = adata
        
        # Check if gene_name is actually a gene or an annotation
        if gene_name in adata.var_names:
            # This is a gene
            if coexpression_mode == 'coexpression' and gene2_name:
                # Use co-expression visualization
                fig = plot_coexpression_embedding(
                    adata=plot_adata,
                    embedding_key=clustering,
                    gene1=gene_name,
                    gene2=gene2_name,
                    x_axis=x_axis,
                    y_axis=y_axis,
                    threshold1=threshold1,
                    threshold2=threshold2,
                    transformation=transformation,
                    color_map=None,  # Use default colors for co-expression
                    marker_size=marker_size,
                    opacity=opacity,
                    legend_show=legend_show,
                    axis_show=axis_show,
                )
            else:
                # Use single gene visualization
                fig = plot_continuous_embedding(
                    adata=plot_adata,
                    embedding_key=clustering,
                    color=gene_name,
                    x_axis=x_axis,
                    y_axis=y_axis,
                    transformation=transformation,
                    order=order,
                    color_map=color_map or 'Viridis',
                    marker_size=marker_size,
                    opacity=opacity,
                    annotation=None,
                    axis_show=axis_show,
                )
        elif is_continuous_annotation(adata, gene_name):
            # This is a continuous annotation
            fig = plot_continuous_annotation(
                adata=plot_adata,
                embedding_key=clustering,
                annotation=gene_name,
                x_axis=x_axis,
                y_axis=y_axis,
                transformation=None,  # Disable transformation for annotation data
                order=order,
                color_map=color_map or 'Viridis',
                marker_size=marker_size,
                opacity=opacity,
                axis_show=axis_show,
            )
        else:
            # This is a categorical annotation
            if discrete_color_map is None:
                discrete_color_map_value = color_config
            else:
                discrete_color_map_value = palette_json["color_palettes"][discrete_color_map]
            
            fig = plot_categorical_embedding_with_fixed_colors(
                adata=plot_adata,
                adata_full=adata,  # Pass full adata for color reference
                gene=None,  # Don't pass gene to avoid unnecessary computation
                embedding_key=clustering,
                color=gene_name,
                x_axis=x_axis,
                y_axis=y_axis,
                color_map=discrete_color_map_value,
                marker_size=marker_size,
                opacity=opacity,
                legend_show=legend_show,
                axis_show=axis_show,
            )
        
        

        fig = apply_relayout(fig, annotation_relayout)

        return fig
    
    # ===== Threshold Slider Update Callback =====
    @app.callback(
        [Output(f'{prefix}-gene1-threshold-slider', 'min'),
         Output(f'{prefix}-gene1-threshold-slider', 'max'),
         Output(f'{prefix}-gene1-threshold-slider', 'value'),
         Output(f'{prefix}-gene2-threshold-slider', 'min'),
         Output(f'{prefix}-gene2-threshold-slider', 'max'),
         Output(f'{prefix}-gene2-threshold-slider', 'value')],
        [Input(f'{prefix}-scatter-gene-selection', 'value'),
         Input(f'{prefix}-scatter-gene2-selection', 'value'),
         Input(f'{prefix}-scatter-log-or-zscore', 'value'),
         Input(f'{prefix}-global-filtered-data', 'data')],
    )
    def update_threshold_ranges(gene1, gene2, transformation, filtered_data):
        """Update threshold slider ranges based on gene expression min/max values"""
        # Get the appropriate dataset
        if (filtered_data and 
            filtered_data.get('cell_indices') is not None and 
            filtered_data.get('n_cells', adata.n_obs) < adata.n_obs):
            plot_adata = adata[filtered_data['cell_indices']]
        else:
            plot_adata = adata
        
        # Default values for sliders (0-1 range)
        default_min, default_max, default_value = 0, 1, 0.5
        
        # Calculate ranges for gene1
        if gene1 and gene1 in plot_adata.var_names:
            from guanaco.pages.matrix.cellplotly.gene_extraction_utils import extract_gene_expression, apply_transformation
            gene1_expr = extract_gene_expression(plot_adata, gene1)
            
            # Apply transformation if specified
            if transformation:
                gene1_expr = apply_transformation(gene1_expr, transformation, copy=True)
            
            # Check if gene is expressed
            if gene1_expr.max() > gene1_expr.min():
                gene1_min = float(gene1_expr.min())
                gene1_max = float(gene1_expr.max())
                gene1_value = (gene1_min + gene1_max) / 2
            else:
                # Gene not expressed, use default range
                gene1_min, gene1_max, gene1_value = default_min, default_max, default_value
        else:
            gene1_min, gene1_max, gene1_value = default_min, default_max, default_value
        
        # Calculate ranges for gene2
        if gene2 and gene2 in plot_adata.var_names:
            from guanaco.pages.matrix.cellplotly.gene_extraction_utils import extract_gene_expression, apply_transformation
            gene2_expr = extract_gene_expression(plot_adata, gene2)
            
            # Apply transformation if specified
            if transformation:
                gene2_expr = apply_transformation(gene2_expr, transformation, copy=True)
            
            # Check if gene is expressed
            if gene2_expr.max() > gene2_expr.min():
                gene2_min = float(gene2_expr.min())
                gene2_max = float(gene2_expr.max())
                gene2_value = (gene2_min + gene2_max) / 2
            else:
                # Gene not expressed, use default range
                gene2_min, gene2_max, gene2_value = default_min, default_max, default_value
        else:
            gene2_min, gene2_max, gene2_value = default_min, default_max, default_value
        
        return gene1_min, gene1_max, gene1_value, gene2_min, gene2_max, gene2_value
    
    # ===== Cell Selection Callback =====
    @app.callback(
        Output(f'{prefix}-selected-cells-store', 'data'),
        Output(f'{prefix}-selection-status', 'children'),
        [Input(f'{prefix}-update-plots-button', 'n_clicks')],
        [State(f'{prefix}-annotation-scatter', 'selectedData'),
         State(f'{prefix}-annotation-dropdown', 'value'),
         State(f'{prefix}-global-filtered-data', 'data')], 
        prevent_initial_call=True
    )
    def store_selected_cells(n_clicks, selected_data, current_annotation, filtered_data):
        """Store indices of selected cells from annotation scatter plot when button is clicked"""
        if n_clicks == 0:
            return None, ""
        
        if (filtered_data and 
            filtered_data.get('cell_indices') is not None and 
            filtered_data.get('n_cells', adata.n_obs) < adata.n_obs):
            plot_adata = adata[filtered_data['cell_indices']]
        else:
            plot_adata = adata
        
        if not selected_data or not selected_data.get('points'):
            all_indices = plot_adata.obs.index.tolist()
            n_cells = len(all_indices)
            status_msg = dbc.Alert(
                f"✓ All {n_cells} cells from scatter plot selected. Other plots updated.",
                color="info",
                dismissable=True,
                duration=4000
            )
            return all_indices, status_msg
        
        selected_points = selected_data['points']
        selected_indices = []
        
        if current_annotation in adata.var_names or is_continuous_annotation(plot_adata, current_annotation):
            for point in selected_points:
                if 'customdata' in point:
                    customdata = point['customdata']
                    if isinstance(customdata, (list, tuple)) and len(customdata) > 1:
                        cell_idx = int(customdata[1])
                    else:
                        # Single value customdata contains the cell index directly
                        cell_idx = int(customdata)
                    selected_indices.append(plot_adata.obs.index[cell_idx])
                else:
                    # Fallback to point number if customdata is not available
                    point_number = point.get('pointNumber', 0)
                    selected_indices.append(plot_adata.obs.index[point_number])
        else:
            for point in selected_points:
                curve_number = point.get('curveNumber', 0)
                point_number = point.get('pointNumber', 0)
                if curve_number == 0:
                    continue
                
                if 'customdata' in point:
                    # Get category from customdata
                    customdata = point['customdata']
                    if isinstance(customdata, (list, tuple)):
                        category = customdata[0]  # Category is first element
                    else:
                        category = customdata
                    
                    # Find all cells with this category
                    category_mask = plot_adata.obs[current_annotation] == category
                    category_indices = plot_adata.obs.index[category_mask].tolist()
                    
                    # The point_number is the index within this category
                    if point_number < len(category_indices):
                        selected_indices.append(category_indices[point_number])
                else:
                    # Fallback to original logic
                    unique_categories = sorted(plot_adata.obs[current_annotation].unique())
                    category_index = curve_number - 1
                    
                    if category_index >= 0 and category_index < len(unique_categories):
                        selected_category = unique_categories[category_index]
                        category_mask = plot_adata.obs[current_annotation] == selected_category
                        category_indices = plot_adata.obs.index[category_mask].tolist()
                        
                        if point_number < len(category_indices):
                            selected_indices.append(category_indices[point_number])
        
        if selected_indices:
            n_selected = len(selected_indices)
            status_msg = dbc.Alert(
                f"✓ {n_selected} cells selected from {current_annotation}. Other plots updated.",
                color="success",
                dismissable=True,
                duration=4000
            )
            
            return selected_indices, status_msg
        else:
            return None, ""
    
    @app.callback(
        Output(f'{prefix}-download-menu', 'disabled'),
        [Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def toggle_download_menu(selected_cells):
        """Enable download menu when cells are selected"""
        return not bool(selected_cells)
    
    @app.callback(
        Output(f'{prefix}-download-cells-data', 'data'),
        [Input(f'{prefix}-download-cellids', 'n_clicks'),
         Input(f'{prefix}-download-adata', 'n_clicks')],
        [State(f'{prefix}-selected-cells-store', 'data')],
        prevent_initial_call=True
    )
    def download_selected_cells(n_clicks_txt, n_clicks_h5ad, selected_cells):
        """Download selected cell IDs or subset AnnData"""
        ctx = callback_context
        if not ctx.triggered or not selected_cells:
            raise PreventUpdate
        
        # Determine which button was clicked
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        if f'{prefix}-download-cellids' in button_id:
            content = "\n".join(selected_cells)
            return dict(
                content=content,
                filename="selected_cells.txt"
            )
        elif f'{prefix}-download-adata' in button_id:
            import tempfile
            import os
            
            subset_adata = adata[selected_cells].copy()
            
            with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp:
                subset_adata.write_h5ad(tmp.name)
                tmp_path = tmp.name
            
            with open(tmp_path, 'rb') as f:
                content = f.read()
            
            os.unlink(tmp_path)
            
            import base64
            return dict(
                content=base64.b64encode(content).decode(),
                filename=f"subset_{len(selected_cells)}_cells.h5ad",
                base64=True
            )
        
        raise PreventUpdate
    
    # ===== Other Plots Callbacks =====
    
    @app.callback(
        Output(f'{prefix}-violin2-group-selection', 'options'),
        Output(f'{prefix}-violin2-group-selection', 'value'),
        [Input(f'{prefix}-meta1-selection', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_group_labels(selected_column, selected_cells):
        # Filter adata if cells are selected
        if selected_cells:
            filtered_adata = adata[selected_cells]
            unique_labels = sorted(filtered_adata.obs[selected_column].dropna().unique(), key=str)
        else:
            unique_labels = sorted(adata.obs[selected_column].dropna().unique(), key=str)
        
        options = [{'label': str(label), 'value': str(label)} for label in unique_labels]
        values = [str(label) for label in unique_labels]
        return options, values
    
    @app.callback(
        Output(f'{prefix}-single-cell-genes-selection', 'options'),
        Input(f'{prefix}-single-cell-genes-selection', 'search_value'),
        State(f'{prefix}-single-cell-genes-selection', 'value')
    )
    def update_genes_dropdown(search_value, value):
        if not search_value:
            raise PreventUpdate
        label_list = adata.var_names.to_list()
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        selected_labels = value if value else []
        all_labels = list(set(selected_labels + matching_labels[:10]))
        return [{'label': label, 'value': label} for label in all_labels]
    
    @app.callback(
        [Output(f'{prefix}-single-cell-label-selection', 'options'),
         Output(f'{prefix}-single-cell-label-selection', 'value')],
        [Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_labels_based_on_annotation(selected_annotation, selected_cells):
        # Filter adata if cells are selected
        if selected_cells:
            filtered_adata = adata[selected_cells]
            unique_labels = filtered_adata.obs[selected_annotation].unique().tolist()
        else:
            unique_labels = adata.obs[selected_annotation].unique().tolist()
        
        label_options = [{'label': label, 'value': label} for label in sorted(unique_labels)]
        label_values = sorted(unique_labels)
        return label_options, label_values
    
    @app.callback(
        Output(f'{prefix}-heatmap', 'figure'),
        [Input(f'{prefix}-single-cell-genes-selection', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-heatmap-transformation', 'value'),
         Input(f'{prefix}-heatmap-colorscale-dropdown', 'value'),
         Input(f'{prefix}-heatmap-label-dropdown', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-heatmap-secondary-colormap-dropdown', 'value'),  # Add secondary colormap input
         Input(f'{prefix}-selected-cells-store', 'data'),
         Input(f'{prefix}-single-cell-tabs', 'value')],
        [State(f'{prefix}-heatmap', 'figure')]  # Keep current figure as state
    )
    def update_heatmap(selected_genes, selected_annotation, selected_labels, transformation, heatmap_color, secondary_annotation, discrete_color_map, secondary_colormap, selected_cells, active_tab, current_figure):
        # Lazy loading: only update if this tab is active
        if active_tab != 'heatmap-tab':
            # Return the current figure if it exists, otherwise return empty
            return current_figure if current_figure else go.Figure()
        
        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
        
        groupby1_label_color_map = None
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            unique_labels1 = sorted(adata.obs[selected_annotation].unique())
            groupby1_label_color_map = {
                label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels1)
            }
        groupby2_label_color_map = None
        if secondary_annotation and secondary_annotation != 'None' and secondary_annotation != selected_annotation:
            unique_labels2 = sorted(adata.obs[secondary_annotation].unique())
            if secondary_colormap:
                secondary_palette = palette_json["color_palettes"][secondary_colormap]
                groupby2_label_color_map = {
                    label: secondary_palette[i % len(secondary_palette)] for i, label in enumerate(unique_labels2)
                }
        return plot_unified_heatmap(
            adata=filtered_adata,
            genes=selected_genes,
            groupby1=selected_annotation,
            groupby2=secondary_annotation if secondary_annotation and secondary_annotation != 'None' and secondary_annotation != selected_annotation else None,
            labels=selected_labels,
            log=(transformation == 'log'),
            z_score=(transformation == 'z_score'),
            boundary=1,
            color_map=heatmap_color,
            groupby1_label_color_map=groupby1_label_color_map,
            groupby2_label_color_map=groupby2_label_color_map,
            transformation=transformation, # For compatibility
            adata_obs=adata.obs,  # Pass original observations for consistent color mapping
        )
    
    # Add store for violin plot cache
    @app.callback(
        Output(f'{prefix}-violin-plot-cache-store', 'data'),
        [Input(f'{prefix}-single-cell-genes-selection', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-violin-log-or-zscore', 'value'),
         Input(f'{prefix}-show-box1', 'value'),
         Input(f'{prefix}-show-scatter1', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')],
        [State(f'{prefix}-violin-plot-cache-store', 'data')]
    )
    def update_violin_cache(selected_genes, selected_annotation, selected_labels,
                           transformation, show_box_plot, show_scatter1, discrete_color_map, 
                           selected_cells, current_cache):
        cache_key = f"{selected_genes}_{selected_annotation}_{selected_labels}_{transformation}_{show_box_plot}_{show_scatter1}_{discrete_color_map}_{selected_cells}"
        
        if current_cache is None:
            current_cache = {}
        
        if 'current_key' in current_cache and current_cache['current_key'] == cache_key:
            return current_cache
        
        color_map = None
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            unique_labels = sorted(adata.obs[selected_annotation].unique())
            color_map = {
                label: discrete_palette[i % len(discrete_palette)] for i, label in enumerate(unique_labels)
            }
        
        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
        
        fig = plot_violin1(
            filtered_adata,
            selected_genes,
            selected_labels,
            groupby=selected_annotation,
            transformation=transformation,
            show_box='show' in show_box_plot if show_box_plot else False,
            show_points='show' in show_scatter1 if show_scatter1 else False,
            groupby_label_color_map=color_map,
            adata_obs=adata.obs  # Pass original observations for consistent color mapping
        )
        
        num_genes = len(selected_genes) if selected_genes else 0
        num_categories = len(selected_labels) if selected_labels else 0
        fig.update_layout(
            height=min(1000, max(300, 80 * num_genes)),
            width=min(500, max(200, 110 * num_categories)),
            margin=dict(l=130, r=10, t=30, b=30)
        )
        
        if len(current_cache) > 10:  
            keys_to_remove = list(current_cache.keys())[:-9]  
            for key in keys_to_remove:
                if key != 'current_key':
                    current_cache.pop(key, None)
        
        current_cache[cache_key] = fig.to_dict()
        current_cache['current_key'] = cache_key
        
        return current_cache

    @app.callback(
        Output(f'{prefix}-violin-plot1', 'figure'),
        [Input(f'{prefix}-violin-plot-cache-store', 'data'),
         Input(f'{prefix}-single-cell-tabs', 'value')],
        [State(f'{prefix}-violin-plot1', 'figure')]
    )
    def display_violin1(cache_data, active_tab, current_figure):
        if active_tab != 'violin-tab':
            return current_figure if current_figure else go.Figure()
        
        if cache_data and 'current_key' in cache_data:
            current_key = cache_data['current_key']
            if current_key in cache_data:
                return go.Figure(cache_data[current_key])
        
        return current_figure if current_figure else go.Figure()
    
    @app.callback(
        Output(f'{prefix}-mode-explanation', 'children'),
        Input(f'{prefix}-mode-selection', 'value')
    )
    def update_mode_explanation(mode):
        explanations = {
            'mode1': "Compare expression across groups in obs1 only. Obs2 will be ignored.",
            'mode2': "Create facets by obs1, compare obs2 groups within each facet.",
            'mode3': "Linear model treating obs2 as a confounder: expression ~ obs1 + obs2",
            'mode4': "Mixed model treating obs2 as random effect: expression ~ meta1 + (1|obs2)"
        }
        return explanations.get(mode, "")
    
    @app.callback(
        Output(f'{prefix}-test-method-selection', 'options'),
        Output(f'{prefix}-test-method-selection', 'value'),
        [Input(f'{prefix}-mode-selection', 'value'),
         Input(f'{prefix}-meta1-selection', 'value'),
         Input(f'{prefix}-meta2-selection', 'value')]
    )
    def update_test_methods(mode, meta1, meta2):
        # Always include auto and none
        base_options = [
            {'label': 'None', 'value': 'none'}
        ]
        
        if mode == 'mode1':
            n_levels = len(adata.obs[meta1].unique()) if meta1 else 0
            if n_levels == 2:
                options = base_options + [
                    {'label': 'Mann-Whitney U', 'value': 'mwu-test'},
                    {'label': 'T-test', 'value': 'ttest'}
                ]
            else:
                options = base_options + [
                    {'label': 'Kruskal-Wallis', 'value': 'kw-test'},
                    {'label': 'ANOVA', 'value': 'anova'}
                ]
        
        elif mode == 'mode2':
            if meta2 and meta2 != 'none':
                n_levels = len(adata.obs[meta2].unique())
                if n_levels == 2:
                    options = base_options + [
                        {'label': 'Mann-Whitney U', 'value': 'mwu-test'},
                        {'label': 'T-test', 'value': 'ttest'}
                    ]
                else:
                    options = base_options + [
                        {'label': 'Kruskal-Wallis', 'value': 'kw-test'},
                        {'label': 'ANOVA', 'value': 'anova'}
                    ]
            else:
                options = base_options
        
        elif mode == 'mode3':
            options = base_options + [
                {'label': 'Linear Model', 'value': 'linear-model'}
            ]
        
        elif mode == 'mode4':
            options = base_options + [
                {'label': 'Mixed Model', 'value': 'mixed-model'}
            ]
        
        return options, 'auto'
    
    @app.callback(
        Output(f'{prefix}-violin2-gene-selection', 'options'),
        Input(f'{prefix}-violin2-gene-selection', 'search_value')
    )
    def update_violin_genes_dropdown(search_value):
        if not search_value:
            raise PreventUpdate
        label_list = adata.var_names.to_list()
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        return [{'label': label, 'value': label} for label in matching_labels[:10]]
    
    @app.callback(
        [Output(f'{prefix}-meta2-selection', 'disabled'),
         Output(f'{prefix}-meta2-selection', 'value')],
        Input(f'{prefix}-mode-selection', 'value')
    )
    def toggle_meta2_dropdown(mode):
        """Disable meta2 dropdown when mode1 is selected."""
        if mode == 'mode1':
            return True, 'none' 
        else:
            return False, dash.no_update  
    
    @app.callback(
        Output(f'{prefix}-violin-plot2', 'figure'),
        [Input(f'{prefix}-violin2-gene-selection', 'value'),
         Input(f'{prefix}-meta1-selection', 'value'),
         Input(f'{prefix}-meta2-selection', 'value'),
         Input(f'{prefix}-mode-selection', 'value'),
         Input(f'{prefix}-test-method-selection', 'value'),
         Input(f'{prefix}-show-box2', 'value'),
         Input(f'{prefix}-show-scatter2', 'value'),
         Input(f'{prefix}-violin2-log-or-zscore', 'value'),
         Input(f'{prefix}-violin2-group-selection', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_violin2(gene_selection, meta1, meta2, mode, test_method,
                       show_box2, show_points, transformation, labels, selected_cells):
        if selected_cells:
            filtered_adata = filter_data(adata, None, None, selected_cells)
        else:
            # If no cells selected but labels are provided, filter by those
            if labels and meta1:
                filtered_adata = filter_data(adata, meta1, labels, None)
            else:
                filtered_adata = adata
        
        if mode == 'mode1':
            meta2 = None
        elif meta2 == 'none':
            meta2 = None
            
        if mode in ['mode2', 'mode3', 'mode4'] and meta2 is None:
            raise PreventUpdate
            
        return plot_violin2_new(
            filtered_adata,
            key=gene_selection,
            meta1=meta1,
            meta2=meta2,
            mode=mode,
            transformation=transformation,
            show_box='show' in show_box2 if show_box2 else False,
            show_points='show' in show_points if show_points else False,
            test_method=test_method,
            labels=labels,
            color_map=None
        )
    
    @app.callback(
        Output(f'{prefix}-dotplot', 'figure'),
        [Input(f'{prefix}-single-cell-genes-selection', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-plot-type-switch', 'value'),
        Input(f'{prefix}-dotplot-log-or-zscore', 'value'),
         Input(f'{prefix}-dotplot-standardization', 'value'),
         Input(f'{prefix}-dotmatrix-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data'),
         Input(f'{prefix}-single-cell-tabs', 'value')],
        [State(f'{prefix}-dotplot', 'figure')]  # Keep current figure
    )
    def update_dotplot(selected_genes, selected_annotation, selected_labels, plot_type,
                       transformation, standardization, color_map, selected_cells, active_tab, current_figure):
        # Lazy loading: only update if this tab is active
        if active_tab != 'dotplot-tab':
            return current_figure if current_figure else go.Figure()

        if adata.n_obs > 10000 or (hasattr(adata, 'isbacked') and adata.isbacked):
            # Pass original adata - the optimized function will handle filtering internally
            return plot_dot_matrix(
                adata,
                selected_genes,
                selected_annotation,
                selected_labels,
                aggregation='mean',  # Default to mean aggregation
                transformation=transformation,
                standardization=standardization,
                color_map=color_map,
                plot_type=plot_type
            )
        else:
            filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
            return plot_dot_matrix(
                filtered_adata,
                selected_genes,
                selected_annotation,
                selected_labels,
                aggregation='mean',  # Default to mean aggregation
                transformation=transformation,
                standardization=standardization,
                color_map=color_map,
                plot_type=plot_type
            )
    
    @app.callback(
        [Output(f'{prefix}-x-axis-draggable-grid', 'children'),
         Output(f'{prefix}-x-axis-groups-state', 'data')],
        [Input(f'{prefix}-x-meta-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-single-cell-tabs', 'value')],
        [State(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         State(f'{prefix}-x-axis-groups-state', 'data')]
    )
    def update_x_axis_groups_grid(x_meta, selected_labels, active_tab, y_meta, current_state):
        if active_tab != 'stacked-bar-tab' or not x_meta:
            return [], {}
        
        x_values = sorted(adata.obs[x_meta].unique())

        if selected_labels and len(selected_labels) < len(x_values):
            # Check if selected labels are a subset of x_values
            if set(selected_labels).issubset(set(x_values)):
                x_values = selected_labels
        
        # Convert to strings for consistency
        x_values = [str(val) for val in x_values]
        
        # Initialize state for new values
        if not current_state:
            current_state = {}
        
        # Update state to include all values (default to enabled)
        new_state = {}
        for val in x_values:
            new_state[val] = current_state.get(val, True)
        
        # Create draggable items
        items = []
        for i, value in enumerate(x_values):
            is_enabled = new_state.get(value, True)
            # Create the item content
            items.append(
                html.Div([
                    html.Div(value, 
                            style={
                                'fontWeight': 'bold',
                                'marginBottom': '5px',
                                'color': '#000' if is_enabled else '#999'
                            }),
                    dbc.Switch(
                        id={'type': f'{prefix}-x-group-switch', 'index': value},
                        value=is_enabled,
                        style={'marginTop': '5px'}
                    )
                ], 
                key=f'x-group-{value}',
                style={
                    'padding': '10px',
                    'backgroundColor': '#fff' if is_enabled else '#f5f5f5',
                    'border': '2px solid #007bff' if is_enabled else '1px solid #ddd',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'cursor': 'move',
                    'height': '100%'
                })
            )
        
        if not items:
            items = [html.Div("No groups available", 
                            key='empty-msg',
                            style={'color': '#6c757d', 'fontStyle': 'italic', 'padding': '20px'})]
        
        return items, new_state
    
    @app.callback(
        Output(f'{prefix}-x-axis-groups-state', 'data', allow_duplicate=True),
        [Input({'type': f'{prefix}-x-group-switch', 'index': ALL}, 'value')],
        [State({'type': f'{prefix}-x-group-switch', 'index': ALL}, 'id'),
         State(f'{prefix}-x-axis-groups-state', 'data')],
        prevent_initial_call=True
    )
    def update_group_state(switch_values, switch_ids, current_state):
        if not switch_ids:
            return current_state
        
        new_state = current_state.copy() if current_state else {}
        
        for i, switch_id in enumerate(switch_ids):
            if 'index' in switch_id:
                group_name = switch_id['index']
                new_state[group_name] = switch_values[i]
        
        return new_state
    
    @app.callback(
        [Output(f'{prefix}-stacked-bar-x-order-grid', 'columnDefs'),
         Output(f'{prefix}-stacked-bar-x-order-grid', 'rowData')],
        [Input(f'{prefix}-stacked-bar-x-axis', 'value'),
         Input(f'{prefix}-single-cell-tabs', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]  # Add selected cells
    )
    def update_x_axis_order_grid(x_axis_meta, active_tab, selected_cells):
        if active_tab != 'stacked-bar-tab' or not x_axis_meta:
            return [], []
        
        if selected_cells:
            filtered_adata = adata[selected_cells]
            x_values = sorted(filtered_adata.obs[x_axis_meta].unique())
        else:
            x_values = sorted(adata.obs[x_axis_meta].unique())
        x_values_str = [str(val) for val in x_values]
        
        # Create column definitions - each x-axis group becomes a column
        column_defs = []
        for val in x_values_str:
            column_defs.append({
                "field": val,
                "headerName": val,
                "width": 150,
                "minWidth": 120,
                "suppressMovable": False,  # Allow dragging
                "headerClass": "ag-header-cell-center",
                "resizable": True
            })
        
        # No rows, only headers
        row_data = []
        
        return column_defs, row_data
    
    # Store to keep track of column order
    @app.callback(
        Output(f'{prefix}-x-axis-column-order-store', 'data'),
        Input(f'{prefix}-stacked-bar-x-order-grid', 'columnState'),
        prevent_initial_call=True
    )
    def update_column_order(column_state):
        if not column_state:
            return []
        
        # Extract column order from column state
        column_order = []
        for col in column_state:
            if 'colId' in col:
                column_order.append(col['colId'])
        
        return column_order
    
    # Stacked bar plot callback
    @app.callback(
        Output(f'{prefix}-stacked-bar-plot', 'figure'),
        [Input(f'{prefix}-norm-box', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data'),
         Input(f'{prefix}-single-cell-tabs', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-stacked-bar-x-axis', 'value'),
         Input(f'{prefix}-x-axis-column-order-store', 'data')],
        [State(f'{prefix}-stacked-bar-plot', 'figure')]
    )
    def update_stacked_bar(norm, discrete_color_map, selected_cells, active_tab, 
                          annotation, selected_labels, x_axis_meta, x_axis_order, current_figure):
        # Lazy loading: only update if this tab is active
        if active_tab != 'stacked-bar-tab':
            return current_figure if current_figure else go.Figure()
        
        if not annotation or not x_axis_meta:
            # Return empty figure with message
            fig = go.Figure()
            fig.add_annotation(
                text="Please select both annotation (for stacking) and x-axis metadata",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False,
                font=dict(size=14),
                xanchor="center",
                yanchor="middle"
            )
            fig.update_layout(
                plot_bgcolor='white',
                paper_bgcolor='white',
                xaxis=dict(visible=False),
                yaxis=dict(visible=False)
            )
            return fig
        
        filtered_adata = filter_data(adata, annotation, selected_labels, selected_cells)
        
        if x_axis_order and len(x_axis_order) > 0:
            # Use the order from dragged columns
            final_x_order = x_axis_order
        elif x_axis_meta:
            # Default to sorted order if no dragging has occurred
            x_values = sorted(adata.obs[x_axis_meta].unique())
            final_x_order = [str(val) for val in x_values]
        else:
            final_x_order = None
        
        # Create fixed color mapping based on annotation categories
        all_categories = sorted(adata.obs[annotation].unique())
        
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            fixed_color_map = {
                cat: discrete_palette[i % len(discrete_palette)] 
                for i, cat in enumerate(all_categories)
            }
        else:
            fixed_color_map = {
                cat: color_config[i % len(color_config)] 
                for i, cat in enumerate(all_categories)
            }
        
        # Use x_axis_meta for x-axis and annotation for stacking (y_meta)
        return plot_stacked_bar(
            x_meta=x_axis_meta,  # From the dropdown in the stacked bar tab
            y_meta=annotation,   # From the left control panel
            norm=norm, 
            adata=filtered_adata, 
            color_map=fixed_color_map,
            y_order=selected_labels,  # Order from Select Labels in left control
            x_order=final_x_order  # Order from the draggable AgGrid or default
        )
    
    
    # Pseudotime plot callback
    @app.callback(
        Output(f'{prefix}-pseudotime-plot', 'figure'),
        [Input(f'{prefix}-single-cell-tabs', 'value'),
         Input(f'{prefix}-single-cell-genes-selection', 'value'),
         Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-single-cell-label-selection', 'value'),
         Input(f'{prefix}-pseudotime-min-expr-slider', 'value'),
         Input(f'{prefix}-pseudotime-transformation', 'value'),
         Input(f'{prefix}-pseudotime-key-dropdown', 'value'),
         Input(f'{prefix}-discrete-color-map-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')],
        [State(f'{prefix}-pseudotime-plot', 'figure')]  # Keep current figure
    )
    def update_pseudotime_plot(selected_tab, selected_genes, selected_annotation, selected_labels,
                               min_expr, transformation, pseudotime_key,
                               discrete_color_map, selected_cells, current_figure):
        # Only process when pseudotime tab is active
        if selected_tab != 'pseudotime-tab':
            return current_figure if current_figure else go.Figure()
        
        if not selected_genes:
            # Return empty figure if no genes selected
            fig = go.Figure()
            fig.add_annotation(
                text="<b>Please select genes to plot</b><br><br>Use the 'Select Variables' dropdown on the left to choose genes",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False,
                font=dict(size=16),
                align="center"
            )
            fig.update_layout(
                plot_bgcolor='#f8f9fa',
                paper_bgcolor='white',
                height=400,
                margin=dict(t=50, b=50, l=50, r=50)
            )
            return fig
        
        # Filter data
        filtered_adata = filter_data(adata, selected_annotation, selected_labels, selected_cells)
        
        # Get color map
        if discrete_color_map:
            discrete_palette = palette_json["color_palettes"][discrete_color_map]
            all_categories = sorted(filtered_adata.obs[selected_annotation].unique())
            color_map = {
                cat: discrete_palette[i % len(discrete_palette)] 
                for i, cat in enumerate(all_categories)
            }
        else:
            all_categories = sorted(filtered_adata.obs[selected_annotation].unique())
            color_map = {
                cat: color_config[i % len(color_config)] 
                for i, cat in enumerate(all_categories)
            }
        
        # Use default pseudotime key if not specified
        if not pseudotime_key:
            # Try to find a pseudotime column
            numeric_cols = []
            for col in filtered_adata.obs.columns:
                if filtered_adata.obs[col].dtype in ['float32', 'float64', 'int32', 'int64']:
                    if filtered_adata.obs[col].nunique() > 50:
                        numeric_cols.append(col)
            
            # Prioritize columns with 'pseudotime' in the name
            pseudotime_cols = [col for col in numeric_cols if 'pseudotime' in col.lower() or 'dpt' in col.lower()]
            
            if pseudotime_cols:
                pseudotime_key = pseudotime_cols[0]
            elif numeric_cols:
                pseudotime_key = numeric_cols[0]
            else:
                # No suitable pseudotime column found
                fig = go.Figure()
                fig.add_annotation(
                    text="<b>No pseudotime column found</b><br><br>Please ensure your data has a numeric column with pseudotime values",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False,
                    font=dict(size=16),
                    align="center"
                )
                fig.update_layout(
                    plot_bgcolor='#f8f9fa',
                    paper_bgcolor='white',
                    height=400,
                    margin=dict(t=50, b=50, l=50, r=50)
                )
                return fig
        
        return plot_genes_in_pseudotime(
            filtered_adata,
            selected_genes,
            pseudotime_key=pseudotime_key,
            groupby=selected_annotation,
            min_expr=min_expr,
            transformation=transformation,
            color_map=color_map
        )
    
    @app.callback(
        [Output(f'{prefix}-pseudotime-key-dropdown', 'options'),
         Output(f'{prefix}-pseudotime-key-dropdown', 'value')],
        Input(f'{prefix}-single-cell-tabs', 'value')
    )
    def update_pseudotime_key_options(active_tab):
        if active_tab == 'pseudotime-tab':
            numeric_cols = []
            for col in adata.obs.columns:
                if adata.obs[col].dtype in ['float32', 'float64', 'int32', 'int64']:
                    if adata.obs[col].nunique() > 50:
                        numeric_cols.append(col)
            pseudotime_cols = [col for col in numeric_cols if 'pseudotime' in col.lower() or 'dpt' in col.lower()]
            other_cols = [col for col in numeric_cols if col not in pseudotime_cols]
            
            options = []
            all_cols = pseudotime_cols + other_cols
            for col in all_cols:
                options.append({'label': col, 'value': col})
            
            # If no pseudotime columns found, return empty with message
            if not options:
                options = [{'label': 'No pseudotime columns found', 'value': None}]
                return options, None
            
            # Set default value to first pseudotime column or first numeric column
            default_value = all_cols[0] if all_cols else None
            
            return options, default_value
        return [], None
    
    # Add callback to update filter status
    @app.callback(
        Output(f'{prefix}-filter-status', 'children'),
        Input(f'{prefix}-selected-cells-store', 'data')
    )
    def update_filter_status(selected_cells):
        if selected_cells:
            n_selected = len(selected_cells)
            n_total = adata.n_obs
            return dbc.Alert(
                f"Showing {n_selected} of {n_total} cells ({n_selected/n_total*100:.1f}%) based on scatter plot selection",
                color="info",
                dismissable=False,
                style={'margin': '0'}
            )
        return None
