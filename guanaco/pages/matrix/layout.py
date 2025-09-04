from pathlib import Path
import json
from dash import dcc, html
import dash_bootstrap_components as dbc
import plotly.express as px
import dash_draggable
import dash_ag_grid as dag
from guanaco.config import common_config


# Load color palettes
cvd_color_path = Path(__file__).parent / "cvd_color.json"
with open(cvd_color_path, "r") as f:
    palette_json = json.load(f)
palette_names = list(palette_json["color_palettes"].keys())

# Get the available color scales
colorscales = px.colors.named_colorscales()

# 1.Heatmap Layyout
def generate_heatmap_layout(adata, prefix):
    label_list = adata.obs.columns.to_list()

    # --- Transformation selection ---
    heatmap_transformation_selection =   html.Div([
        html.Label(
            'Transformation:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-heatmap-transformation',
            options=[
                {'label': 'None', 'value': 'None'},
                {'label': 'Log', 'value': 'log'},
                {'label': 'Z-score (across cell)', 'value': 'z_score'}
            ],
            value='None',
            inline=True,
            style={'marginBottom': '10px'}
        )
    ])

    # --- Secondary dropdown ---
    heatmap_secondary_dropdown =  html.Div([
        html.Label('Secondary Annotation:', style={'fontWeight': 'bold', 'marginTop': '10px'}),
        dcc.Dropdown(
            id=f'{prefix}-heatmap-label-dropdown',
            options=[{'label': 'None', 'value': 'None'}] + [{'label': label, 'value': label} for label in label_list],
            value='None',
            style={'width': '200px'}
        )
    ], style={'marginBottom': '10px'})

    # --- Color map selection ---
    dotmatrix_color_map_dropdown = html.Div([
        html.Label(
            'Color Map:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dcc.Dropdown(
            id=f'{prefix}-heatmap-colorscale-dropdown',
            options=[
                {'label': scale, 'value': scale} for scale in colorscales
            ],
            value='viridis', 
            clearable=False,
            style={'width': '200px', 'marginBottom': '10px'}
        )
    ])

    # --- Secondary Annotation Color Map ---
    secondary_annotation_colormap_dropdown = html.Div([
        html.Label(
            'Secondary Annotation ColorMap:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dcc.Dropdown(
            id=f'{prefix}-heatmap-secondary-colormap-dropdown',
            options=[{'label': name, 'value': name} for name in palette_names],
            value='tab20',  # Default to tab20 as you requested
            placeholder="Select colormap for secondary annotation",
            clearable=True,
            style={'width': '200px', 'marginBottom': '10px'}
        )
    ])

    # --- Draggable container with plot ---
    draggable_container = dash_draggable.GridLayout(
        id=f'{prefix}-draggable-heatmap',
        className='grid-layout-no-border',
        children=[
            html.Div(children=[
               dcc.Graph(
                  id=f'{prefix}-heatmap',
                  config=common_config,
                  responsive=True,
                  style={
                        "min-height":"0",
                        "flex-grow":"1"
                  })],
               style={
                  "height":'100%',
                  "width":'100%',
                  "display":"flex",
                  "flex-direction":"column",
                  "flex-grow":"0"
               }),
      ])
    
    # --- Final layout container ---
    heatmap_layout = html.Div([
        heatmap_transformation_selection,
        heatmap_secondary_dropdown,
        dotmatrix_color_map_dropdown,
        secondary_annotation_colormap_dropdown,
        draggable_container
    ], style={'padding': '20px', 'marginBottom': '15px'})

    return heatmap_layout

# 2.violin layout
def generate_violin_layout(default_gene_markers,discrete_label_list,prefix):
    # Pre-create some common dropdowns and toggles
    violin1_transformation_selection = dbc.RadioItems(
        id=f'{prefix}-violin-log-or-zscore',
        options=[
            {'label': 'None', 'value': 'False'},
            {'label': 'Log', 'value': 'log'}
        ],
        value='False',
        inline=True,
        style={'marginLeft': '10px'}
    )

    violin2_transformation_selection = dbc.RadioItems(
        id=f'{prefix}-violin2-log-or-zscore',
        options=[
            {'label': 'None', 'value': 'False'},
            {'label': 'Log', 'value': 'log'}
        ],
        value='False',
        inline=False
    )

    violin_show_box1 = dbc.Checklist(
        id=f'{prefix}-show-box1',
        options=[{'label': 'Show Box Plot', 'value': 'show'}],
        value=[],
        switch=True
    )

    violin_show_scatter1 = dbc.Checklist(
        id=f'{prefix}-show-scatter1',
        options=[{'label': 'Show Scatter Points', 'value': 'show'}],
        value=[],
        switch=True
    )

    violin_show_box2 = dbc.Checklist(
        id=f'{prefix}-show-box2',
        options=[{'label': 'Show Box Plot', 'value': 'show'}],
        value=[],
        switch=True
    )

    violin_show_scatter2 = dbc.Checklist(
        id=f'{prefix}-show-scatter2',
        options=[{'label': 'Show Scatter Points', 'value': 'show'}],
        value=[],
        switch=True
    )

    # Meta1 dropdown (primary metadata)
    meta1_selection = dcc.Dropdown(
        id=f'{prefix}-meta1-selection',
        options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
        value=discrete_label_list[0],
        clearable=False,
        placeholder="Select primary metadata"
    )

    # Meta2 dropdown (secondary metadata - optional)
    meta2_selection = dcc.Dropdown(
        id=f'{prefix}-meta2-selection',
        options=[{'label': 'None', 'value': 'none'}] + [{'label': meta, 'value': meta} for meta in discrete_label_list],
        value='none',
        clearable=False,
        placeholder="Select secondary metadata (optional)"
    )

    # Mode selection dropdown
    mode_selection = dcc.Dropdown(
        id=f'{prefix}-mode-selection',
        options=[
            {'label': 'Mode 1: One metadata only', 'value': 'mode1'},
            {'label': 'Mode 2: Facet by meta1, compare meta2', 'value': 'mode2'},
            {'label': 'Mode 3: Linear model (meta1 + confounder)', 'value': 'mode3'},
            {'label': 'Mode 4: Mixed model (meta1 + random effect)', 'value': 'mode4'},
        ],
        value='mode1',
        clearable=False
    )

    # Test method dropdown
    test_method_selection = dcc.Dropdown(
        id=f'{prefix}-test-method-selection',
        options=[
            {'label': 'None', 'value': 'none'},
            {'label': 'Mann-Whitney U', 'value': 'mwu-test'},
            {'label': 'T-test', 'value': 'ttest'},
            {'label': 'Kruskal-Wallis', 'value': 'kw-test'},
            {'label': 'ANOVA', 'value': 'anova'},
            {'label': 'Linear Model', 'value': 'linear-model'},
            {'label': 'Mixed Model', 'value': 'mixed-model'}
        ],
        value='none',
        clearable=False
    )

    violin2_gene_selection = dcc.Dropdown(
        id=f'{prefix}-violin2-gene-selection',
        options=[{'label': gene, 'value': gene} for gene in default_gene_markers],
        value=default_gene_markers[0] if default_gene_markers else None
    )

    group_selection = dcc.Dropdown(
        id=f'{prefix}-violin2-group-selection',
        multi=True,
        style={'font-size': '12px'}
    )

    # Final layout
    violin_layout = html.Div([
        # Violin Plot 1
        html.Div([
            html.Div([
                html.Label("Transformation:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                violin1_transformation_selection
            ], style={'marginBottom': '15px'}),
            html.Div([
                html.Label("More Options:", style={'fontWeight': 'bold'}),
                html.Div([violin_show_box1, violin_show_scatter1],
                        style={'display': 'flex', 'gap': '20px'})
            ], style={'marginBottom': '15px'}),
            # Add cache store for violin plot optimization
            dcc.Store(id=f'{prefix}-violin-plot-cache-store'),
            dcc.Loading(
                id="loading-violin1",
                type="circle",
                children=[
                    dcc.Graph(
                        id=f'{prefix}-violin-plot1',
                        config=common_config,
                        style={'width': '100%', 'minHeight': '400px'}
                    )
                ]
            ),
        ], style={'marginBottom': '30px', 'padding': '10px'}),

        # Violin Plot 2
        html.Div([
            html.H4("Split Violin/Grouped Violin", style={'textAlign': 'center', 'margin': '10px 0', 'fontWeight': 'bold'}),
            html.Label("Select Gene", style={'fontWeight': 'bold'}),
            violin2_gene_selection,
            html.Div([
                html.Label("Transformation:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                violin2_transformation_selection
            ], style={'marginBottom': '10px'}),
            html.Div([
                html.Div([
                    html.Label("Obs1 (Primary):", style={'fontWeight': 'bold'}),
                    meta1_selection
                ], style={'flex': '1'}),
                html.Div([
                    html.Label("Obs2 (Secondary):", style={'fontWeight': 'bold'}),
                    meta2_selection
                ], style={'flex': '1'}),
            ], style={'display': 'flex', 'marginBottom': '10px', 'gap': '10px'}),

            html.Div([
                html.Div([
                    html.Label("Analysis Mode:", style={'fontWeight': 'bold'}),
                    mode_selection
                ], style={'flex': '1'}),
                html.Div([
                    html.Label("Statistical Test:", style={'fontWeight': 'bold'}),
                    test_method_selection
                ], style={'flex': '1'}),
            ], style={'display': 'flex', 'marginBottom': '10px', 'gap': '10px'}),
            
            # Mode explanation helper text
            html.Div(
                id=f'{prefix}-mode-explanation',
                style={'fontSize': '12px', 'color': 'gray', 'marginBottom': '10px', 'fontStyle': 'italic'}
            ),

            html.Label("More Options:", style={'fontWeight': 'bold'}),
            violin_show_box2,
            violin_show_scatter2,

            # Draggable Graph for Violin 2
            dcc.Loading(
                id="loading-violin2",
                type="circle",
                children=[
                    html.Div([
                        dash_draggable.GridLayout(
                            id=f'{prefix}-draggable-violin2',
                            className='grid-layout-no-border',
                            children=[
                                html.Div(
                                    children=dcc.Graph(
                                        id=f'{prefix}-violin-plot2',
                                        config=common_config,
                                        responsive=True,
                                        style={
                                            "min-height": "0",
                                            "flex-grow": "1",
                                            "border": "none",
                                            "box-shadow": "none"
                                        }
                                    ),
                                    style={
                                        "height": '100%',
                                        "width": '100%',
                                        "display": "flex",
                                        "flex-direction": "column",
                                        "flex-grow": "0",
                                        "border": "none",
                                        "box-shadow": "none"
                                    }
                                )
                            ],
                            isResizable=True,
                            isDraggable=True,
                            height=30,
                            gridCols=12,
                            style={
                                "min-width": "300px",   # optional: prevents collapse
                                "min-height": "300px",  # optional: ensures initial size
                                "background": "transparent",
                                "padding": "0px"
                            }
                        )
                    ])
                ],
                style={
                    'padding': '10px',
                    'display': 'flex',
                    'flex-direction': 'column',
                    'height': '100%'
                }
            ),

            html.Label("Labels in Group:", style={'fontWeight': 'bold', 'marginTop': '10px'}),
            group_selection,
        ], style={'padding': '10px'})
    ], style={'width': '100%'})



    return violin_layout
# 3. Dotplot Layout
def generate_dotplot_layout(prefix):

    # --- Plot type switch ---
    plot_type_switch = html.Div([
        html.Label(
            'Plot Type:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-plot-type-switch',
            options=[
                {'label': 'Dotplot', 'value': 'dotplot'},
                {'label': 'Matrixplot', 'value': 'matrixplot'}
            ],
            value='dotplot',
            inline=True,
            style={'marginBottom': '15px'}
        )
    ])


    # --- Transformation selection ---
    dotplot_transformation_selection = html.Div([
        html.Label(
            'Transformation:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-dotplot-log-or-zscore',
            options=[
                {'label': 'None', 'value': 'None'},
                {'label': 'Log', 'value': 'log'}
            ],
            value='None',
            inline=True,
            style={'marginBottom': '10px'}
        )
    ])

    dotplot_standardization_selection = html.Div([
        html.Label(
            'Standardization:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dbc.RadioItems(
            id=f'{prefix}-dotplot-standardization',
            options=[
                {'label': 'None', 'value': 'None'},
                {'label': 'By variable', 'value': 'var'},
                {'label': 'By group', 'value': 'group'}
            ],
            value='None',
            inline=True,
            style={'marginBottom': '10px'}
        )
    ])
    dotmatrix_color_map_dropdown = html.Div([
        html.Label(
            'Color Map:',
            style={'fontWeight': 'bold', 'marginBottom': '5px'}
        ),
        dcc.Dropdown(
            id=f'{prefix}-dotmatrix-color-map-dropdown',
            options=[
                {'label': scale, 'value': scale} for scale in colorscales
            ],
            value='viridis',
            clearable=False,
            style={'width': '200px', 'marginBottom': '10px'}
        )
    ])

    draggable_container = dash_draggable.GridLayout(
    id=f'{prefix}-draggable-dotplot',
    className='grid-layout-no-border',
    children=[
        html.Div([
            dcc.Graph(
                id=f'{prefix}-dotplot',
                config=common_config,
                responsive=True,
                style={
                    'min-height': '0',
                    'flex-grow': '1',
                    'backgroundColor': 'transparent'
                }
            )
        ],
            style={
                'backgroundColor': 'transparent', 
                'border': 'none', 
                'boxShadow': 'none',
                'height': '100%',
                'width': '100%',
                'display': 'flex',
                'flex-direction': 'column',
                'flex-grow': '0'
            }) 
        ],
        isResizable=True,
        isDraggable=True,
        height=30,
        gridCols=12,  
        style={
            'backgroundColor': 'transparent',
            'padding': '0px',
            'border': 'none',
            'boxShadow': 'none'
        }
    )
    # --- Final layout container ---
    dotplot_layout = html.Div([
        dotplot_transformation_selection,
        dotplot_standardization_selection,
        dotmatrix_color_map_dropdown,
        plot_type_switch,
        draggable_container
    ], style={'padding': '20px', 'marginBottom': '15px'})

    return dotplot_layout

# 4. Stacked Bar Layout
def generate_stacked_bar_layout(discrete_label_list, prefix):

    # Dropdown to select x-axis metadata
    x_axis_dropdown = html.Div([
        html.Label("Cell info (x-axis):", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dcc.Dropdown(
            id=f'{prefix}-stacked-bar-x-axis',
            options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
            value=discrete_label_list[0] if discrete_label_list else None,
            clearable=False,
            placeholder="Select metadata for x-axis",
            style={'marginBottom': '15px'}
        )
    ], style={'marginBottom': '15px'})
    
    norm_box = html.Div([
        html.Label("Plot value:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dbc.RadioItems(
            id=f'{prefix}-norm-box',
            options=[
                {'label': 'Proportion', 'value': 'prop'},
                {'label': 'Count', 'value': 'count'},
            ],
            value='prop',
            inline=True
        )
    ], style={'marginBottom': '15px'})
    
    draggable_bar = dash_draggable.GridLayout(
        id=f'{prefix}-draggable',
        className='grid-layout-no-border',
        children=[
            html.Div(children=[
                dcc.Graph(
                    id=f'{prefix}-stacked-bar-plot',
                    config=common_config,
                    responsive=True,
                    style={
                        "min-height": "400px",
                        "flex-grow": "1",
                    }
                )
            ],
            style={
                "height": '100%',
                "width": '100%',
                "display": "flex",
                "flex-direction": "column",
                "flex-grow": "0",
            }),
        ],
        isResizable=True,
        isDraggable=True,
        height=30,
        gridCols=12,
    )
    
    info_text = html.Div([
        html.P([
            html.I(className="fas fa-info-circle", style={'marginRight': '5px'}),
            "The stacked layers come from 'Select Annotation' and 'Select Labels' in the left control panel"
        ], style={'color': '#6c757d', 'fontSize': '14px', 'marginBottom': '15px'})
    ])
    

    x_axis_order_component = html.Div([
        html.Label("X-axis group order:", 
                   style={'fontWeight': 'bold', 'marginBottom': '10px'}),
        dag.AgGrid(
            id=f'{prefix}-stacked-bar-x-order-grid',
            rowData=[],
            columnDefs=[],
            defaultColDef={
                "sortable": False,
                "filter": False,
                "resizable": True,
                "suppressMenu": True,
                "headerHeight": 40,
                "minWidth": 120,
                "width": 150,
                "headerClass": "ag-header-cell-center"
            },
            dashGridOptions={
                "headerHeight": 40,
                "rowHeight": 0,
                "suppressRowClickSelection": True,
                "suppressCellSelection": True,
                "suppressMovableColumns": False, 
                "animateRows": False,
                "suppressHorizontalScroll": False,
                "onColumnMoved": True,
                "suppressLoadingOverlay": True,
                "suppressNoRowsOverlay": True,
                "suppressDisplayTotal": True
            },
            style={"height": "40px", "marginBottom": "10px", "overflow": "hidden"},
            className="ag-theme-alpine"
        ),
        html.P("Drag column headers to reorder x-axis groups.", 
               style={'fontSize': '12px', 'color': '#6c757d', 'marginTop': '5px', 'marginBottom': '15px'})
    ])
    
    column_order_store = dcc.Store(id=f'{prefix}-x-axis-column-order-store', data=[])
    
    bar_layout = html.Div([
        column_order_store,
        info_text,
        x_axis_dropdown,
        norm_box,
        draggable_bar,
        x_axis_order_component 
    ], style={
        'padding': '20px',
        'marginBottom': '15px',
    })

    return bar_layout

# 5. Pseudotime Layout
def generate_pseudotime_layout(prefix):
    """Generate the layout for the pseudotime plot tab."""
    
    layout = html.Div([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.Label('Minimum Expression Threshold:', 
                              className="control-label"),
                    dcc.Slider(
                        id=f'{prefix}-pseudotime-min-expr-slider',
                        min=0,
                        max=5,
                        step=0.1,
                        value=0.5,
                        marks=None,
                        tooltip={"placement": "bottom", "always_visible": True},
                        className="dbc-slider"
                    )
                ], style={'marginBottom': '20px'}),
                
                html.Div([
                    html.Label('Transformation:', 
                              style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    dbc.RadioItems(
                        id=f'{prefix}-pseudotime-transformation',
                        options=[
                            {'label': 'None', 'value': 'none'},
                            {'label': 'Log', 'value': 'log'},
                            {'label': 'Z-score', 'value': 'z_score'}
                        ],
                        value='none',
                        inline=True,
                        style={'fontSize': '14px'}
                    )
                ], style={'marginBottom': '20px'}),
                
                html.Div([
                    html.Label('Pseudotime Column:', 
                              style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    dcc.Dropdown(
                        id=f'{prefix}-pseudotime-key-dropdown',
                        options=[], 
                        value=None, 
                        placeholder='Select pseudotime column',
                        style={'marginBottom': '15px'},
                        clearable=False
                    )
                ], style={'marginBottom': '20px'}),
                
            ], width=12)
        ], style={'marginBottom': '20px', 'padding': '15px', 'backgroundColor': '#f8f9fa', 'borderRadius': '5px'}),
        
        dbc.Row([
            dbc.Col([
                dcc.Loading(
                    id=f"{prefix}-pseudotime-loading",
                    type="default",
                    children=[
                        dcc.Graph(
                            id=f'{prefix}-pseudotime-plot',
                            style={'height': 'auto', 'minHeight': '600px'},
                            config={
                                'displayModeBar': True,
                                'displaylogo': False,
                                'modeBarButtonsToRemove': ['lasso2d', 'select2d'],
                                'toImageButtonOptions': {
                                    'format': 'png',
                                    'filename': 'pseudotime_plot',
                                    'height': 800,
                                    'width': 1200,
                                    'scale': 2
                                }
                            }
                        )
                    ]
                )
            ], width=12)
        ])
    ])
    
    return layout
