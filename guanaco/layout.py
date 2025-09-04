import dash_bootstrap_components as dbc
from dash import html, dcc
from guanaco.pages.track.gene_browser import gene_browser_layout
from guanaco.pages.matrix.callbacks import generate_embedding_plots, generate_other_plots
import muon as mu

# tip
def resize_tip_toast():
    return dbc.Toast(
        children=[
            html.P(
                "Drag the draggable corner at the bottom right of any chart to resize it",
                className="mb-1",
            ),
            dbc.Button("Got it", id="close-tip", size="sm", color="primary"),
        ],
        id="tip-modal",
        header=html.Div(
            [
                html.Img(
                    src="/assets/lamp_guanaco.png",
                    style={"height": "1.25rem", "marginRight": "0.5rem"},
                ),
                html.Span("Tip!"),
            ],
            style={"display": "flex", "alignItems": "center"},
        ),
        icon=None,
        duration=30000,
        is_open=False,
        style={
            "position": "fixed",
            "bottom": 20,
            "right": 20,
            "width": 320,
            "zIndex": 1000,
        },
    )

# Footer and footprint
footprint = html.Div(
    style={
        'backgroundImage': 'url("/assets/footprint.png")',
        'backgroundRepeat': 'repeat-x',
        'backgroundPosition': 'left',
        'backgroundSize': 'contain',
        'width': '100%',
        'height': '10px',
        'margin': '20px 0'
    }
)

guanaco_footer = html.Footer(
    html.Div([
        "This webpage was made using ",
        html.A("GUANACO", href="https://github.com/Systems-Immunometabolism-Lab/guanaco-viz", target="_blank"),
        "."
    ],
    style={
        "textAlign": "center",
        "fontSize": "14px",
        "padding": "10px",
        "color": "#666"
    })
)

# Navbar layout
def navbar(datasets):
    return html.Div(
        dbc.Navbar(
            dbc.Container(
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/assets/logo.png", height="70px"), width="auto", align="center"),
                        dbc.Col(
                            dbc.NavLink(
                                "GUANACO",
                                href="/",
                                style={
                                    'fontSize': '36px',
                                    'fontWeight': 'bold',
                                    'color': 'white',
                                    'textDecoration': 'none'
                                }
                            ),
                            width="auto",
                            align="center"
                        ),
                        dbc.Col(width=True),
                        dbc.Col(
                            dbc.Tabs(
                                id="tabs-dataset",
                                active_tab=list(datasets.keys())[0],
                                children=[
                                    dbc.Tab(label=dataset.title, tab_id=name)
                                    for name, dataset in datasets.items()
                                ],
                                className="dataset-tabs"
                            ),
                            width="auto",
                            className="tabs-align-bottom"
                        ),
                    ],
                    align="center",
                    justify="between",
                    className="w-100"
                )
            ),
            color="grey",
            dark=True,
            className="px-3 custom-navbar"
        ),
        style={
            "position": "fixed",
            "top": 0,
            "width": "100%",
            "zIndex": 1000
        }
    )

def description_layout(dataset):
    summary_items = []
    
    # Add description
    summary_items.append(html.P(dataset.description))
    
    # Add AnnData information if available
    if dataset.adata is not None:
        adata = dataset.adata
        # N. Cells (same for all modalities)
        n_cells = adata.n_obs
        summary_items.append(
            html.P([html.B("N. Cells: "), f"{n_cells:,}"], className="summary-item")
        )
        meta_list = []
        # N. Variables per modality
        if isinstance(adata, mu.MuData):
            for mod_name, mod_adata in adata.mod.items():
                meta_list.extend([k for k in mod_adata.obs.keys()])
                summary_items.append(
                    html.P([
                        html.B(f"{mod_name.upper()} - N. Variables: "),
                        f"{mod_adata.n_vars:,}"
                    ], className="summary-item")
                )
        else:
            meta_list.extend([k for k in adata.obs.keys()])
            summary_items.append(
                html.P([
                    html.B("N. Variables: "),
                    f"{adata.n_vars:,}"
                ], className="summary-item")
            )
        meta_list = list(set(meta_list))  # Remove duplicates
        summary_items.append(
            html.P([
                html.B("Metadata: "),
                ", ".join(meta_list)
            ], className="summary-item")
        )
    
    # Add genome browser information if available
    if dataset.genome_tracks is not None:
        track_count = sum(len(tracks) for tracks in dataset.genome_tracks.values())
        summary_items.append(
            html.P([html.B("Genome Browser Tracks: "), f"{track_count}"], className="summary-item")
        )

    # Final layout
    return html.Div(
        summary_items,
        className="content-container"
    )



def create_modality_tabs(dataset, tab):
    if dataset.adata is None:
        return html.Div()  # Return empty div if no adata
    
    modalities = dataset.adata.mod.keys() if isinstance(dataset.adata, mu.MuData) else ["rna"]
    tabs = [
        dbc.Tab(label=mod.upper(), tab_id=mod) for mod in modalities
    ]

    # Wrap the tabs in a Card to match style
    return dbc.Container(
        fluid=True,
        children=[
                dbc.CardHeader(
                    dbc.Tabs(
                        id={"type": "modality-tabs", "index": tab},
                        active_tab=list(modalities)[0],
                        children=tabs,
                        class_name="modality-tabs"
                    ),
                    style={'padding': '0px'}  # remove default padding
                ),
            
        ]
    )

# AnnData layout (scatter and other plots)
def anndata_layout(adata, default_gene_markers, discrete_label_list, prefix):
    return dbc.Container(
        fluid=True,
        children=[
            dbc.Card(
                [
                    html.Div(generate_embedding_plots(adata, prefix), style={'padding': '20px'}),
                    html.Hr(style={'border': '1px solid #ddd'}),
                    html.Div(
                        generate_other_plots(adata, default_gene_markers, discrete_label_list, prefix),
                        className="dbc",
                        style={'padding': '20px'}
                    ),
                    # html.Hr(style={'border': '1px solid #ddd'}),
                ],
                style={'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)'}
            )
        ]
    )

# IGV layout
def igv_layout(session_names, prefix):
    if not session_names:  # If genome tracks not ready yet
        return dbc.Container(
            fluid=True,
            children=[
                dbc.Card(
                    html.Div(
                        html.P("Genome browser is loading...", 
                               style={'textAlign': 'center', 'color': '#868e96', 'padding': '40px'}),
                        style={'padding': '20px'}
                    ),
                    style={'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)', 'marginTop': '20px'}
                )
            ]
        )
    
    return dbc.Container(
        fluid=True,
        children=[
            dbc.Card(
            html.Div(
                gene_browser_layout(prefix, session_names),
                style={'padding': '20px'}
            ),
            style={'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)', 'marginTop': '20px'}
            )
        ]
    )

# Modality selector
def create_modality_selector(dataset, tab):
    if dataset.adata is None:
        return html.Div()  # Return empty div if no adata
    
    return html.Div([
        html.Label("Select modality:"),
        dcc.Dropdown(
            id={"type": "modality-selector", "index": tab},
            options=[
                {"label": mod, "value": mod} for mod in dataset.adata.mod.keys()
            ] if isinstance(dataset.adata, mu.MuData) else [{"label": "rna", "value": "rna"}],
            value=list(dataset.adata.mod.keys())[0] if isinstance(dataset.adata, mu.MuData) else "rna",
            clearable=False
        )
    ], style={"padding": "10px"})

# Entire tab content
def tab_content(dataset, tab):
    return html.Div([
        html.Div(id={"type": "description-layout-div", "index": tab}),
        html.Div([
            create_modality_tabs(dataset, tab),
            html.Div(id={"type": "ann-layout-div", "index": tab})
        ]),
        # Wrap IGV in loading component with delay
        dcc.Loading(
            id=f"loading-igv-{tab}",
            type="circle",
            children=[
                html.Div(id={"type": "igv-layout-div", "index": tab}),
            ],
            style={"marginTop": "20px"},
            delay_show=1000,  # Wait 1 second before showing loading spinner
        )
    ])

