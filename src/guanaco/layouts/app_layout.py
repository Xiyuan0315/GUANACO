import dash_bootstrap_components as dbc
from dash import html
from guanaco.pages.matrix.layouts.embedding_layout import generate_embedding_plots
from guanaco.pages.visualizations import generate_visualization_sections
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
                    alt="Tip",
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
    role="presentation",
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
        "color": "#6c757d"
    })
)

# Navbar layout
def navbar(datasets):
    return html.Div(
        dbc.Navbar(
            dbc.Container(
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/assets/logo.png", alt="GUANACO logo", height="70px"), width="auto", align="center"),
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
    modality_track_sets = [
        cfg.get("genome_tracks")
        for cfg in (dataset.modality_configs or {}).values()
        if cfg.get("genome_tracks")
    ]
    configured_track_sets = modality_track_sets or (
        [dataset.genome_tracks] if dataset.genome_tracks else []
    )
    if configured_track_sets:
        track_count = sum(
            len(tracks)
            for genome_tracks in configured_track_sets
            for tracks in genome_tracks.values()
        )
        summary_items.append(
            html.P(
                [html.B("Genome Browser Tracks: "), f"{track_count}"],
                className="summary-item",
            )
        )


    # Final layout
    return html.Div(
        summary_items,
        className="content-container"
    )



def create_modality_tabs(dataset, tab, multiomics_source=None):
    if dataset.adata is None:
        track_modalities = [
            name
            for name, cfg in (dataset.modality_configs or {}).items()
            if cfg.get("genome_tracks")
        ]
        if track_modalities:
            modalities = track_modalities
        elif dataset.genome_tracks:
            modalities = ["genome"]
        else:
            return html.Div()
    else:
        modalities = (
            dataset.adata.mod.keys()
            if isinstance(dataset.adata, mu.MuData)
            else ["rna"]
        )
    modalities = list(modalities)
    tabs = [dbc.Tab(label=mod.upper(), tab_id=mod) for mod in modalities]
    if multiomics_source is not None:
        tabs.append(
            dbc.Tab(label=multiomics_source.label, tab_id="multiomics")
        )

    # Wrap the tabs in a Card to match style
    return dbc.Container(
        fluid=True,
        children=[
            dbc.CardHeader(
                dbc.Tabs(
                    id={"type": "modality-tabs", "index": tab},
                    active_tab=list(modalities)[0],
                    children=tabs,
                    class_name="modality-tabs",
                ),
                style={"padding": "0px"},  # remove default padding
            ),
        ],
    )


# AnnData layout (scatter and other plots)
def anndata_layout(
    adata,
    default_gene_markers,
    discrete_label_list,
    prefix,
    optional_plot_components=None,
    scatter_defaults=None,
    gene_annotation_path=None,
    genome_tracks=None,
    ref_track=None,
    multiomics_source=None,
    modality_name=None,
    organism="human",
):
    is_unpaired_multiomics = (
        multiomics_source is not None and not multiomics_source.is_paired
    )
    sections = generate_visualization_sections(
        adata,
        default_gene_markers,
        discrete_label_list,
        prefix,
        optional_plot_components=optional_plot_components,
        gene_annotation_path=gene_annotation_path,
        genome_tracks=genome_tracks,
        ref_track=ref_track,
        multiomics_source=multiomics_source,
        modality_name=modality_name,
        organism=organism,
    )
    children = []
    if adata is not None or is_unpaired_multiomics:
        children.append(
            dbc.Card(
                html.Div(
                    generate_embedding_plots(
                        adata,
                        prefix,
                        scatter_defaults=scatter_defaults,
                        multiomics_source=multiomics_source,
                    ),
                    className="plot-section",
                ),
                className="card-elevated",
            )
        )
    children.extend(sections)
    return dbc.Container(
        fluid=True,
        children=children,
    )


# Entire tab content
def tab_content(dataset, tab, multiomics_source=None):
    return html.Div(
        [
            html.Div(id={"type": "description-layout-div", "index": tab}),
            html.Div(
                [
                    create_modality_tabs(
                        dataset, tab, multiomics_source=multiomics_source
                    ),
                    html.Div(id={"type": "ann-layout-div", "index": tab}),
                ]
            ),
        ]
    )
