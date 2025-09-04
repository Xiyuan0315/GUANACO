from pyjaspar import jaspardb
import dash_bio as dashbio
from dash import dcc, html,Input, Output, State
from guanaco.pages.track.utils import plot_motif

def gene_browser_callbacks(app, genome_tracks, ref_track, prefix):
    """
    Register genome browser callbacks for a specific dataset.
    Arguments:
        app: Dash app
        genome_tracks: dict of genome tracks (from DatasetBundle)
        ref_track: reference genome (from DatasetBundle)
        prefix: dataset name (used to ensure unique IDs)
    """
    
    # Safety check - should not happen with updated main.py, but good to have
    if genome_tracks is None or ref_track is None:
        return

    @app.callback(
        Output(f'{prefix}-igv-container', 'children'),
        Input(f'{prefix}-igv-genome-select', 'value')
    )
    def return_igv(selected_atac):
        # If no session is selected, show a prompt message
        if selected_atac is None:
            return html.Div(
                [
                    html.P(
                        "Please select an IGV session from the dropdown above to view the genome browser.",
                        style={
                            'textAlign': 'center',
                            'color': '#868e96',
                            'fontSize': '16px',
                            'padding': '40px',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '8px',
                            'margin': '20px 0'
                        }
                    )
                ]
            )
        
        if genome_tracks.get(selected_atac) is None:
            raise Exception(f"No tracks configured for ATAC {selected_atac}")

        return html.Div([
            dashbio.Igv(
                id=f'igv-{prefix}',
                genome=ref_track['label'],
                locus='chr1:1-10000000',  # default or dataset-specific
                tracks=genome_tracks[selected_atac]
            )
        ])

    @app.callback(
        Output(f'{prefix}-search-results', 'children'),
        Input(f'{prefix}-search-button', 'n_clicks'),
        State(f'{prefix}-search-input', 'value')
    )
    def handle_search(n_clicks, search_value):
        if n_clicks is None or n_clicks == 0 or not search_value:
            return html.Div("Please enter a valid motif ID and click Search.", style={"color": "gray"})

        try:
            # Fetch motif data from JASPAR
            jdb_obj = jaspardb(release='JASPAR2024')
            motif = jdb_obj.fetch_motif_by_id(search_value)
            motif_info, img_data = plot_motif(motif)

            return html.Div(
                children=[
                    html.Div("Motif Information", style={"fontWeight": "bold", "marginBottom": "10px"}),
                    html.Table(
                        [
                            *[
                                html.Tr(
                                    [html.Th(f'{key}: ', style={"backgroundColor": "#f2f2f2"}),
                                     html.Td(value, style={"backgroundColor": "#f9f9f9"})]
                                )
                                for key, value in zip(
                                    ["TF Name", "Matrix ID", "Collection", "TF Class", "TF Family", "Data Type", "Medline"],
                                    motif_info
                                )
                            ],
                        ],
                        style={"width": "100%", "border": "2px solid #6699CC", "borderCollapse": "collapse", "textAlign": "left"}
                    ),
                    html.Div("Sequence Logo", style={"fontWeight": "bold", "marginTop": "20px"}),
                    html.Img(src=f"data:image/png;base64,{img_data}", style={"maxWidth": "100%"})
                ]
            )
        except Exception:
            return html.Div(f"Incorrect Motif ID", style={"color": "red"})


def gene_browser_layout(prefix, hosted_genome_dict):
    """
    Generate genome browser layout for a given dataset (prefix).
    Arguments:
        prefix: dataset name, e.g. 'Dataset1'
        hosted_genome_dict: list of ATAC track options for this dataset (like ['Overall', 'Splitted'])
    """
    return html.Div(
        style={'display': 'flex', 'flexDirection': 'column'},
        children=[
            html.Div(
                style={'display': 'flex', 'flexDirection': 'row'},
                children=[
                    # Left column: IGV container
                    html.Div(
                        style={'flex': '2', 'padding': '10px'},
                        children=[
                            html.P('Select the IGV session to display below:'),
                            dcc.Dropdown(
                                id=f'{prefix}-igv-genome-select',
                                options=[
                                    {'label': s, 'value': s}
                                    for s in hosted_genome_dict
                                ],
                                value=None,  # No default selection
                                placeholder="Select an IGV session..."
                            ),
                            html.Hr(),
                            dcc.Loading(id=f'{prefix}-igv-container'),
                        ]
                    ),
                    # Right column: Search box
                    html.Div(
                        style={
                            'flex': '1',
                            'padding': '10px',
                            'borderLeft': '1px solid #ccc'
                        },
                        children=[
                            html.H4('Motif Search Box'),
                            html.P('Search with motif id, from JASPAR database'),
                            html.Div(
                                style={
                                    'display': 'flex',
                                    'alignItems': 'center',
                                    'gap': '10px'
                                },
                                children=[
                                    dcc.Input(
                                        id=f'{prefix}-search-input',
                                        type='text',
                                        placeholder='Enter a motif id (e.g., MA1972.1)',
                                        style={'flex': '1', 'marginBottom': '10px'}
                                    ),
                                    html.Button(
                                        'Search',
                                        id=f'{prefix}-search-button',
                                        n_clicks=0,
                                        style={
                                            'padding': '10px 20px',
                                            'whiteSpace': 'nowrap'
                                        }
                                    )
                                ]
                            ),
                            html.Div(
                                id=f'{prefix}-search-results',
                                style={'marginTop': '20px'},
                                children=[]
                            )
                        ]
                    )
                ]
            ),
        ]
    )
