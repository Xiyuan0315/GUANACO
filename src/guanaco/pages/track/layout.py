from dash import dcc, html


def gene_browser_layout(prefix, hosted_genome_dict):
    """
    Generate genome browser layout for a given dataset.
    Arguments:
        prefix: dataset name, e.g. 'Dataset1'
        hosted_genome_dict: list of ATAC track options for this dataset
    """
    return html.Div(
        style={'display': 'flex', 'flexDirection': 'column'},
        children=[
            html.Div(
                style={'display': 'flex', 'flexDirection': 'row'},
                children=[
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
                                value=None,
                                placeholder="Select an IGV session..."
                            ),
                            html.Hr(),
                            dcc.Loading(id=f'{prefix}-igv-container'),
                        ]
                    ),
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
