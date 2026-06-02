from pyjaspar import jaspardb
import dash_bio as dashbio
from dash import dcc, html, Input, Output, State

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
    if genome_tracks is None or ref_track is None:
        return

    @app.callback(
        Output(f'{prefix}-igv-container', 'children'),
        Input(f'{prefix}-igv-genome-select', 'value')
    )
    def return_igv(selected_atac):
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
                locus='chr1:1-10000000',
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
            return html.Div("Incorrect Motif ID", style={"color": "red"})
