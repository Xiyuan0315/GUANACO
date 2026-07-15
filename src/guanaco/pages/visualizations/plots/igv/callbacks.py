"""Callback registration for the IGV exploratory plot."""

from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate

from guanaco.utils.colors import (
    DEFAULT_DISCRETE_COLORMAP,
    resolve_discrete_palette,
)

DEFAULT_IGV_LOCUS = "chr1:1-10000000"
JASPAR_RELEASE = "JASPAR2024"
_EXPLORATION_WORKSPACE = "dataset-exploration"
_IGV_TAB = "igv-tab"


def _tracks_with_discrete_palette(tracks, palette_name):
    """Copy IGV tracks and apply the shared discrete palette."""
    fallback = [track.get("color") for track in tracks if track.get("color")]
    palette = resolve_discrete_palette(
        palette_name or DEFAULT_DISCRETE_COLORMAP,
        len(tracks),
        default=fallback,
    )
    if not palette:
        return [dict(track) for track in tracks]
    return [
        {**track, "color": palette[index % len(palette)]}
        for index, track in enumerate(tracks)
    ]


def _default_igv_factory(**kwargs):
    import dash_bio as dashbio

    return dashbio.Igv(**kwargs)


def _default_motif_lookup(motif_id):
    from pyjaspar import jaspardb

    return jaspardb(release=JASPAR_RELEASE).fetch_motif_by_id(motif_id)


def register_igv_callbacks(
    app,
    genome_tracks,
    ref_track,
    prefix,
    *,
    discrete_color_prefix=None,
    igv_factory=None,
    motif_lookup=None,
):
    """Register IGV as a modality-scoped exploratory plot."""
    if not genome_tracks or not ref_track:
        return

    igv_factory = igv_factory or _default_igv_factory
    motif_lookup = motif_lookup or _default_motif_lookup
    igv_inputs = [
        Input(f"{prefix}-igv-genome-select", "value"),
        Input(f"{prefix}-visualization-workspace-tabs", "value"),
        Input(f"{prefix}-exploratory-tabs", "value"),
    ]
    if discrete_color_prefix:
        igv_inputs.append(
            Input(
                f"{discrete_color_prefix}-discrete-color-map-dropdown",
                "value",
            )
        )

    @app.callback(Output(f"{prefix}-igv-container", "children"), igv_inputs)
    def render_igv(
        selected_session,
        active_workspace,
        active_tab,
        discrete_color_map=DEFAULT_DISCRETE_COLORMAP,
    ):
        if active_workspace != _EXPLORATION_WORKSPACE or active_tab != _IGV_TAB:
            raise PreventUpdate
        if selected_session is None:
            return html.Div(
                html.P(
                    "Please select an IGV session from the dropdown above to view "
                    "the genome browser.",
                    style={
                        "textAlign": "center",
                        "color": "#868e96",
                        "fontSize": "16px",
                        "padding": "40px",
                        "backgroundColor": "#f8f9fa",
                        "borderRadius": "8px",
                        "margin": "20px 0",
                    },
                )
            )

        tracks = genome_tracks.get(selected_session)
        if tracks is None:
            raise ValueError(f"No tracks configured for IGV session {selected_session}")
        return html.Div(
            igv_factory(
                id=f"igv-{prefix}",
                genome=ref_track["label"],
                locus=DEFAULT_IGV_LOCUS,
                tracks=_tracks_with_discrete_palette(
                    tracks,
                    discrete_color_map,
                ),
            )
        )

    @app.callback(
        Output(f"{prefix}-search-results", "children"),
        Input(f"{prefix}-search-button", "n_clicks"),
        State(f"{prefix}-search-input", "value"),
    )
    def search_motif(n_clicks, search_value):
        if not n_clicks or not search_value:
            return html.Div(
                "Please enter a valid motif ID and click Search.",
                style={"color": "gray"},
            )

        try:
            from .motif import MOTIF_INFO_LABELS, render_motif

            motif_info, image_data = render_motif(motif_lookup(search_value))
            return html.Div(
                [
                    html.Div(
                        "Motif Information",
                        style={
                            "fontWeight": "bold",
                            "marginBottom": "10px",
                        },
                    ),
                    html.Table(
                        [
                            html.Tr(
                                [
                                    html.Th(
                                        f"{key}: ",
                                        style={"backgroundColor": "#f2f2f2"},
                                    ),
                                    html.Td(
                                        value,
                                        style={"backgroundColor": "#f9f9f9"},
                                    ),
                                ]
                            )
                            for key, value in zip(
                                MOTIF_INFO_LABELS,
                                motif_info,
                                strict=True,
                            )
                        ],
                        style={
                            "width": "100%",
                            "border": "2px solid #6699CC",
                            "borderCollapse": "collapse",
                            "textAlign": "left",
                        },
                    ),
                    html.Div(
                        "Sequence Logo",
                        style={
                            "fontWeight": "bold",
                            "marginTop": "20px",
                        },
                    ),
                    html.Img(
                        src=f"data:image/png;base64,{image_data}",
                        style={"maxWidth": "100%"},
                    ),
                ]
            )
        except Exception:
            return html.Div("Incorrect Motif ID", style={"color": "red"})

