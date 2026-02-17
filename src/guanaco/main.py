import warnings
warnings.filterwarnings("ignore", category=FutureWarning, message=".*__version__.*deprecated.*")

from dash import dcc, html, Output, Input, MATCH, State
from guanaco.app import app
from guanaco.layout import (
    navbar, tab_content, footprint, guanaco_footer, description_layout,
    anndata_layout, igv_layout, resize_tip_toast
)
from guanaco.pages.track.gene_browser import gene_browser_callbacks
from guanaco.pages.matrix.callbacks import matrix_callbacks
from guanaco.data_loader import datasets, get_discrete_labels
import muon as mu

mu.set_options(pull_on_update=False)

_label_cache: dict[tuple[int, int, int], list[str]] = {}


def _cached_discrete_labels(adata):
    cache_key = (id(adata), adata.n_obs, adata.n_vars)
    labels = _label_cache.get(cache_key)
    if labels is None:
        labels = get_discrete_labels(adata)
        _label_cache[cache_key] = labels
    return labels


app.layout = html.Div([
    dcc.Location(id="url", refresh=False),
    dcc.Store(id="tip-store", storage_type="session", data={"shown": False}),
    navbar(datasets),
    resize_tip_toast(),
    html.Div(id="tabs-content", style={"paddingTop": "70px"}),
    footprint,
    guanaco_footer,
])

# Register callbacks for scatter and other plots for each dataset
for name, dataset in datasets.items():
    dataset_adata = dataset.adata

    # Register AnnData callbacks if adata exists
    if dataset_adata is not None:
        if isinstance(dataset_adata, mu.MuData):
            for mod in dataset_adata.mod.keys():
                mod_adata = dataset_adata.mod[mod]
                prefix = f"{name}-{mod}"
                matrix_callbacks(
                    app,
                    mod_adata,
                    prefix,
                    embedding_render_backend=dataset.embedding_render_backend,
                )
        else:
            prefix = name
            matrix_callbacks(
                app,
                dataset_adata,
                prefix,
                embedding_render_backend=dataset.embedding_render_backend,
            )

    if dataset.genome_tracks is not None and dataset.ref_track is not None:
        gene_browser_callbacks(app, dataset.genome_tracks, dataset.ref_track, dataset.title)


@app.callback(
    Output("tabs-content", "children"),
    Input("tabs-dataset", "active_tab")
)
def update_tab_content(tab):
    dataset = datasets[tab]
    return tab_content(dataset, tab)

# Update description layout
@app.callback(
    Output({"type": "description-layout-div", "index": MATCH}, "children"),
    Input("tabs-dataset", "active_tab")
)
def update_description_layout(active_tab):
    dataset = datasets[active_tab]

    return description_layout(dataset)


# Update AnnData layout
@app.callback(
    Output({"type": "ann-layout-div", "index": MATCH}, "children"),
    Input({"type": "modality-tabs", "index": MATCH}, "active_tab"),
    Input("tabs-dataset", "active_tab")
)
def update_anndata_layout(selected_modality, active_tab):
    dataset = datasets[active_tab]
    dataset_adata = dataset.adata
    if dataset_adata is None:
        return html.Div("No AnnData available for this dataset", style={"padding": "20px"})
    
    is_multimodal = isinstance(dataset_adata, mu.MuData)
    adata = dataset_adata.mod[selected_modality] if is_multimodal else dataset_adata
    label_list = _cached_discrete_labels(adata)
    prefix = f"{active_tab}-{selected_modality}" if is_multimodal else active_tab
    
    if selected_modality == 'rna' and dataset.gene_markers is not None:
        modality_markers = dataset.gene_markers
    else:
        modality_markers = adata.var_names[:6].tolist() if adata else []
    
    return anndata_layout(
        adata,
        modality_markers,
        label_list,
        prefix,
        optional_plot_components=dataset.optional_plot_components,
    )

@app.callback(
    [Output("tip-modal", "is_open"), Output("tip-store", "data")],
    [Input("url", "pathname"), Input("close-tip", "n_clicks")],
    [State("tip-modal", "is_open"), State("tip-store", "data")],
)
def toggle_tip(pathname, n_clicks, is_open, store):
    if store is None:
        store = {"shown": False}
    
    if n_clicks:
        return False, {"shown": True}
    
    if not store.get("shown", False):
        return True, {"shown": True}

    return False, store

@app.callback(
    Output({"type": "igv-layout-div", "index": MATCH}, "children"),
    Input("tabs-dataset", "active_tab")
)
def update_igv_layout(active_tab):
    dataset = datasets.get(active_tab)
    if dataset is None:
        return html.Div("Dataset not found.", style={"padding": "20px"})

    genome_tracks = dataset.genome_tracks or {}
    session_names = list(genome_tracks.keys())
    return igv_layout(session_names, prefix=active_tab)

server = app.server

if __name__ == "__main__":
    app.run_server(host='0.0.0.0', debug=True, port=4399)
