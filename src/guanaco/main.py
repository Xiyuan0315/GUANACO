import warnings
warnings.filterwarnings("ignore", category=FutureWarning, message=".*__version__.*deprecated.*")

from dash import dcc, html, Output, Input, MATCH, State
from guanaco.dash_app import app
from guanaco.layouts import (
    navbar, tab_content, footprint, guanaco_footer, description_layout,
    anndata_layout, igv_layout, resize_tip_toast
)
from guanaco.pages.matrix.callbacks import matrix_callbacks
from guanaco.data.loader import get_discrete_labels
from guanaco.data.registry import datasets, embedding_render_backend
from guanaco.utils.gene_extraction_utils import pin_genes
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

# Helper to look up a per-modality config value, falling back to dataset-level.
def _modality_config(dataset, mod: str) -> dict:
    """Return the effective per-modality config dict for *mod*."""
    mc = (dataset.modality_configs or {}).get(mod, {})
    return {
        "gene_markers": mc.get("gene_markers", dataset.gene_markers),
        "scatter_defaults": mc.get("scatter_defaults", dataset.scatter_defaults),
        "optional_plot_components": mc.get("optional_plot_components", dataset.optional_plot_components),
        "gene_annotation_path": mc.get("gene_annotation_path", dataset.gene_annotation_path),
    }

# Register callbacks for scatter and other plots for each dataset
for name, dataset in datasets.items():
    dataset_adata = dataset.adata

    # Register AnnData callbacks if adata exists
    if dataset_adata is not None:
        if isinstance(dataset_adata, mu.MuData):
            for mod in dataset_adata.mod.keys():
                mod_adata = dataset_adata.mod[mod]
                prefix = f"{name}-{mod}"
                mod_cfg = _modality_config(dataset, mod)
                # Pre-load config marker genes into memory for backed RNA data.
                if dataset.backed_mode and mod == "rna" and mod_cfg["gene_markers"]:
                    pin_genes(mod_adata, mod_cfg["gene_markers"])
                matrix_callbacks(
                    app,
                    mod_adata,
                    prefix,
                    embedding_render_backend=embedding_render_backend,
                    color_config=dataset.color_config,
                    gene_annotation_path=mod_cfg["gene_annotation_path"],
                )
        else:
            prefix = name
            mod_cfg = _modality_config(dataset, "rna")
            # Pre-load config marker genes into memory so the first access is instant.
            if dataset.backed_mode and mod_cfg["gene_markers"]:
                pin_genes(dataset_adata, mod_cfg["gene_markers"])
            matrix_callbacks(
                app,
                dataset_adata,
                prefix,
                embedding_render_backend=embedding_render_backend,
                color_config=dataset.color_config,
                gene_annotation_path=mod_cfg["gene_annotation_path"],
            )

    if dataset.genome_tracks is not None and dataset.ref_track is not None:
        # Imported lazily so the genome-browser deps (pyjaspar, logomaker, dash-bio)
        # are only required when a dataset actually configures genome tracks --
        # installed via the `guanaco-viz[tracks]` extra.
        from guanaco.pages.track.callbacks import gene_browser_callbacks
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

    # Per-modality config (falls back to dataset-level for single-modality / legacy configs).
    mod_cfg = _modality_config(dataset, selected_modality)

    if mod_cfg["gene_markers"] is not None:
        modality_markers = mod_cfg["gene_markers"]
    else:
        modality_markers = adata.var_names[:6].tolist() if adata else []

    return anndata_layout(
        adata,
        modality_markers,
        label_list,
        prefix,
        optional_plot_components=mod_cfg["optional_plot_components"],
        scatter_defaults=mod_cfg["scatter_defaults"],
        gene_annotation_path=mod_cfg["gene_annotation_path"],
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
    app.run_server(host = '127.0.0.1', debug=True, port=4399)
