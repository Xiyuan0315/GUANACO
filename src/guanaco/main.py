import warnings
from dash import dcc, html, Output, Input, MATCH, State
from guanaco.dash_app import app
from guanaco.layouts import (
    navbar,
    tab_content,
    footprint,
    guanaco_footer,
    description_layout,
    anndata_layout,
    resize_tip_toast,
)
from guanaco.pages.visualizations.callbacks import register_visualization_callbacks
from guanaco.data.loader import get_discrete_labels
from guanaco.data.multiomics import (
    JOINT_TAB_ID,
    expose_global_obs_to_modalities,
    try_build_multiomics_source,
)
from guanaco.data.registry import datasets, embedding_render_backend
from guanaco.utils.gene_extraction_utils import pin_genes
import muon as mu

warnings.filterwarnings(
    "ignore", category=FutureWarning, message=".*__version__.*deprecated.*"
)
mu.set_options(pull_on_update=False)

_label_cache: dict[tuple[int, int, int], list[str]] = {}
_multiomics_sources = {}


def _cached_discrete_labels(adata):
    cache_key = (id(adata), adata.n_obs, adata.n_vars)
    labels = _label_cache.get(cache_key)
    if labels is None:
        labels = get_discrete_labels(adata)
        _label_cache[cache_key] = labels
    return labels


app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="tip-store", storage_type="session", data={"shown": False}),
        navbar(datasets),
        resize_tip_toast(),
        html.Div(id="tabs-content", style={"paddingTop": "70px"}),
        footprint,
        guanaco_footer,
    ]
)


# Helper to look up a per-modality config value, falling back to dataset-level.
def _modality_config(dataset, mod: str) -> dict:
    """Return the effective per-modality config dict for *mod*."""
    mc = (dataset.modality_configs or {}).get(mod, {})
    genome_tracks = mc.get("genome_tracks")
    ref_track = mc.get("ref_track")

    # Backward compatibility: legacy bucket_urls are dataset-level. Attach them
    # to the ATAC modality when present, otherwise to the sole/first modality.
    if genome_tracks is None and dataset.genome_tracks is not None:
        adata = dataset.adata
        if isinstance(adata, mu.MuData):
            modality_names = list(adata.mod.keys())
            track_modality = "atac" if "atac" in modality_names else modality_names[0]
            if mod == track_modality:
                genome_tracks = dataset.genome_tracks
                ref_track = dataset.ref_track
        else:
            genome_tracks = dataset.genome_tracks
            ref_track = dataset.ref_track

    return {
        "gene_markers": mc.get("gene_markers", dataset.gene_markers),
        "scatter_defaults": mc.get("scatter_defaults", dataset.scatter_defaults),
        "optional_plot_components": mc.get(
            "optional_plot_components", dataset.optional_plot_components
        ),
        "gene_annotation_path": mc.get(
            "gene_annotation_path", dataset.gene_annotation_path
        ),
        "genome_tracks": genome_tracks,
        "ref_track": ref_track,
    }


# Register callbacks for scatter and other plots for each dataset
for name, dataset in datasets.items():
    dataset_adata = dataset.adata

    # Register AnnData callbacks if adata exists
    if dataset_adata is not None:
        if isinstance(dataset_adata, mu.MuData):
            modality_names = list(dataset_adata.mod.keys())
            multiomics_source, unavailable_reason = try_build_multiomics_source(
                dataset_adata
            )
            # Shared metadata stays stored once in mdata.obs. The existing
            # single-modality callbacks receive aligned, in-memory views of it.
            expose_global_obs_to_modalities(dataset_adata)
            for mod in modality_names:
                mod_adata = dataset_adata.mod[mod]
                prefix = f"{name}-{mod}"
                mod_cfg = _modality_config(dataset, mod)
                # Pre-load config marker genes into memory for backed RNA data.
                if dataset.backed_mode and mod == "rna" and mod_cfg["gene_markers"]:
                    pin_genes(mod_adata, mod_cfg["gene_markers"])
                register_visualization_callbacks(
                    app,
                    mod_adata,
                    prefix,
                    embedding_render_backend=embedding_render_backend,
                    color_config=dataset.color_config,
                    gene_annotation_path=mod_cfg["gene_annotation_path"],
                    optional_plot_components=mod_cfg["optional_plot_components"],
                    genome_tracks=mod_cfg["genome_tracks"],
                    ref_track=mod_cfg["ref_track"],
                    modality_name=mod,
                    organism=dataset.organism,
                    default_features=mod_cfg["gene_markers"],
                )
            if multiomics_source is not None:
                _multiomics_sources[name] = multiomics_source
                register_visualization_callbacks(
                    app,
                    multiomics_source.base_adata,
                    f"{name}-{JOINT_TAB_ID}",
                    embedding_render_backend=embedding_render_backend,
                    color_config=dataset.color_config,
                    optional_plot_components=None,
                    multiomics_source=multiomics_source,
                    modality_name=JOINT_TAB_ID,
                    organism=dataset.organism,
                )
            elif unavailable_reason:
                print(
                    f"[guanaco] Multi-omics view unavailable for {name}: "
                    f"{unavailable_reason}"
                )
        else:
            prefix = name
            mod_cfg = _modality_config(dataset, "rna")
            # Pre-load config marker genes into memory so the first access is instant.
            if dataset.backed_mode and mod_cfg["gene_markers"]:
                pin_genes(dataset_adata, mod_cfg["gene_markers"])
            register_visualization_callbacks(
                app,
                dataset_adata,
                prefix,
                embedding_render_backend=embedding_render_backend,
                color_config=dataset.color_config,
                gene_annotation_path=mod_cfg["gene_annotation_path"],
                optional_plot_components=mod_cfg["optional_plot_components"],
                genome_tracks=mod_cfg["genome_tracks"],
                ref_track=mod_cfg["ref_track"],
                modality_name="rna",
                organism=dataset.organism,
                default_features=mod_cfg["gene_markers"],
            )
    else:
        # Legacy tracks-only datasets use a synthetic genome modality. Explicit
        # modality track configs keep their own prefixes.
        track_modalities = [
            mod
            for mod, cfg in (dataset.modality_configs or {}).items()
            if cfg.get("genome_tracks")
        ] or ["genome"]
        for mod in track_modalities:
            mod_cfg = _modality_config(dataset, mod)
            prefix = name if mod == "genome" else f"{name}-{mod}"
            register_visualization_callbacks(
                app,
                None,
                prefix,
                optional_plot_components=mod_cfg["optional_plot_components"],
                genome_tracks=mod_cfg["genome_tracks"],
                ref_track=mod_cfg["ref_track"],
                has_palette_control=False,
                modality_name=mod,
                organism=dataset.organism,
            )


@app.callback(Output("tabs-content", "children"), Input("tabs-dataset", "active_tab"))
def update_tab_content(tab):
    dataset = datasets[tab]
    return tab_content(dataset, tab, multiomics_source=_multiomics_sources.get(tab))


# Update description layout
@app.callback(
    Output({"type": "description-layout-div", "index": MATCH}, "children"),
    Input("tabs-dataset", "active_tab"),
)
def update_description_layout(active_tab):
    dataset = datasets[active_tab]

    return description_layout(dataset)


# Update AnnData layout
@app.callback(
    Output({"type": "ann-layout-div", "index": MATCH}, "children"),
    Input({"type": "modality-tabs", "index": MATCH}, "active_tab"),
    Input("tabs-dataset", "active_tab"),
)
def update_anndata_layout(selected_modality, active_tab):
    dataset = datasets[active_tab]
    dataset_adata = dataset.adata
    is_multimodal = isinstance(dataset_adata, mu.MuData)
    multiomics_source = (
        _multiomics_sources.get(active_tab)
        if selected_modality == JOINT_TAB_ID
        else None
    )
    if multiomics_source is not None:
        adata = multiomics_source.base_adata
    else:
        adata = dataset_adata.mod[selected_modality] if is_multimodal else dataset_adata
    label_list = _cached_discrete_labels(adata) if adata is not None else []
    has_explicit_track_modality = (
        dataset_adata is None and selected_modality != "genome"
    )
    prefix = (
        f"{active_tab}-{selected_modality}"
        if is_multimodal or has_explicit_track_modality
        else active_tab
    )

    # Per-modality config (falls back to dataset-level for single-modality / legacy configs).
    mod_cfg = _modality_config(dataset, selected_modality)

    if multiomics_source is not None and multiomics_source.is_paired:
        markers_by_modality = {
            modality: _modality_config(dataset, modality).get("gene_markers")
            for modality in multiomics_source.modalities
        }
        modality_markers = multiomics_source.default_features(markers_by_modality)
        label_list = multiomics_source.discrete_obs_names
        left_embedding = multiomics_source.preferred_embedding(
            multiomics_source.modalities[0]
        )
        right_embedding = multiomics_source.preferred_embedding(
            multiomics_source.modalities[1]
        )
        default_color = label_list[0] if label_list else None
        mod_cfg = {
            **mod_cfg,
            "scatter_defaults": {
                "embedding_left": left_embedding,
                "embedding_right": right_embedding,
                "color_left": default_color,
                "color_right": default_color,
            },
            "optional_plot_components": None,
            "gene_annotation_path": None,
            "genome_tracks": None,
            "ref_track": None,
        }
    elif multiomics_source is not None:
        embedding_modalities = [
            modality
            for modality in multiomics_source.modalities
            if multiomics_source.modality_embeddings(modality)
        ]
        left_modality, right_modality = embedding_modalities[:2]
        left_cfg = _modality_config(dataset, left_modality)
        right_cfg = _modality_config(dataset, right_modality)

        def resolve_panel_color(modality, configured):
            if not configured:
                return None
            modality_adata = multiomics_source.modality_adata(modality)
            if configured in modality_adata.obs.columns:
                return configured
            return multiomics_source.resolve_text_feature(
                f"{modality}::{configured}"
            )

        modality_markers = []
        label_list = []
        mod_cfg = {
            **mod_cfg,
            "scatter_defaults": {
                "embedding_left": multiomics_source.preferred_embedding(
                    left_modality
                ),
                "embedding_right": multiomics_source.preferred_embedding(
                    right_modality
                ),
                "color_left": resolve_panel_color(
                    left_modality,
                    (left_cfg.get("scatter_defaults") or {}).get("color_left"),
                ),
                "color_right": resolve_panel_color(
                    right_modality,
                    (right_cfg.get("scatter_defaults") or {}).get("color_right"),
                ),
            },
            "optional_plot_components": None,
            "gene_annotation_path": None,
            "genome_tracks": None,
            "ref_track": None,
        }
    elif mod_cfg["gene_markers"] is not None:
        modality_markers = mod_cfg["gene_markers"]
    else:
        modality_markers = adata.var_names[:6].tolist() if adata is not None else []

    return anndata_layout(
        adata,
        modality_markers,
        label_list,
        prefix,
        optional_plot_components=mod_cfg["optional_plot_components"],
        scatter_defaults=mod_cfg["scatter_defaults"],
        gene_annotation_path=mod_cfg["gene_annotation_path"],
        genome_tracks=mod_cfg["genome_tracks"],
        ref_track=mod_cfg["ref_track"],
        multiomics_source=multiomics_source,
        modality_name=selected_modality,
        organism=dataset.organism,
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


server = app.server

if __name__ == "__main__":
    app.run(host="127.0.0.1", debug=True, port=4399)
