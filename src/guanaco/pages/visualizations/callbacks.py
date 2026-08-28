"""One callback-registration entry point for every visualization plot."""

from guanaco.pages.matrix.callbacks import matrix_callbacks
from guanaco.pages.matrix.callbacks.unpaired_multiomics_callbacks import (
    register_unpaired_multiomics_callbacks,
)
from guanaco.pages.matrix.callbacks.cross_modal_concordance_callbacks import (
    register_cross_modal_concordance_callbacks,
)
from guanaco.pages.matrix.callbacks.multiomics_composition_callbacks import (
    register_multiomics_composition_callbacks,
)
from guanaco.pages.visualizations.registry import multiomics_plot_components
from guanaco.pages.visualizations.registry import resolve_plot_components


def register_visualization_callbacks(
    app,
    adata,
    prefix,
    *,
    embedding_render_backend="scattergl",
    color_config=None,
    gene_annotation_path=None,
    optional_plot_components=None,
    genome_tracks=None,
    ref_track=None,
    has_palette_control=None,
    multiomics_source=None,
    modality_name=None,
    organism="human",
):
    """Register all available plots for one modality through one lifecycle."""
    has_igv = bool(genome_tracks) and bool(ref_track)
    enabled = resolve_plot_components(
        adata,
        optional_plot_components,
        has_igv=has_igv,
        modality_name=modality_name,
        feature_data_available=bool(multiomics_source.feature_names)
        if multiomics_source is not None
        else None,
        discrete_data_available=bool(multiomics_source.discrete_obs_names)
        if multiomics_source is not None
        else None,
    )
    if multiomics_source is not None:
        enabled = multiomics_plot_components(enabled)

    if multiomics_source is not None and "multiomics-composition" in enabled:
        register_multiomics_composition_callbacks(
            app,
            multiomics_source,
            prefix,
            color_config=color_config,
        )

    if multiomics_source is not None and not multiomics_source.is_paired:
        if multiomics_source.supports_embedding_view:
            register_unpaired_multiomics_callbacks(
                app,
                multiomics_source,
                prefix,
                embedding_render_backend=embedding_render_backend,
                color_config=color_config,
            )
        register_cross_modal_concordance_callbacks(
            app,
            multiomics_source,
            prefix,
        )
        return enabled

    if adata is not None:
        matrix_enabled = tuple(
            key
            for key in enabled
            if key not in {"igv", "multiomics-composition"}
        )
        matrix_callbacks(
            app,
            adata,
            prefix,
            enabled_components=matrix_enabled,
            embedding_render_backend=embedding_render_backend,
            color_config=color_config,
            gene_annotation_path=gene_annotation_path,
            multiomics_source=multiomics_source,
            organism=organism,
        )

    if "igv" in enabled:
        # Keep optional Dash Bio / JASPAR dependencies out of matrix-only startup.
        from guanaco.pages.visualizations.plots.igv.callbacks import (
            register_igv_callbacks,
        )

        palette_control = (
            adata is not None
            if has_palette_control is None
            else has_palette_control
        )
        register_igv_callbacks(
            app,
            genome_tracks,
            ref_track,
            prefix,
            discrete_color_prefix=prefix if palette_control else None,
        )

    return enabled
