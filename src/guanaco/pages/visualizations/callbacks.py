"""One callback-registration entry point for every visualization plot."""

from guanaco.pages.matrix.callbacks import matrix_callbacks
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
):
    """Register all available plots for one modality through one lifecycle."""
    has_igv = bool(genome_tracks) and bool(ref_track)
    enabled = resolve_plot_components(
        adata,
        optional_plot_components,
        has_igv=has_igv,
    )

    if adata is not None:
        matrix_enabled = tuple(key for key in enabled if key != "igv")
        matrix_callbacks(
            app,
            adata,
            prefix,
            enabled_components=matrix_enabled,
            embedding_render_backend=embedding_render_backend,
            color_config=color_config,
            gene_annotation_path=gene_annotation_path,
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
