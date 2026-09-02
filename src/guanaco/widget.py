"""Scanpy-style plotting API for GUANACO, for interactive use in notebooks.

Designed to feel familiar to scanpy users::

    import guanaco as gc

    gc.pl.umap(adata, color="leiden")
    gc.pl.embedding(adata, basis="X_pca", color="CD8A")
    gc.pl.violin(adata, keys=["CD8A", "CD4"], groupby="cell_type")
    gc.pl.heatmap(adata, var_names=markers, groupby="cell_type")
    gc.pl.dotplot(adata, var_names=markers, groupby="cell_type")
    gc.pl.stacked_bar(adata, x="sample", color="cell_type")
    gc.pl.pseudotime(adata, genes=["CD8A"], pseudotime_key="dpt_pseudotime")
    gc.pl.volcano(adata, group="group1")

Most functions return an interactive Plotly figure. Following scanpy's
convention, each takes ``show`` (default ``True``: render the figure inline) and
``return_fig`` (default ``False``: return the ``plotly.graph_objects.Figure``
instead of showing it).

PAGA uses Cytoscape, so ``gc.pl.paga`` returns an ``ipycytoscape`` widget
instead of a Plotly figure (use ``return_widget`` rather than ``return_fig``).
It needs ``pip install ipycytoscape`` and renders in
Jupyter/JupyterLab; classic ipywidgets are not currently supported in marimo.
"""

from __future__ import annotations

from typing import Literal, Mapping, Sequence

import plotly.graph_objects as go

from guanaco.utils.plot_style import GUANACO_QUALITATIVE, categorical_color_map


def _obs_col(obs, col):
    s = obs[col]
    return s.to_series() if hasattr(s, "to_series") else s


# --------------------------------------------------------------------------- #
# Internal helpers
# --------------------------------------------------------------------------- #
def _render(
    fig: go.Figure,
    show: bool,
    return_fig: bool,
    *,
    title: str | None = None,
    height: str | int | None = None,
):
    """Apply notebook presentation, then return or display the figure."""
    from guanaco.utils.plot_style import apply_guanaco_figure_style

    if isinstance(height, int):
        fig.update_layout(height=height)
    elif isinstance(height, str) and height.endswith("px"):
        try:
            fig.update_layout(height=int(height[:-2]))
        except ValueError:
            pass
    apply_guanaco_figure_style(fig, title=title)
    if return_fig:
        return fig
    if show:
        fig.show()
        return None
    return fig


def _show_or_return(obj, *, show: bool, return_widget: bool):
    """Notebook output convention for cytoscape renders (widget / HTML objects).

    The analogue of :func:`_render` for things that aren't Plotly figures: return
    the object, or display it inline via IPython (falling back to returning it
    when no display is available).
    """
    if return_widget:
        return obj
    if show:
        try:
            from IPython.display import display

            display(obj)
            return None
        except Exception:
            return obj
    return obj


def _as_list(value) -> list:
    """Accept a single name or a sequence; always return a list."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return list(value)


def _resolve_basis(adata, basis: str) -> str:
    """Map a scanpy-style basis ('umap', 'X_umap', 'pca'...) to an obsm key."""
    if basis in adata.obsm:
        return basis
    if f"X_{basis}" in adata.obsm:
        return f"X_{basis}"
    return basis  # let the plot function raise a clear error if it's invalid


def _palette_list(palette, n_colors: int):
    """Resolve a palette through GUANACO's shared categorical color contract."""
    return list(categorical_color_map(range(n_colors), palette).values())


def _label_color_map(adata, groupby: str, palette):
    """Build a stable category-to-color map for ``groupby``."""
    labels = sorted(_obs_col(adata.obs, groupby).unique().tolist())
    return categorical_color_map(labels, palette)


def _all_labels(adata, groupby: str) -> list:
    return sorted(_obs_col(adata.obs, groupby).unique().tolist())


_VIEW_CONTROL = frozenset({"adata", "id", "data", "title", "show", "return_fig"})


def _linked_spec(plot: str, values: Mapping):
    """Turn a plotting call into a ViewSpec when the user supplied ``id=``."""

    view_id = values.get("id")
    source = values.get("data")
    adata = values.get("adata")
    if view_id is None:
        if source is not None:
            raise TypeError("`data=` is only used with the linked-view `id=` form.")
        return None
    if adata is not None:
        raise TypeError(
            "Pass data objects to `linked_view()`; linked plot `data=` selects a "
            "named source."
        )
    from guanaco.linking.model import view as build_view

    normalized = {
        name: value for name, value in values.items() if name not in _VIEW_CONTROL
    }
    normalized.update(normalized.pop("kwargs", {}))
    if normalized.get("height") is None:
        normalized.pop("height", None)
    return build_view(
        plot,
        id=view_id,
        data=source,
        title=values.get("title"),
        **normalized,
    )


# --------------------------------------------------------------------------- #
# Embeddings / scatter
# --------------------------------------------------------------------------- #
def embedding(
    adata=None,
    basis: str | None = None,
    color: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    color_mode: str = "auto",
    transformation: str | None = None,
    layer: str | None = None,
    order: str | None = None,
    color_map: str = "Viridis",
    palette=None,
    size: int = 5,
    opacity: float = 1.0,
    legend_loc: str = "right margin",
    axis_show: bool = True,
    render_backend: str = "scattergl",
    img_key: str | None = None,
    library_id: str | None = None,
    img_alpha: float = 1.0,
    color_range: Sequence[float] | None = None,
    color_midpoint: float | None = None,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Scatter of an embedding colored by a gene or an ``obs`` column.

    Parameters mirror scanpy's ``sc.pl.embedding``: ``basis`` is the embedding
    (e.g. ``"X_umap"`` or ``"umap"``), ``color`` a gene or ``obs`` column,
    ``color_map`` the continuous colormap and ``palette`` the categorical one
    (palette name or list of colors). Continuous vs. categorical is detected
    automatically.

    ``render_backend="datashader"`` rasterizes the points server-side (fast for
    very large datasets); it needs the ``datashader`` package and falls back to
    ``"scattergl"`` if it is missing. Note: a datashader raster has no per-point
    hover/lasso, so use ``"scattergl"`` when you need box/lasso selection.
    """
    if (linked := _linked_spec("embedding", locals())) is not None:
        return linked
    if adata is None or basis is None or color is None:
        raise TypeError(
            "`adata`, `basis`, and `color` are required for a standalone plot."
        )
    from guanaco.pages.matrix.plots.embedding import plot_embedding

    n_cat = adata.obs[color].nunique() if color in adata.obs.columns else 50
    fig = plot_embedding(
        adata=adata,
        adata_full=adata,
        source_adata=adata,
        embedding_key=_resolve_basis(adata, basis),
        color=color,
        mode=color_mode,
        transformation=transformation,
        layer=layer,
        order=order,
        continuous_color_map=color_map,
        discrete_color_map=_palette_list(palette, n_cat),
        marker_size=size,
        opacity=opacity,
        legend_show="on data" if legend_loc == "on data" else "on legend",
        axis_show=axis_show,
        render_backend=render_backend,
        img_key=img_key,
        library_id=library_id,
        img_alpha=img_alpha,
    )
    if color_range is not None or color_midpoint is not None:
        for trace in fig.data:
            marker = getattr(trace, "marker", None)
            if marker is None or getattr(marker, "colorscale", None) is None:
                continue
            if color_range is not None:
                marker.cmin, marker.cmax = color_range
            if color_midpoint is not None:
                marker.cmid = color_midpoint
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


def _named_embedding(plot, basis, adata, color, id, data, kwargs):
    """Implement the three named embedding aliases once."""

    values = {"adata": adata, "color": color, "id": id, "data": data, **kwargs}
    if (linked := _linked_spec(plot, values)) is not None:
        return linked
    if adata is None or color is None:
        raise TypeError(f"`adata` and `color` are required for a standalone {plot}.")
    return embedding(adata, basis=basis, color=color, **kwargs)


def umap(adata=None, color: str | None = None, *, id=None, data=None, **kwargs):
    """``embedding`` with ``basis="X_umap"`` (scanpy ``sc.pl.umap``)."""
    return _named_embedding("umap", "X_umap", adata, color, id, data, kwargs)


def pca(adata=None, color: str | None = None, *, id=None, data=None, **kwargs):
    """``embedding`` with ``basis="X_pca"`` (scanpy ``sc.pl.pca``)."""
    return _named_embedding("pca", "X_pca", adata, color, id, data, kwargs)


def tsne(adata=None, color: str | None = None, *, id=None, data=None, **kwargs):
    """``embedding`` with ``basis="X_tsne"`` (scanpy ``sc.pl.tsne``)."""
    return _named_embedding("tsne", "X_tsne", adata, color, id, data, kwargs)


def spatial(adata=None, color: str | None = None, *, id=None, data=None, **kwargs):
    """``embedding`` with ``basis="spatial"`` and AnnData tissue-image support."""
    return _named_embedding("spatial", "spatial", adata, color, id, data, kwargs)


def coexpression(
    adata=None,
    basis: str | None = None,
    gene1: str | None = None,
    gene2: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    threshold1: float = 0.5,
    threshold2: float = 0.5,
    transformation: str | None = None,
    layer: str | None = None,
    size: int = 5,
    opacity: float = 1.0,
    legend_loc: str = "right margin",
    axis_show: bool = True,
    img_key: str | None = None,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Two-gene co-expression on an embedding (cells binned into 4 groups)."""
    if (linked := _linked_spec("coexpression", locals())) is not None:
        return linked
    if adata is None or basis is None or gene1 is None or gene2 is None:
        raise TypeError(
            "`adata`, `basis`, `gene1`, and `gene2` are required for a standalone plot."
        )
    from guanaco.pages.matrix.plots.embedding import plot_coexpression_embedding

    fig = plot_coexpression_embedding(
        adata=adata,
        source_adata=adata,
        embedding_key=_resolve_basis(adata, basis),
        gene1=gene1,
        gene2=gene2,
        threshold1=threshold1,
        threshold2=threshold2,
        transformation=transformation,
        layer=layer,
        marker_size=size,
        opacity=opacity,
        legend_show="on data" if legend_loc == "on data" else "right",
        axis_show=axis_show,
        img_key=img_key,
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


# --------------------------------------------------------------------------- #
# Violin (stacked: genes as rows, grouped by a category)
# --------------------------------------------------------------------------- #
def violin(
    adata=None,
    keys=None,
    groupby: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    labels: Sequence | None = None,
    transformation: str | None = None,
    layer: str | None = None,
    show_box: bool = False,
    bandwidth: float | None = None,
    palette=None,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Stacked violins of ``keys`` (genes) split by ``groupby`` (scanpy ``sc.pl.violin``).

    ``labels`` restricts/orders the ``groupby`` categories (default: all).
    ``show_box`` overlays a neutral Median + IQR box.
    """
    if (linked := _linked_spec("violin", locals())) is not None:
        return linked
    if adata is None or keys is None or groupby is None:
        raise TypeError(
            "`adata`, `keys`, and `groupby` are required for a standalone violin."
        )
    from guanaco.pages.matrix.plots.violin1 import plot_violin1

    used_labels = list(labels) if labels is not None else _all_labels(adata, groupby)
    fig = plot_violin1(
        adata=adata,
        genes=_as_list(keys),
        groupby=groupby,
        labels=used_labels,
        transformation=transformation,
        layer=layer,
        show_box=show_box,
        groupby_label_color_map=_label_color_map(adata, groupby, palette),
        adata_obs=adata.obs,
        bandwidth=bandwidth,
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


def stacked_violin(adata=None, var_names=None, groupby: str | None = None, **kwargs):
    """Alias of :func:`violin` (scanpy ``sc.pl.stacked_violin``)."""
    return violin(adata, var_names, groupby, **kwargs)


def violin_grouped(
    adata=None,
    key: str | None = None,
    groupby: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    groupby2: str | None = None,
    mode: str = "mode1",
    test_method: str = "none",
    transformation: str | None = None,
    layer: str | None = None,
    show_box: bool = True,
    show_points: bool = False,
    bandwidth: float | None = None,
    labels: Sequence | None = None,
    palette=None,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Grouped violin of a single gene across one/two metadata, with optional stats.

    This is GUANACO's "Violin Plot" comparison view (``plot_violin2_new``); it
    mirrors the web app when ``mode`` is passed.

    ``mode``:
      - ``"mode1"`` -- one metadata (``groupby`` only)
      - ``"mode2"`` -- facet by ``groupby``, compare ``groupby2``
      - ``"mode3"`` -- linear model (``groupby`` + ``groupby2``)
      - ``"mode4"`` -- mixed model (``groupby`` + (1|``groupby2``))

    Modes 2-4 require ``groupby2``.

    ``test_method`` (significance annotations): ``"none"``, ``"mwu-test"``,
    ``"ttest"``, ``"kw-test"``, ``"anova"``, ``"linear-model"``,
    ``"linear-model-interaction"``, ``"mixed-model"``.
    """
    if (linked := _linked_spec("violin_grouped", locals())) is not None:
        return linked
    if adata is None or key is None or groupby is None:
        raise TypeError(
            "`adata`, `key`, and `groupby` are required for a standalone grouped violin."
        )
    from guanaco.pages.matrix.plots.violin2 import plot_violin2_new

    meta2 = None if mode == "mode1" or groupby2 == "none" else groupby2
    if mode in ("mode2", "mode3", "mode4") and meta2 is None:
        raise ValueError(f"mode '{mode}' requires `groupby2`.")
    n_colors = _obs_col(adata.obs, groupby).nunique()
    if meta2:
        n_colors = max(n_colors, _obs_col(adata.obs, meta2).nunique())
    fig = plot_violin2_new(
        adata,
        key=key,
        meta1=groupby,
        meta2=meta2,
        mode=mode,
        transformation=transformation,
        layer=layer,
        show_box=show_box,
        show_points=show_points,
        test_method=test_method,
        labels=list(labels) if labels is not None else None,
        bandwidth=bandwidth,
        color_map=None,
        palette=_palette_list(palette, n_colors),
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


# --------------------------------------------------------------------------- #
# Heatmap
# --------------------------------------------------------------------------- #
def heatmap(
    adata=None,
    var_names=None,
    groupby: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    groupby2: str | None = None,
    labels: Sequence | None = None,
    log: bool = False,
    z_score: bool = False,
    color_map: str = "Viridis",
    transformation: str | None = None,
    standardization: str | None = None,
    layer: str | None = None,
    max_cells: int = 10000,
    n_bins: int = 4000,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Expression heatmap of ``var_names`` grouped by ``groupby`` (scanpy ``sc.pl.heatmap``)."""
    if (linked := _linked_spec("heatmap", locals())) is not None:
        return linked
    if adata is None or var_names is None or groupby is None:
        raise TypeError(
            "`adata`, `var_names`, and `groupby` are required for a standalone heatmap."
        )
    from guanaco.pages.matrix.plots.heatmap import plot_unified_heatmap

    fig = plot_unified_heatmap(
        adata=adata,
        genes=_as_list(var_names),
        groupby1=groupby,
        groupby2=groupby2,
        labels=list(labels) if labels is not None else None,
        log=log,
        z_score=z_score,
        color_map=color_map,
        transformation=transformation,
        standardization=standardization,
        layer=layer,
        max_cells=max_cells,
        n_bins=n_bins,
        adata_obs=adata.obs,
        color_config=list(GUANACO_QUALITATIVE),
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


# --------------------------------------------------------------------------- #
# Dotplot / matrixplot (both from dotmatrix.plot_dot_matrix)
# --------------------------------------------------------------------------- #
def _dot_matrix(
    adata,
    var_names,
    groupby,
    *,
    plot_type,
    title,
    labels,
    color_map,
    transformation,
    standardization,
    layer,
    cluster,
    transpose,
    height,
    show,
    return_fig,
):
    from guanaco.pages.matrix.plots.dotmatrix import plot_dot_matrix

    used_labels = list(labels) if labels is not None else _all_labels(adata, groupby)
    fig = plot_dot_matrix(
        adata=adata,
        genes=_as_list(var_names),
        groupby=groupby,
        selected_labels=used_labels,
        transformation=transformation,
        standardization=standardization,
        layer=layer,
        color_map=color_map,
        plot_type=plot_type,
        cluster=cluster,
        transpose=transpose,
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


def dotplot(
    adata=None,
    var_names=None,
    groupby: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    labels: Sequence | None = None,
    color_map: str = "Viridis",
    transformation: str | None = None,
    standardization: str | None = None,
    layer: str | None = None,
    cluster: str = "none",
    transpose: bool = False,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Dotplot of ``var_names`` grouped by ``groupby`` (scanpy ``sc.pl.dotplot``).

    Dot color = mean expression, dot size = fraction of cells expressing.
    """
    if (linked := _linked_spec("dotplot", locals())) is not None:
        return linked
    if adata is None or var_names is None or groupby is None:
        raise TypeError(
            "`adata`, `var_names`, and `groupby` are required for a standalone dotplot."
        )
    return _dot_matrix(
        adata,
        var_names,
        groupby,
        plot_type="dotplot",
        title=title,
        labels=labels,
        color_map=color_map,
        transformation=transformation,
        standardization=standardization,
        layer=layer,
        cluster=cluster,
        transpose=transpose,
        height=height,
        show=show,
        return_fig=return_fig,
    )


def matrixplot(
    adata=None,
    var_names=None,
    groupby: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    labels: Sequence | None = None,
    color_map: str = "Viridis",
    transformation: str | None = None,
    standardization: str | None = None,
    layer: str | None = None,
    cluster: str = "none",
    transpose: bool = False,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Matrixplot: a mean-expression heatmap of ``var_names`` x ``groupby``
    (scanpy ``sc.pl.matrixplot``).

    Same data as :func:`dotplot` but drawn as a colored grid (no dot sizing).
    ``cluster`` adds dendrograms: ``"none"`` / ``"row"`` / ``"col"`` / ``"both"``.
    ``standardization`` scales the means: ``None`` / ``"var"`` (per gene) /
    ``"group"`` (per group). ``transpose`` swaps the gene/group axes.
    """
    if (linked := _linked_spec("matrixplot", locals())) is not None:
        return linked
    if adata is None or var_names is None or groupby is None:
        raise TypeError(
            "`adata`, `var_names`, and `groupby` are required for a standalone matrixplot."
        )
    return _dot_matrix(
        adata,
        var_names,
        groupby,
        plot_type="matrixplot",
        title=title,
        labels=labels,
        color_map=color_map,
        transformation=transformation,
        standardization=standardization,
        layer=layer,
        cluster=cluster,
        transpose=transpose,
        height=height,
        show=show,
        return_fig=return_fig,
    )


# --------------------------------------------------------------------------- #
# Stacked bar (composition)
# --------------------------------------------------------------------------- #
def stacked_bar(
    adata=None,
    x: str | None = None,
    color: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    normalize: str | bool = "proportion",
    x_order: Sequence | None = None,
    palette=None,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Stacked composition bars: ``x`` on the x-axis, stacked/colored by ``color``.

    ``normalize`` -> proportion ("proportion"/"prop"/True) or raw counts.
    """
    if (linked := _linked_spec("stacked_bar", locals())) is not None:
        return linked
    if adata is None or x is None or color is None:
        raise TypeError(
            "`adata`, `x`, and `color` are required for a standalone stacked bar."
        )
    from guanaco.pages.matrix.plots.stacked_bar import plot_stacked_bar

    norm = "prop" if normalize in ("proportion", "prop", True) else "count"
    fig = plot_stacked_bar(
        x_meta=x,
        y_meta=color,
        norm=norm,
        adata=adata,
        color_map=_label_color_map(adata, color, palette),
        x_order=list(x_order) if x_order is not None else None,
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


# --------------------------------------------------------------------------- #
# Expression trend along pseudotime
# --------------------------------------------------------------------------- #
def pseudotime(
    adata=None,
    genes=None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    pseudotime_key: str = "pseudotime",
    groupby: str | None = None,
    min_expr: float = 0.5,
    transformation: str = "none",
    layer: str | None = None,
    size: int = 3,
    opacity: float = 0.6,
    palette=None,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Smoothed gene-expression trend along a pseudotime in ``obs``."""
    if (linked := _linked_spec("pseudotime", locals())) is not None:
        return linked
    if adata is None or genes is None:
        raise TypeError(
            "`adata` and `genes` are required for a standalone pseudotime plot."
        )
    from guanaco.pages.matrix.plots.pseudotime import plot_genes_in_pseudotime

    color_map = _label_color_map(adata, groupby, palette) if groupby else None
    fig = plot_genes_in_pseudotime(
        adata=adata,
        genes=_as_list(genes),
        pseudotime_key=pseudotime_key,
        groupby=groupby,
        min_expr=min_expr,
        transformation=transformation,
        layer=layer,
        color_map=color_map,
        marker_size=size,
        opacity=opacity,
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


# --------------------------------------------------------------------------- #
# ATAC Peak Browser (genome browser over peak-like var_names)
# --------------------------------------------------------------------------- #
def _resolve_peak_region(adata, region, gene_index):
    """Turn ``region`` into a ``{"chrom","start","end"}`` dict.

    ``region`` may be ``None`` (a populated default window), a locus string, a gene
    name (resolved via ``gene_index``), or a coordinate dict — the same "locus or
    gene" search the web app's box accepts.
    """
    from guanaco.pages.matrix.plots.atac_browser import default_region, parse_locus

    if region is None:
        return default_region(adata)
    if not isinstance(region, str):
        return dict(region)

    parsed = parse_locus(region)
    if parsed is not None:
        chrom, start, end = parsed
        return {"chrom": chrom, "start": start, "end": end}

    # Not a locus -> treat it as a gene name (needs the annotation).
    if gene_index is not None:
        from guanaco.pages.matrix.plots.gene_annotation import find_gene_region

        gene_region = find_gene_region(gene_index, region)
        if gene_region is not None:
            return gene_region
        raise ValueError(
            f"Gene {region!r} not found in the annotation. Pass a locus like "
            "'chr1:1,000,000-2,000,000', or a gene present in gene_annotation."
        )
    raise ValueError(
        f"Could not parse {region!r} as a locus. To search by gene name, also pass "
        "gene_annotation='hg38' (or a path to a GTF/GFF3 file)."
    )


def peak_browser(
    adata,
    region: "str | dict | None" = None,
    *,
    title: str | None = None,
    groupby: str | None = None,
    labels=None,
    metric: str = "mean",
    max_peaks: int = 400,
    gene_annotation: str | None = None,
    selected_cells=None,
    palette=None,
    y_mode: str = "shared",
    show: bool = True,
    return_fig: bool = False,
):
    """ATAC accessibility genome browser over peak-like features (GUANACO's "Peak Browser").

    Needs genomic peaks: either ``var_names`` like ``"chr1:10000-10500"`` or
    ``adata.var`` columns ``["chrom", "start", "end"]``. Draws one accessibility bar
    track per group (``groupby``).

    ``region`` accepts a **gene name** (``"CD8A"``) or a locus string
    (``"chr1:1,000,000-2,000,000"``) — the same "locus or gene" search the web app's
    box accepts — or a ``{"chrom", "start", "end"}`` dict; it defaults to a populated
    window. Gene-name search requires ``gene_annotation`` to be set. ``metric`` is
    ``"mean"`` accessibility or ``"detection"`` fraction. ``y_mode`` is ``"shared"``
    (default — every track uses the same y-range so heights are directly comparable)
    or ``"auto"`` (each track scales to its own peak). Pass ``gene_annotation``
    (a genome id like ``"hg38"`` / ``"mm10"``, or a path to a GTF/GFF3 file) to add a
    gene-model track above the signal tracks.
    """
    from guanaco.pages.matrix.plots.atac_browser import (
        compute_atac_signal,
        has_genomic_peak_features,
        plot_atac_browser,
    )

    if not has_genomic_peak_features(adata):
        raise ValueError(
            "peak_browser needs genomic peak features: var_names like "
            "'chr1:10000-10500', or var columns ['chrom', 'start', 'end']."
        )

    # Load the gene annotation up front: it powers both the optional gene track and
    # gene-name search in ``region`` (mirrors the web app's single search box).
    gene_index = None
    if gene_annotation:
        from guanaco.pages.matrix.plots.gene_annotation import (
            load_gene_annotation,
            resolve_annotation_source,
        )

        gene_index = load_gene_annotation(resolve_annotation_source(gene_annotation))

    region_dict = _resolve_peak_region(adata, region, gene_index)

    group_order = None
    color_map = None
    if groupby and groupby in adata.obs.columns:
        group_order = _all_labels(adata, groupby)
        color_map = _label_color_map(adata, groupby, palette)

    payload = compute_atac_signal(
        adata,
        region_dict,
        selected_cells=selected_cells,
        groupby=groupby,
        labels=_as_list(labels) if labels is not None else None,
        group_order=group_order,
        metric=metric,
        max_peaks=max_peaks,
    )

    gene_models = None
    if gene_index is not None:
        from guanaco.pages.matrix.plots.gene_annotation import query_gene_models

        r = payload["region"]
        gene_models = query_gene_models(
            gene_index, str(r["chrom"]), int(r["start"]), int(r["end"])
        )

    fig = plot_atac_browser(
        payload,
        gene_models=gene_models,
        color_map=color_map,
        y_mode=y_mode,
    )
    return _render(fig, show, return_fig, title=title)


# --------------------------------------------------------------------------- #
# Volcano (from precomputed DE in adata.uns)
# --------------------------------------------------------------------------- #
def volcano(
    adata=None,
    group: str | None = None,
    *,
    id: str | None = None,
    data=None,
    title: str | None = None,
    x_field: str = "logfoldchange",
    padj_threshold: float = 0.05,
    x_threshold: float = 1.0,
    top_n: int = 12,
    height: str | int | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Volcano plot from precomputed DE in ``adata.uns`` (``rank_genes_groups`` or ``volcano``).

    ``group`` selects the comparison/entry; defaults to the first available.
    """
    if (linked := _linked_spec("volcano", locals())) is not None:
        return linked
    if adata is None:
        raise TypeError("`adata` is required for a standalone volcano plot.")
    from guanaco.pages.matrix.plots.volcano import load_volcano_payload, plot_volcano

    entries = load_volcano_payload(adata)["entries"]
    if group is None:
        group = next(iter(entries))
    if group not in entries:
        raise KeyError(f"Volcano group '{group}' not found. Available: {list(entries)}")
    fig = plot_volcano(
        entry_name=group,
        entry=entries[group],
        x_field=x_field,
        padj_threshold=padj_threshold,
        x_threshold=x_threshold,
        top_n=top_n,
    )
    return _render(fig, title=title, height=height, show=show, return_fig=return_fig)


# --------------------------------------------------------------------------- #
# PAGA Cytoscape rendering -- interactive ipycytoscape widgets (Jupyter)
# --------------------------------------------------------------------------- #
def _cytoscape_legend_html(legend) -> str:
    """A small HTML legend (color swatches or a continuous ramp) for the widget."""
    if not legend:
        return ""
    title = legend.get("title") or ""
    if legend.get("kind") == "continuous":
        vmin, vmax, stops = (
            legend.get("vmin"),
            legend.get("vmax"),
            legend.get("stops") or [],
        )
        if vmin is None or not stops:
            body = "<div style='color:#6B7280'>no finite values</div>"
        else:
            gradient = ", ".join(stops)
            body = (
                "<div style='display:flex;align-items:stretch;gap:8px;height:150px;'>"
                f"<div style='width:18px;border:1px solid #AEB6C2;border-radius:2px;"
                f"background:linear-gradient(to top, {gradient});'></div>"
                "<div style='display:flex;flex-direction:column;justify-content:space-between;"
                "font-size:11px;color:#2F3E46;'>"
                f"<span>{vmax:.4g}</span><span>{vmin:.4g}</span></div></div>"
            )
    else:  # categorical
        body = "".join(
            "<div style='display:flex;align-items:center;margin:5px 0;'>"
            f"<span style='width:14px;height:14px;border-radius:50%;background:{e['color']};"
            "display:inline-block;margin-right:8px;border:1px solid rgba(0,0,0,.2);"
            "flex:0 0 auto;'></span>"
            f"<span>{e['label']}</span></div>"
            for e in legend.get("entries", [])
        )
    return (
        "<div style='font:13px/1.4 -apple-system,Segoe UI,Helvetica,sans-serif;color:#2F3E46;"
        "padding:10px 14px;max-height:560px;overflow-y:auto;'>"
        f"<div style='font-weight:700;margin-bottom:8px;'>{title}</div>{body}</div>"
    )


_CYTOSCAPE_CDN = "https://cdn.jsdelivr.net/npm/cytoscape@3.30.2/dist/cytoscape.min.js"
_CYTOSCAPE_JS_CACHE = None


def _get_cytoscape_js() -> str:
    """Fetch & cache cytoscape.min.js so it can be inlined (offline HTML export).

    Returns the JS source, or ``""`` if it can't be fetched -- the caller then
    falls back to a CDN ``<script src>`` (which needs internet to view).
    """
    global _CYTOSCAPE_JS_CACHE
    if _CYTOSCAPE_JS_CACHE is None:
        try:
            import urllib.request

            with urllib.request.urlopen(_CYTOSCAPE_CDN, timeout=20) as resp:
                _CYTOSCAPE_JS_CACHE = resp.read().decode("utf-8")
        except Exception:
            _CYTOSCAPE_JS_CACHE = ""
    return _CYTOSCAPE_JS_CACHE


def _cytoscape_html(graph, *, height: int = 560) -> str:
    """Self-contained cytoscape.js network (in a sandboxed iframe).

    Renders identically to the ipycytoscape widget (same elements, stylesheet,
    layout, pie nodes and legend) but as plain HTML+JS, so it survives a static
    notebook -> HTML export and works offline (cytoscape.js is inlined). The
    iframe ``srcdoc`` also bypasses Jupyter's output-script sanitizer. Returns an
    HTML string for ``IPython.display.HTML``.
    """
    import json as _json

    def _js(obj):
        # Embed JSON inside a <script>; guard against a literal </script>.
        return _json.dumps(obj).replace("</", "<\\/")

    legend_html = _cytoscape_legend_html(graph.get("legend"))
    cyto_js = _get_cytoscape_js()
    if cyto_js:
        head_js = "<script>" + cyto_js.replace("</script", "<\\/script") + "</script>"
    else:
        head_js = "<script src='" + _CYTOSCAPE_CDN + "'></script>"

    doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        + head_js
        + "<style>html,body{margin:0;height:100%}"
        "#wrap{display:flex;height:" + str(height) + "px;"
        "font-family:-apple-system,Segoe UI,Helvetica,sans-serif}"
        "#cy{flex:1 1 auto;min-width:0;height:100%}"
        "#legend{flex:0 0 auto;overflow:auto;max-height:100%}</style></head>"
        "<body><div id='wrap'><div id='cy'></div>"
        "<div id='legend'>" + legend_html + "</div></div><script>"
        "var cy=cytoscape({container:document.getElementById('cy'),"
        "elements:"
        + _js(graph["elements"])
        + ",style:"
        + _js(graph["stylesheet"])
        + ",layout:"
        + _js(graph.get("layout", {"name": "cose"}))
        + ",minZoom:0.2,maxZoom:3});</script></body></html>"
    )
    srcdoc = doc.replace("&", "&amp;").replace('"', "&quot;")
    return (
        "<iframe sandbox='allow-scripts' srcdoc=\"" + srcdoc + '" '
        "style='width:100%;height:" + str(height + 24) + "px;border:1px solid #e5e7eb;"
        "border-radius:6px;'></iframe>"
    )


def _render_cytoscape(graph, *, show, return_widget, renderer="widget"):
    """Render a backend-agnostic cytoscape.js graph dict for a notebook.

    ``renderer="widget"`` -> an interactive ipycytoscape widget (live Jupyter
    only). ``renderer="html"`` -> a self-contained cytoscape.js iframe that also
    renders in a static, offline HTML export. Both show the network beside an
    HTML legend. These return an ipywidget / HTML object (not a Plotly figure),
    so they use ``return_widget`` instead of ``return_fig``.
    """
    if renderer == "html":
        from IPython.display import HTML

        return _show_or_return(
            HTML(_cytoscape_html(graph)), show=show, return_widget=return_widget
        )

    try:
        import ipycytoscape
    except ImportError as exc:  # optional dependency
        raise ImportError(
            "PAGA in notebooks needs ipycytoscape. "
            "Install it with `pip install ipycytoscape` (Jupyter/JupyterLab)."
        ) from exc

    elements = graph["elements"]
    # Edges carry data.source/data.target; nodes don't -- split for ipycytoscape.
    nodes = [e for e in elements if "source" not in e.get("data", {})]
    edges = [e for e in elements if "source" in e.get("data", {})]

    widget = ipycytoscape.CytoscapeWidget()
    widget.graph.add_graph_from_json(
        {"nodes": nodes, "edges": edges},
        directed=graph.get("directed", False),
    )
    # cytoscape.js accepts the same stylesheet (selectors + style, incl. pie and
    # data() mappers) we feed the Dash component, so it transfers unchanged.
    widget.set_style(graph["stylesheet"])
    layout = graph.get("layout", {"name": "cose"})
    widget.set_layout(
        name=layout.get("name", "cose"), padding=layout.get("padding", 40)
    )

    # Hover a node to see its composition / value (the Dash app's side panel).
    if any("hover_text" in n.get("data", {}) for n in nodes):
        try:
            widget.set_tooltip_source("hover_text")
        except Exception:
            pass

    # Pair the network with the legend so colors are labeled, like the web app.
    legend_html = _cytoscape_legend_html(graph.get("legend"))
    rendered = widget
    if legend_html:
        try:
            import ipywidgets as W

            widget.layout.width = "100%"
            widget.layout.min_height = "540px"
            net = W.Box([widget], layout=W.Layout(flex="1 1 auto", min_width="0"))
            side = W.HTML(value=legend_html, layout=W.Layout(flex="0 0 auto"))
            rendered = W.HBox(
                [net, side], layout=W.Layout(align_items="stretch", width="100%")
            )
        except Exception:
            rendered = widget

    return _show_or_return(rendered, show=show, return_widget=return_widget)


def paga(
    adata,
    *,
    color: str | None = None,
    gene: str | None = None,
    color_map: str = "Viridis",
    palette=None,
    edge_threshold: float = 0.03,
    node_font_size: int = 18,
    renderer: str = "widget",
    show: bool = True,
    return_widget: bool = False,
):
    """PAGA graph (pie-chart nodes) as an interactive cytoscape network.

    Reads ``adata.uns['paga']`` (run ``sc.tl.paga`` first). Needs
    ``pip install ipycytoscape``.

    ``renderer="widget"`` (default) is an interactive ipycytoscape widget (live
    Jupyter/JupyterLab). ``renderer="html"`` is a self-contained cytoscape.js
    iframe that also shows up in a **static, offline HTML export** of the notebook
    (use this when you `nbconvert --to html`).

    Color the nodes by an ``obs`` column (``color=...``: categorical -> pie nodes
    showing composition, numeric -> a solid continuous color) or by a ``gene``
    (mean expression per group). With neither, nodes are colored by the column
    PAGA was computed on. ``palette`` overrides categorical colors (default: the
    dataset's stored colors, else GUANACO's palette); ``edge_threshold`` hides
    weak connectivities; ``node_font_size`` sets the cluster-label size.

    The network is shown next to a **legend** mapping pie colors to categories
    (and, in widget mode, hovering a node reveals its exact composition).
    """
    from guanaco.pages.matrix.plots.paga import paga_graph

    graph = paga_graph(
        adata,
        color_mode="gene" if gene else "obs",
        obs_key=color,
        gene=gene,
        continuous_color_map=color_map,
        discrete_palette=_palette_list(palette, 256) if palette is not None else None,
        edge_threshold=edge_threshold,
        node_font_size=node_font_size,
    )
    return _render_cytoscape(
        graph, show=show, return_widget=return_widget, renderer=renderer
    )


# --------------------------------------------------------------------------- #
# Declarative linked views
# --------------------------------------------------------------------------- #
def view(
    plot: str,
    id: str,
    data=None,
    title: str | None = None,
    **options,
):
    """Describe one plot in a :func:`linked_view` workspace."""

    from guanaco.linking import view as build_view

    return build_view(plot, id=id, data=data, title=title, **options)


def link(
    source: str,
    target: str,
    *,
    by: Literal["cell", "feature"] | None = None,
    action: Literal["highlight", "filter"] | None = None,
    key: str | None = None,
):
    """Connect two views by their DataFrame, cell, or feature index.

    Omit ``by`` for AnnData-to-AnnData cell links and DataFrame-to-DataFrame
    row links. Cross-container links declare ``by="cell"`` or
    ``by="feature"``. Cell links highlight by default; row links filter;
    feature links update the target's feature context. ``key`` names a shared
    table column containing the logical ID; for mixed table/AnnData links it
    maps table rows to native ``obs_names`` or ``var_names``. The DataFrame
    index remains the unique atomic row ID.
    """

    from guanaco.linking import link as build_link

    return build_link(source, target, by=by, action=action, key=key)


def linked_view(
    data,
    *,
    views: Sequence,
    links: Sequence,
    layout: Literal["overview-detail", "grid", "column"] = "overview-detail",
    registry=None,
    prefix: str | None = None,
    title: str = "GUANACO Linked View",
):
    """Compile user-selected plots and biological links into a Dash workspace."""

    from guanaco.linking import linked_view as build_linked_view

    return build_linked_view(
        data,
        views=views,
        links=links,
        layout=layout,
        registry=registry,
        prefix=prefix,
        title=title,
    )


def linked_plot_types(plot: str | None = None):
    """Describe built-in linked plots, optionally returning just one contract."""

    from guanaco.linking import default_plot_registry

    registry = default_plot_registry()
    return registry.contract(plot) if plot is not None else registry.contracts()


__all__ = [
    "embedding",
    "umap",
    "pca",
    "tsne",
    "coexpression",
    "violin",
    "stacked_violin",
    "violin_grouped",
    "heatmap",
    "dotplot",
    "matrixplot",
    "stacked_bar",
    "pseudotime",
    "peak_browser",
    "volcano",
    "paga",
    "view",
    "link",
    "linked_view",
    "linked_plot_types",
]
