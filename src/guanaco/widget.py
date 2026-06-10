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

PAGA and GRN are Cytoscape networks, so ``gc.pl.paga`` / ``gc.pl.grn`` return an
``ipycytoscape`` widget instead of a Plotly figure (use ``return_widget`` rather
than ``return_fig``). They need ``pip install ipycytoscape`` and render in
Jupyter/JupyterLab; classic ipywidgets are not currently supported in marimo.
"""

from __future__ import annotations

from typing import Any, Sequence

import plotly.graph_objects as go


# --------------------------------------------------------------------------- #
# Internal helpers
# --------------------------------------------------------------------------- #
def _render(fig: go.Figure, show: bool, return_fig: bool):
    """scanpy-style output: return the figure, or show it inline."""
    if return_fig:
        return fig
    if show:
        fig.show()
        return None
    return fig


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
    """Turn a palette name or list into a list of colors (or None)."""
    if palette is None:
        return None
    if isinstance(palette, str):
        from guanaco.utils.colors import resolve_discrete_palette
        return resolve_discrete_palette(palette, n_colors)
    return list(palette)


def _label_color_map(adata, groupby: str, palette):
    """Build a stable {category: color} map for ``groupby`` (or None)."""
    if palette is None:
        return None
    labels = sorted(adata.obs[groupby].unique().tolist())
    colors = _palette_list(palette, len(labels))
    if not colors:
        return None
    return {lab: colors[i % len(colors)] for i, lab in enumerate(labels)}


def _all_labels(adata, groupby: str) -> list:
    return sorted(adata.obs[groupby].unique().tolist())


def _default_palette() -> list:
    """A config-independent default categorical palette.

    Some plot functions otherwise fall back to ``guanaco.data.registry``, which
    loads a config file at import time -- absent in a notebook. Passing this
    keeps the notebook API self-contained.
    """
    try:
        from guanaco.data.loader import DEFAULT_COLORS
        return list(DEFAULT_COLORS)
    except Exception:
        import plotly.express as px
        return list(px.colors.qualitative.Plotly)


def _effective_palette(palette):
    """The user's palette, or GUANACO's shared default (the violin/Okabe-Ito set).

    Routing every categorical plot through this keeps the default colors
    consistent across the notebook API.
    """
    return palette if palette is not None else _default_palette()


# --------------------------------------------------------------------------- #
# Embeddings / scatter
# --------------------------------------------------------------------------- #
def embedding(
    adata,
    basis: str,
    color: str,
    *,
    transformation: str | None = None,
    order: str | None = None,
    color_map: str = "Viridis",
    palette=None,
    size: int = 5,
    opacity: float = 1.0,
    legend_loc: str = "right margin",
    axis_show: bool = True,
    render_backend: str = "scattergl",
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
    from guanaco.pages.matrix.plots.embedding import plot_embedding

    embedding_key = _resolve_basis(adata, basis)
    n_cat = adata.obs[color].nunique() if color in adata.obs.columns else 50
    legend_show = "on data" if legend_loc == "on data" else "on legend"
    fig = plot_embedding(
        adata=adata,
        adata_full=adata,
        source_adata=adata,
        embedding_key=embedding_key,
        color=color,
        mode="auto",
        transformation=transformation,
        order=order,
        continuous_color_map=color_map,
        discrete_color_map=_palette_list(_effective_palette(palette), n_cat),
        marker_size=size,
        opacity=opacity,
        legend_show=legend_show,
        axis_show=axis_show,
        render_backend=render_backend,
    )
    return _render(fig, show, return_fig)


def umap(adata, color: str, **kwargs):
    """``embedding`` with ``basis="X_umap"`` (scanpy ``sc.pl.umap``)."""
    return embedding(adata, "X_umap", color, **kwargs)


def pca(adata, color: str, **kwargs):
    """``embedding`` with ``basis="X_pca"`` (scanpy ``sc.pl.pca``)."""
    return embedding(adata, "X_pca", color, **kwargs)


def tsne(adata, color: str, **kwargs):
    """``embedding`` with ``basis="X_tsne"`` (scanpy ``sc.pl.tsne``)."""
    return embedding(adata, "X_tsne", color, **kwargs)


def coexpression(
    adata,
    basis: str,
    gene1: str,
    gene2: str,
    *,
    threshold1: float = 0.5,
    threshold2: float = 0.5,
    transformation: str | None = None,
    size: int = 5,
    opacity: float = 1.0,
    show: bool = True,
    return_fig: bool = False,
):
    """Two-gene co-expression on an embedding (cells binned into 4 groups)."""
    from guanaco.pages.matrix.plots.embedding import plot_coexpression_embedding

    fig = plot_coexpression_embedding(
        adata=adata,
        embedding_key=_resolve_basis(adata, basis),
        gene1=gene1,
        gene2=gene2,
        threshold1=threshold1,
        threshold2=threshold2,
        transformation=transformation,
        marker_size=size,
        opacity=opacity,
    )
    return _render(fig, show, return_fig)


# --------------------------------------------------------------------------- #
# Violin (stacked: genes as rows, grouped by a category)
# --------------------------------------------------------------------------- #
def violin(
    adata,
    keys,
    groupby: str,
    *,
    labels: Sequence | None = None,
    transformation: str | None = None,
    show_box: bool = False,
    palette=None,
    show: bool = True,
    return_fig: bool = False,
):
    """Stacked violins of ``keys`` (genes) split by ``groupby`` (scanpy ``sc.pl.violin``).

    ``labels`` restricts/orders the ``groupby`` categories (default: all).
    ``show_box`` overlays a neutral Median + IQR box.
    """
    from guanaco.pages.matrix.plots.violin1 import plot_violin1

    genes = _as_list(keys)
    used_labels = list(labels) if labels is not None else _all_labels(adata, groupby)
    fig = plot_violin1(
        adata=adata,
        genes=genes,
        groupby=groupby,
        labels=used_labels,
        transformation=transformation,
        show_box=show_box,
        groupby_label_color_map=_label_color_map(adata, groupby, _effective_palette(palette)),
        adata_obs=adata.obs,
    )
    return _render(fig, show, return_fig)


def stacked_violin(adata, var_names, groupby: str, **kwargs):
    """Alias of :func:`violin` (scanpy ``sc.pl.stacked_violin``)."""
    return violin(adata, var_names, groupby, **kwargs)


def violin_grouped(
    adata,
    key: str,
    groupby: str,
    *,
    groupby2: str | None = None,
    mode: str = "mode1",
    test_method: str = "none",
    transformation: str | None = None,
    show_box: bool = True,
    show_points: bool = False,
    labels: Sequence | None = None,
    palette=None,
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
    from guanaco.pages.matrix.plots.violin2 import plot_violin2_new

    # Same meta2/mode handling as the web callback.
    meta2 = None if mode == "mode1" else groupby2
    if meta2 == "none":
        meta2 = None
    if mode in ("mode2", "mode3", "mode4") and meta2 is None:
        raise ValueError(f"mode '{mode}' requires `groupby2`.")

    n_colors = adata.obs[groupby].nunique()
    if meta2:
        n_colors = max(n_colors, adata.obs[meta2].nunique())

    fig = plot_violin2_new(
        adata,
        key=key,
        meta1=groupby,
        meta2=meta2,
        mode=mode,
        transformation=transformation,
        show_box=show_box,
        show_points=show_points,
        test_method=test_method,
        labels=list(labels) if labels is not None else None,
        color_map=None,
        palette=_palette_list(_effective_palette(palette), n_colors),
    )
    return _render(fig, show, return_fig)


# --------------------------------------------------------------------------- #
# Heatmap
# --------------------------------------------------------------------------- #
def heatmap(
    adata,
    var_names,
    groupby: str,
    *,
    groupby2: str | None = None,
    labels: Sequence | None = None,
    log: bool = False,
    z_score: bool = False,
    color_map: str = "Viridis",
    transformation: str | None = None,
    show: bool = True,
    return_fig: bool = False,
):
    """Expression heatmap of ``var_names`` grouped by ``groupby`` (scanpy ``sc.pl.heatmap``)."""
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
        adata_obs=adata.obs,
        color_config=_default_palette(),
    )
    return _render(fig, show, return_fig)


# --------------------------------------------------------------------------- #
# Dotplot / matrixplot (both from dotmatrix.plot_dot_matrix)
# --------------------------------------------------------------------------- #
def _dot_matrix(
    adata,
    var_names,
    groupby,
    *,
    plot_type,
    labels,
    color_map,
    transformation,
    standardization,
    cluster,
    transpose,
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
        color_map=color_map,
        plot_type=plot_type,
        cluster=cluster,
        transpose=transpose,
    )
    return _render(fig, show, return_fig)


def dotplot(
    adata,
    var_names,
    groupby: str,
    *,
    labels: Sequence | None = None,
    color_map: str = "Viridis",
    transformation: str | None = None,
    standardization: str | None = None,
    cluster: str = "none",
    transpose: bool = False,
    show: bool = True,
    return_fig: bool = False,
):
    """Dotplot of ``var_names`` grouped by ``groupby`` (scanpy ``sc.pl.dotplot``).

    Dot color = mean expression, dot size = fraction of cells expressing.
    """
    return _dot_matrix(
        adata, var_names, groupby, plot_type="dotplot",
        labels=labels, color_map=color_map, transformation=transformation,
        standardization=standardization, cluster=cluster, transpose=transpose,
        show=show, return_fig=return_fig,
    )


def matrixplot(
    adata,
    var_names,
    groupby: str,
    *,
    labels: Sequence | None = None,
    color_map: str = "Viridis",
    transformation: str | None = None,
    standardization: str | None = None,
    cluster: str = "none",
    transpose: bool = False,
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
    return _dot_matrix(
        adata, var_names, groupby, plot_type="matrixplot",
        labels=labels, color_map=color_map, transformation=transformation,
        standardization=standardization, cluster=cluster, transpose=transpose,
        show=show, return_fig=return_fig,
    )


# --------------------------------------------------------------------------- #
# Stacked bar (composition)
# --------------------------------------------------------------------------- #
def stacked_bar(
    adata,
    x: str,
    color: str,
    *,
    normalize: str | bool = "proportion",
    x_order: Sequence | None = None,
    palette=None,
    show: bool = True,
    return_fig: bool = False,
):
    """Stacked composition bars: ``x`` on the x-axis, stacked/colored by ``color``.

    ``normalize`` -> proportion ("proportion"/"prop"/True) or raw counts.
    """
    from guanaco.pages.matrix.plots.stacked_bar import plot_stacked_bar

    norm = "prop" if normalize in ("proportion", "prop", True) else "count"
    fig = plot_stacked_bar(
        x_meta=x,
        y_meta=color,
        norm=norm,
        adata=adata,
        color_map=_label_color_map(adata, color, _effective_palette(palette)),
        x_order=list(x_order) if x_order is not None else None,
    )
    return _render(fig, show, return_fig)


# --------------------------------------------------------------------------- #
# Expression trend along pseudotime
# --------------------------------------------------------------------------- #
def pseudotime(
    adata,
    genes,
    *,
    pseudotime_key: str = "pseudotime",
    groupby: str | None = None,
    min_expr: float = 0.5,
    transformation: str = "none",
    palette=None,
    show: bool = True,
    return_fig: bool = False,
):
    """Smoothed gene-expression trend along a pseudotime in ``obs``."""
    from guanaco.pages.matrix.plots.pseudotime import plot_genes_in_pseudotime

    color_map = _label_color_map(adata, groupby, _effective_palette(palette)) if groupby else None
    fig = plot_genes_in_pseudotime(
        adata=adata,
        genes=_as_list(genes),
        pseudotime_key=pseudotime_key,
        groupby=groupby,
        min_expr=min_expr,
        transformation=transformation,
        color_map=color_map,
    )
    return _render(fig, show, return_fig)


# --------------------------------------------------------------------------- #
# Volcano (from precomputed DE in adata.uns)
# --------------------------------------------------------------------------- #
def volcano(
    adata,
    group: str | None = None,
    *,
    x_field: str = "logfoldchange",
    padj_threshold: float = 0.05,
    x_threshold: float = 1.0,
    top_n: int = 12,
    show: bool = True,
    return_fig: bool = False,
):
    """Volcano plot from precomputed DE in ``adata.uns`` (``rank_genes_groups`` or ``volcano``).

    ``group`` selects the comparison/entry; defaults to the first available.
    """
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
    return _render(fig, show, return_fig)


# --------------------------------------------------------------------------- #
# Cytoscape networks (PAGA, GRN) -- interactive ipycytoscape widgets (Jupyter)
# --------------------------------------------------------------------------- #
def _cytoscape_legend_html(legend) -> str:
    """A small HTML legend (color swatches or a continuous ramp) for the widget."""
    if not legend:
        return ""
    title = legend.get("title") or ""
    if legend.get("kind") == "continuous":
        vmin, vmax, stops = legend.get("vmin"), legend.get("vmax"), legend.get("stops") or []
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
        "elements:" + _js(graph["elements"]) + ",style:" + _js(graph["stylesheet"])
        + ",layout:" + _js(graph.get("layout", {"name": "cose"}))
        + ",minZoom:0.2,maxZoom:3});</script></body></html>"
    )
    srcdoc = doc.replace("&", "&amp;").replace('"', "&quot;")
    return (
        "<iframe sandbox='allow-scripts' srcdoc=\"" + srcdoc + "\" "
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

        obj = HTML(_cytoscape_html(graph))
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

    try:
        import ipycytoscape
    except ImportError as exc:  # optional dependency
        raise ImportError(
            "PAGA/GRN networks in notebooks need ipycytoscape. "
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
    widget.set_layout(name=layout.get("name", "cose"), padding=layout.get("padding", 40))

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
            rendered = W.HBox([net, side], layout=W.Layout(align_items="stretch", width="100%"))
        except Exception:
            rendered = widget

    if return_widget:
        return rendered
    if show:
        try:
            from IPython.display import display

            display(rendered)
            return None
        except Exception:
            return rendered
    return rendered


def grn(
    adata,
    *,
    context: str = "All",
    edge_threshold: float | None = None,
    layout: str = "cose",
    node_font_size: int = 16,
    renderer: str = "widget",
    show: bool = True,
    return_widget: bool = False,
):
    """Gene-regulatory-network graph as an interactive cytoscape network.

    Reads ``adata.uns['grn']`` (columns ``source``/``target``/``regulation``,
    optional ``weight`` and a context column). Needs ``pip install ipycytoscape``.

    ``renderer="widget"`` (default) is an interactive ipycytoscape widget (live
    Jupyter/JupyterLab). ``renderer="html"`` is a self-contained cytoscape.js
    iframe that also shows up in a **static, offline HTML export** of the notebook
    (use this when you `nbconvert --to html`). ``context`` filters to one
    cell-type/condition ("All" = every context); ``edge_threshold`` hides weak
    edges; ``layout`` is any cytoscape.js layout ("cose", "circle", "concentric",
    "breadthfirst"); ``node_font_size`` sets the gene-label size.
    """
    from guanaco.pages.matrix.plots.grn_demo import grn_graph

    graph = grn_graph(
        adata,
        selected_context=context,
        edge_threshold=edge_threshold,
        layout_name=layout,
        node_font_size=node_font_size,
    )
    return _render_cytoscape(graph, show=show, return_widget=return_widget, renderer=renderer)


def paga(
    adata,
    *,
    color: str | None = None,
    gene: str | None = None,
    color_map: str = "Viridis",
    palette=None,
    edge_threshold: float = 0.03,
    annotation: str | None = None,
    labels: Sequence | None = None,
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
    weak connectivities; ``node_font_size`` sets the cluster-label size;
    ``annotation``/``labels`` optionally restrict the cells.

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
        selected_annotation=annotation,
        selected_labels=list(labels) if labels is not None else None,
        node_font_size=node_font_size,
    )
    return _render_cytoscape(graph, show=show, return_widget=return_widget, renderer=renderer)


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
    "volcano",
    "grn",
    "paga",
]
