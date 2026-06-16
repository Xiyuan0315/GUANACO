"""Interactive control-bar panels for GUANACO plots in `marimo <https://marimo.io>`_.

Each function here builds a *control bar* -- a row of dropdowns / sliders wired to
a GUANACO plot -- and returns a ``Panel``. A ``Panel`` is a real marimo UI element
(a subclass of ``mo.ui.dictionary``) with one extra thing: a ``.plot`` property that
renders the figure from the current control values.

Because marimo's reactivity is cell-based -- a cell cannot both *create* a UI
element and *reactively read* it -- use the two-cell idiom::

    import guanaco.marimo as gm

    dp = gm.dotplot(adata)   # cell 1: build the control bar
    dp.hstack()              #         (display it; `dp` alone works too)

    dp.plot                  # cell 2: the reactive figure

Touching a control in cell 1 re-runs cell 2 and redraws the figure. ``gm.show(dp)``
is an alias for ``dp.plot`` if you prefer a function call.

This wraps the same plotting code as the ``gc.pl`` notebook API; it only adds the
marimo controls. ``marimo`` itself is an optional dependency (``pip install marimo``).
"""

from __future__ import annotations

from typing import Any, Callable

__all__ = [
    "embedding",
    "violin",
    "dotplot",
    "heatmap",
    "stacked_bar",
    "pseudotime",
    "peak_browser",
    "volcano",
    "show",
]


# --------------------------------------------------------------------------- #
# marimo + Panel plumbing
# --------------------------------------------------------------------------- #
def _mo():
    try:
        import marimo as mo
    except ImportError as exc:  # pragma: no cover - trivial guard
        raise ImportError(
            "guanaco.marimo requires marimo. Install it with `pip install marimo`."
        ) from exc
    return mo


_PANEL_CLASS = None


def _panel_class():
    """Lazily build the Panel class (subclass of mo.ui.dictionary).

    Defined lazily so importing ``guanaco.marimo`` never imports marimo itself.
    """
    global _PANEL_CLASS
    if _PANEL_CLASS is None:
        mo = _mo()
        base = type(mo.ui.dictionary({}))

        class Panel(base):  # type: ignore[valid-type, misc]
            """A control bar wired to a GUANACO plot.

            It *is* an ``mo.ui.dictionary`` (so marimo reacts to it), plus a
            ``.plot`` property that draws the figure from the current values.
            """

            @property
            def plot(self):
                build: Callable[[dict], Any] | None = getattr(self, "_gx_build", None)
                if build is None:
                    return None
                return build(self.value)

            def bar(self, **kwargs):
                """The controls laid out in a single horizontal row."""
                return self.hstack(**kwargs)

        _PANEL_CLASS = Panel
    return _PANEL_CLASS


def _panel(controls: dict, build: Callable[[dict], Any]):
    panel = _panel_class()(controls)
    panel._gx_build = build
    return panel


def show(panel):
    """Render a panel's figure from its current control values (alias of ``panel.plot``)."""
    return panel.plot


# --------------------------------------------------------------------------- #
# Dataset introspection helpers (build the option lists)
# --------------------------------------------------------------------------- #
def _genes(adata) -> list:
    return list(adata.var_names)


def _obs_cols(adata) -> list:
    return list(adata.obs.columns)


# Above this many distinct values, a numeric obs column is treated as continuous
# rather than a grouping key.
_MAX_DISCRETE_CARDINALITY = 50


def _discrete_obs(adata) -> list:
    """obs columns suitable for grouping (categorical / low-cardinality)."""
    cols = []
    for c in adata.obs.columns:
        s = adata.obs[c]
        if str(s.dtype) in ("category", "object", "bool") or s.nunique(dropna=True) <= _MAX_DISCRETE_CARDINALITY:
            cols.append(c)
    return cols or list(adata.obs.columns)


def _numeric_obs(adata) -> list:
    return list(adata.obs.select_dtypes("number").columns)


def _embeddings(adata) -> list:
    return list(adata.obsm.keys())


def _default_basis(adata):
    embs = _embeddings(adata)
    if not embs:
        return None
    for pref in ("X_umap", "umap", "X_tsne", "X_pca"):
        for e in embs:
            if pref.lower() in e.lower():
                return e
    return embs[0]


_COMMON_MARKERS = [
    "MS4A1", "IL7R", "CD8A", "FCER1A", "CD14", "MS4A7", "NKG7", "GNLY",
    "CD3D", "CD79A", "LYZ", "PPBP", "FCGR3A", "CD4", "Gata1", "Gata2",
]


def _default_markers(adata, n: int = 8) -> list:
    genes = _genes(adata)
    present = [g for g in _COMMON_MARKERS if g in set(genes)]
    return (present or genes)[:n]


def _coerce(value, options, fallback=None):
    """Keep ``value`` only if it's a valid option, else fall back."""
    if value in options:
        return value
    if fallback is not None and fallback in options:
        return fallback
    return options[0] if options else None


# --------------------------------------------------------------------------- #
# Panels
# --------------------------------------------------------------------------- #
# Shown by a gene-driven panel's `.plot` when no gene is selected.
_NO_GENES_MSG = "Select at least one gene."


def embedding(adata, *, color: str | None = None, basis: str | None = None, size: int = 5):
    """Control bar for an embedding scatter (``gc.pl.embedding``).

    Controls: colour (obs column or gene), embedding, marker size, datashader.
    """
    mo = _mo()
    obs_cols, genes, embs = _obs_cols(adata), _genes(adata), _embeddings(adata)
    color_opts = obs_cols + genes
    color = _coerce(color, color_opts, obs_cols[0] if obs_cols else (genes[0] if genes else None))
    basis = _coerce(basis, embs, _default_basis(adata))

    controls = {
        "color": mo.ui.dropdown(color_opts, value=color, label="Color"),
        "basis": mo.ui.dropdown(embs, value=basis, label="Embedding"),
        "size": mo.ui.slider(1, 15, value=size, label="Marker size"),
        "datashader": mo.ui.checkbox(value=False, label="Datashader"),
    }

    def build(v):
        from guanaco import widget as pl
        return pl.embedding(
            adata, basis=v["basis"], color=v["color"], size=v["size"],
            render_backend="datashader" if v["datashader"] else "scattergl",
            return_fig=True,
        )

    return _panel(controls, build)


def violin(adata, *, keys=None, groupby: str | None = None):
    """Control bar for stacked violins (``gc.pl.violin``).

    Controls: genes (multi), group-by, box overlay.
    """
    mo = _mo()
    genes, disc = _genes(adata), _discrete_obs(adata)
    default_genes = [g for g in (keys or _default_markers(adata, 7)) if g in set(genes)]
    groupby = _coerce(groupby, disc)

    controls = {
        "genes": mo.ui.multiselect(genes, value=default_genes, max_selections=15, label="Genes"),
        "groupby": mo.ui.dropdown(disc, value=groupby, label="Group by"),
        "show_box": mo.ui.checkbox(value=True, label="Median + IQR box"),
    }

    def build(v):
        from guanaco import widget as pl
        if not v["genes"]:
            return _NO_GENES_MSG
        return pl.violin(
            adata, keys=v["genes"], groupby=v["groupby"], show_box=v["show_box"],
            return_fig=True,
        )

    return _panel(controls, build)


def dotplot(adata, *, var_names=None, groupby: str | None = None):
    """Control bar for a dotplot (``gc.pl.dotplot``).

    Controls: genes (multi), group-by, standardization (None / Min-max / Z-score).
    """
    mo = _mo()
    genes, disc = _genes(adata), _discrete_obs(adata)
    default_genes = [g for g in (var_names or _default_markers(adata, 8)) if g in set(genes)]
    groupby = _coerce(groupby, disc)

    controls = {
        "genes": mo.ui.multiselect(genes, value=default_genes, max_selections=30, label="Genes"),
        "groupby": mo.ui.dropdown(disc, value=groupby, label="Group by"),
        "standardize": mo.ui.dropdown(["none", "minmax", "zscore"], value="minmax", label="Standardize"),
    }

    def build(v):
        from guanaco import widget as pl
        if not v["genes"]:
            return _NO_GENES_MSG
        std = None if v["standardize"] == "none" else v["standardize"]
        return pl.dotplot(
            adata, var_names=v["genes"], groupby=v["groupby"],
            standardization=std, return_fig=True,
        )

    return _panel(controls, build)


def heatmap(adata, *, var_names=None, groupby: str | None = None):
    """Control bar for an expression heatmap (GUANACO's unified heatmap).

    Controls: genes (multi), group-by, optional second group, standardization.
    """
    mo = _mo()
    genes, disc = _genes(adata), _discrete_obs(adata)
    default_genes = [g for g in (var_names or _default_markers(adata, 10)) if g in set(genes)]
    groupby = _coerce(groupby, disc)

    controls = {
        "genes": mo.ui.multiselect(genes, value=default_genes, max_selections=50, label="Genes"),
        "groupby": mo.ui.dropdown(disc, value=groupby, label="Group by"),
        "groupby2": mo.ui.dropdown(["(none)"] + disc, value="(none)", label="Second group"),
        "standardize": mo.ui.dropdown(["none", "minmax", "zscore"], value="minmax", label="Standardize"),
    }

    def build(v):
        from guanaco.pages.matrix.plots.heatmap import plot_unified_heatmap
        from guanaco.widget import _default_palette
        if not v["genes"]:
            return _NO_GENES_MSG
        g2 = None if v["groupby2"] == "(none)" else v["groupby2"]
        std = None if v["standardize"] == "none" else v["standardize"]
        return plot_unified_heatmap(
            adata=adata, genes=v["genes"], groupby1=v["groupby"], groupby2=g2,
            standardization=std, adata_obs=adata.obs, color_config=_default_palette(),
        )

    return _panel(controls, build)


def stacked_bar(adata, *, x: str | None = None, color: str | None = None):
    """Control bar for a stacked composition bar chart (``gc.pl.stacked_bar``).

    Controls: x-axis group, stack-by colour, normalization.
    """
    mo = _mo()
    disc = _discrete_obs(adata)
    x = _coerce(x, disc)
    color = _coerce(color, disc, disc[1] if len(disc) > 1 else (disc[0] if disc else None))

    controls = {
        "x": mo.ui.dropdown(disc, value=x, label="X-axis (group)"),
        "color": mo.ui.dropdown(disc, value=color, label="Stack by (color)"),
        "normalize": mo.ui.dropdown(["proportion", "count"], value="proportion", label="Normalize"),
    }

    def build(v):
        from guanaco import widget as pl
        fig = pl.stacked_bar(
            adata, x=v["x"], color=v["color"], normalize=v["normalize"], return_fig=True,
        )
        # marimo's inline Plotly can fall back to grouped bars for multi-series bar
        # charts; render via Plotly's own HTML so barmode="stack" is always honored.
        return mo.iframe(fig.to_html(include_plotlyjs="cdn", default_height="520px"), height="560px")

    return _panel(controls, build)


def pseudotime(adata, *, genes=None, pseudotime_key: str | None = None, groupby: str | None = None):
    """Control bar for a gene-trend-along-pseudotime plot (``gc.pl.pseudotime``).

    Controls: genes (multi), pseudotime obs column, optional group-by.
    """
    mo = _mo()
    all_genes, disc, numeric = _genes(adata), _discrete_obs(adata), _numeric_obs(adata)
    default_genes = [g for g in (genes or _default_markers(adata, 4)) if g in set(all_genes)]
    pt = _coerce(pseudotime_key, numeric, "dpt_pseudotime")

    controls = {
        "genes": mo.ui.multiselect(all_genes, value=default_genes, max_selections=8, label="Genes"),
        "pseudotime_key": mo.ui.dropdown(numeric, value=pt, label="Pseudotime (obs)"),
        "groupby": mo.ui.dropdown(["(none)"] + disc, value="(none)", label="Group by"),
    }

    def build(v):
        from guanaco import widget as pl
        if not v["genes"]:
            return _NO_GENES_MSG
        gb = None if v["groupby"] == "(none)" else v["groupby"]
        return pl.pseudotime(
            adata, genes=v["genes"], pseudotime_key=v["pseudotime_key"],
            groupby=gb, return_fig=True,
        )

    return _panel(controls, build)


def peak_browser(adata, *, region: str | None = None, groupby: str | None = None, gene_annotation: str | None = None):
    """Control bar for the ATAC Peak Browser (``gc.pl.peak_browser``).

    Controls: search box (a **gene name** or a locus), group-by, metric (mean /
    detection), and y-axis scale (``shared`` — same range on every track for easy
    comparison, the default — or ``auto``). Needs peak-like ``var_names``
    (``"chr1:10000-10500"``) or ``var[['chrom','start','end']]``. Pass
    ``gene_annotation`` (a genome id like ``"hg38"`` or a GTF/GFF path) to add a
    gene-model track and enable gene search.
    """
    mo = _mo()
    disc = _discrete_obs(adata)
    groupby = _coerce(groupby, disc) if groupby else (disc[0] if disc else None)

    # Seed the locus box with a populated default window.
    region_text = region
    if region_text is None:
        try:
            from guanaco.pages.matrix.plots.atac_browser import default_region, format_locus

            r = default_region(adata)
            region_text = format_locus(str(r["chrom"]), int(r["start"]), int(r["end"]))
        except Exception:
            region_text = "chr1:1,000,000-2,000,000"

    controls = {
        "region": mo.ui.text(value=region_text, label="Gene or locus"),
        "groupby": mo.ui.dropdown(["(none)"] + disc, value=groupby or "(none)", label="Group by"),
        "metric": mo.ui.dropdown(["mean", "detection"], value="mean", label="Metric"),
        "y_mode": mo.ui.dropdown(["shared", "auto"], value="shared", label="Y-axis"),
    }

    def build(v):
        from guanaco import widget as pl

        gb = None if v["groupby"] == "(none)" else v["groupby"]
        try:
            return pl.peak_browser(
                adata,
                region=v["region"],
                groupby=gb,
                metric=v["metric"],
                y_mode=v["y_mode"],
                gene_annotation=gene_annotation,
                return_fig=True,
            )
        except ValueError as exc:
            return str(exc)

    return _panel(controls, build)


def volcano(adata, *, group: str | None = None):
    """Control bar for a volcano plot from precomputed DE in ``adata.uns`` (``gc.pl.volcano``).

    Controls: comparison, padj threshold, effect-size threshold, top-N labels.
    """
    mo = _mo()
    try:
        from guanaco.pages.matrix.plots.volcano import load_volcano_payload
        groups = list(load_volcano_payload(adata)["entries"].keys())
    except Exception:
        groups = []
    group = _coerce(group, groups) if groups else group

    controls = {
        "group": mo.ui.dropdown(groups, value=group, label="Comparison"),
        "padj_threshold": mo.ui.slider(0.0, 0.1, value=0.05, step=0.005, label="padj threshold"),
        "x_threshold": mo.ui.slider(0.0, 10.0, value=1.0, step=0.5, label="effect-size threshold"),
        "top_n": mo.ui.slider(0, 50, value=12, step=1, label="Top N labels"),
    }

    def build(v):
        from guanaco import widget as pl
        if not v["group"]:
            return "No differential-expression results found in adata.uns."
        return pl.volcano(
            adata, group=v["group"], padj_threshold=v["padj_threshold"],
            x_threshold=v["x_threshold"], top_n=int(v["top_n"]), return_fig=True,
        )

    return _panel(controls, build)
