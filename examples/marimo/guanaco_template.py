import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import scanpy as sc
    import guanaco as gc

    return gc, mo, sc


@app.cell
def _(mo):
    mo.md(
        """
        # GUANACO × marimo — single-cell explorer

        A **reactive** dashboard built on GUANACO's `gc.pl` API: change any control
        and the plot beneath it updates automatically — no re-running cells.

        **Use it on your own data:** edit `DATA_PATH` in the next cell, then run all
        cells. Every control rebuilds itself from your `obs` columns and genes.

        **Run & share:**
        ```bash
        pip install marimo guanaco-sc
        marimo edit guanaco_template.py                          # interactive editor
        marimo run  guanaco_template.py                          # serve as a web app
        marimo export html-wasm guanaco_template.py -o app.html  # self-contained file
        ```
        """
    )
    return


@app.cell
def _(sc):
    from pathlib import Path

    # ─────────────────────────────────────────────────────────────────────────
    #  EDIT THIS — path to your .h5ad file.
    #  (Defaults to the bundled demo so the notebook runs as-is.)
    # ─────────────────────────────────────────────────────────────────────────
    DATA_PATH = Path(__file__).resolve().parents[3] / "data" / "volcano_test.h5ad"
    # e.g. DATA_PATH = "/path/to/your_dataset.h5ad"

    adata = sc.read_h5ad(DATA_PATH)
    adata
    return DATA_PATH, adata


@app.cell
def _(adata):
    # Option lists, derived once from the dataset so every control adapts to it.
    obs_cols = list(adata.obs.columns)
    genes = list(adata.var_names)
    bases = [k for k in adata.obsm if k.lower().startswith("x_")] or list(adata.obsm)
    color_options = obs_cols + genes
    return bases, color_options, genes, obs_cols


@app.cell
def _(mo):
    mo.md("""## 1. Embedding""")
    return


@app.cell
def _(bases, color_options, mo):
    color = mo.ui.dropdown(color_options, value=color_options[0], label="Color by (obs column or gene)")
    basis = mo.ui.dropdown(bases, value=bases[0], label="Embedding")
    size = mo.ui.slider(1, 15, value=5, label="Marker size")
    # Datashader rasterizes points server-side -> fast for very large datasets.
    # (Needs the `datashader` package; falls back to scattergl if it's missing.)
    ds_mode = mo.ui.checkbox(value=False, label="Datashader (fast for large data)")
    mo.hstack([color, basis, size, ds_mode])
    return basis, color, ds_mode, size


@app.cell
def _(adata, basis, color, ds_mode, gc, size):
    gc.pl.embedding(
        adata,
        basis=basis.value,
        color=color.value,
        size=size.value,
        render_backend="datashader" if ds_mode.value else "scattergl",
        return_fig=True,
    )
    return


@app.cell
def _(mo):
    mo.md("""## 2. Co-expression  ·  *two genes on the embedding*""")
    return


@app.cell
def _(bases, genes, mo):
    ce_g1 = mo.ui.dropdown(genes, value=genes[0], label="Gene 1")
    ce_g2 = mo.ui.dropdown(genes, value=genes[1] if len(genes) > 1 else genes[0], label="Gene 2")
    ce_basis = mo.ui.dropdown(bases, value=bases[0], label="Embedding")
    ce_t1 = mo.ui.slider(0.0, 5.0, value=1.0, step=0.5, label="Gene 1 threshold")
    ce_t2 = mo.ui.slider(0.0, 5.0, value=1.0, step=0.5, label="Gene 2 threshold")
    mo.hstack([ce_g1, ce_g2, ce_basis, ce_t1, ce_t2])
    return ce_basis, ce_g1, ce_g2, ce_t1, ce_t2


@app.cell
def _(adata, ce_basis, ce_g1, ce_g2, ce_t1, ce_t2, gc):
    gc.pl.coexpression(
        adata, basis=ce_basis.value, gene1=ce_g1.value, gene2=ce_g2.value,
        threshold1=ce_t1.value, threshold2=ce_t2.value, return_fig=True,
    )
    return


@app.cell
def _(mo):
    mo.md("""## 3. Violin  ·  *several genes by group*""")
    return


@app.cell
def _(genes, mo, obs_cols):
    # max_selections caps the count and removes the "select all" option.
    v_genes = mo.ui.multiselect(genes, value=genes[:3], max_selections=15, label="Genes")
    v_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
    v_box = mo.ui.checkbox(value=True, label="Median + IQR box")
    mo.hstack([v_genes, v_group, v_box])
    return v_box, v_genes, v_group


@app.cell
def _(adata, gc, v_box, v_genes, v_group):
    (
        gc.pl.violin(
            adata, keys=v_genes.value, groupby=v_group.value,
            show_box=v_box.value, return_fig=True,
        )
        if v_genes.value
        else "Select at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md("""## 4. Dotplot""")
    return


@app.cell
def _(genes, mo, obs_cols):
    d_genes = mo.ui.multiselect(genes, value=genes[:5], max_selections=30, label="Genes")
    d_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
    d_std = mo.ui.dropdown(["none", "var", "group"], value="var", label="Standardize")
    mo.hstack([d_genes, d_group, d_std])
    return d_genes, d_group, d_std


@app.cell
def _(adata, d_genes, d_group, d_std, gc):
    (
        gc.pl.dotplot(
            adata, var_names=d_genes.value, groupby=d_group.value,
            standardization=d_std.value, return_fig=True,
        )
        if d_genes.value
        else "Select at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md("""## 5. Heatmap""")
    return


@app.cell
def _(genes, mo, obs_cols):
    h_genes = mo.ui.multiselect(genes, value=genes[:8], max_selections=50, label="Genes")
    h_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
    h_group2 = mo.ui.dropdown(["(none)"] + obs_cols, value="(none)", label="Second group (optional)")
    h_zscore = mo.ui.checkbox(value=True, label="z-score per gene")
    mo.hstack([h_genes, h_group, h_group2, h_zscore])
    return h_genes, h_group, h_group2, h_zscore


@app.cell
def _(adata, gc, h_genes, h_group, h_group2, h_zscore):
    (
        gc.pl.heatmap(
            adata, var_names=h_genes.value, groupby=h_group.value,
            groupby2=None if h_group2.value == "(none)" else h_group2.value,
            z_score=h_zscore.value, return_fig=True,
        )
        if h_genes.value
        else "Select at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md("""## 6. Stacked bar  ·  *composition*""")
    return


@app.cell
def _(mo, obs_cols):
    sb_x = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="X-axis (group)")
    sb_color = mo.ui.dropdown(
        obs_cols, value=obs_cols[1] if len(obs_cols) > 1 else obs_cols[0], label="Stack by (color)"
    )
    sb_norm = mo.ui.dropdown(["proportion", "count"], value="proportion", label="Normalize")
    mo.hstack([sb_x, sb_color, sb_norm])
    return sb_color, sb_norm, sb_x


@app.cell
def _(adata, gc, mo, sb_color, sb_norm, sb_x):
    _fig = gc.pl.stacked_bar(
        adata, x=sb_x.value, color=sb_color.value, normalize=sb_norm.value, return_fig=True
    )
    # Render through Plotly's own HTML (in an iframe) so the stacked layout
    # (barmode="stack") is always honored. marimo's inline Plotly display can
    # otherwise fall back to grouped bars for multi-series bar charts.
    mo.iframe(_fig.to_html(include_plotlyjs="cdn", default_height="520px"), height="560px")
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 7. Gene trend along pseudotime

        Needs a **numeric** pseudotime column in `obs` (e.g. `dpt_pseudotime`).
        """
    )
    return


@app.cell
def _(adata, genes, mo, obs_cols):
    pt_numeric = list(adata.obs.select_dtypes("number").columns)
    if pt_numeric:
        _def_pt = next(
            (c for c in pt_numeric if "pseudotime" in c.lower() or "dpt" in c.lower()),
            pt_numeric[0],
        )
        pt_genes = mo.ui.multiselect(genes, value=genes[:2], max_selections=8, label="Genes")
        pt_key = mo.ui.dropdown(pt_numeric, value=_def_pt, label="Pseudotime (numeric obs)")
        pt_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
        _controls = mo.hstack([pt_genes, pt_key, pt_group])
    else:
        pt_genes = pt_key = pt_group = None
        _controls = mo.md("*No numeric `obs` column found — add e.g. `dpt_pseudotime` to use this panel.*")
    _controls
    return pt_genes, pt_group, pt_key


@app.cell
def _(adata, gc, pt_genes, pt_group, pt_key):
    (
        gc.pl.pseudotime(
            adata, genes=pt_genes.value, pseudotime_key=pt_key.value,
            groupby=pt_group.value, return_fig=True,
        )
        if (pt_genes is not None and pt_genes.value)
        else "Pseudotime needs a numeric obs column and at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 8. Volcano (differential expression)

        Needs **precomputed** DE in `adata.uns` (run `sc.tl.rank_genes_groups` first).
        """
    )
    return


@app.cell
def _(adata):
    # Volcano reads precomputed DE (rank_genes_groups or a 'volcano' entry).
    try:
        from guanaco.pages.matrix.plots.volcano import load_volcano_payload

        vol_groups = list(load_volcano_payload(adata)["entries"].keys())
    except Exception:
        vol_groups = []
    return (vol_groups,)


@app.cell
def _(mo, vol_groups):
    if vol_groups:
        vol_group = mo.ui.dropdown(vol_groups, value=vol_groups[0], label="Comparison")
        vol_padj = mo.ui.slider(0.0, 0.1, value=0.05, step=0.005, label="padj threshold")
        vol_x = mo.ui.slider(0.0, 10.0, value=1.0, step=0.5, label="effect-size threshold")
        vol_topn = mo.ui.slider(0, 50, value=12, step=1, label="Top N labels")
        _controls = mo.hstack([vol_group, vol_padj, vol_x, vol_topn])
    else:
        vol_group = vol_padj = vol_x = vol_topn = None
        _controls = mo.md("*No DE found in `adata.uns` — run `sc.tl.rank_genes_groups` first.*")
    _controls
    return vol_group, vol_padj, vol_topn, vol_x


@app.cell
def _(adata, gc, vol_group, vol_padj, vol_topn, vol_x):
    (
        gc.pl.volcano(
            adata, group=vol_group.value, padj_threshold=vol_padj.value,
            x_threshold=vol_x.value, top_n=vol_topn.value, return_fig=True,
        )
        if vol_group is not None
        else "Volcano needs precomputed DE in adata.uns."
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---
        **Notes**
        - Same `gc.pl` API as the Jupyter notebook; see `docs/notebook_api.md` for every
          function and parameter, and `docs/marimo_guide.md` for this app.
        - PAGA and GRN are Cytoscape networks (not Plotly), so they live in the full
          web app (`guanaco -c config.json`), not here.
        """
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
