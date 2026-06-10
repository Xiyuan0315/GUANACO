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
    mo.md("""
    # GUANACO × marimo — interactive single-cell explorer

    A **reactive** panel built on GUANACO's `gc.pl` API
    """)
    return


@app.cell
def _(mo):
    upload = mo.ui.file(
        filetypes=[".h5ad"],
        kind="area",
        max_size=2_000_000_000,  # 2 GB (marimo's default cap is 100 MB)
        label="Optional: drop an .h5ad here (otherwise the bundled demo is used)",
    )
    upload
    return (upload,)


@app.cell
def _(sc, upload):
    from pathlib import Path

    # Bundled demo data, resolved relative to this file so it works from any cwd.
    data_dir = Path(__file__).resolve().parents[3] / "data"

    if upload.value:
        # Read the uploaded .h5ad from memory (write to a temp file for h5py).
        import tempfile
        import anndata as ad

        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as _tmp:
            _tmp.write(upload.value[0].contents)
            _uploaded_path = _tmp.name
        adata = ad.read_h5ad(_uploaded_path)
    else:
        adata = sc.read_h5ad(data_dir / "volcano_test.h5ad")
    return adata, data_dir


@app.cell
def _(adata):
    # Option lists, built once from the dataset.
    obs_cols = list(adata.obs.columns)
    genes = list(adata.var_names)
    bases = [k for k in adata.obsm.keys()] or list(adata.obsm.keys())
    color_options = obs_cols + genes
    return bases, color_options, genes, obs_cols


@app.cell
def _(mo):
    mo.md("""
    ## 1. Embedding explorer  ·  *box/lasso-select cells to filter the plots below*
    """)
    return


@app.cell
def _(bases, color_options, mo):
    color = mo.ui.dropdown(color_options, value=color_options[0], label="Color by (obs column or gene)")
    basis = mo.ui.dropdown(bases, value=bases[1], label="Embedding")
    size = mo.ui.slider(1, 15, value=5, label="Marker size")
    # Datashader rasterizes points server-side -> fast for very large datasets,
    # but a raster has no per-point lasso, so selection is disabled while it's on.
    ds_mode = mo.ui.checkbox(value=False, label="Datashader")
    mo.hstack([color, basis, size, ds_mode])
    return basis, color, ds_mode, size


@app.cell
def _(adata, basis, color, ds_mode, gc, mo, size):
    fig = gc.pl.embedding(
        adata,
        basis=basis.value,
        color=color.value,
        size=size.value,
        render_backend="datashader" if ds_mode.value else "scattergl",
        return_fig=True,
    )
    # Default to lasso so you can just drag to select; switch to box-select (or
    # zoom/pan) in the toolbar that appears at the figure's top-right on hover.
    fig.update_layout(dragmode="lasso")
    # mo.ui.plotly captures the box/lasso selection into umap.value. In datashader
    # mode the plot is a single raster image, so there are no points to lasso.
    umap = mo.ui.plotly(fig)
    umap
    return fig, umap


@app.cell
def _(adata, fig, mo, umap):
    # marimo may return the selected points as a list, or wrapped in a dict.
    _sel = umap.value
    _points = (_sel.get("points") or _sel.get("selection") or []) if isinstance(_sel, dict) else (_sel or [])

    _names = set(adata.obs_names)

    def _cell_id(p):
        if not isinstance(p, dict):
            return None
        # 1) customdata carried on the point itself; 2) else look it up on the
        # figure trace by (curveNumber, pointNumber). gc.pl encodes identity as
        # an obs_name (categorical color) or a cell index (continuous color).
        cd = p.get("customdata")
        if cd is None:
            try:
                cd = fig.data[p["curveNumber"]].customdata[p.get("pointNumber", p.get("pointIndex"))]
            except Exception:
                return None
        first = cd[0] if (hasattr(cd, "__len__") and not isinstance(cd, str)) else cd
        if isinstance(first, str):
            return first if first in _names else None
        try:
            return adata.obs_names[int(first)]
        except (TypeError, ValueError, IndexError):
            return None

    _ids = [i for i in (_cell_id(p) for p in _points) if i is not None]
    sub = adata[_ids] if _ids else adata  # no/empty selection -> all cells

    # Diagnostic: if points came through but none mapped, show their keys so the
    # extractor can be adjusted to this marimo version.
    _debug = ""
    if _points and not _ids:
        _p0 = _points[0]
        _debug = f"  ·  ⚠ couldn't map — sample point: {sorted(_p0) if isinstance(_p0, dict) else type(_p0).__name__}"
    elif not _points and _sel not in (None, [], {}):
        _debug = f"  ·  ⚠ no points (value type: {type(_sel).__name__})"

    mo.md(
        f"**{sub.n_obs:,}** of **{adata.n_obs:,}** cells selected{_debug} — "
        f"lasso the UMAP to focus the plots below; clear it to use all cells."
    )
    return (sub,)


@app.cell
def _(mo):
    mo.md("""
    ## 2. Co-expression  ·  *two genes on the embedding (all cells)*
    """)
    return


@app.cell
def _(bases, genes, mo):
    ce_g1 = mo.ui.dropdown(genes, value='CD14', label="Gene 1")
    ce_g2 = mo.ui.dropdown(genes, value='FCGR3A' if len(genes) > 1 else genes[0], label="Gene 2")
    ce_basis = mo.ui.dropdown(bases, value=bases[1], label="Embedding")
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
    mo.md("""
    ## 3. Violin explorer  ·  *of the selected cells*
    """)
    return


@app.cell
def _():
    markers = ['MS4A1', 'IL7R', 'CD8A', 'FCER1A', 'CD14', 'MS4A7', 'NKG7']
    return (markers,)


@app.cell
def _(genes, markers, mo, obs_cols):
    # max_selections caps the count and removes the "select all" option
    # (gene lists are often thousands long).
    v_genes = mo.ui.multiselect(genes, value=markers, max_selections=15, label="Genes")
    v_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
    v_box = mo.ui.checkbox(value=True, label="Median + IQR box")
    mo.hstack([v_genes, v_group, v_box])
    return v_box, v_genes, v_group


@app.cell
def _(gc, sub, v_box, v_genes, v_group):
    # `sub` = the lassoed cells (or all cells when nothing is selected).
    (
        gc.pl.violin(
            sub, keys=v_genes.value, groupby=v_group.value,
            show_box=v_box.value, return_fig=True, transformation= 'log'
        )
        if v_genes.value
        else "Select at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## 4. Dotplot explorer  ·  *of the selected cells*
    """)
    return


@app.cell
def _(genes, markers, mo, obs_cols):
    d_genes = mo.ui.multiselect(genes, value=markers, max_selections=30, label="Genes")
    d_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
    d_std = mo.ui.dropdown(["none", "var", "group"], value="var", label="Standardize")
    mo.hstack([d_genes, d_group, d_std])
    return d_genes, d_group, d_std


@app.cell
def _(d_genes, d_group, d_std, gc, sub):
    # `sub` = the lassoed cells (or all cells when nothing is selected).
    (
        gc.pl.dotplot(
            sub, var_names=d_genes.value, groupby=d_group.value,
            standardization=d_std.value, return_fig=True,
        )
        if d_genes.value
        else "Select at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## 5. Heatmap  ·  *of the selected cells*
    """)
    return


@app.cell
def _(genes, markers, mo, obs_cols):
    h_genes = mo.ui.multiselect(genes, value=markers, max_selections=50, label="Genes")
    h_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
    h_group2 = mo.ui.dropdown(["(none)"] + obs_cols, value=obs_cols[1], label="Second group (optional)")
    h_zscore = mo.ui.checkbox(value=True, label="z-score per gene")
    mo.hstack([h_genes, h_group, h_group2, h_zscore])
    return h_genes, h_group, h_group2, h_zscore


@app.cell
def _(gc, h_genes, h_group, h_group2, h_zscore, sub):
    # `sub` = the lassoed cells (or all cells when nothing is selected).
    (
        gc.pl.heatmap(
            sub, var_names=h_genes.value, groupby=h_group.value,
            groupby2=None if h_group2.value == "(none)" else h_group2.value,
            z_score=h_zscore.value, return_fig=True
        )
        if h_genes.value
        else "Select at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## 6. Stacked bar  ·  *cell-type composition (all cells)*
    """)
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
    mo.md("""
    ## 7. Gene trend along pseudotime

    Uses a hematopoiesis trajectory dataset with a pseudotime column in `obs`.
    """)
    return


@app.cell
def _(data_dir, sc):
    adata_traj = sc.read_h5ad(data_dir / "paul15_processed_v2.h5ad")
    return (adata_traj,)


@app.cell
def _(adata_traj):
    adata_traj
    return


@app.cell
def _(adata_traj):
    traj_genes = list(adata_traj.var_names)
    traj_obs = list(adata_traj.obs.columns)

    # Pseudotime must be a numeric obs column.
    traj_pt_keys = list(adata_traj.obs.select_dtypes("number").columns)
    traj_bases = [k for k in adata_traj.obsm.keys()] or list(adata_traj.obsm.keys())
    traj_color_options = traj_genes + traj_obs
    return traj_bases, traj_color_options, traj_genes, traj_obs, traj_pt_keys


@app.cell
def _(mo, traj_bases, traj_color_options, traj_genes, traj_obs, traj_pt_keys):
    pt_genes = mo.ui.multiselect(
        traj_genes,
        value=[g for g in ["Gata1", "Gata2"] if g in traj_genes] or traj_genes[:2],
        max_selections=8,
        label="Genes",
    )
    pt_key = mo.ui.dropdown(
        traj_pt_keys,
        value="dpt_pseudotime" if "dpt_pseudotime" in traj_pt_keys else traj_pt_keys[0],
        label="Pseudotime (obs)",
    )
    pt_group = mo.ui.dropdown(
        traj_obs,
        value="cluster" if "cluster" in traj_obs else traj_obs[0],
        label="Group by",
    )
    t_color = mo.ui.dropdown(traj_color_options, value='cluster', label="Color by (obs column or gene)")
    t_basis = mo.ui.dropdown(traj_bases, value=traj_bases[1], label="Embedding")
    t_size = mo.ui.slider(1, 15, value=5, label="Marker size")
    t_ds_mode = mo.ui.checkbox(value=False, label="Datashader")
    mo.hstack([t_color, t_basis, t_size, t_ds_mode])
    return pt_genes, pt_group, pt_key, t_basis, t_color, t_ds_mode, t_size


@app.cell
def _(adata_traj, gc, t_basis, t_color, t_ds_mode, t_size):
    gc.pl.embedding(
        adata_traj,
        basis=t_basis.value,
        color=t_color.value,
        size=t_size.value,
        render_backend="datashader" if t_ds_mode.value else "scattergl",
        return_fig=True,
    )

    return


@app.cell
def _(mo, pt_genes, pt_group, pt_key):
    mo.hstack([pt_genes, pt_key, pt_group])
    return


@app.cell
def _(adata_traj, gc, pt_genes, pt_group, pt_key):
    (
        gc.pl.pseudotime(
            adata_traj, genes=pt_genes.value, pseudotime_key=pt_key.value,
            groupby=pt_group.value, return_fig=True,
        )
        if pt_genes.value
        else "Select at least one gene."
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## 8. Volcano (differential expression)

    Uses a bulk RNA-seq dataset with precomputed DE in `adata.uns`.
    """)
    return


@app.cell
def _(data_dir, sc):
    bulk = sc.read_h5ad(data_dir / "kich_demo.h5ad")
    return (bulk,)


@app.cell
def _(bulk):
    from guanaco.pages.matrix.plots.volcano import load_volcano_payload

    # Available comparisons for the volcano dropdown.
    vol_groups = list(load_volcano_payload(bulk)["entries"].keys())
    return (vol_groups,)


@app.cell
def _(mo, vol_groups):
    vol_group = mo.ui.dropdown(
        vol_groups,
        value="KICH: Tumor vs Normal" if "KICH: Tumor vs Normal" in vol_groups else vol_groups[0],
        label="Comparison",
    )
    vol_padj = mo.ui.slider(0.0, 0.1, value=0.05, step=0.005, label="padj threshold")
    vol_x = mo.ui.slider(0.0, 10.0, value=1.0, step=0.5, label="effect-size threshold")
    vol_topn = mo.ui.slider(0, 50, value=12, step=1, label="Top N labels")
    mo.hstack([vol_group, vol_padj, vol_x, vol_topn])
    return vol_group, vol_padj, vol_topn, vol_x


@app.cell
def _(bulk, gc, vol_group, vol_padj, vol_topn, vol_x):
    gc.pl.volcano(
        bulk,
        group=vol_group.value,
        padj_threshold=vol_padj.value,
        x_threshold=vol_x.value,
        top_n=vol_topn.value,
        return_fig=True,
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    **Notes**
    - This drives the same `gc.pl` API documented in `docs/notebook_api.md`; add more
      panels with `gc.pl.heatmap`, `gc.pl.matrixplot`, `gc.pl.stacked_bar`, etc.
    - `marimo run` serves this as a shareable web app; `marimo export html-wasm`
      produces a self-contained interactive app that runs in the browser.
    """)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
