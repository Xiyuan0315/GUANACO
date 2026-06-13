import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md("""
    # GUANACO × marimo — interactive plot panels

    `guanaco.marimo` gives you a **plot + its control bar in one line**. Each panel
    builds dropdowns / sliders from your `adata`, wired to a GUANACO plot.

    **Install & run**

    ```bash
    pip install "guanaco-viz" marimo
    marimo edit guanaco_panels.py     # interactive editor
    # or:  marimo run guanaco_panels.py   (shareable read-only app)
    ```

    **The pattern.** marimo reactivity is *cell-based*, so each plot is **two tiny
    cells** — build the control bar in one, render the figure in the next:

    ```python
    dp = gm.dotplot(adata)   # cell 1 — the control bar
    dp.hstack()
    dp.plot                  # cell 2 — the figure (redraws when you change a control)
    ```

    Touch a control above and the figure below updates automatically.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import scanpy as sc

    import guanaco.marimo as gm

    return gm, mo, sc


@app.cell
def _(sc):
    from pathlib import Path

    # Bundled demo data, resolved relative to this file.
    data_dir = Path(__file__).resolve().parents[3] / "data"
    adata = sc.read_h5ad(data_dir / "volcano_test.h5ad")
    adata
    return adata, data_dir


@app.cell
def _(mo):
    mo.md("""
    ## 1 · Embedding

    Colour by any `obs` column or gene; switch embedding, marker size, or turn on
    datashader for very large datasets.
    """)
    return


@app.cell
def _(adata, gm):
    emb = gm.embedding(adata)
    emb.hstack()
    return (emb,)


@app.cell
def _(emb):
    emb.plot
    return


@app.cell
def _(mo):
    mo.md("""
    ## 2 · Dotplot

    Dot colour = mean expression, dot size = fraction of cells expressing.
    **Standardize** = None / Min–max (0–1) / Z-score.
    """)
    return


@app.cell
def _(adata, gm):
    dot = gm.dotplot(adata)
    dot.hstack()
    return (dot,)


@app.cell
def _(dot):
    dot.plot
    return


@app.cell
def _(mo):
    mo.md("""
    ## 3 · Violin

    Stacked violins, one row per gene. Toggle the median/IQR box.
    """)
    return


@app.cell
def _(adata, gm):
    vio = gm.violin(adata)
    vio.hstack()
    return (vio,)


@app.cell
def _(vio):
    vio.plot
    return


@app.cell
def _(mo):
    mo.md("""
    ## 4 · Heatmap

    Optionally add a second grouping; the same None / Min–max / Z-score
    standardization as the dotplot.
    """)
    return


@app.cell
def _(adata, gm):
    hm = gm.heatmap(adata)
    hm.hstack()
    return (hm,)


@app.cell
def _(hm):
    hm.plot
    return


@app.cell
def _(mo):
    mo.md("""
    ## 5 · Stacked bar (composition)

    Cell-type composition across a grouping; proportion or raw counts.
    """)
    return


@app.cell
def _(adata, gm):
    sb = gm.stacked_bar(adata)
    sb.hstack()
    return (sb,)


@app.cell
def _(sb):
    sb.plot
    return


@app.cell
def _(mo):
    mo.md("""
    ## 6 · Pseudotime trend

    Smoothed gene-expression along a pseudotime `obs` column (uses a separate
    trajectory dataset).
    """)
    return


@app.cell
def _(data_dir, sc):
    adata_traj = sc.read_h5ad(data_dir / "paul15_processed_v2.h5ad")
    return (adata_traj,)


@app.cell
def _(adata_traj, gm):
    pt = gm.pseudotime(adata_traj)
    pt.hstack()
    return (pt,)


@app.cell
def _(pt):
    pt.plot
    return


@app.cell
def _(mo):
    mo.md("""
    ## 7 · Volcano

    From precomputed differential expression in `adata.uns` (uses a bulk dataset).
    """)
    return


@app.cell
def _(data_dir, sc):
    bulk = sc.read_h5ad(data_dir / "kich_demo.h5ad")
    return (bulk,)


@app.cell
def _(bulk, gm):
    vol = gm.volcano(bulk)
    vol.hstack()
    return (vol,)


@app.cell
def _(vol):
    vol.plot
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ### Tips

    - **Seed the defaults:** every panel takes optional arguments, e.g.
      `gm.dotplot(adata, groupby="cell_type")` or
      `gm.embedding(adata, color="CD8A", basis="X_umap")`.
    - **Function form:** `gm.show(panel)` is the same as `panel.plot`.
    - **Read the values:** `panel.value` is a plain dict of the current control
      settings — handy if you want to branch on them.
    - A `panel` is a real `mo.ui` element: display it as `panel` (labelled column)
      or `panel.hstack()` (single row).
    """)
    return


if __name__ == "__main__":
    app.run()
