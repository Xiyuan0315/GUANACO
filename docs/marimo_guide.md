# Using GUANACO in marimo

[marimo](https://marimo.io) is a reactive Python notebook: when you change a
control (a dropdown, a slider), every cell that depends on it re-runs
automatically. Combined with GUANACO's scanpy-style plotting API (`gc.pl.*`),
that gives you a live single-cell dashboard you can build in a few lines, run as
a web app, and share as a single file — no server, no callbacks to wire up.

This guide covers the two ready-made notebooks in `examples/marimo/` and how to
adapt them to your own data.

---

## 1. Install & run

```bash
pip install marimo guanaco-sc          # marimo + GUANACO
pip install datashader                 # optional: fast rendering for large data
```

Three ways to use a notebook (e.g. `guanaco_template.py`):

```bash
marimo edit guanaco_template.py                          # interactive editor
marimo run  guanaco_template.py                          # serve as a read-only web app
marimo export html-wasm guanaco_template.py -o app.html  # one self-contained file
```

- **`marimo edit`** — develop and explore; you can see and change the code.
- **`marimo run`** — share a clean app (controls + plots, code hidden) over a URL.
- **`export html-wasm`** — produces a single `.html` that runs entirely in the
  browser (via WebAssembly). Host it anywhere static (GitHub Pages, a shared
  drive) and send your supervisor a link — no Python install on their side.

> Unlike Jupyter, a marimo notebook **is** a plain `.py` file (no stored
> outputs), so it diffs cleanly in git and there's no stale-output confusion.

---

## 2. The two notebooks

| File | What it's for |
|------|---------------|
| `examples/marimo/guanaco_template.py` | **Start here.** A clean template — point it at your own `.h5ad` and go. One data path, every control adapts to your dataset. |
| `examples/marimo/guanaco_explorer.py` | A fuller demo with **file upload** (drag-and-drop an `.h5ad`, up to 2 GB) and **lasso cross-filtering** (select cells on the UMAP to filter the violin/dotplot/heatmap below). |

Both drive the exact same `gc.pl` API documented in
[`notebook_api.md`](notebook_api.md).

---

## 3. Adapt the template to your data

Open `guanaco_template.py` and edit one line:

```python
DATA_PATH = Path(__file__).resolve().parents[3] / "data" / "volcano_test.h5ad"
# →
DATA_PATH = "/path/to/your_dataset.h5ad"
```

Run all cells. The notebook reads your `AnnData` once and builds every dropdown
from it:

```python
obs_cols = list(adata.obs.columns)                       # group-by / color options
genes    = list(adata.var_names)                          # gene options
bases    = [k for k in adata.obsm if k.startswith("X_")]  # embeddings (UMAP, PCA, …)
```

That's the whole pattern: **derive option lists from `adata`, feed them to
`mo.ui` controls, pass `control.value` to `gc.pl`.** To add a panel, copy a
controls cell + a plot cell and swap the `gc.pl` call.

### Panels included

| # | Panel | `gc.pl` call | Needs |
|---|-------|--------------|-------|
| 1 | Embedding | `gc.pl.embedding` | an `obsm` embedding |
| 2 | Co-expression | `gc.pl.coexpression` | two genes |
| 3 | Violin | `gc.pl.violin` | genes + a group column |
| 4 | Dotplot | `gc.pl.dotplot` | genes + a group column |
| 5 | Heatmap | `gc.pl.heatmap` | genes + a group column |
| 6 | Stacked bar | `gc.pl.stacked_bar` | two `obs` columns |
| 7 | Pseudotime | `gc.pl.pseudotime` | a numeric `obs` (e.g. `dpt_pseudotime`) |
| 8 | Volcano | `gc.pl.volcano` | precomputed DE in `adata.uns` |
| 9 | Peak browser | `gc.pl.peak_browser` | peak-like `var_names` (`chr1:10000-10500`) — ATAC data |

Panels 7, 8 and 9 **adapt to your data**: if your dataset has no numeric pseudotime
column, no DE in `uns`, or no genomic peaks, the panel shows a short note instead of
erroring, so the rest of the dashboard still works.

---

## 4. How a reactive panel works

A panel is just two cells — controls, then a plot that reads `control.value`:

```python
# controls
v_genes = mo.ui.multiselect(genes, value=genes[:3], max_selections=15, label="Genes")
v_group = mo.ui.dropdown(obs_cols, value=obs_cols[0], label="Group by")
mo.hstack([v_genes, v_group])

# plot — re-runs automatically whenever a control above changes
gc.pl.violin(adata, keys=v_genes.value, groupby=v_group.value, return_fig=True)
```

Notes that matter in practice:

- **Always pass `return_fig=True`.** It returns the Plotly figure as the cell's
  last expression, which marimo renders. (Without it, `gc.pl` would call
  `fig.show()`, meant for Jupyter.)
- **`max_selections`** on gene multiselects removes the "select all" option —
  important when there are thousands of genes.
- **Stacked bar uses an iframe.** marimo's inline Plotly display can render a
  multi-series bar chart as *grouped* instead of *stacked*, so the template
  renders it via Plotly's own HTML:
  ```python
  _fig = gc.pl.stacked_bar(adata, x=..., color=..., return_fig=True)
  mo.iframe(_fig.to_html(include_plotlyjs="cdn", default_height="520px"), height="560px")
  ```
  This needs internet (it loads `plotly.js` from a CDN); use
  `include_plotlyjs=True` to inline it for fully offline use.

---

## 5. Large datasets — datashader

The embedding panel has a **Datashader** checkbox. Datashader rasterizes the
points server-side, so the figure stays fast and small no matter how many cells
you have:

```python
gc.pl.embedding(adata, basis=..., color=..., render_backend="datashader", return_fig=True)
```

- Needs `pip install datashader`; if it's missing, GUANACO falls back to the
  normal WebGL scatter automatically.
- A datashader plot is a single rasterized image, so it has **no per-point hover
  or lasso** — turn it off when you need the `guanaco_explorer.py` lasso
  selection.

**Color is auto-detected.** A gene or a numeric `obs` column with many distinct
values (e.g. `n_counts`, `n_genes`) is drawn with a continuous colormap; a
categorical column (or a low-cardinality integer like cluster ids) gets discrete
colors. This is what keeps a continuous covariate from being (slowly) drawn as
thousands of categories.

---

## 6. File upload (explorer notebook)

`guanaco_explorer.py` lets users drag-and-drop their own `.h5ad`:

```python
upload = mo.ui.file(filetypes=[".h5ad"], kind="area", max_size=2_000_000_000)  # 2 GB
```

marimo's default cap is **100 MB**; the explorer raises it to **2 GB**. Uploads
travel through the kernel in memory (base64 over the websocket), so very large
files are slow this way — for big data, prefer the template's `DATA_PATH`
(read directly from disk).

---

## 7. Sharing

- **A link, app only:** `marimo run guanaco_template.py` and share the URL (code
  hidden, controls live). Needs a running Python process.
- **A single file, no server:** `marimo export html-wasm guanaco_template.py -o app.html`
  and host the `.html` statically. It runs in the browser via WebAssembly — ideal
  for sending to a supervisor or attaching to a paper.

---

## 8. PAGA & GRN

PAGA and GRN are interactive **Cytoscape** networks, not Plotly figures. They are
available as `gc.pl.paga(adata)` / `gc.pl.grn(adata)` in **Jupyter/JupyterLab**
(via `pip install ipycytoscape`), and of course in the full web app
(`guanaco -c config.json`).

They are **not** currently available in marimo — they're classic ipywidgets with
their own JS extension, which marimo doesn't render. So in a marimo dashboard,
use the web app for PAGA/GRN, or open those two in a Jupyter notebook.

Everything else — embeddings, co-expression, violin, dotplot, heatmap, stacked
bar, pseudotime, volcano — works identically in marimo, Jupyter, and the web app,
because they all call the same `gc.pl` functions.

See also: [`notebook_api.md`](notebook_api.md) for the full `gc.pl` reference.
