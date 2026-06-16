# GUANACO notebook plotting API (`gc.pl`)

GUANACO ships a scanpy-style plotting API for interactive use in Jupyter
notebooks. If you know `scanpy`'s `sc.pl.*`, you already know this:

```python
import guanaco as gc

gc.pl.umap(adata, color="cell_type")
```

Every function takes an [AnnData](https://anndata.readthedocs.io) object and
returns an **interactive Plotly figure** (zoom, pan, hover) — except `gc.pl.paga`
and `gc.pl.grn`, which return interactive Cytoscape (`ipycytoscape`) networks
(see [PAGA & GRN](#paga--grn-cytoscape-networks)).

---

## Conventions (shared by every function)

| Argument | Default | Meaning |
|----------|---------|---------|
| `show` | `True` | Render the figure inline (calls `fig.show()`) and return `None`. |
| `return_fig` | `False` | Return the `plotly.graph_objects.Figure` instead of showing it. Use this to save (`fig.write_html(...)` / `fig.write_image(...)`) or post-process. |

Other recurring arguments:

- **`color` / `keys` / `var_names`** — a gene name (in `adata.var_names`) or an
  `adata.obs` column. `keys`/`var_names` also accept a **list**; `color` is a single value.
- **`groupby`** — an `adata.obs` categorical column to group cells by.
- **`color_map`** — a *continuous* colormap name (e.g. `"Viridis"`, `"magma"`, `"plasma"`).
- **`palette`** — a *categorical* palette: a palette **name** (e.g. `"tab10"`) or a **list of colors** (`["#E69F00", "#56B4E9", ...]`).
- **`transformation`** — `None` (raw), `"log"` (log1p), or `"z_score"` (z-scored, clipped). *(Per-plot support noted below.)*

> **Tip:** to save instead of display:
> ```python
> fig = gc.pl.umap(adata, color="CD8A", return_fig=True)
> fig.write_html("umap.html")     # interactive
> fig.write_image("umap.png")     # static (needs `kaleido`)
> ```

---

## Embedding / scatter

### `gc.pl.embedding(adata, basis, color, ...)`

Scatter of any embedding, colored by a gene or an `obs` column. Continuous vs.
categorical coloring is detected automatically.

```python
gc.pl.embedding(adata, basis="X_umap", color="CD8A")          # gene (continuous)
gc.pl.embedding(adata, basis="umap",   color="cell_type",     # category
                palette="tab10", legend_loc="on data")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `adata` | — | AnnData object. |
| `basis` | — | Embedding key. Accepts `"X_umap"` or `"umap"` (an `X_` prefix is added if needed). |
| `color` | — | Gene or `obs` column to color by. |
| `transformation` | `None` | `None` / `"log"` / `"z_score"` (applied when `color` is a gene). |
| `order` | `None` | Point draw order for continuous color: `None`, `"max"` (high on top), `"min"`, `"random"`. |
| `color_map` | `"Viridis"` | Continuous colormap. |
| `palette` | `None` | Categorical palette (name or list). |
| `size` | `5` | Marker size. |
| `opacity` | `1.0` | Marker opacity (0–1). |
| `legend_loc` | `"right margin"` | `"right margin"` or `"on data"` (labels drawn on the clusters). |
| `axis_show` | `True` | Show axis ticks/labels. |

### `gc.pl.umap` / `gc.pl.pca` / `gc.pl.tsne`

Shortcuts for `embedding` with `basis` set to `"X_umap"`, `"X_pca"`, `"X_tsne"`.
All other arguments are forwarded to `embedding`.

```python
gc.pl.umap(adata, color="leiden")
gc.pl.pca(adata,  color="CD8A", order="max")
```

### `gc.pl.coexpression(adata, basis, gene1, gene2, ...)`

Two-gene co-expression on an embedding: cells are binned into four groups
(neither / gene1 only / gene2 only / both) by expression thresholds.

```python
gc.pl.coexpression(adata, "X_umap", "CD8A", "CD4", threshold1=0.5, threshold2=0.5)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `basis` | — | Embedding key. |
| `gene1`, `gene2` | — | The two genes. |
| `threshold1`, `threshold2` | `0.5` | Expression cutoffs that define "positive" for each gene. |
| `transformation` | `None` | `None` / `"log"` / `"z_score"`. |
| `size` | `5` | Marker size. |
| `opacity` | `1.0` | Marker opacity. |

---

## Violin

### `gc.pl.violin(adata, keys, groupby, ...)`

Stacked violins — one row per gene in `keys`, split along the x-axis by
`groupby` categories. Large groups are randomly subsampled (to 5,000 points) so
the distribution shape is preserved while staying fast.

```python
gc.pl.violin(adata, keys=["CD8A", "CD4", "MS4A1"], groupby="cell_type",
             show_box=True, transformation="log")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `keys` | — | Gene(s); a single name or a list (one violin row each). |
| `groupby` | — | `obs` category for the x-axis groups. |
| `labels` | `None` | Subset/order of `groupby` categories (default: all). |
| `transformation` | `None` | `None` / `"log"` / `"z_score"`. |
| `show_box` | `False` | Overlay a neutral **Median + IQR** box. |
| `palette` | `None` | Category palette (name or list). |

### `gc.pl.stacked_violin(adata, var_names, groupby, ...)`

Alias of `gc.pl.violin` (same arguments), named to match scanpy's
`sc.pl.stacked_violin`.

### `gc.pl.violin_grouped(adata, key, groupby, ...)`

Grouped violin of a **single** gene across one or two metadata, with optional
statistical-significance annotations. This is GUANACO's "Violin Plot" comparison
view (the same as the web app when `mode` is passed).

```python
# one metadata
gc.pl.violin_grouped(adata, key="CD8A", groupby="cell_type", test_method="kw-test")

# compare a condition within each cell type
gc.pl.violin_grouped(adata, key="CD8A", groupby="cell_type",
                     groupby2="condition", mode="mode2", test_method="mwu-test")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `key` | — | A **single** gene. |
| `groupby` | — | Primary metadata (`meta1`). |
| `groupby2` | `None` | Secondary metadata (`meta2`); required for modes 2–4. |
| `mode` | `"mode1"` | `"mode1"` one metadata · `"mode2"` facet by `groupby`, compare `groupby2` · `"mode3"` linear model (`groupby` + `groupby2`) · `"mode4"` mixed model (`groupby` + (1\|`groupby2`)). |
| `test_method` | `"none"` | `"none"`, `"mwu-test"`, `"ttest"`, `"kw-test"`, `"anova"`, `"linear-model"`, `"linear-model-interaction"`, `"mixed-model"`. |
| `transformation` | `None` | `None` / `"log"` / `"z_score"`. |
| `show_box` | `True` | Overlay a box plot. |
| `show_points` | `False` | Overlay individual points. |
| `labels` | `None` | Subset/order of categories. |
| `palette` | `None` | Category palette (name or list). |

---

## Heatmap

### `gc.pl.heatmap(adata, var_names, groupby, ...)`

Expression heatmap of `var_names` (rows) across cells grouped by `groupby`.

```python
gc.pl.heatmap(adata, var_names=markers, groupby="cell_type", z_score=True)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `var_names` | — | Gene(s) to show. |
| `groupby` | — | Primary `obs` grouping. |
| `groupby2` | `None` | Optional secondary grouping (adds a second annotation bar). |
| `labels` | `None` | Subset/order of `groupby` categories (default: all). |
| `log` | `False` | Apply log1p to expression. |
| `z_score` | `False` | Z-score each gene across cells. |
| `color_map` | `"Viridis"` | Continuous colormap. |
| `transformation` | `None` | Extra transformation (`None` / `"log"` / `"z_score"`); `log`/`z_score` flags above are the common path. |

---

## Dotplot & matrixplot

Both come from the same engine (`dotmatrix.plot_dot_matrix`) and share parameters.

### `gc.pl.dotplot(adata, var_names, groupby, ...)`

Dotplot: dot **color** = mean expression, dot **size** = fraction of cells
expressing, per gene × group.

```python
gc.pl.dotplot(adata, var_names=markers, groupby="cell_type", standardization="var")
```

### `gc.pl.matrixplot(adata, var_names, groupby, ...)`

Matrixplot: the same mean-expression values drawn as a colored **heatmap grid**
(no dot sizing) — scanpy's `sc.pl.matrixplot`.

```python
gc.pl.matrixplot(adata, var_names=markers, groupby="cell_type",
                 standardization="var", cluster="both")   # with dendrograms
```

Shared parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `var_names` | — | Gene(s). |
| `groupby` | — | `obs` grouping. |
| `labels` | `None` | Subset/order of groups (default: all). |
| `color_map` | `"Viridis"` | Continuous colormap for mean expression. |
| `transformation` | `None` | `None` / `"log"` / `"z_score"`. |
| `standardization` | `None` | Scale mean expression: `None`, `"var"` (per gene), or `"group"` (per group). |
| `cluster` | `"none"` | Dendrograms: `"none"` / `"row"` / `"col"` / `"both"`. |
| `transpose` | `False` | Swap the gene and group axes. |

---

## Stacked bar (composition)

### `gc.pl.stacked_bar(adata, x, color, ...)`

Composition bar chart: each bar is an `x` category; the bar is split/colored by
`color`.

```python
# Fraction of each cell type within each sample:
gc.pl.stacked_bar(adata, x="sample", color="cell_type", normalize="proportion")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `x` | — | `obs` column for the x-axis groups. |
| `color` | — | `obs` column for the stacked layers (the colored segments). |
| `normalize` | `"proportion"` | `"proportion"` (fractions per bar) or `"count"` (raw counts). |
| `x_order` | `None` | Explicit order of the x-axis groups. |
| `palette` | `None` | Palette for the stacked layers (name or list). |

---

## Expression trend along pseudotime

### `gc.pl.pseudotime(adata, genes, ...)`

Smoothed expression trend of `genes` along a pseudotime stored in `adata.obs`.

```python
gc.pl.pseudotime(adata, genes=["CD8A", "GZMB"],
                 pseudotime_key="dpt_pseudotime", groupby="cell_type")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `genes` | — | Gene(s) to plot. |
| `pseudotime_key` | `"pseudotime"` | Name of the `obs` column holding the pseudotime values. |
| `groupby` | `None` | Optional `obs` category to color/group the trend lines. |
| `min_expr` | `0.5` | Drop cells below this expression when fitting the trend. |
| `transformation` | `"none"` | `"none"` / `"log"` / `"z_score"`. |

---

## Volcano (differential expression)

### `gc.pl.volcano(adata, group=None, ...)`

Volcano plot from **precomputed** DE results in `adata.uns` — either
`adata.uns["rank_genes_groups"]` (scanpy's `sc.tl.rank_genes_groups`) or
`adata.uns["volcano"]`.

```python
sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")   # produces the DE
gc.pl.volcano(adata, group="T cell", padj_threshold=0.05, x_threshold=1.0, top_n=15)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `group` | `None` | Which comparison/entry to plot; defaults to the first available. |
| `x_field` | `"logfoldchange"` | Statistic on the x-axis. |
| `padj_threshold` | `0.05` | Adjusted-p significance cutoff (horizontal line). |
| `x_threshold` | `1.0` | Effect-size cutoff (vertical lines). |
| `top_n` | `12` | Number of top genes to label. |

> Requires DE to be computed first. If neither `rank_genes_groups` nor `volcano`
> is in `adata.uns`, the call raises an error.

---

## Peak browser (ATAC)

### `gc.pl.peak_browser(adata, region=None, ...)`

A genome browser over **peak-like features**: one accessibility bar track per group.
It needs genomic peaks — either `var_names` like `"chr1:10000-10500"`, or `adata.var`
columns `["chrom", "start", "end"]` — so it is typically used on an ATAC (or the ATAC
modality of a multiome) dataset.

```python
# Search by gene NAME (resolved via gene_annotation) — like the web app's box:
gc.pl.peak_browser(adata_atac, region="CD8A", groupby="cell_type",
                   gene_annotation="hg38")

# ...or by locus. gene_annotation also draws the gene-model track.
gc.pl.peak_browser(adata_atac, region="chr1:1,000,000-3,000,000",
                   groupby="cell_type", metric="mean", gene_annotation="hg38")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `region` | `None` | A **gene name** (needs `gene_annotation`), a locus `"chr:start-end"`, or a `{"chrom","start","end"}` dict; defaults to a populated window. |
| `groupby` | `None` | `obs` category — one accessibility track per group. |
| `labels` | `None` | Restrict to these groups (subset / order of `groupby`). |
| `metric` | `"mean"` | `"mean"` accessibility or `"detection"` fraction (fraction of cells with a peak). |
| `y_mode` | `"shared"` | `"shared"` = same y-range on every track (heights comparable); `"auto"` = each track scales to its own peak. |
| `max_peaks` | `400` | Cap on peaks drawn in the window (evenly downsampled above this). |
| `gene_annotation` | `None` | Genome id (`"hg38"`, `"mm10"` …) or a GTF/GFF3 path → adds a gene-model track. |
| `selected_cells` | `None` | Restrict the signal to a subset of cell barcodes. |

> Raises an error if `adata` has no peak-like features.

---

## Data requirements at a glance

| Function | Needs in `adata` |
|----------|------------------|
| `embedding` / `umap` / `pca` / `tsne` / `coexpression` | the embedding in `adata.obsm` (e.g. `X_umap`) |
| `violin` / `stacked_violin` / `heatmap` / `dotplot` | genes in `var_names`, a categorical `obs` column |
| `stacked_bar` | two categorical `obs` columns |
| `pseudotime` | a pseudotime column in `obs` |
| `peak_browser` | peak-like `var_names` (`chr1:10000-10500`) or `var[["chrom","start","end"]]` |
| `volcano` | `adata.uns["rank_genes_groups"]` or `adata.uns["volcano"]` |
| `paga` | `adata.uns["paga"]` (run `sc.tl.paga`) |
| `grn` | `adata.uns["grn"]` (a DataFrame: `source`, `target`, `regulation`, optional `weight`/context) |

---

## PAGA & GRN (Cytoscape networks)

PAGA and GRN are interactive **Cytoscape** networks rather than Plotly figures, so
`gc.pl.paga` / `gc.pl.grn` return an [`ipycytoscape`](https://ipycytoscape.readthedocs.io)
widget. They:

- need `pip install ipycytoscape`;
- render in **Jupyter / JupyterLab** (and the full web app `guanaco -c config.json`);
- are **not** available in marimo (classic ipywidgets aren't rendered there);
- use `return_widget` (default `False`) instead of `return_fig`, since they return
  a widget, not a figure.

### `gc.pl.paga(adata, color=None, gene=None, ...)`

Pie-chart PAGA graph (each node = a group; node positions come from the PAGA
layout). Run `sc.tl.paga(adata, groups="cell_type")` first.

```python
gc.pl.paga(adata)                       # nodes colored by the PAGA groupby
gc.pl.paga(adata, color="phase")        # categorical obs -> pie nodes of composition
gc.pl.paga(adata, gene="CD8A")          # mean gene expression per node
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `color` | `None` | An `obs` column: categorical → pie nodes, numeric → solid continuous color. Defaults to the column PAGA was computed on. |
| `gene` | `None` | Color nodes by mean expression of this gene (overrides `color`). |
| `color_map` | `"Viridis"` | Continuous colormap (numeric `color` / `gene`). |
| `palette` | `None` | Categorical palette; default uses the dataset's stored colors, else GUANACO's palette. |
| `edge_threshold` | `0.03` | Hide connectivities below this weight. |
| `annotation` / `labels` | `None` | Optionally restrict to a subset of cells. |

### `gc.pl.grn(adata, context="All", ...)`

Gene-regulatory-network graph from `adata.uns["grn"]` (source/target nodes,
activation/repression edges).

```python
gc.pl.grn(adata)
gc.pl.grn(adata, context="Monocyte", edge_threshold=0.5, layout="circle")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `context` | `"All"` | Filter edges to one cell-type/condition (`"All"` = every context). |
| `edge_threshold` | `None` | Hide edges with `weight` below this. |
| `layout` | `"cose"` | Any cytoscape.js layout: `"cose"`, `"circle"`, `"concentric"`, `"breadthfirst"`. |
