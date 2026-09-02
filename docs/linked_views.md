# Linked views

GUANACO links plots through stable IDs:

```text
selected mark  ->  row / cell / feature IDs  ->  detail plot
```

Ligand–receptor, spatial-neighborhood, volcano, and cross-omics views are
different data arrangements built on this one rule. They do not have separate
linking APIs.

## Quick start

```python
import guanaco as gc

demo = gc.pl.linked_view(
    adata,
    views=[
        gc.pl.umap(id="cells", color="cell_type"),
        gc.pl.violin(
            id="expression",
            keys=["CD4"],
            groupby="cell_type",
        ),
    ],
    links=[gc.pl.link("cells", "expression")],
)

demo.show_jupyter(port=8060, height=620)
```

Both plots use the same `AnnData`, so the omitted `by` means shared
`obs_names`. Lassoing cells in the UMAP updates the violin as a
Selected-versus-Others comparison.

Named plotting functions work in both contexts:

```python
gc.pl.umap(adata, color="cell_type")        # standalone figure
gc.pl.umap(id="cells", color="cell_type")  # linked-view specification
```

Use `gc.pl.view()` for registered types without a named helper:

```python
gc.pl.view(
    "plotly.scatter",
    id="scores",
    data="scores",
    x="attribute_1",
    y="attribute_2",
)
```

Both forms use the same GUANACO plot builders and style.

## The complete link API

```python
gc.pl.link(source, target, *, by=None, action=None, key=None)
```

Defaults cover the common cases:

| Data behind both views | Default identity | Default action |
|---|---|---|
| `AnnData` → `AnnData` | cells (`obs_names`) | `highlight` |
| `DataFrame` → `DataFrame` | table index | `filter` |
| mixed `DataFrame` / `AnnData` | declare `by="cell"` or `by="feature"`; optional `key` maps table rows | inferred |

`action="highlight"` keeps the complete target and distinguishes selected
items. `action="filter"` keeps only selected items. Feature links have no
action: the selected feature IDs become the target's feature context.

```python
gc.pl.link("umap", "violin")
gc.pl.link("umap", "heatmap", action="filter")
gc.pl.link("volcano", "expression", by="feature")
gc.pl.link("network", "pairs")
# one pair row can update many child rows in both detail tables
gc.pl.link("neighborhoods", "locations", key="pair_id")
# table membership rows can update native AnnData observations
gc.pl.link("neighborhoods", "spatial", by="cell", key="cell_id", action="filter")
```

Click versus lasso is inferred from the source plot. A click may emit one ID,
a lasso many IDs, and an aggregated mark many member IDs. One-to-one,
one-to-many, many-to-one, and many-to-many therefore need no separate syntax.

Links are directed from overview to terminal detail. One overview may update
several details, and several overviews may update one detail. A detail cannot
become another source; this prevents redraws from creating ambiguous reverse
selections. Detail plots are zoom-only.

## The four data cases

### 1. Cells inside AnnData or MuData

Both plots refer to observations with matching `obs_names`:

```python
gc.pl.link("rna_umap", "protein_umap")
```

The default is `by="cell", action="highlight"`. Use
`action="filter"` when unselected cells should disappear from the detail.

### 2. Features inside AnnData

Feature-oriented marks emit `var_names`:

```python
demo = gc.pl.linked_view(
    adata,
    views=[
        gc.pl.matrixplot(
            id="summary",
            var_names=["CD4", "IL7R", "NKG7"],
            groupby="cell_type",
        ),
        gc.pl.umap(id="expression", color="CD4"),
    ],
    links=[gc.pl.link("summary", "expression", by="feature")],
)
```

Clicking a matrix tile changes the feature shown by the UMAP. The same link
works if the target is a violin, ridge, heatmap, dotplot, or matrixplot.
Targets that display one feature, such as an expression-colored UMAP, use the
first selected feature; multi-feature targets retain the complete selection.

### 3. An external table with AnnData

The table index must match the declared AnnData identity.

Cell-indexed table:

```python
# cell_scores.index matches adata.obs_names
demo = gc.pl.linked_view(
    {"scores": cell_scores, "cells": adata},
    views=[
        gc.pl.view(
            "plotly.scatter", id="scores", data="scores",
            x="attribute_1", y="attribute_2",
        ),
        gc.pl.umap(id="umap", data="cells", color="cell_type"),
    ],
    links=[gc.pl.link("scores", "umap", by="cell")],
)
```

Feature-indexed table:

```python
# differential.index matches adata.var_names
demo = gc.pl.linked_view(
    {"de": differential, "cells": adata},
    views=[
        gc.pl.view(
            "plotly.scatter", id="volcano", data="de",
            x="log2_fold_change", y="minus_log10_padj",
        ),
        gc.pl.violin(
            id="expression", data="cells",
            keys=["CD4"], groupby="cell_type",
        ),
    ],
    links=[gc.pl.link("volcano", "expression", by="feature")],
)
```

### 4. External table to external table

External views use the DataFrame index as their atomic row identity. Usually
both views should use one atomic long-form table:

```python
# one row is one ligand–receptor interaction
demo = gc.pl.linked_view(
    significant_interactions,
    views=[
        gc.pl.view(
            "network", id="network",
            source="source", target="target",
        ),
        gc.pl.view(
            "plotly.scatter", id="pairs",
            x="ligand", y="receptor", color="magnitude_rank",
        ),
    ],
    links=[gc.pl.link("network", "pairs")],
)
```

The network groups rows by `source` and `target`, but each arrow retains the
indices of its contributing interactions. Clicking the arrow filters the
detail to exactly those rows. The same contract powers neighborhood heatmap →
spatial locations, pathway → genes, and clone summary → cells.

When the two tables have different child granularities, add a shared logical
parent column and pass it as `key`. The source mark emits its atomic row ID,
which is converted to that key; every matching target row is then selected.
This is still an ordinary table link, not a spatial- or network-specific mode:

```python
gc.pl.link("overview", "spots", key="pair_id")
gc.pl.link("overview", "curves", key="pair_id")
```

Here `overview.index` contains one row per pair, while `spots.index` and
`curves.index` remain unique atomic IDs and their `pair_id` columns may repeat.

For a mixed table/AnnData link, the table key is translated to the AnnData
identity domain: `by="cell"` matches `key` values to `obs_names`, and
`by="feature"` matches them to `var_names`. This lets an aggregated table mark
drive the native `gc.pl.spatial(...)`, `gc.pl.umap(...)`, or feature plots
without materializing a second coordinates table.

## External table contract

Use a `pandas.DataFrame` with:

1. a stable, unique, non-null index;
2. one row for the smallest record the user may retrieve;
3. ordinary columns for coordinates, groups, labels, scores, and weights;
4. an optional repeated parent column such as `pair_id` for one-to-many links.

For example:

| index | source | target | ligand | receptor | score |
|---|---|---|---|---|---:|
| `lr_001` | NK | B | CD4 | HLA-DRA | 0.93 |
| `lr_002` | NK | B | CD4 | HLA-DRB | 0.88 |

Do not reset the index after filtering and do not pre-aggregate away the rows
needed by the detail plot. This is ordinary long/tidy data; "atomic" is the
important property.

`linked_view()` accepts `AnnData`, `MuData`, a DataFrame, or a mapping of names
to those objects. Set `data=` on each view when the input has several sources.
Here `data=` is the source name; data objects belong only in `linked_view()`.
Construction validates unique IDs, plot capabilities, and source/target index
overlap before starting the app.

## Notebook display and selection retrieval

```python
demo.show_jupyter(port=8060, height=680)
demo.show_marimo(height=680)

selection = demo.get_selection("cells")
if selection is not None:
    print(selection.by)   # "cell"
    print(selection.ids)  # obs_names tuple
```

Rerun the reading cell after interacting. If one source emits both cells and
features, pass `by="cell"` or `by="feature"` to `get_selection()`.

Use `demo.describe()` to inspect resolved links and
`gc.pl.linked_plot_types()` to inspect registered plot contracts.

## Adding a plot type

Most new biological applications require only a correctly indexed table, not
new linking code. Add a `PlotAdapter` only for a new visual primitive. It has
three responsibilities:

1. validate one `ViewSpec`;
2. render one `ViewState`;
3. decode a browser mark into `MarkMembers`.

The adapter never implements callbacks between plots. Register it once and
pass the registry to `linked_view()`:

```python
from guanaco.linking import PlotAdapter, default_plot_registry

registry = default_plot_registry()
registry.register("my_plot", MyPlotAdapter())

demo = gc.pl.linked_view(
    data,
    views=[...],
    links=[...],
    registry=registry,
)
```

See the runnable examples in
[`examples/notebooks/Linked_views_demo.ipynb`](../examples/notebooks/Linked_views_demo.ipynb).
