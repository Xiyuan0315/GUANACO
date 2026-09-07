# Linked-view examples

## Interactive showcase

Five runnable examples in one Dash application: cell selection, RNA/protein,
spatial relationships, pathways/genes, and gene-to-ATAC navigation. They reuse
GUANACO's existing linking runtime; no additional web framework is required.

From the repository root, using the project's Python environment:

```bash
python examples/linked_views/app.py
```

Open http://127.0.0.1:8050. Set `PORT` if that port is already occupied.
Select an example in the sidebar. Each page includes instructions, live plots,
reset, the builder's Python code, and data notes. Individual examples can be
linked directly, e.g. `/?demo=atac` or `/?demo=spatial`.

The spatial example uses the **original Visium mouse-brain data** from the
[Squidpy H&E tutorial](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html)
(10x Genomics). The bundled spatial-only file preserves all 2,688 spots, 15 groups,
original H&E images and scale factors, the complete 15 × 15 neighborhood-enrichment
matrix (including diagonal and off-diagonal values), and precomputed co-occurrence
curves. No spatial statistics are recomputed. Unused gene expression and unrelated
metadata are omitted; the matrix in this excerpt is a zero-valued placeholder,
not expression data. Visium spots can contain multiple cells.

The other four examples use **deterministic synthetic data**, including invented
genomic coordinates in the ATAC example. No other local research files are
discovered or downloaded. These examples are not scale benchmarks.

To reproduce the spatial-only export from the original notebook dataset:

```bash
python examples/linked_views/prepare_spatial.py /path/to/visium_hne_spatial.h5ad
```

This writes `data/visium_hne_spatial.h5ad` within the example directory and refuses
to overwrite an existing file. Use `--output /new/path.h5ad` for another export.
The original file is opened read-only. Missing gallery data raises an error rather
than silently substituting a synthetic tissue.

### Publish to Plotly Cloud

Build a clean upload folder, including a wheel of **this checkout**, so Cloud does
not accidentally install an older published GUANACO without the linking API:

```bash
python examples/linked_views/prepare_cloud.py
```

The command prints a new temporary folder. Upload that folder's contents at
https://cloud.plotly.com, choose **app.py** as the main file and **Python 3.12**,
then publish. No environment variables or external data credentials are needed.
In Sharing, set **Anyone with the link → Can view** when ready for public access.
Check your account's resource limits and plan before publishing. Upload only the
prepared folder, not the repository or your environment.

See [Plotly Cloud publishing](https://dash.plotly.com/plotly-cloud/publish) and
[access controls](https://dash.plotly.com/plotly-cloud/share).

The preparation command packages code, synthetic gene annotation, and only the
bundled public spatial excerpt; it does not publish anything. Include the `data`
folder when uploading. Review the folder before upload. Dependency ranges come from
GUANACO's wheel metadata; this is not a fully locked deployment environment.

### Reuse or extend

`cases.py` contains one builder per example and the small `CASES` catalogue.
From a notebook with this directory on `sys.path`:

```python
from cases import rna_protein

demo = rna_protein()
demo.show_jupyter()
```

Builders also accept prepared datasets, so real public data can replace the
synthetic examples without duplicating the plotting code. Update visible data
notes and provenance when changing data. Never silently substitute private files.
For the gene-to-ATAC builder, supply the data dictionary, gene list and annotation.

The showcase constructs five small demos and registers callbacks once per
process. Only the selected example's components are sent to the browser; this is
not lazy loading of datasets. Navigation and reset remount a pristine layout, and
each browser owns its selection in `dcc.Store`. The runtime's notebook selection
mirror is not used to initialize subsequent visitors. Keep this property if
changing routing. Large real datasets need a separate memory/concurrency review.

### Tests

```bash
python -m pytest tests/test_linked_showcase.py
```

Tests use Dash's HTTP endpoints for navigation, all five linking cases, target
redraws, clear/reset responses, and fresh visitor state. They do not replace a
visual browser check or a hosted smoke test.

## Full notebook walkthrough

Open the original 11-example walkthrough:

- [`Linked_views_demo.ipynb`](../notebooks/Linked_views_demo.ipynb)

The notebook covers additional linking patterns and can use real local datasets.
[`demo_data.py`](demo_data.py) provides shared data preparation. The
public API and external-table contract are documented in
[`docs/linked_views.md`](../../docs/linked_views.md).
