# Volcano Plot in GUANACO

GUANACO supports Volcano Plot as an optional plot in the lower matrix plot
section, alongside heatmap, dotplot, violin plot, stacked bar, pseudotime, and
PAGA.

Enable it in a dataset config with:

```json
"optional_plot_components": ["dotplot", "heatmap", "volcano"]
```

The Volcano Plot reads differential expression results from `adata.uns` and
renders one contrast or group at a time. The plot includes:

- group or contrast dropdown
- x-axis dropdown
- p-value or adjusted p-value threshold
- x-axis threshold
- top-gene label count
- DEG summary after threshold
- DEG CSV download

The graph uses the same Plotly graph configuration as the other lower-section
plots, including the shared mode bar settings.

## Supported DEG Sources

GUANACO currently supports two DEG result sources for Volcano Plot.

### Scanpy

Scanpy DEG results are read directly from:

```python
adata.uns["rank_genes_groups"]
```

This is the standard output from:

```python
sc.tl.rank_genes_groups(...)
```

GUANACO extracts:

- `names`
- `logfoldchanges`
- `pvals_adj` if present, otherwise `pvals`
- `scores` if present

No conversion helper is needed for standard Scanpy output.

### PyDESeq2

PyDESeq2 results are usually stored as a pandas DataFrame or CSV file. GUANACO
does not read arbitrary CSV files directly during app rendering. Instead, use
the helper below to save the result into `adata.uns["volcano"]`.

```python
from guanaco.utils import save_pydeseq_results_to_adata_uns

save_pydeseq_results_to_adata_uns(
    adata,
    results_df_or_csv,
    entry_name="treated_vs_control",
)
```

`results_df_or_csv` can be either:

- a pandas DataFrame
- a CSV file path

The default expected PyDESeq2 columns are:

```text
gene_name
log2FoldChange
padj
stat
```

Common extra columns are preserved for DEG export:

```text
baseMean
lfcSE
pvalue
gene_id / Ensembl_gene_id_v / DataFrame index
```

The helper writes:

```python
adata.uns["volcano"]
```

After that, GUANACO can use the result in the Volcano Plot tab.

## Generic GUANACO Volcano Payload

The generic payload has this shape:

```python
adata.uns["volcano"] = {
    "entries": {
        "treated_vs_control": {
            "gene": [...],
            "logfoldchange": [...],
            "padj": [...],
            "score": [...],
        }
    },
    "default_entry": "treated_vs_control",
}
```

Required fields per entry:

```text
gene
logfoldchange
padj
```

Optional fields:

```text
score
any extra per-gene columns
```

Extra per-gene columns are included in downloaded DEG CSV files.

## DEG CSV Download

The downloaded DEG CSV contains genes passing both active thresholds:

```text
p or padj <= p-value threshold
abs(selected x-axis value) >= x-axis threshold
```

The filename includes the contrast, p-value type and threshold, x-axis type and
threshold:

```text
treated_vs_control_padj005_log2fc05_deg.csv
```

In this example:

- `padj005` means adjusted p-value threshold `0.05`
- `log2fc05` means log2 fold-change threshold `0.5`

