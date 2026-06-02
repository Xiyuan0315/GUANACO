# Expression Trend Plot

The Expression Trend plot shows how selected gene expression changes along a continuous observation variable from `adata.obs`.

Although this plot was originally built for pseudotime, the x-axis can be any numeric per-cell observation, such as:

- pseudotime
- diffusion pseudotime
- latent time
- trajectory score
- cell-cycle score
- spatial coordinate
- disease severity score
- quality-control metric
- any other continuous `adata.obs` column

Cells are colored by the selected annotation in the left control panel. If that annotation is categorical, the plot uses the same discrete colormap selected in the embedding plot controls. The point size and opacity also follow the embedding plot marker size and opacity controls.

## What The Plot Shows

For each selected gene, the plot displays:

- One scatter panel per gene.
- X-axis: the selected continuous `adata.obs` variable.
- Y-axis: expression of that gene.
- Point color: the selected annotation or group.
- Black curve: a smoothed expression trend across the continuous variable.

This is useful for checking whether genes increase, decrease, peak, or change differently across cell groups along a biological or technical gradient.

## Data Requirements

Your AnnData object must contain:

- Gene expression in `adata.X`.
- Gene names in `adata.var_names`.
- A numeric continuous column in `adata.obs`.
- Optional categorical columns in `adata.obs` for coloring/grouping cells.

The continuous column must be numeric. Columns with strings such as `"early"`, `"middle"`, and `"late"` will not work unless converted to numeric values.

Example:

```python
adata.obs["latent_time"] = adata.obs["latent_time"].astype(float)
```

The app automatically prioritizes numeric columns whose names contain `pseudotime` or `dpt`, but any numeric `adata.obs` column with enough unique values can be selected from the Continuous Variable dropdown.

## Example AnnData Preparation

```python
import scanpy as sc

adata = sc.read_h5ad("my_dataset.h5ad")

# Required: a continuous per-cell value
adata.obs["latent_time"] = adata.obs["latent_time"].astype(float)

# Optional: a categorical label for coloring cells
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

adata.write_h5ad("my_dataset_with_expression_trend.h5ad")
```

## Example Config

Use `expression-trend` to enable this plot:

```yaml
optional_plot_components:
  - heatmap
  - violin
  - dotplot
  - stacked-bar
  - expression-trend
```

The legacy config value `pseudotime` is still supported for compatibility.

## Notes

- Infinite values and missing values in the continuous variable are ignored.
- Cells can be filtered through the embedding plot selection workflow before updating this plot.
- The Minimum Expression Threshold keeps cells where at least one selected gene passes the threshold.
- Log and z-score transformations affect the displayed expression values and smoothed trend.
