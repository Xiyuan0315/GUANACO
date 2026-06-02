# GRN Plot in GUANACO

GUANACO supports a data-driven gene regulatory network (GRN) plot as an
optional plot in the lower matrix plot section.

Enable it in a dataset config with:

```json
"optional_plot_components": ["grn"]
```

The older key `"grn-demo"` is still accepted as an alias, but `"grn"` is the
preferred name.

## Required Input

The GRN plot reads a pandas DataFrame from:

```python
adata.uns["grn"]
```

The DataFrame must contain these required columns:

```text
source
target
regulation
```

Column meanings:

| Column | Required | Meaning |
| --- | --- | --- |
| `source` | yes | Regulator/source node, usually a transcription factor |
| `target` | yes | Target gene/node |
| `regulation` | yes | Direction of regulation |
| `weight` | no | Edge strength; used for edge thickness and threshold filtering |
| `context` | no | Grouping context such as cell type, condition, treatment, or time point |

The `regulation` column supports:

```text
+           activation
-           repression
activation  activation
repression  repression
```

## Context Column

The optional context column controls the dropdown shown above the GRN plot.
It can represent different cell types, conditions, treatments, time points, or
other analysis groups.

Preferred column name:

```text
context
```

Supported aliases:

```text
cell_type
condition
group
category
```

If no context-like column is present, the dropdown shows only:

```text
All
```

## Weight Column

If the DataFrame contains a `weight` column, GUANACO uses it to:

- set edge thickness
- show the edge-threshold control
- filter edges with `weight >= threshold`

If `weight` is not present, all edges use the same thickness and the
edge-threshold control is hidden.

The threshold input range is derived from the actual weight range in the data.
This supports both normalized weights, such as `0.0` to `1.0`, and larger
importance scores, such as `3.2` to `14.1`.

## Example Input

```python
import pandas as pd

grn = pd.DataFrame(
    {
        "source": ["STAT1", "STAT1", "NFKB1", "RELA"],
        "target": ["IRF1", "CXCL10", "IL1B", "CCL2"],
        "regulation": ["+", "+", "+", "-"],
        "weight": [0.92, 0.89, 0.87, 0.81],
        "cell_type": ["Monocyte", "Monocyte", "Macrophage", "Macrophage"],
    }
)

adata.uns["grn"] = grn
```

Equivalent CSV:

```csv
source,target,regulation,weight,cell_type
STAT1,IRF1,+,0.92,Monocyte
STAT1,CXCL10,+,0.89,Monocyte
NFKB1,IL1B,+,0.87,Macrophage
RELA,CCL2,-,0.81,Macrophage
```

Load it with:

```python
adata.uns["grn"] = pd.read_csv("demo_grn_cytoscape_edges.csv")
```

## Adapting Common GRN Outputs

If your GRN result uses columns like this:

```text
TF
target
importance
rho
regulation
direction
```

convert it before launching GUANACO:

```python
grn = adata.uns["grn"].copy()

grn = grn.rename(
    columns={
        "TF": "source",
        "importance": "weight",
    }
)

grn["regulation"] = grn["direction"].map(
    {
        "activating": "+",
        "repressing": "-",
        "activation": "+",
        "repression": "-",
    }
).fillna(grn["regulation"])

adata.uns["grn"] = grn
```

After conversion, the DataFrame has the required fields:

```text
source
target
regulation
```

and optional fields such as:

```text
weight
rho
direction
```

Extra columns are allowed. GUANACO currently uses only `source`, `target`,
`regulation`, optional `weight`, and optional context columns.

## Plot Behavior

The GRN plot renders:

- source nodes as rounded rectangles
- target-only nodes as circles
- activation edges in green with arrowheads
- repression edges in red with tee-shaped heads
- edge thickness from `weight`, if provided

The controls are:

- context dropdown, if a context-like column is present
- graph layout dropdown
- edge threshold, only when `weight` is present

