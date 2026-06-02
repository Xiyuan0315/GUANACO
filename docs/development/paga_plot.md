# PAGA Plot in GUANACO

GUANACO renders PAGA as an interactive Cytoscape network. Each node represents
one PAGA group, each edge represents PAGA connectivity between groups, and the
node color summarizes either categorical cell composition or mean gene
expression.

## Data Requirements

The PAGA panel is available when the input `AnnData` object contains:

```python
adata.uns["paga"]["connectivities"]
```

GUANACO uses the PAGA grouping stored in:

```python
adata.uns["paga"]["groups"]
```

If `adata.uns["paga"]["groups"]` is not present, GUANACO tries to infer the PAGA
grouping from a categorical `adata.obs` column whose category count matches the
number of PAGA nodes.

Optional PAGA positions can be provided with:

```python
adata.uns["paga"]["pos"]
```

If positions are missing, GUANACO computes a graph layout from the connectivity
matrix.

## How to Use

Open the `PAGA` tab from the matrix visualization page.

### Color By obs

Choose `Color By: obs` to show categorical composition inside each PAGA node.

In this mode:

- `obs Column` only shows categorical `adata.obs` columns.
- Each node is drawn as a pie chart.
- Pie slices show the fraction of selected cells in that PAGA group.
- The right panel shows the category color legend.
- Hovering over a node shows the group name, cell count, category counts, and
  category percentages.

The discrete colors come from the shared embedding discrete colormap control.
There is no separate PAGA discrete colormap dropdown. The default shared
discrete colormap is:

```text
cc/glasbey
```

### Color By gene

Choose `Color By: gene` and select a gene to show expression on the PAGA graph.

In this mode:

- Each node is drawn as a solid circle.
- The node color is the mean expression of the selected gene in that PAGA group.
- The right panel shows a continuous colorbar.
- Hovering over a node shows the exact mean expression value and the number of
  cells used for that group.

The continuous colors come from the shared embedding continuous colormap control.
There is no separate PAGA continuous colormap dropdown.

### Connectivity Threshold

Use `Connectivity Threshold` to hide weak edges. Edges with connectivity below
the threshold are not shown. Visible edge widths are scaled by connectivity
strength.

### Selection Behavior

PAGA responds to the active cell selection used by the matrix page.

When cells are filtered or selected:

- Node composition in `obs` mode is recalculated from the selected cells.
- Mean expression in `gene` mode is recalculated from the selected cells.
- Node size reflects the number of selected cells in each PAGA group.
- The graph structure still comes from `adata.uns["paga"]["connectivities"]`.

## Visual Rules

PAGA is Cytoscape-only in GUANACO.

Current behavior:

- No renderer selector is shown.
- Node labels are fixed to the PAGA group labels.
- Gene coloring changes node color, colorbar, and hover values, but does not
  change the PAGA group labels.
- `obs` mode and `gene` mode use the same node diameter scale.
- Hover details are shown in the right panel, not below the graph, so long text
  does not resize the plot.

## Implementation

The PAGA implementation is split across layout, callback, and plot files.

```text
src/guanaco/pages/matrix/layouts/paga_layout.py
src/guanaco/pages/matrix/callbacks/paga_callbacks.py
src/guanaco/pages/matrix/plots/paga.py
```

### Layout

`generate_paga_layout()` builds the PAGA controls and graph container.

Important layout choices:

- The renderer option is removed.
- The only renderer is Cytoscape.
- The local PAGA colormap dropdowns are removed.
- `obs Column` is restricted to non-numeric `adata.obs` columns.
- The gene dropdown is search-based, so large gene lists do not need to be
  rendered all at once.

### Callbacks

`register_paga_callbacks()` wires the PAGA controls to the Cytoscape plot.

Important callback inputs:

```text
paga color mode
paga obs column
paga gene
shared continuous colormap dropdown
shared discrete colormap dropdown
connectivity threshold
current annotation selection
current label selection
selected cells
active matrix tab
```

The PAGA callback only updates when the active tab is `paga-tab`.

For categorical `obs` coloring, the callback resolves the selected shared
discrete colormap with:

```python
resolve_discrete_palette(...)
```

For gene coloring, the callback passes the selected shared continuous colormap
to the plot builder.

Hover text is stored on each Cytoscape node as `hover_text`. A small callback
reads `mouseoverNodeData` and renders the hover details in the right panel.

### Plot Builder

`build_paga_cytoscape()` is the main rendering function.

It performs these steps:

1. Read PAGA connectivity, groups, positions, and current cell selection.
2. Compute node sizes from selected-cell counts.
3. Filter edges by connectivity threshold.
4. Build Cytoscape node and edge elements.
5. Build either categorical pie node data or continuous solid-node data.
6. Add a right-side legend panel.
7. Return a Dash Cytoscape component wrapped in a responsive layout.

Key helper functions in `paga.py`:

```text
_prepare_paga_context()
_categorical_group_composition()
_resolve_obs_color_mapping()
_values_to_colors()
_continuous_colorbar_legend()
_format_pie_hover_text()
_format_continuous_hover_text()
```

### obs Mode Internals

For categorical `obs` coloring, GUANACO computes the category composition inside
each PAGA group after applying the active selection.

The result is encoded into Cytoscape pie fields:

```text
pie_1_size, pie_1_color
pie_2_size, pie_2_color
...
```

At most 16 pie slices are shown. If a column has more categories than this,
small categories are collapsed into an `Other` slice.

When a user-selected discrete colormap is active, it is preferred over existing
`adata.uns["<obs>_colors"]` colors so that PAGA matches the shared embedding
color control.

### Gene Mode Internals

For gene coloring, GUANACO extracts the selected gene expression, computes the
mean expression per PAGA group, and maps those values through the shared
continuous colormap.

The node stores:

```text
background_color
hover_text
```

The right panel shows a generated continuous colorbar for the same value range.

To prevent stale pie styling when switching from `obs` mode to `gene` mode, all
pie fields are explicitly reset to empty values before mode-specific node data
is assigned.

