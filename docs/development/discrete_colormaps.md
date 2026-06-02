# Discrete Colormaps in GUANACO

GUANACO uses discrete palettes for categorical annotations, cluster labels,
groups, composition plots, and categorical PAGA views. All discrete palette
choices are defined through:

```python
guanaco.utils.colors
```

## What Is Included

| Source | Dropdown options | Naming | Notes |
| --- | ---: | --- | --- |
| Custom/self-defined | 8 | `<name>` | Curated colorblind-friendly palettes in `cvd_color.json` |
| Plotly qualitative | 19 | `plotly/<name>` | Original qualitative palettes only; reversed `_r` palettes are removed |
| Colorcet Glasbey | 13 | `cc/<name>` | Fixed Glasbey palettes from Colorcet |
| Crameri categorical | 16 | `cmc/<name>` | Scientific categorical `S` palettes such as `cmc/batlowS` |
| Dynamic Glasbey | 2 | `glasbey/<mode>` | Generated on demand for the requested number of categories |
| **Total** | **58** |  |  |

## Naming Rule

Programmatic sources use a prefix:

```text
plotly/Plotly
cc/glasbey
cmc/batlowS
glasbey/default
glasbey/safe
```

Custom palettes keep plain names:

```text
okabe_ito
tol_bright
tol_vibrant
tol_muted
polychrome_colorsafe
krzywinski_12
krzywinski_15
krzywinski_24
```

## Inclusion Rules

1. Custom palettes are curated for colorblind-aware categorical display.
2. Plotly qualitative palettes include originals only; reversed `_r` variants
   are not shown for discrete palettes.
3. Colorcet contributes fixed Glasbey palettes only.
4. Crameri contributes categorical maps ending in `S`.
5. `glasbey/default` and `glasbey/safe` are dynamic options. They generate a
   palette sized to the number of categories in the plot and cache it for the
   current Python session.

## Custom Palettes

The custom palettes are grouped by practical cluster count:

For up to 8 clusters:

```text
okabe_ito
tol_bright
tol_vibrant
```

For 9-12 clusters:

```text
tol_muted
polychrome_colorsafe
krzywinski_12
```

For 13-24 clusters:

```text
krzywinski_15
krzywinski_24
```

For larger or unknown category counts, use:

```text
glasbey/default
glasbey/safe
```

## Inspect From Python

Get the grouped overview:

```python
from guanaco.utils.colors import discrete_palette_overview

overview = discrete_palette_overview()
custom = overview["custom"]
plotly = overview["plotly"]
colorcet = overview["cc"]
crameri = overview["cmc"]
glasbey = overview["glasbey"]
```

Each entry contains:

```python
{
    "source": "custom",
    "name": "okabe_ito",
    "label": "okabe_ito",
    "value": "okabe_ito",
    "n_colors": 8,
    "dynamic": False,
    "palette": ["#000000", "..."],
}
```

Get dropdown options:

```python
from guanaco.utils.colors import discrete_palette_options

options = discrete_palette_options()
```

Resolve one palette:

```python
from guanaco.utils.colors import resolve_discrete_palette

palette = resolve_discrete_palette("okabe_ito")
palette = resolve_discrete_palette("glasbey/safe", n_colors=37)
```

## Preview Palettes

In a notebook:

```python
from guanaco.utils.colors import plot_discrete_palette_overview

fig = plot_discrete_palette_overview()
```

Preview one source:

```python
fig = plot_discrete_palette_overview(source="custom")
fig = plot_discrete_palette_overview(source="plotly")
fig = plot_discrete_palette_overview(source="cc")
fig = plot_discrete_palette_overview(source="cmc")
fig = plot_discrete_palette_overview(source="glasbey")
```

Preview dynamic Glasbey with a larger category count:

```python
fig = plot_discrete_palette_overview(source="glasbey", dynamic_n_colors=80)
```

In iTerm or a terminal, save and open the preview:

```python
fig = plot_discrete_palette_overview()
fig.savefig("discrete_colormaps.png", dpi=200, bbox_inches="tight")
```

```bash
open discrete_colormaps.png
```
