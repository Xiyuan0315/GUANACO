# Continuous Colormaps in GUANACO

GUANACO exposes continuous colormaps for expression values, continuous
annotations, heatmaps, dotplots, and gene-colored PAGA views. All continuous
colormap choices are defined in one place:

```python
guanaco.utils.colors
```

## What Is Included

GUANACO includes three continuous colormap sources:

| Source | Original maps | Dropdown maps | Notes |
| --- | ---: | ---: | --- |
| Plotly-delivered scientific maps | 23 | 46 | Curated perceptually uniform maps from Matplotlib viridis family, cmocean, and ColorBrewer |
| Colorcet | 56 | 112 | Linear + diverging maps only; rainbow maps removed |
| Crameri | 31 | 62 | Curated scientific linear + diverging maps |
| **Total** | **110** | **220** | Dropdown includes each original map plus its reversed `_r` version |

## Naming Rule

Dropdown labels use:

```text
source/type_name
```

Examples:

```text
plotly/linear_viridis
plotly/linear_viridis_r
cc/diverging_bwr_40_95_c42
cmc/diverging_vik
cmc/diverging_vik_r
```

Stored internal values are source-specific:

```text
plotly/linear_viridis -> viridis
plotly/linear_viridis_r -> viridis_r
cc/diverging_bwr_40_95_c42 -> cc:diverging_bwr_40_95_c42
cmc/diverging_vik -> cmc:vik
cmc/diverging_vik_r -> cmc:vik_r
```

## Inclusion Rules

GUANACO includes only continuous colormaps that are appropriate for scientific
data display:

1. Include linear/sequential maps for ordered low-to-high values.
2. Include diverging maps for values centered around a meaningful midpoint.
3. Add reversed `_r` variants for every dropdown map.
4. Do not include cyclic maps.
5. Do not include rainbow maps.
6. Do not include multi-sequential maps such as `oleron`, `bukavu`, and `fes`.
7. Do not expose arbitrary Plotly continuous colormaps; only the curated
   perceptually uniform list is shown.

## Plotly-Delivered Curated Maps

Linear:

```text
viridis, plasma, inferno, magma, cividis,
thermal, haline, solar, ice, deep, dense, matter, speed, amp, tempo
```

Diverging:

```text
balance, curl, delta, rdbu, brbg, puor, piyg, prgn
```

## Colorcet Rule

Colorcet maps are selected with:

```python
linear = cc.all_original_names(group="linear")
diverging = cc.all_original_names(group="diverging")
names = sorted(set(linear + diverging))
names = [name for name in names if "rainbow" not in name.lower()]
```

This gives 56 original Colorcet maps. The Colorcet map
`diverging_rainbow_bgymr_45_85_c67` is excluded because it is a rainbow map.

## Crameri Curated Maps

Linear:

```text
batlow, batlowK, batlowW, acton, bamako, bilbao, buda, davos, devon,
grayC, hawaii, imola, lajolla, lapaz, lipari, navia, nuuk, oslo,
tokyo, turku, glasgow
```

Diverging:

```text
vik, roma, broc, cork, bam, berlin, lisbon, tofino, vanimo, managua
```

## Inspect From Python

Get original colormaps only:

```python
from guanaco.utils.colors import continuous_colormap_overview

overview = continuous_colormap_overview()
linear = overview["linear"]
diverging = overview["diverging"]
```

Get dropdown options, including reversed `_r` variants:

```python
from guanaco.utils.colors import continuous_colormap_options

options = continuous_colormap_options()
```

## Preview Colormaps

In a notebook:

```python
from guanaco.utils.colors import plot_continuous_colormap_overview

fig = plot_continuous_colormap_overview()
```

Preview only one type:

```python
fig = plot_continuous_colormap_overview(colormap_type="linear")
fig = plot_continuous_colormap_overview(colormap_type="diverging")
```

Preview only one source:

```python
fig = plot_continuous_colormap_overview(source="plotly")
fig = plot_continuous_colormap_overview(source="cc")
fig = plot_continuous_colormap_overview(source="cmc")
```

In iTerm or a terminal, save the preview and open it:

```python
fig = plot_continuous_colormap_overview()
fig.savefig("continuous_colormaps.png", dpi=200, bbox_inches="tight")
```

```bash
open continuous_colormaps.png
```

The preview intentionally shows original colormaps only. Reversed `_r` variants
are available in the GUANACO dropdown, but are omitted from the overview plot to
avoid duplicate rows.
