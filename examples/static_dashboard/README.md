# Offline dashboard demo

Build from the repository root, in your GUANACO Python environment:

```bash
python examples/static_dashboard/build_demo.py
```

Open `demo.html` directly in a modern browser. Disconnecting from the internet or
closing Python does not affect it. Each chart independently offers **Cell type**
and **Condition**, both computed during export. No linking is defined or inferred.
Use `--overwrite` only to intentionally replace an existing output.

The data are the existing gallery's deterministic simulation, not patient data.
The demo reuses `gc.pl.dotplot`, `gc.pl.violin` and `gc.pl.heatmap`:

- Dot means and detection fractions are calculated by GUANACO during export.
- The heatmap contains prepared cell-level/binned expression, ordered by group;
  grouping is not the same as aggregating to one mean per group.
- `freeze_violin` converts the existing violin samples into precomputed Gaussian
  density polygons, using the supplied bandwidth and a finite evaluation grid.
  This is an approximation of the interactive violin outline, not an image.
  Browser-side KDE is avoided. Sample means are marked; raw violin values are not
  included by that adapter. Existing GUANACO violin sampling still applies.

## Use your own figures

```python
from guanaco.static_export import StaticPanel, freeze_violin, write_dashboard

write_dashboard("results.html", [
    StaticPanel("expression", "Gene expression", {
        "Cell type": cell_type_dotplot,
        "Sample": sample_dotplot,
    }),
    StaticPanel("distribution", "Expression distributions", {
        "Cell type": freeze_violin(violin_figure),
    }),
])
```

Inputs are ordinary Plotly figures from notebooks, scripts or dashboard code.
This prototype supports scatter, bar and heatmap traces. The violin adapter
supports GUANACO's vertical, width-scaled violins without points or box overlays.
Other statistical traces fail explicitly rather than silently recomputing.

The viewer uses embedded Plotly.js, CSS, JavaScript and figure data; there are no
CDN dependencies or backend callbacks. It blocks network connections via CSP.
Figure rendering and hover/zoom still require normal browser computation. New
groupings, new genes, new statistical results and cell-selection summaries require
re-exporting. Only the two supplied presets are available in this demo.

This is a minimal reusable export surface, **not** a full wizard integration or
automatic export of every GUANACO plot. Stable panel IDs reserve a clear place for
a future optional, user-configured linking component; none is implemented here.

The HTML contains plotted data (including the heatmap values). It is not encrypted
and cannot enforce recipient permissions. More presets increase file size. This
small demo is not a million-cell browser-memory benchmark.
