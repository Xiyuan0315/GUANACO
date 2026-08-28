# GUANACO dash-draggable fork

This directory vendors the `dash-draggable` 0.1.2 responsive-grid behaviour
for Dash 4. The original MIT license and README are retained as `LICENSE` and
`README.upstream.md`.

The compatibility fork replaces Dash's removed `_dashprivate_layout` access
with the public Dash 4 `componentPath` / `dash_component_api.getLayout` API.
It otherwise keeps the original grid settings, resize behaviour, and browser
layout persistence.

After changing the React source, rebuild the committed runtime from the
GUANACO repository root:

```bash
scripts/build_dash_draggable.sh
```
