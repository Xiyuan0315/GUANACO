"""The bare Dash application object.

This module *only* constructs the ``dash.Dash`` instance (and its WSGI ``server``)
and sets the title. It deliberately does **not** assign ``app.layout`` or register
any callbacks -- that happens in :mod:`guanaco.main`. Keeping the raw object here
(rather than in :mod:`guanaco.app`) lets ``guanaco.main`` import it without a
circular dependency, while :mod:`guanaco.app` re-exports the *fully initialized*
app from ``guanaco.main`` for cloud builders that import ``guanaco.app`` directly.
"""

import json
import os
from pathlib import Path

import dash
import dash_bootstrap_components as dbc


def load_config(json_path):
    if not json_path.exists():
        raise FileNotFoundError(f"Config file not found: {json_path}")
    return json.loads(json_path.read_text())


# Get config path from environment variable or use default.
JSON_PATH = Path(os.environ.get("GUANACO_CONFIG", "guanaco.json"))

# Initialize the app.
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.LUX, "/assets/scientific_style.css"],
    suppress_callback_exceptions=True,
)

# Enable dash-cytoscape's extra layouts and SVG image export, used by the
# PAGA / GRN "Download SVG" buttons. No-op if dash-cytoscape isn't installed.
try:
    import dash_cytoscape as _cyto

    _cyto.load_extra_layouts()
except Exception:
    pass

config = load_config(JSON_PATH)
app.title = config.get("title", "GUANACO")
server = app.server
