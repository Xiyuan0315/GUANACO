"""Dash 4-compatible GUANACO fork of dash-draggable."""

import json
from pathlib import Path

from ._imports_ import *  # noqa: F403
from ._imports_ import __all__


_package_path = Path(__file__).with_name("package-info.json")
with _package_path.open() as package_file:
    _package = json.load(package_file)

__version__ = _package["version"]
_namespace = _package["name"].replace(" ", "_").replace("-", "_")

_js_dist = [
    {
        "relative_package_path": "dash_draggable.min.js",
        "namespace": _namespace,
    },
    {
        "relative_package_path": "dash_draggable.min.js.map",
        "namespace": _namespace,
        "dynamic": True,
    },
]
_css_dist = []

for _component in __all__:
    setattr(locals()[_component], "_js_dist", _js_dist)
    setattr(locals()[_component], "_css_dist", _css_dist)
