"""Compatibility deployment entrypoint.

The raw Dash object lives in :mod:`guanaco.dash_app`, but the app is not usable
until :mod:`guanaco.main` assigns the layout and registers callbacks. Some cloud
builders discover/import ``guanaco.app`` directly, so this module intentionally
exports the fully initialized app from ``guanaco.main`` -- importing it always
yields an app whose ``layout`` is set (avoiding ``NoLayoutException``).
"""

from guanaco.main import app, server  # noqa: F401
