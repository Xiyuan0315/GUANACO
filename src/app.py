"""Plotly Cloud deployment entrypoint.

Plotly Cloud (dash-app-builder) flattens this file's directory to the app root
and serves it with gunicorn as ``app:server``. We import the *fully initialized*
app from ``guanaco.main`` -- which assigns ``app.layout`` and registers every
callback -- so the WSGI ``server`` is ready to serve (no ``NoLayoutException``).

Set the platform's "Main file" to ``app.py``.
"""

from guanaco.main import app, server  # noqa: F401

if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=4399)
