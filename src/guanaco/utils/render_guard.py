"""Helpers to avoid recomputing/redrawing a plot when nothing relevant changed.

Tab-gated callbacks take ``{prefix}-single-cell-tabs`` as an Input, so they fire
on every tab switch even when the plot's real parameters are unchanged. Pairing
``rendered_key_store`` (a per-session dcc.Store) with ``signature`` lets a
callback short-circuit with ``dash.no_update`` in that case — skipping both the
expensive recompute and the browser redraw. Per-session (not server-global) so
it stays correct for the multi-user browser deployment.
"""

import hashlib
import json

from dash import dcc


def signature(*parts) -> str:
    """Stable, compact hash of the parameters that define a figure."""
    payload = json.dumps(parts, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def rendered_key_store(prefix: str, plot: str) -> dcc.Store:
    """Per-session store holding the signature of the figure currently shown."""
    return dcc.Store(id=f"{prefix}-{plot}-rendered-key")
