"""Small shared helpers for interactive Cytoscape plots."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any


def _selector_value(value: Any) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def adjacent_highlight_stylesheet(
    base: list[dict[str, Any]],
    edges: Iterable[dict[str, Any]],
    selected_node: str,
) -> list[dict[str, Any]]:
    """Dim a graph except for one node, its neighbors, and connecting edges."""
    neighbors = {selected_node}
    selected_edges = []
    for edge in edges:
        if edge.get("source") == selected_node or edge.get("target") == selected_node:
            neighbors.update((str(edge.get("source")), str(edge.get("target"))))
            selected_edges.append(str(edge.get("id")))

    overlay = [
        {"selector": "node", "style": {"opacity": 0.12}},
        {"selector": "edge", "style": {"opacity": 0.05}},
    ]
    overlay.extend(
        {
            "selector": f'node[id = "{_selector_value(node_id)}"]',
            "style": {"opacity": 1, "z-index": 9999},
        }
        for node_id in neighbors
    )
    overlay.extend(
        {
            "selector": f'edge[id = "{_selector_value(edge_id)}"]',
            "style": {"opacity": 0.95, "z-index": 9999},
        }
        for edge_id in selected_edges
    )
    return base + overlay
