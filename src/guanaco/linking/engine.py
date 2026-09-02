"""Pure reducer that routes source identity snapshots to target views."""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from typing import Any

from .model import LinkSpec, MarkMembers, ViewState


def reduce_members(
    current: Mapping[str, Any] | None,
    *,
    source_view: str,
    members: MarkMembers | None,
) -> dict[str, Any]:
    """Store one source event as the workspace's sole selection state."""

    sources = dict((current or {}).get("sources") or {})
    if members is None:
        sources.pop(source_view, None)
    else:
        updated = 1 + max(
            (int(item.get("updated", 0)) for item in sources.values()),
            default=0,
        )
        sources[source_view] = {
            "members": members.to_dict(),
            "updated": updated,
        }
    return {"sources": sources}


def state_for_target(
    state: Mapping[str, Any] | None,
    target: str,
    links: Sequence[LinkSpec],
    table_sources: Collection[str] = (),
    row_keys: Mapping[tuple[str, str | None], Mapping[str, str]] | None = None,
) -> ViewState:
    """Resolve the latest selection and action on each target identity axis.

    ``row_keys`` optionally maps each table's atomic row ID to a logical row
    key.  This keeps the table index unique while allowing one overview row to
    update many child rows (for example, one interaction pair to many spots or
    distance bins).
    """

    sources = (state or {}).get("sources") or {}
    latest: dict[str, tuple[tuple[int, int], tuple[str, ...], str | None]] = {}
    for order, link in enumerate(links):
        if link.target != target:
            continue
        snapshot = sources.get(link.source)
        if not isinstance(snapshot, Mapping):
            continue
        members = MarkMembers.from_dict(snapshot.get("members") or {})
        source_map = (row_keys or {}).get((link.source, link.key)) or None
        target_map = (row_keys or {}).get((link.target, link.key)) or None
        if link.source in table_sources:
            selected_rows = members.rows
            selected_ids = tuple(
                dict.fromkeys(
                    source_map.get(row, row) if source_map is not None else row
                    for row in selected_rows
                )
            )
        else:
            selected_ids = members.project(link.selection_by)

        # A table endpoint may contain atomic rows (for example one row per
        # spot/pair membership) while its link key represents the identities
        # consumed by an AnnData view. Expand the selected logical IDs back to
        # all matching table rows; for an AnnData target the logical IDs are
        # already its native obs_names/var_names.
        if target_map is not None:
            selected = set(selected_ids)
            ids = tuple(row for row, key in target_map.items() if key in selected)
        else:
            ids = selected_ids
        if not ids:
            continue
        rank = (int(snapshot.get("updated", 0)), order)
        previous = latest.get(link.selection_by)
        if previous is None or rank >= previous[0]:
            latest[link.selection_by] = (rank, ids, link.resolved_action)

    values: dict[str, tuple[str, ...] | None] = {
        "rows": None,
        "cells": None,
        "features": None,
    }
    highlights: list[str] = []
    for by, (_rank, ids, action) in latest.items():
        values[f"{by}s"] = ids
        if action == "highlight":
            highlights.append(by)
    return ViewState(**values, highlight_axes=tuple(highlights))
