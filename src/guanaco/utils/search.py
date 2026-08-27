"""Small helpers for predictable dropdown search ordering."""

from __future__ import annotations

from collections.abc import Iterable


def ranked_substring_matches(
    values: Iterable[str],
    query: str | None,
    *,
    limit: int = 10,
    match_values: Iterable[str] | None = None,
) -> list[str]:
    """Return exact, prefix, then substring matches while preserving source order.

    ``match_values`` allows a display label such as ``RNA · CD4`` to be ranked
    by its underlying feature name (``CD4``).
    """
    if limit <= 0:
        return []

    labels = list(values)
    candidates = labels if match_values is None else list(match_values)
    if len(labels) != len(candidates):
        raise ValueError("values and match_values must have the same length")

    needle = str(query or "").strip().casefold()
    if not needle:
        return labels[:limit]

    exact: list[str] = []
    prefixes: list[str] = []
    substrings: list[str] = []
    for label, candidate in zip(labels, candidates):
        normalized = str(candidate).casefold()
        if normalized == needle:
            bucket = exact
        elif normalized.startswith(needle):
            bucket = prefixes
        elif needle in normalized:
            bucket = substrings
        else:
            continue
        if len(bucket) < limit:
            bucket.append(label)

    return (exact + prefixes + substrings)[:limit]
