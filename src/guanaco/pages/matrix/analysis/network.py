"""Provider-independent construction of expression-annotated network payloads."""

from __future__ import annotations

import re
from typing import Any, Iterable

import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from guanaco.data.network_sources import NetworkEdge
from guanaco.utils.gene_extraction_utils import (
    GENE_SYMBOL_COLUMNS,
    extract_gene_expression,
    gene_symbol_lookup,
    prewarm_gene_cache,
)


MAX_QUERY_GENES = 100
MAX_NETWORK_NODES = 200
MAX_NETWORK_EDGES = 800
NETWORK_VIEW_MODES = {"first-order", "input-only", "minimum"}
MIRNA_REGULATOR_FILTERS = {
    "shared-20": (20, 2),
    "shared-50": (50, 2),
    "shared-100": (100, 2),
    "all": (None, 1),
}


def parse_gene_list(value: str | Iterable[str] | None) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        tokens = re.split(r"[,;\s]+", value)
    else:
        tokens = [str(item) for item in value]
    genes = []
    seen = set()
    for token in tokens:
        gene = str(token).strip()
        folded = gene.casefold()
        if gene and folded not in seen:
            genes.append(gene)
            seen.add(folded)
    return genes


def validate_gene_list(value: str | Iterable[str] | None) -> list[str]:
    genes = parse_gene_list(value)
    if not genes:
        raise ValueError("Enter at least one gene.")
    if len(genes) > MAX_QUERY_GENES:
        raise ValueError(f"Enter no more than {MAX_QUERY_GENES} genes.")
    return genes


def _background_gene_symbols(adata) -> set[str]:
    var_names = adata.var_names.astype(str)
    symbol_column = next(
        (column for column in GENE_SYMBOL_COLUMNS if column in adata.var.columns),
        None,
    )
    if symbol_column is None:
        return {name.casefold() for name in var_names if name}

    background = set()
    for var_name, symbol in zip(var_names, adata.var[symbol_column], strict=False):
        text = "" if symbol is None else str(symbol).strip()
        if not text or text.lower() == "nan":
            text = var_name
        if text:
            background.add(text.casefold())
    return background


def trrust_key_regulator_enrichment(
    adata,
    query_genes: Iterable[str],
    trrust_edges: Iterable[NetworkEdge],
) -> dict[str, Any]:
    """Prioritize TRRUST TFs by one-sided target-set over-representation."""
    query_genes = validate_gene_list(query_genes)
    background = _background_gene_symbols(adata)
    query_labels = {gene.casefold(): gene for gene in query_genes}
    valid_query = set(query_labels).intersection(background)
    result = {
        "background_size": len(background),
        "input_query_count": len(query_genes),
        "valid_query_count": len(valid_query),
        "results": [],
        "reason": "",
    }
    if len(valid_query) < 5:
        result["reason"] = "Enter at least 5 genes measured in this RNA matrix."
        return result

    targets_by_tf: dict[str, set[str]] = {}
    target_labels: dict[str, str] = {}
    tf_labels: dict[str, str] = {}
    for edge in trrust_edges:
        tf = edge.source.casefold()
        target = edge.target.casefold()
        if target not in background:
            continue
        tf_labels.setdefault(tf, edge.source)
        target_labels.setdefault(target, edge.target)
        targets_by_tf.setdefault(tf, set()).add(target)

    candidates = []
    p_values = []
    for tf, targets in targets_by_tf.items():
        overlap = targets.intersection(valid_query)
        if len(targets) < 5 or len(overlap) < 2:
            continue
        p_value = float(
            hypergeom.sf(
                len(overlap) - 1,
                len(background),
                len(targets),
                len(valid_query),
            )
        )
        candidates.append(
            {
                "tf": tf_labels[tf],
                "overlap_count": len(overlap),
                "target_count": len(targets),
                "p_value": p_value,
                "overlap_genes": [target_labels[target] for target in sorted(overlap)],
            }
        )
        p_values.append(p_value)

    adjusted = multipletests(p_values, method="fdr_bh")[1] if p_values else []
    for candidate, fdr in zip(candidates, adjusted, strict=False):
        candidate["fdr"] = float(fdr)
    candidates.sort(
        key=lambda item: (
            item["fdr"],
            item["p_value"],
            -item["overlap_count"],
            item["tf"].casefold(),
        )
    )
    result["results"] = candidates
    if not candidates:
        result["reason"] = (
            "No TF has at least 5 measurable targets and 2 overlapping query genes."
        )
    return result


def mean_expression_by_gene(
    adata,
    labels: Iterable[str],
) -> tuple[dict[str, float], list[str]]:
    """Return mean expression keyed by requested gene label."""
    labels = list(dict.fromkeys(str(label) for label in labels))
    lookup = gene_symbol_lookup(adata)
    resolved = {label: lookup.get(label.casefold()) for label in labels}
    var_names = list(
        dict.fromkeys(value for value in resolved.values() if value is not None)
    )
    # h5py/backed sparse datasets require increasing fancy-index positions. The
    # shared prewarmer reads all columns in one slice, so meet that constraint
    # without changing the user-facing node order.
    positions = adata.var_names.get_indexer(var_names)
    var_names = [name for _, name in sorted(zip(positions, var_names, strict=False))]
    prewarm_gene_cache(adata, var_names, dtype=np.float32)
    means_by_var = {}
    for var_name in var_names:
        values = extract_gene_expression(
            adata, var_name, use_cache=True, dtype=np.float32
        )
        finite = values[np.isfinite(values)]
        means_by_var[var_name] = float(finite.mean()) if finite.size else 0.0
    measured_means = {
        label: means_by_var[var_name]
        for label, var_name in resolved.items()
        if var_name is not None
    }
    missing = [label for label, var_name in resolved.items() if var_name is None]
    return measured_means, missing


def _deduplicate_edges(edges: Iterable[NetworkEdge]) -> list[NetworkEdge]:
    best: dict[tuple, NetworkEdge] = {}
    for edge in edges:
        if not edge.source or not edge.target or edge.source == edge.target:
            continue
        endpoints = (edge.source, edge.target)
        if not edge.directed:
            endpoints = tuple(sorted(endpoints, key=str.casefold))
        key = (*endpoints, edge.directed, edge.effect)
        current = best.get(key)
        current_score = (
            float("-inf") if current is None or current.score is None else current.score
        )
        candidate_score = float("-inf") if edge.score is None else edge.score
        if current is None or candidate_score > current_score:
            best[key] = edge
    return list(best.values())


def _limit_graph(
    edges: list[NetworkEdge],
    query_genes: list[str],
    *,
    max_nodes: int,
    max_edges: int,
) -> tuple[list[NetworkEdge], bool]:
    query_folded = {gene.casefold() for gene in query_genes}
    node_scores: dict[str, float] = {}
    node_degrees: dict[str, int] = {}
    seed_neighbors: dict[str, set[str]] = {}
    canonical: dict[str, str] = {}
    for edge in edges:
        score = edge.score if edge.score is not None else 0.0
        endpoints = (edge.source, edge.target)
        for label in endpoints:
            folded = label.casefold()
            canonical.setdefault(folded, label)
            node_scores[folded] = max(score, node_scores.get(folded, float("-inf")))
            node_degrees[folded] = node_degrees.get(folded, 0) + 1
        source, target = (value.casefold() for value in endpoints)
        if source not in query_folded and target in query_folded:
            seed_neighbors.setdefault(source, set()).add(target)
        if target not in query_folded and source in query_folded:
            seed_neighbors.setdefault(target, set()).add(source)

    # Query genes belong to the graph even if no provider edge resolves them.
    # Reserve their node budget before ranking database-added neighbors.
    for gene in query_genes:
        folded = gene.casefold()
        canonical.setdefault(folded, gene)
        node_scores.setdefault(folded, float("inf"))

    retained = sorted(query_folded)
    max_nodes = max(max_nodes, len(retained))
    external = sorted(
        (folded for folded in canonical if folded not in query_folded),
        key=lambda folded: (
            -len(seed_neighbors.get(folded, set())),
            -node_degrees.get(folded, 0),
            -node_scores.get(folded, 0.0),
            canonical[folded].casefold(),
        ),
    )
    retained.extend(external[: max(0, max_nodes - len(retained))])
    retained_set = set(retained)
    limited = [
        edge
        for edge in edges
        if edge.source.casefold() in retained_set
        and edge.target.casefold() in retained_set
    ]
    limited.sort(
        key=lambda edge: (
            -(edge.score if edge.score is not None else 0.0),
            edge.source.casefold(),
            edge.target.casefold(),
        )
    )
    truncated = len(retained_set) < len(canonical) or len(limited) > max_edges
    return limited[:max_edges], truncated


def _minimum_network_edges(
    edges: list[NetworkEdge],
    query_genes: list[str],
) -> list[NetworkEdge]:
    """Approximate a seed-connecting network with a union of shortest paths."""
    adjacency: dict[str, list[tuple[str, int]]] = {}
    for index, edge in enumerate(edges):
        source = edge.source.casefold()
        target = edge.target.casefold()
        adjacency.setdefault(source, []).append((target, index))
        adjacency.setdefault(target, []).append((source, index))

    remaining = {gene.casefold() for gene in query_genes}
    retained_edges: set[int] = set()
    while remaining:
        root = min(remaining)
        remaining.remove(root)
        if root not in adjacency:
            continue

        parents: dict[str, tuple[str, int] | None] = {root: None}
        queue = [root]
        for node in queue:
            for neighbor, edge_index in adjacency.get(node, []):
                if neighbor in parents:
                    continue
                parents[neighbor] = (node, edge_index)
                queue.append(neighbor)

        connected_seeds = sorted(remaining.intersection(parents))
        for seed in connected_seeds:
            cursor = seed
            while cursor != root:
                parent = parents[cursor]
                if parent is None:
                    break
                cursor, edge_index = parent
                retained_edges.add(edge_index)
        remaining.difference_update(connected_seeds)

    return [edge for index, edge in enumerate(edges) if index in retained_edges]


def _filter_network_view(
    edges: list[NetworkEdge],
    query_genes: list[str],
    view_mode: str,
) -> list[NetworkEdge]:
    if view_mode not in NETWORK_VIEW_MODES:
        raise ValueError(f"Unsupported network view: {view_mode}")
    if view_mode == "first-order":
        return edges
    if view_mode == "minimum":
        return _minimum_network_edges(edges, query_genes)
    query = {gene.casefold() for gene in query_genes}
    return [
        edge
        for edge in edges
        if edge.source.casefold() in query and edge.target.casefold() in query
    ]


def _query_connected_edges(
    edges: list[NetworkEdge],
    query_genes: list[str],
) -> list[NetworkEdge]:
    """Remove database-only components disconnected from every input gene."""
    adjacency: dict[str, set[str]] = {}
    for edge in edges:
        source = edge.source.casefold()
        target = edge.target.casefold()
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set()).add(source)

    reachable = {
        gene.casefold() for gene in query_genes if gene.casefold() in adjacency
    }
    stack = list(reachable)
    while stack:
        node = stack.pop()
        for neighbor in adjacency.get(node, ()):
            if neighbor not in reachable:
                reachable.add(neighbor)
                stack.append(neighbor)

    return [
        edge
        for edge in edges
        if edge.source.casefold() in reachable and edge.target.casefold() in reachable
    ]


def _filter_mirna_regulators(
    edges: list[NetworkEdge],
    query_genes: list[str],
    regulator_filter: str,
) -> tuple[list[NetworkEdge], dict[str, int | None]]:
    limit, requested_minimum = MIRNA_REGULATOR_FILTERS.get(
        regulator_filter,
        MIRNA_REGULATOR_FILTERS["shared-20"],
    )
    minimum_targets = 1 if len(query_genes) == 1 else requested_minimum
    query = {gene.casefold() for gene in query_genes}
    targets_by_regulator: dict[str, set[str]] = {}
    evidence_by_regulator: dict[str, float] = {}
    canonical: dict[str, str] = {}
    for edge in edges:
        if edge.target.casefold() not in query:
            continue
        regulator = edge.source.casefold()
        canonical.setdefault(regulator, edge.source)
        targets_by_regulator.setdefault(regulator, set()).add(edge.target.casefold())
        evidence_by_regulator[regulator] = evidence_by_regulator.get(
            regulator, 0.0
        ) + float(edge.score or 0.0)

    eligible = [
        regulator
        for regulator, targets in targets_by_regulator.items()
        if len(targets) >= minimum_targets
    ]
    eligible.sort(
        key=lambda regulator: (
            -len(targets_by_regulator[regulator]),
            -evidence_by_regulator[regulator],
            canonical[regulator].casefold(),
        )
    )
    selected = set(eligible if limit is None else eligible[:limit])
    return (
        [edge for edge in edges if edge.source.casefold() in selected],
        {
            "eligible_regulator_count": len(eligible),
            "retained_regulator_count": len(selected),
            "mirna_min_targets": minimum_targets,
        },
    )


def build_network_payload(
    adata,
    query_genes: Iterable[str],
    network_type: str,
    edges: Iterable[NetworkEdge],
    *,
    view_mode: str = "first-order",
    max_nodes: int = MAX_NETWORK_NODES,
    max_edges: int = MAX_NETWORK_EDGES,
    mirna_regulator_filter: str = "shared-20",
) -> dict[str, Any]:
    query_genes = validate_gene_list(query_genes)
    deduplicated = _deduplicate_edges(edges)
    mirna_stats: dict[str, int | None] = {
        "eligible_regulator_count": 0,
        "retained_regulator_count": 0,
        "mirna_min_targets": 0,
    }
    if network_type == "mirna":
        deduplicated, mirna_stats = _filter_mirna_regulators(
            deduplicated,
            query_genes,
            mirna_regulator_filter,
        )
    view_edges = _query_connected_edges(
        _filter_network_view(deduplicated, query_genes, view_mode),
        query_genes,
    )
    available_nodes = {
        label.casefold() for edge in view_edges for label in (edge.source, edge.target)
    }
    available_nodes.update(gene.casefold() for gene in query_genes)
    limited_edges, truncated = _limit_graph(
        view_edges,
        query_genes,
        max_nodes=max_nodes,
        max_edges=max_edges,
    )
    degrees: dict[str, int] = {}
    for edge in limited_edges:
        for label in (edge.source, edge.target):
            folded = label.casefold()
            degrees[folded] = degrees.get(folded, 0) + 1
    query_lookup = {gene.casefold() for gene in query_genes}
    node_types: dict[str, str] = {}
    canonical: dict[str, str] = {}
    for edge in limited_edges:
        source_type = (
            "tf"
            if network_type == "tf-gene" and edge.source_type in {"gene", "protein"}
            else edge.source_type
        )
        endpoint_types = (
            (edge.source, source_type),
            (edge.target, edge.target_type),
        )
        for label, entity_type in endpoint_types:
            folded = label.casefold()
            canonical.setdefault(folded, label)
            existing = node_types.get(folded)
            priority = {
                "gene": 0,
                "protein": 0,
                "tf": 1,
                "mirna": 2,
                "microrna": 2,
                "metabolite": 2,
            }
            if existing is None or priority.get(entity_type, 0) > priority.get(
                existing, 0
            ):
                node_types[folded] = entity_type
    # Preserve query genes as isolated nodes when the provider returns no edge for
    # one of them. This makes unresolved queries visible instead of silently lost.
    for gene in query_genes:
        folded = gene.casefold()
        canonical.setdefault(folded, gene)
        node_types.setdefault(folded, "gene")

    expression_labels = [
        canonical[folded]
        for folded, entity_type in node_types.items()
        if entity_type
        not in {"mirna", "microrna", "metabolite", "small_molecule", "drug", "lipid"}
    ]
    expression, missing_expression = mean_expression_by_gene(adata, expression_labels)
    measured_values = np.asarray(list(expression.values()), dtype=float)
    if measured_values.size:
        expression_min, expression_max = np.nanpercentile(measured_values, [5, 95])
        if not np.isfinite(expression_min) or not np.isfinite(expression_max):
            expression_min, expression_max = 0.0, 1.0
        if expression_min == expression_max:
            expression_max = expression_min + 1.0
    else:
        expression_min, expression_max = 0.0, 1.0

    nodes = []
    for folded, label in sorted(canonical.items(), key=lambda item: item[1].casefold()):
        entity_type = node_types.get(folded, "gene")
        value = expression.get(label)
        degree = degrees.get(folded, 0)
        nodes.append(
            {
                "id": label,
                "label": label,
                "entity_type": entity_type,
                "is_query": folded in query_lookup,
                "measured": value is not None,
                "expression": value,
                "node_size": 10.0,
                "degree": degree,
            }
        )

    maximum_degree = max((int(node["degree"]) for node in nodes), default=0)
    for node in nodes:
        degree = int(node["degree"])
        if degree <= 0:
            node["node_size"] = 10.0
        elif maximum_degree <= 1:
            node["node_size"] = 18.0
        else:
            node["node_size"] = float(
                10.0 + 34.0 * np.log1p(degree) / np.log1p(maximum_degree)
            )

    edge_payload = []
    for index, edge in enumerate(limited_edges):
        edge_payload.append(
            {
                "id": f"edge-{index}",
                **edge.to_dict(),
            }
        )

    return {
        "network_type": network_type,
        "nodes": nodes,
        "edges": edge_payload,
        "expression_min": float(expression_min),
        "expression_max": float(expression_max),
        "missing_expression": missing_expression,
        "truncated": truncated,
        "available_edge_count": len(view_edges),
        "available_node_count": len(available_nodes),
        **mirna_stats,
    }
