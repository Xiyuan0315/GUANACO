"""External interaction sources used by the exploratory Network plot.

The rest of the application consumes one small, provider-independent edge model.
HTTP details, provider field names, timeouts, and response caching stay here so a
database outage cannot leak into the Cytoscape renderer or layout code.
"""

from __future__ import annotations

from collections import OrderedDict
import csv
from dataclasses import asdict, dataclass
import io
import json
import threading
import time
from typing import Any, Callable, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import requests


STRING_API = "https://version-12-0.string-db.org/api/json/network"
KEGG_REST_API = "https://rest.kegg.jp"
MIRTARBASE_QUERY_API = "https://awi.cuhk.edu.cn/miRTarBase/search/results/download/"
TRRUST_URLS = {
    "human": "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    "mouse": "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv",
}

ORGANISM_TAXON_IDS = {
    "human": 9606,
    "mouse": 10090,
    "rat": 10116,
}
KEGG_ORGANISM_CODES = {"human": "hsa", "mouse": "mmu", "rat": "rno"}
MIRTARBASE_SPECIES_NAMES = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "rat": "Rattus norvegicus",
}

NETWORK_TYPES = {"ppi", "tf-gene", "metabolite", "mirna"}


class NetworkSourceError(RuntimeError):
    """A remote interaction source could not satisfy a network query."""


@dataclass(frozen=True)
class NetworkEdge:
    source: str
    target: str
    source_type: str = "gene"
    target_type: str = "gene"
    directed: bool = False
    effect: str = "unknown"
    score: float | None = None
    database: str = ""
    references: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class _TTLCache:
    def __init__(self, max_items: int = 96, ttl_seconds: float = 24 * 60 * 60):
        self.max_items = int(max_items)
        self.ttl_seconds = float(ttl_seconds)
        self._items: OrderedDict[tuple, tuple[float, tuple[NetworkEdge, ...]]] = (
            OrderedDict()
        )
        self._lock = threading.Lock()

    def get(self, key: tuple) -> tuple[NetworkEdge, ...] | None:
        with self._lock:
            item = self._items.get(key)
            if item is None:
                return None
            created_at, value = item
            if time.time() - created_at >= self.ttl_seconds:
                self._items.pop(key, None)
                return None
            self._items.move_to_end(key)
            return value

    def set(self, key: tuple, value: Iterable[NetworkEdge]) -> tuple[NetworkEdge, ...]:
        stored = tuple(value)
        with self._lock:
            self._items[key] = (time.time(), stored)
            self._items.move_to_end(key)
            while len(self._items) > self.max_items:
                self._items.popitem(last=False)
        return stored

    def clear(self) -> None:
        with self._lock:
            self._items.clear()


_QUERY_CACHE = _TTLCache()
_TRRUST_TABLE_CACHE = _TTLCache(max_items=2)
_MIRTARBASE_GENE_CACHE = _TTLCache(max_items=512)
_KEGG_GENE_MAP_CACHE: dict[
    str, tuple[float, dict[str, str], dict[str, tuple[str, ...]]]
] = {}
_KEGG_GENE_MAP_LOCK = threading.Lock()


def normalize_organism(value: str | int | None) -> tuple[str, int]:
    """Return the canonical name and NCBI taxonomy id for a supported organism."""
    if isinstance(value, int) or (isinstance(value, str) and value.strip().isdigit()):
        taxon_id = int(value)
        for name, candidate in ORGANISM_TAXON_IDS.items():
            if candidate == taxon_id:
                return name, taxon_id
    normalized = str(value or "human").strip().lower().replace("_", " ")
    aliases = {
        "homo sapiens": "human",
        "hsa": "human",
        "hg19": "human",
        "hg38": "human",
        "mus musculus": "mouse",
        "mm10": "mouse",
        "mm39": "mouse",
        "rattus norvegicus": "rat",
        "rn6": "rat",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in ORGANISM_TAXON_IDS:
        supported = ", ".join(ORGANISM_TAXON_IDS)
        raise ValueError(f"Unsupported organism '{value}'. Choose one of: {supported}.")
    return normalized, ORGANISM_TAXON_IDS[normalized]


def _request_json(
    url: str,
    params: dict[str, Any],
    *,
    method: str = "GET",
    timeout: float = 10,
) -> list[dict[str, Any]]:
    encoded = urlencode(params).encode("utf-8")
    if method == "POST":
        request = Request(
            url,
            data=encoded,
            headers={
                "Accept": "application/json",
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": "guanaco-viz/1.0",
            },
            method="POST",
        )
    else:
        request = Request(
            f"{url}?{encoded.decode('utf-8')}",
            headers={"Accept": "application/json", "User-Agent": "guanaco-viz/1.0"},
        )
    try:
        with urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        if exc.code == 404:
            return []
        raise NetworkSourceError(
            f"Interaction database returned HTTP {exc.code}."
        ) from exc
    except (URLError, TimeoutError) as exc:
        raise NetworkSourceError(
            "Interaction database is unavailable or timed out."
        ) from exc
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise NetworkSourceError(
            "Interaction database returned an invalid response."
        ) from exc
    if isinstance(payload, dict):
        payload = payload.get("data", payload.get("results", []))
    if not isinstance(payload, list):
        raise NetworkSourceError(
            "Interaction database returned an unexpected response shape."
        )
    return [item for item in payload if isinstance(item, dict)]


def _request_text(url: str, *, timeout: float = 15) -> str:
    request = Request(
        url,
        headers={
            "Accept": "text/tab-separated-values",
            "User-Agent": "guanaco-viz/1.0",
        },
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        raise NetworkSourceError(
            f"Interaction database returned HTTP {exc.code}."
        ) from exc
    except (URLError, TimeoutError) as exc:
        raise NetworkSourceError(
            "Interaction database is unavailable or timed out."
        ) from exc
    except UnicodeDecodeError as exc:
        raise NetworkSourceError(
            "Interaction database returned an invalid response."
        ) from exc


def _request_mirtarbase_text(url: str) -> str:
    try:
        response = requests.get(
            url,
            timeout=(15, 45),
            headers={"Accept": "text/csv", "User-Agent": "guanaco-viz/1.0"},
        )
        response.raise_for_status()
        return response.content.decode("utf-8-sig")
    except (requests.RequestException, UnicodeDecodeError) as exc:
        raise NetworkSourceError("miRTarBase is unavailable or timed out.") from exc


def _as_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result


def _string_edges(
    genes: tuple[str, ...],
    taxon_id: int,
    request_json: Callable[..., list[dict[str, Any]]],
    ppi_network_type: str = "functional",
) -> list[NetworkEdge]:
    resolved_network_type = str(ppi_network_type or "functional").strip().lower()
    if resolved_network_type not in {"functional", "physical"}:
        raise ValueError("STRING network type must be functional or physical.")
    rows = request_json(
        STRING_API,
        {
            "identifiers": "\r".join(genes),
            "species": taxon_id,
            "network_type": resolved_network_type,
            # Cache medium-confidence and stronger interactions so the UI can
            # change its display threshold without another STRING request.
            "required_score": 400,
            "add_nodes": 10,
            "show_query_node_labels": 1,
            "caller_identity": "guanaco-viz",
        },
        method="POST",
    )
    edges = []
    for row in rows:
        source = row.get("preferredName_A") or row.get("stringId_A")
        target = row.get("preferredName_B") or row.get("stringId_B")
        if not source or not target:
            continue
        edges.append(
            NetworkEdge(
                source=str(source),
                target=str(target),
                score=_as_float(row.get("score")),
                database=f"STRING v12 {resolved_network_type}",
            )
        )
    return edges


def _trrust_table(
    organism: str,
    request_text: Callable[..., str],
    *,
    use_cache: bool,
) -> list[NetworkEdge]:
    if organism not in TRRUST_URLS:
        raise ValueError(
            "TRRUST networks currently support human and mouse datasets only."
        )

    cache_key = ("trrust-v2", organism)
    if use_cache and request_text is _request_text:
        cached = _TRRUST_TABLE_CACHE.get(cache_key)
        if cached is not None:
            return list(cached)

    rows = request_text(TRRUST_URLS[organism])
    effects = {"activation": "activation", "repression": "inhibition"}
    edges = []
    for line in rows.splitlines():
        fields = line.strip().split("\t")
        if len(fields) != 4:
            continue
        source, target, effect, reference = fields
        if not source or not target:
            continue
        edges.append(
            NetworkEdge(
                source=source,
                target=target,
                directed=True,
                effect=effects.get(effect.strip().lower(), "unknown"),
                database="TRRUST v2",
                references=reference,
            )
        )
    if not edges:
        raise NetworkSourceError("TRRUST returned no valid TF–gene interactions.")

    if use_cache and request_text is _request_text:
        return list(_TRRUST_TABLE_CACHE.set(cache_key, edges))
    return edges


def _trrust_edges(
    genes: tuple[str, ...],
    organism: str,
    request_text: Callable[..., str],
    *,
    use_cache: bool,
) -> list[NetworkEdge]:
    query = {gene.casefold() for gene in genes}
    return [
        edge
        for edge in _trrust_table(organism, request_text, use_cache=use_cache)
        if edge.source.casefold() in query or edge.target.casefold() in query
    ]


def fetch_trrust_table(
    organism: str | int | None = "human",
    *,
    request_text: Callable[..., str] = _request_text,
    use_cache: bool = True,
) -> list[NetworkEdge]:
    """Return the complete cached TRRUST v2 table for enrichment analysis."""
    organism_name, _ = normalize_organism(organism)
    return _trrust_table(organism_name, request_text, use_cache=use_cache)


def _kegg_gene_maps(
    organism_code: str,
    request_text: Callable[..., str],
    *,
    use_cache: bool,
) -> tuple[dict[str, str], dict[str, tuple[str, ...]]]:
    if use_cache and request_text is _request_text:
        with _KEGG_GENE_MAP_LOCK:
            cached = _KEGG_GENE_MAP_CACHE.get(organism_code)
            if cached is not None and time.time() - cached[0] < 24 * 60 * 60:
                return cached[1], cached[2]

    primary: dict[str, str] = {}
    aliases: dict[str, set[str]] = {}
    for line in request_text(f"{KEGG_REST_API}/list/{organism_code}").splitlines():
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        gene_id, description = fields[0], fields[-1]
        symbols = description.split(";", 1)[0]
        names = [name.strip() for name in symbols.split(",") if name.strip()]
        if not gene_id or not names:
            continue
        primary.setdefault(names[0].casefold(), gene_id)
        for name in names:
            aliases.setdefault(name.casefold(), set()).add(gene_id)

    frozen_aliases = {name: tuple(sorted(ids)) for name, ids in aliases.items()}
    if use_cache and request_text is _request_text:
        with _KEGG_GENE_MAP_LOCK:
            _KEGG_GENE_MAP_CACHE[organism_code] = (time.time(), primary, frozen_aliases)
    return primary, frozen_aliases


def _kegg_link(
    target: str,
    source_ids: Iterable[str],
    request_text: Callable[..., str],
) -> dict[str, set[str]]:
    source_ids = list(dict.fromkeys(source_ids))
    links: dict[str, set[str]] = {}
    for start in range(0, len(source_ids), 10):
        identifiers = "+".join(source_ids[start : start + 10])
        rows = request_text(f"{KEGG_REST_API}/link/{target}/{identifiers}")
        for line in rows.splitlines():
            fields = line.split("\t", 1)
            if len(fields) == 2:
                links.setdefault(fields[0], set()).add(fields[1])
    return links


def _kegg_compound_names(
    compound_ids: Iterable[str],
    request_text: Callable[..., str],
) -> dict[str, str]:
    compound_ids = list(dict.fromkeys(compound_ids))
    names = {}
    for start in range(0, len(compound_ids), 10):
        identifiers = "+".join(compound_ids[start : start + 10])
        rows = request_text(f"{KEGG_REST_API}/list/{identifiers}")
        for line in rows.splitlines():
            fields = line.split("\t", 1)
            if len(fields) != 2:
                continue
            compound_id = fields[0].removeprefix("cpd:")
            names[compound_id] = fields[1].split(";", 1)[0].strip()
    return names


def _kegg_metabolite_edges(
    genes: tuple[str, ...],
    organism: str,
    request_text: Callable[..., str],
    *,
    use_cache: bool,
) -> list[NetworkEdge]:
    organism_code = KEGG_ORGANISM_CODES[organism]
    primary, aliases = _kegg_gene_maps(
        organism_code,
        request_text,
        use_cache=use_cache,
    )
    gene_labels: dict[str, str] = {}
    for gene in genes:
        folded = gene.casefold()
        gene_id = primary.get(folded)
        if gene_id is None:
            candidates = aliases.get(folded, ())
            gene_id = candidates[0] if len(candidates) == 1 else None
        if gene_id is not None:
            gene_labels[gene_id] = gene

    gene_to_ec = _kegg_link("ec", gene_labels, request_text)
    ec_to_reaction = _kegg_link(
        "reaction",
        sorted({ec for values in gene_to_ec.values() for ec in values}),
        request_text,
    )
    reaction_to_compound = _kegg_link(
        "compound",
        sorted(
            {
                reaction
                for ec_values in gene_to_ec.values()
                for ec in ec_values
                for reaction in ec_to_reaction.get(ec, ())
            }
        ),
        request_text,
    )

    pair_reactions: dict[tuple[str, str], set[str]] = {}
    for gene_id, ec_values in gene_to_ec.items():
        gene = gene_labels[gene_id]
        for ec in ec_values:
            for reaction in ec_to_reaction.get(ec, ()):
                for compound in reaction_to_compound.get(reaction, ()):
                    pair_reactions.setdefault((compound, gene), set()).add(reaction)

    compound_names = _kegg_compound_names(
        sorted({compound for compound, _ in pair_reactions}),
        request_text,
    )
    edges = []
    for (compound, gene), reactions in sorted(pair_reactions.items()):
        compound_id = compound.removeprefix("cpd:")
        name = compound_names.get(compound_id, compound_id)
        edges.append(
            NetworkEdge(
                source=f"{name} ({compound_id})",
                target=gene,
                source_type="metabolite",
                target_type="gene",
                directed=False,
                database="KEGG reactions",
                references=", ".join(
                    reaction.removeprefix("rn:") for reaction in sorted(reactions)
                ),
            )
        )
    return edges


def _mirtarbase_gene_edges(
    gene: str,
    organism: str,
    request_text: Callable[[str], str],
    *,
    use_cache: bool,
) -> list[NetworkEdge]:
    cache_key = ("mirtarbase-query", organism, gene.casefold())
    if use_cache and request_text is _request_mirtarbase_text:
        cached = _MIRTARBASE_GENE_CACHE.get(cache_key)
        if cached is not None:
            return list(cached)

    query = urlencode(
        {
            "mode": "target",
            "keyword": gene,
            "species": KEGG_ORGANISM_CODES[organism],
        }
    )
    text = request_text(f"{MIRTARBASE_QUERY_API}?{query}")
    evidence: dict[tuple[str, str], dict[str, Any]] = {}
    expected_species = MIRTARBASE_SPECIES_NAMES[organism]
    for row in csv.DictReader(io.StringIO(text.lstrip("\ufeff"))):
        if (
            str(row.get("Species (miRNA)") or "").strip() != expected_species
            or str(row.get("Species (Target)") or "").strip() != expected_species
        ):
            continue
        mirna = str(row.get("miRNA") or "").strip()
        target = str(row.get("Target Gene") or "").strip()
        if not mirna or target.casefold() != gene.casefold():
            continue
        target = gene
        item = evidence.setdefault(
            (mirna, target),
            {"score": 0.0, "references": set()},
        )
        tier = str(row.get("Score") or "").strip()
        tier_score = {"A": 4.0, "B": 3.0, "C": 2.0, "D": 1.0}.get(
            tier.removeprefix("Tier ")[:1],
            1.0,
        )
        item["score"] = max(item["score"], tier_score)
        mti_id = str(row.get("MTI ID") or "").strip()
        if mti_id:
            item["references"].add(mti_id)

    edges = [
        NetworkEdge(
            source=mirna,
            target=target,
            source_type="mirna",
            target_type="gene",
            directed=True,
            score=item["score"],
            database="miRTarBase",
            references=", ".join(sorted(item["references"])),
        )
        for (mirna, target), item in evidence.items()
    ]
    if use_cache and request_text is _request_mirtarbase_text:
        return list(_MIRTARBASE_GENE_CACHE.set(cache_key, edges))
    return edges


def _mirtarbase_edges(
    genes: tuple[str, ...],
    organism: str,
    request_text: Callable[[str], str],
    *,
    use_cache: bool,
) -> list[NetworkEdge]:
    return [
        edge
        for gene in genes
        for edge in _mirtarbase_gene_edges(
            gene,
            organism,
            request_text,
            use_cache=use_cache,
        )
    ]


def fetch_network_edges(
    network_type: str,
    genes: Iterable[str],
    organism: str | int | None = "human",
    *,
    request_json: Callable[..., list[dict[str, Any]]] = _request_json,
    request_text: Callable[..., str] = _request_text,
    mirtarbase_request: Callable[[str], str] = _request_mirtarbase_text,
    ppi_network_type: str = "functional",
    use_cache: bool = True,
) -> list[NetworkEdge]:
    """Fetch one-hop interactions and normalize them to :class:`NetworkEdge`."""
    normalized_type = str(network_type).strip().lower()
    if normalized_type not in NETWORK_TYPES:
        raise ValueError(f"Unsupported network type: {network_type}")
    normalized_genes = tuple(
        dict.fromkeys(str(gene).strip() for gene in genes if str(gene).strip())
    )
    if not normalized_genes:
        raise ValueError("Enter at least one gene.")
    organism_name, taxon_id = normalize_organism(organism)
    if normalized_type == "tf-gene" and organism_name not in TRRUST_URLS:
        raise ValueError(
            "TRRUST networks currently support human and mouse datasets only."
        )

    resolved_ppi_type = str(ppi_network_type or "functional").strip().lower()
    if normalized_type == "ppi" and resolved_ppi_type not in {"functional", "physical"}:
        raise ValueError("STRING network type must be functional or physical.")

    cache_key = (
        normalized_type,
        organism_name,
        resolved_ppi_type if normalized_type == "ppi" else None,
        tuple(sorted(g.casefold() for g in normalized_genes)),
    )
    if normalized_type == "mirna":
        default_transport = mirtarbase_request is _request_mirtarbase_text
    elif normalized_type in {"tf-gene", "metabolite"}:
        default_transport = request_text is _request_text
    else:
        default_transport = request_json is _request_json
    if use_cache and default_transport:
        cached = _QUERY_CACHE.get(cache_key)
        if cached is not None:
            return list(cached)

    if normalized_type == "ppi":
        edges = _string_edges(
            normalized_genes,
            taxon_id,
            request_json,
            resolved_ppi_type,
        )
    elif normalized_type == "tf-gene":
        edges = _trrust_edges(
            normalized_genes,
            organism_name,
            request_text,
            use_cache=use_cache,
        )
    elif normalized_type == "metabolite":
        edges = _kegg_metabolite_edges(
            normalized_genes,
            organism_name,
            request_text,
            use_cache=use_cache,
        )
    else:
        edges = _mirtarbase_edges(
            normalized_genes,
            organism_name,
            mirtarbase_request,
            use_cache=use_cache,
        )

    if use_cache and default_transport:
        return list(_QUERY_CACHE.set(cache_key, edges))
    return edges


def clear_network_source_cache() -> None:
    _QUERY_CACHE.clear()
    _TRRUST_TABLE_CACHE.clear()
    _MIRTARBASE_GENE_CACHE.clear()
    with _KEGG_GENE_MAP_LOCK:
        _KEGG_GENE_MAP_CACHE.clear()
