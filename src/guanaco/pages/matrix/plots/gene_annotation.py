"""GTF/GFF3 gene annotation indexing for the ATAC browser."""

from __future__ import annotations

import gzip
import hashlib
import os
import pickle
import re
import shutil
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np


_GTF_ATTR_RE = re.compile(r'\s*([^";\s]+)\s+"([^"]*)"')

# IGV-style "type a genome id and we fetch the reference for you". Each id maps to
# a GENCODE *basic* annotation (small, ~50-100 MB gzipped); the file is downloaded
# once into a shared cache and then read by the normal local GTF path. Aliases let
# users write either the UCSC-style id (hg38) or the assembly name (GRCh38).
_GENCODE = "https://ftp.ebi.ac.uk/pub/databases/gencode"
GENOME_ANNOTATION_REGISTRY: dict[str, str] = {
    "hg38": f"{_GENCODE}/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz",
    "hg19": f"{_GENCODE}/Gencode_human/release_50/GRCh37_mapping/gencode.v50lift37.basic.annotation.gtf.gz",
    "mm39": f"{_GENCODE}/Gencode_mouse/release_M39/gencode.vM39.basic.annotation.gtf.gz",
    "mm10": f"{_GENCODE}/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz",
}
_GENOME_ALIASES: dict[str, str] = {
    "hg38": "hg38", "grch38": "hg38",
    "hg19": "hg19", "grch37": "hg19",
    "mm39": "mm39", "grcm39": "mm39",
    "mm10": "mm10", "grcm38": "mm10",
}


@dataclass
class GeneAnnotationIndex:
    path: str
    chroms: tuple[str, ...]
    starts: dict[str, np.ndarray]
    ends: dict[str, np.ndarray]
    gene_ids: dict[str, np.ndarray]
    genes: dict[str, dict]
    transcripts_by_gene: dict[str, list[dict]]
    exons_by_transcript: dict[str, list[tuple[int, int]]]
    gene_name_to_ids: dict[str, list[str]]


_annotation_cache: dict[tuple[str, float, int], GeneAnnotationIndex] = {}
_region_cache: dict[tuple[str, str, int, int, int, int], list[dict]] = {}


def _annotation_cache_dir() -> Path:
    """A writable, machine-wide cache for downloaded references (downloaded once)."""
    candidates: list[Path] = []
    env_dir = os.environ.get("GUANACO_CACHE_DIR")
    if env_dir:
        candidates.append(Path(env_dir) / "annotations")
    candidates.append(Path.home() / ".cache" / "guanaco" / "annotations")
    candidates.append(Path(tempfile.gettempdir()) / "guanaco_gene_annotation_cache")
    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            return candidate
        except Exception:
            continue
    fallback = Path(tempfile.gettempdir()) / "guanaco_gene_annotation_cache"
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def _download_to_cache(url: str) -> Path:
    """Fetch ``url`` into the cache once; reuse the local copy on later calls."""
    cache_dir = _annotation_cache_dir()
    digest = hashlib.sha1(url.encode("utf-8")).hexdigest()[:12]
    basename = url.split("?")[0].rstrip("/").split("/")[-1] or "annotation.gtf.gz"
    dest = cache_dir / f"{digest}-{basename}"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_name(dest.name + ".partial")
    try:
        with urllib.request.urlopen(url, timeout=120) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
        tmp.replace(dest)
    except Exception:
        try:
            tmp.unlink()
        except OSError:
            pass
        raise
    return dest


def is_known_genome_id(spec: str | None) -> bool:
    """True if ``spec`` is a registry genome id (hg38/GRCh38/mm10 ...)."""
    return bool(spec) and str(spec).strip().lower() in _GENOME_ALIASES


def resolve_annotation_source(spec: str | None) -> str | None:
    """Turn a config value into a local GTF/GFF3 path.

    Accepts three forms:
      * a known genome id (``hg38``/``GRCh38``/``mm10`` ...) -> download the matching
        GENCODE basic annotation into the cache and return its local path;
      * an ``http(s)``/``ftp`` URL -> download into the cache and return its path;
      * anything else -> treated as a local filesystem path (unchanged behaviour).
    """
    if not spec:
        return None
    text = str(spec).strip()
    alias = _GENOME_ALIASES.get(text.lower())
    if alias is not None:
        return str(_download_to_cache(GENOME_ANNOTATION_REGISTRY[alias]))
    if text.lower().startswith(("http://", "https://", "ftp://")):
        return str(_download_to_cache(text))
    return text


def _open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _parse_attrs(raw: str, *, is_gff3: bool = False) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if is_gff3:
        for field in raw.strip().split(";"):
            if not field or "=" not in field:
                continue
            key, value = field.split("=", 1)
            attrs[key.strip()] = value.strip()
        return attrs

    for field in raw.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        match = _GTF_ATTR_RE.match(field)
        if match:
            attrs[match.group(1)] = match.group(2)
            continue
        parts = field.split(maxsplit=1)
        if len(parts) == 2:
            attrs[parts[0]] = parts[1].strip('"')
    return attrs


def _chrom_sort_key(chrom: str):
    normalized = chrom.removeprefix("chr")
    if normalized.isdigit():
        return (0, int(normalized))
    if normalized == "X":
        return (0, 23)
    if normalized == "Y":
        return (0, 24)
    if normalized in {"M", "MT"}:
        return (0, 25)
    return (1, normalized)


def _cache_candidates(path: Path, stat) -> list[Path]:
    resolved = str(path.resolve())
    digest = hashlib.sha1(f"{resolved}:{stat.st_mtime_ns}:{stat.st_size}".encode("utf-8")).hexdigest()[:20]
    filename = f"{path.name}.{digest}.guanaco-gene-index.pkl"
    return [
        path.parent / ".guanaco_cache" / filename,
        Path(tempfile.gettempdir()) / "guanaco_gene_annotation_cache" / filename,
    ]


def _read_disk_cache(candidates: list[Path]) -> GeneAnnotationIndex | None:
    for candidate in candidates:
        try:
            if candidate.exists():
                with candidate.open("rb") as handle:
                    return pickle.load(handle)
        except Exception:
            continue
    return None


def _write_disk_cache(candidates: list[Path], index: GeneAnnotationIndex) -> None:
    for candidate in candidates:
        try:
            candidate.parent.mkdir(parents=True, exist_ok=True)
            with candidate.open("wb") as handle:
                pickle.dump(index, handle, protocol=pickle.HIGHEST_PROTOCOL)
            return
        except Exception:
            continue


def load_gene_annotation(path: str | None) -> GeneAnnotationIndex | None:
    if not path:
        return None
    resolved_source = resolve_annotation_source(path)
    if resolved_source is None:
        return None
    if str(resolved_source).lower().startswith(("s3://", "gs://", "gcs://")):
        raise ValueError(
            "Object-store paths are not supported; use a local GTF/GFF3 file, an http(s) URL, "
            "or a genome id like 'hg38'."
        )
    expanded = Path(str(resolved_source)).expanduser()
    if not expanded.exists():
        raise FileNotFoundError(f"Gene annotation file not found: {expanded}")
    resolved = str(expanded.resolve())
    stat = expanded.stat()
    cache_key = (resolved, stat.st_mtime, stat.st_size)
    cached = _annotation_cache.get(cache_key)
    if cached is not None:
        return cached
    cache_candidates = _cache_candidates(expanded, stat)
    disk_cached = _read_disk_cache(cache_candidates)
    if disk_cached is not None:
        _annotation_cache[cache_key] = disk_cached
        return disk_cached

    genes: dict[str, dict] = {}
    transcripts_by_gene: dict[str, list[dict]] = {}
    exons_by_transcript: dict[str, list[tuple[int, int]]] = {}
    gene_name_to_ids: dict[str, list[str]] = {}

    lower_name = expanded.name.lower()
    is_gff3 = lower_name.endswith((".gff3", ".gff3.gz", ".gff", ".gff.gz"))
    with _open_text(str(expanded)) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start_raw, end_raw, score, strand, phase, attr_raw = parts[:9]
            feature = feature.lower()
            if feature not in {"gene", "transcript", "mrna", "exon"}:
                continue
            try:
                start = int(start_raw)
                end = int(end_raw)
            except ValueError:
                continue
            attrs = _parse_attrs(attr_raw, is_gff3=is_gff3)

            gene_id = (
                attrs.get("gene_id")
                or attrs.get("gene")
                or attrs.get("Parent")
                or attrs.get("ID")
            )
            transcript_id = (
                attrs.get("transcript_id")
                or attrs.get("transcript")
                or attrs.get("ID")
            )
            if feature in {"transcript", "mrna", "exon"} and attrs.get("Parent") and not attrs.get("gene_id"):
                parent = attrs.get("Parent", "")
                if feature == "exon":
                    transcript_id = transcript_id or parent
                else:
                    gene_id = parent.split(",")[0]

            gene_name = attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene") or gene_id
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype") or attrs.get("biotype") or ""
            transcript_name = attrs.get("transcript_name") or attrs.get("Name") or transcript_id
            transcript_type = attrs.get("transcript_type") or attrs.get("transcript_biotype") or gene_type
            tags = attrs.get("tag", "")

            if feature == "gene":
                if not gene_id:
                    continue
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "source": source,
                }
                gene_name_to_ids.setdefault(str(gene_name).lower(), []).append(gene_id)
            elif feature in {"transcript", "mrna"}:
                if not gene_id or not transcript_id:
                    continue
                transcripts_by_gene.setdefault(gene_id, []).append(
                    {
                        "transcript_id": transcript_id,
                        "transcript_name": transcript_name,
                        "transcript_type": transcript_type,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "is_basic": "basic" in tags,
                    }
                )
            elif feature == "exon":
                if not transcript_id:
                    continue
                exons_by_transcript.setdefault(transcript_id, []).append((start, end))

    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for gene_id, gene in genes.items():
        by_chrom.setdefault(gene["chrom"], []).append((int(gene["start"]), int(gene["end"]), gene_id))

    starts: dict[str, np.ndarray] = {}
    ends: dict[str, np.ndarray] = {}
    gene_ids: dict[str, np.ndarray] = {}
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda row: (row[0], row[1]))
        starts[chrom] = np.asarray([r[0] for r in rows], dtype=np.int64)
        ends[chrom] = np.asarray([r[1] for r in rows], dtype=np.int64)
        gene_ids[chrom] = np.asarray([r[2] for r in rows], dtype=object)

    for exons in exons_by_transcript.values():
        exons.sort()
    for txs in transcripts_by_gene.values():
        txs.sort(key=lambda tx: (not tx["is_basic"], tx["start"], tx["end"], tx["transcript_id"]))

    index = GeneAnnotationIndex(
        path=resolved,
        chroms=tuple(sorted(by_chrom.keys(), key=_chrom_sort_key)),
        starts=starts,
        ends=ends,
        gene_ids=gene_ids,
        genes=genes,
        transcripts_by_gene=transcripts_by_gene,
        exons_by_transcript=exons_by_transcript,
        gene_name_to_ids=gene_name_to_ids,
    )
    _annotation_cache[cache_key] = index
    _write_disk_cache(cache_candidates, index)
    return index


def query_gene_models(
    index: GeneAnnotationIndex | None,
    chrom: str,
    start: int,
    end: int,
    *,
    max_genes: int = 120,
    max_transcripts_per_gene: int = 3,
) -> list[dict]:
    if index is None or chrom not in index.starts:
        return []
    cache_key = (index.path, chrom, int(start), int(end), max_genes, max_transcripts_per_gene)
    cached = _region_cache.get(cache_key)
    if cached is not None:
        return cached

    starts = index.starts[chrom]
    stop = np.searchsorted(starts, end, side="right")
    candidates = np.arange(stop, dtype=np.int64)
    if candidates.size:
        candidates = candidates[index.ends[chrom][candidates] >= start]
    if candidates.size > max_genes:
        # Sample evenly across the window rather than taking the first N by
        # coordinate -- otherwise a zoomed-out view only shows genes clustered at
        # the left edge and looks broken.
        take = np.linspace(0, candidates.size - 1, max_genes).round().astype(np.int64)
        candidates = candidates[np.unique(take)]

    models = []
    for pos in candidates:
        gene_id = str(index.gene_ids[chrom][pos])
        gene = dict(index.genes[gene_id])
        tx_models = []
        for tx in index.transcripts_by_gene.get(gene_id, [])[:max_transcripts_per_gene]:
            tx_copy = dict(tx)
            exons = index.exons_by_transcript.get(tx_copy["transcript_id"], [])
            tx_copy["exons"] = [(s, e) for s, e in exons if e >= start and s <= end]
            tx_models.append(tx_copy)
        gene["transcripts"] = tx_models
        models.append(gene)

    _region_cache[cache_key] = models
    return models


def find_gene_region(index: GeneAnnotationIndex | None, query: str, *, flank: int = 500) -> dict | None:
    if index is None or not query:
        return None
    key = str(query).strip().lower()
    gene_ids = index.gene_name_to_ids.get(key)
    if not gene_ids:
        bare_key = key.split(".")[0]
        gene_ids = [
            gene_id for gene_id, gene in index.genes.items()
            if gene_id.lower() == key or gene_id.split(".")[0].lower() == bare_key
        ]
    if not gene_ids:
        return None
    gene = index.genes[gene_ids[0]]
    # Frame the gene by the extent of the transcripts we actually index/draw, not by
    # the gene record's start/end. In a GENCODE *basic* file the gene line keeps the
    # comprehensive boundaries, so a single long non-basic isoform can make the gene
    # record span far wider than any transcript present (CD96: 388 kb record vs a
    # 110 kb basic-transcript extent), which framed the search far too zoomed-out.
    transcripts = index.transcripts_by_gene.get(gene_ids[0], [])
    if transcripts:
        start = min(int(tx["start"]) for tx in transcripts)
        end = max(int(tx["end"]) for tx in transcripts)
    else:
        start = int(gene["start"])
        end = int(gene["end"])
    return {
        "chrom": gene["chrom"],
        "start": max(0, start - flank),
        "end": end + flank,
    }
