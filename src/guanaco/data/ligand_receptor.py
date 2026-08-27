"""Adapters for precomputed ligand–receptor interaction results."""

from __future__ import annotations

import base64
import json
from collections.abc import Mapping, Sequence
from io import BytesIO
from pathlib import Path
from typing import Any
import weakref

import numpy as np
import pandas as pd


IDENTITY_COLUMNS = ("source", "target", "ligand", "receptor")

_IDENTITY_ALIASES = {
    "source": ("source", "sender", "source_group", "source_cell", "celltype_a"),
    "target": ("target", "receiver", "target_group", "target_cell", "celltype_b"),
    "ligand": (
        "ligand_complex",
        "ligand",
        "ligand_name",
        "gene_a",
        "partner_a",
    ),
    "receptor": (
        "receptor_complex",
        "receptor",
        "receptor_name",
        "gene_b",
        "partner_b",
    ),
}
_OPTIONAL_TEXT_ALIASES = {
    "sample": ("sample", "sample_id", "dataset"),
    "condition": ("condition", "group", "contrast"),
    "method": ("method", "resource", "database"),
}
_MAGNITUDE_PRIORITY = (
    "magnitude_rank",
    "lr_means",
    "prob",
    "communication_probability",
    "interaction_score",
    "score",
    "strength",
    "expr_prod",
    "scaled_weight",
    "lrscore",
    "mean",
)
_SPECIFICITY_PRIORITY = (
    "specificity_rank",
    "cellphone_pvals",
    "pvalue",
    "pvalues",
    "p_val",
    "pval",
    "fdr",
    "qvalue",
    "spec_weight",
)


class LigandReceptorResultError(ValueError):
    """Raised when a supplied result cannot be converted to interaction rows."""


_DEFAULT_RESULT_CACHE: dict[
    int,
    tuple[weakref.ReferenceType, tuple[tuple[str, int], ...], dict[str, Any] | None],
] = {}


def _result_fingerprint(adata) -> tuple[tuple[str, int], ...]:
    return tuple(
        (str(key), id(value))
        for key, value in adata.uns.items()
        if _can_adapt(value)
    )


def _cached_default_result(adata, fingerprint):
    entry = _DEFAULT_RESULT_CACHE.get(id(adata))
    if entry is None:
        return False, None
    reference, cached_fingerprint, payload = entry
    if reference() is adata and cached_fingerprint == fingerprint:
        return True, payload
    _DEFAULT_RESULT_CACHE.pop(id(adata), None)
    return False, None


def _store_default_result(adata, fingerprint, payload) -> None:
    key = id(adata)

    def discard(reference):
        entry = _DEFAULT_RESULT_CACHE.get(key)
        if entry is not None and entry[0] is reference:
            _DEFAULT_RESULT_CACHE.pop(key, None)

    _DEFAULT_RESULT_CACHE[key] = (weakref.ref(adata, discard), fingerprint, payload)


def _normalized_name(value: Any) -> str:
    return (
        str(value)
        .strip()
        .lower()
        .replace(" ", "_")
        .replace(".", "_")
        .replace("-", "_")
    )


def _column_lookup(frame: pd.DataFrame) -> dict[str, Any]:
    return {_normalized_name(column): column for column in frame.columns}


def _find_column(frame: pd.DataFrame, aliases: Sequence[str]):
    lookup = _column_lookup(frame)
    return next((lookup[alias] for alias in aliases if alias in lookup), None)


def _metric_direction(name: str) -> str:
    normalized = _normalized_name(name)
    lower_tokens = ("pval", "p_value", "pvalues", "fdr", "qvalue", "rank")
    return "lower" if any(token in normalized for token in lower_tokens) else "higher"


def _default_metric(metrics: list[str], priorities: Sequence[str]) -> str | None:
    lookup = {_normalized_name(metric): metric for metric in metrics}
    return next((lookup[name] for name in priorities if name in lookup), None)


def _canonical_long_table(
    frame: pd.DataFrame,
    *,
    name: str,
    source_format: str,
) -> dict[str, Any]:
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise LigandReceptorResultError("The ligand–receptor result table is empty.")

    identity_sources = {
        name: _find_column(frame, aliases)
        for name, aliases in _IDENTITY_ALIASES.items()
    }
    missing = [name for name, column in identity_sources.items() if column is None]
    if missing:
        raise LigandReceptorResultError(
            "The interaction table is missing required columns: "
            + ", ".join(missing)
            + "."
        )

    canonical = pd.DataFrame(index=frame.index)
    for name, column in identity_sources.items():
        canonical[name] = frame[column].astype("string").str.strip()

    optional_sources = {
        name: _find_column(frame, aliases)
        for name, aliases in _OPTIONAL_TEXT_ALIASES.items()
    }
    for name, column in optional_sources.items():
        if column is not None:
            canonical[name] = frame[column].astype("string").str.strip()

    used_columns = {
        column
        for column in (*identity_sources.values(), *optional_sources.values())
        if column is not None
    }
    metrics = []
    for column in frame.columns:
        if column in used_columns:
            continue
        numeric = pd.to_numeric(frame[column], errors="coerce")
        if numeric.notna().any():
            metric = str(column)
            canonical[metric] = numeric
            metrics.append(metric)

    valid = np.ones(len(canonical), dtype=bool)
    for column in IDENTITY_COLUMNS:
        values = canonical[column]
        valid &= values.notna().to_numpy()
        valid &= values.ne("").to_numpy()
    canonical = canonical.loc[valid].reset_index(drop=True)
    if canonical.empty:
        raise LigandReceptorResultError(
            "No complete source–target ligand–receptor rows were found."
        )
    if not metrics:
        raise LigandReceptorResultError(
            "The interaction table must contain at least one numeric score."
        )

    default_magnitude = _default_metric(metrics, _MAGNITUDE_PRIORITY) or metrics[0]
    default_specificity = _default_metric(metrics, _SPECIFICITY_PRIORITY)
    if default_specificity == default_magnitude:
        default_specificity = None

    records = json.loads(canonical.to_json(orient="records"))
    return {
        "name": name,
        "format": source_format,
        "records": records,
        "metrics": [
            {
                "label": metric,
                "value": metric,
                "direction": _metric_direction(metric),
            }
            for metric in metrics
        ],
        "default_magnitude": default_magnitude,
        "default_specificity": default_specificity,
    }


def _is_squidpy_mapping(value: Mapping[str, Any]) -> bool:
    means = value.get("means")
    return isinstance(means, pd.DataFrame) and (
        isinstance(means.index, pd.MultiIndex)
        or isinstance(means.columns, pd.MultiIndex)
    )


def _adapt_squidpy(
    value: Mapping[str, Any],
    *,
    name: str,
) -> dict[str, Any]:
    means = value.get("means")
    pvalues = value.get("pvalues")
    if not isinstance(means, pd.DataFrame):
        raise LigandReceptorResultError(
            "A Squidpy result must contain a 'means' DataFrame."
        )
    if not isinstance(pvalues, pd.DataFrame):
        pvalues = None

    rows = []
    for row_index, interaction in enumerate(means.index):
        interaction_parts = (
            tuple(interaction)
            if isinstance(interaction, tuple)
            else tuple(str(interaction).split("|", maxsplit=1))
        )
        if len(interaction_parts) < 2:
            continue
        for column_index, group_pair in enumerate(means.columns):
            group_parts = (
                tuple(group_pair)
                if isinstance(group_pair, tuple)
                else tuple(str(group_pair).split("|", maxsplit=1))
            )
            if len(group_parts) < 2:
                continue
            magnitude = means.iloc[row_index, column_index]
            specificity = (
                pvalues.iloc[row_index, column_index]
                if pvalues is not None
                and row_index < pvalues.shape[0]
                and column_index < pvalues.shape[1]
                else np.nan
            )
            if pd.isna(magnitude) and pd.isna(specificity):
                continue
            rows.append(
                {
                    "source": group_parts[0],
                    "target": group_parts[1],
                    "ligand": interaction_parts[0],
                    "receptor": interaction_parts[1],
                    "mean": magnitude,
                    "pvalue": specificity,
                }
            )
    return _canonical_long_table(
        pd.DataFrame(rows),
        name=name,
        source_format="Squidpy",
    )


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def _adapt_cellphonedb(
    value: Mapping[str, Any],
    *,
    name: str,
) -> dict[str, Any]:
    means = value.get("means")
    if not isinstance(means, pd.DataFrame):
        means = value.get("significant_means")
    pvalues = value.get("pvalues")
    if not isinstance(means, pd.DataFrame):
        raise LigandReceptorResultError(
            "A CellPhoneDB result must contain a 'means' DataFrame."
        )
    if not isinstance(pvalues, pd.DataFrame):
        pvalues = None

    lookup = _column_lookup(means)
    gene_a = lookup.get("gene_a") or lookup.get("partner_a")
    gene_b = lookup.get("gene_b") or lookup.get("partner_b")
    interacting_pair = lookup.get("interacting_pair")
    receptor_a = lookup.get("receptor_a")
    receptor_b = lookup.get("receptor_b")
    pair_columns = [
        column
        for column in means.columns
        if "|" in str(column)
        and (pvalues is None or column in pvalues.columns)
    ]
    if not pair_columns:
        raise LigandReceptorResultError(
            "No CellPhoneDB cell-group pair columns were found."
        )

    rows = []
    for row_index, row in means.iterrows():
        ligand = row.get(gene_a) if gene_a is not None else None
        receptor = row.get(gene_b) if gene_b is not None else None
        if (pd.isna(ligand) or pd.isna(receptor)) and interacting_pair is not None:
            parts = str(row.get(interacting_pair, "")).split("|", maxsplit=1)
            if len(parts) == 2:
                ligand, receptor = parts
        reverse = (
            receptor_a is not None
            and _truthy(row.get(receptor_a))
            and not (receptor_b is not None and _truthy(row.get(receptor_b)))
        )
        for pair_column in pair_columns:
            groups = str(pair_column).split("|", maxsplit=1)
            if len(groups) != 2:
                continue
            source, target = groups
            row_ligand, row_receptor = ligand, receptor
            if reverse:
                source, target = target, source
                row_ligand, row_receptor = receptor, ligand
            rows.append(
                {
                    "source": source,
                    "target": target,
                    "ligand": row_ligand,
                    "receptor": row_receptor,
                    "mean": row[pair_column],
                    "pvalue": (
                        pvalues.loc[row_index, pair_column]
                        if pvalues is not None
                        and row_index in pvalues.index
                        and pair_column in pvalues.columns
                        else np.nan
                    ),
                }
            )
    return _canonical_long_table(
        pd.DataFrame(rows),
        name=name,
        source_format="CellPhoneDB",
    )


def _mapping_table(value: Mapping[str, Any], *, name: str) -> dict[str, Any]:
    if _is_squidpy_mapping(value):
        return _adapt_squidpy(value, name=name)
    return _adapt_cellphonedb(value, name=name)


def _can_adapt(value: Any) -> str | None:
    if isinstance(value, pd.DataFrame):
        columns = _column_lookup(value)
        if all(
            any(alias in columns for alias in aliases)
            for aliases in _IDENTITY_ALIASES.values()
        ):
            return "interaction table"
    if isinstance(value, Mapping):
        means = value.get("means")
        if not isinstance(means, pd.DataFrame):
            means = value.get("significant_means")
        if isinstance(means, pd.DataFrame):
            return "Squidpy" if _is_squidpy_mapping(value) else "CellPhoneDB"
    return None


def discover_ligand_receptor_results(adata) -> list[dict[str, str]]:
    """Return compatible precomputed results stored in ``adata.uns``."""
    if adata is None:
        return []
    options = []
    for key, value in adata.uns.items():
        source_format = _can_adapt(value)
        if source_format:
            options.append(
                {
                    "label": f"{key} · {source_format}",
                    "value": str(key),
                }
            )
    return sorted(options, key=lambda option: option["label"].casefold())


def load_ligand_receptor_result(adata, key: str) -> dict[str, Any]:
    """Normalize one embedded precomputed result."""
    if adata is None or key not in adata.uns:
        raise LigandReceptorResultError(f"No result named {key!r} was found.")
    value = adata.uns[key]
    if isinstance(value, pd.DataFrame):
        return _canonical_long_table(
            value,
            name=str(key),
            source_format="interaction table",
        )
    if isinstance(value, Mapping):
        return _mapping_table(value, name=str(key))
    raise LigandReceptorResultError(
        f"The result {key!r} is not a supported interaction table."
    )


def load_default_ligand_receptor_result(adata) -> dict[str, Any] | None:
    """Load ``liana_res`` when available, otherwise the first compatible result."""
    if adata is None:
        return None
    fingerprint = _result_fingerprint(adata)
    found, cached = _cached_default_result(adata, fingerprint)
    if found:
        return cached
    keys = [option["value"] for option in discover_ligand_receptor_results(adata)]
    if "liana_res" in keys:
        keys.remove("liana_res")
        keys.insert(0, "liana_res")
    for key in keys:
        try:
            payload = load_ligand_receptor_result(adata, key)
            _store_default_result(adata, fingerprint, payload)
            return payload
        except LigandReceptorResultError:
            continue
    _store_default_result(adata, fingerprint, None)
    return None


def _read_upload(contents: str, filename: str) -> pd.DataFrame:
    try:
        encoded = contents.split(",", maxsplit=1)[1]
        raw = base64.b64decode(encoded)
    except (IndexError, ValueError) as exc:
        raise LigandReceptorResultError(
            f"Could not decode uploaded file {filename!r}."
        ) from exc

    suffix = Path(filename).suffix.lower()
    try:
        if suffix == ".parquet":
            return pd.read_parquet(BytesIO(raw))
        if suffix in {".tsv", ".txt"}:
            return pd.read_csv(BytesIO(raw), sep="\t")
        if suffix == ".csv":
            return pd.read_csv(BytesIO(raw))
    except Exception as exc:
        raise LigandReceptorResultError(
            f"Could not read uploaded file {filename!r}: {exc}"
        ) from exc
    raise LigandReceptorResultError(
        "Upload a CSV, TSV, TXT, or Parquet interaction result."
    )


def parse_uploaded_ligand_receptor_results(
    contents: str | Sequence[str],
    filenames: str | Sequence[str],
) -> dict[str, Any]:
    """Normalize one long table or a paired CellPhoneDB upload."""
    content_list = [contents] if isinstance(contents, str) else list(contents or [])
    filename_list = (
        [filenames] if isinstance(filenames, str) else list(filenames or [])
    )
    if not content_list or len(content_list) != len(filename_list):
        raise LigandReceptorResultError("Select one or more result files to upload.")

    tables = {
        filename: _read_upload(content, filename)
        for content, filename in zip(content_list, filename_list, strict=True)
    }
    if len(tables) == 1:
        filename, table = next(iter(tables.items()))
        try:
            return _canonical_long_table(
                table,
                name=filename,
                source_format="uploaded interaction table",
            )
        except LigandReceptorResultError:
            return _adapt_cellphonedb(
                {"means": table},
                name=filename,
            )

    cellphonedb_tables: dict[str, pd.DataFrame] = {}
    for filename, table in tables.items():
        normalized = _normalized_name(Path(filename).stem)
        if "pvalue" in normalized:
            cellphonedb_tables["pvalues"] = table
        elif "significant_mean" in normalized:
            cellphonedb_tables["significant_means"] = table
        elif "mean" in normalized:
            cellphonedb_tables["means"] = table
    if "means" not in cellphonedb_tables and "significant_means" not in cellphonedb_tables:
        raise LigandReceptorResultError(
            "Multiple-file upload currently expects CellPhoneDB means and pvalues files."
        )
    return _adapt_cellphonedb(
        cellphonedb_tables,
        name="CellPhoneDB upload",
    )


def ligand_receptor_frame(payload: Mapping[str, Any]) -> pd.DataFrame:
    """Restore a normalized payload as a DataFrame."""
    records = payload.get("records") if isinstance(payload, Mapping) else None
    if not isinstance(records, list):
        raise LigandReceptorResultError("No normalized interaction data are loaded.")
    return pd.DataFrame.from_records(records)


def metric_direction(payload: Mapping[str, Any], metric: str | None) -> str:
    for descriptor in payload.get("metrics", []):
        if descriptor.get("value") == metric:
            return str(descriptor.get("direction", "higher"))
    return "higher"
