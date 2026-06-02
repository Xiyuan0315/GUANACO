from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def save_pydeseq_results_to_adata_uns(
    adata: Any,
    results: pd.DataFrame | str | Path,
    entry_name: str,
    *,
    key: str = "volcano",
    gene_name_col: str = "gene_name",
    gene_id_col: str | None = None,
    logfoldchange_col: str = "log2FoldChange",
    pvalue_col: str = "padj",
    score_col: str = "stat",
    default_entry: bool = True,
) -> dict[str, Any]:
    """Store PyDESeq2-style results in ``adata.uns`` for GUANACO volcano plots.

    Parameters
    ----------
    adata
        AnnData-like object with an ``uns`` mapping.
    results
        PyDESeq2 results as a pandas DataFrame or CSV path.
    entry_name
        Contrast/group name shown in GUANACO's volcano dropdown.
    key
        Destination key in ``adata.uns``. GUANACO reads ``adata.uns["volcano"]``.
    gene_name_col
        Preferred gene label column. Missing labels fall back to gene IDs or row index.
    gene_id_col
        Optional gene ID column. If omitted, a non-default DataFrame index is used.
    logfoldchange_col
        Column containing log fold change values, usually ``log2FoldChange``.
    pvalue_col
        Column containing adjusted p-values or p-values, usually ``padj``.
    score_col
        Optional score/statistic column, usually ``stat``.
    default_entry
        Whether to set this entry as the default volcano dropdown selection.

    Returns
    -------
    dict
        The payload stored at ``adata.uns[key]``.
    """
    if not hasattr(adata, "uns"):
        raise TypeError("adata must have an 'uns' mapping.")

    df = _coerce_results_dataframe(results)
    _require_columns(df, [logfoldchange_col, pvalue_col])

    reserved_columns = {gene_name_col, logfoldchange_col, pvalue_col, score_col}
    gene_ids, gene_id_source = _resolve_gene_ids(df, gene_id_col, reserved_columns)
    genes = _resolve_gene_labels(df, gene_name_col, gene_ids)
    logfoldchanges = pd.to_numeric(df[logfoldchange_col], errors="coerce")
    pvalues = pd.to_numeric(df[pvalue_col], errors="coerce")
    scores = (
        pd.to_numeric(df[score_col], errors="coerce")
        if score_col in df.columns
        else pd.Series(np.nan, index=df.index)
    )
    if gene_id_source is not None:
        reserved_columns.add(gene_id_source)
    extras = _extra_columns(df, gene_ids, reserved=reserved_columns)

    entry: dict[str, Any] = {
        "gene": genes,
        "logfoldchange": _series_to_clean_list(logfoldchanges),
        "score": _series_to_clean_list(scores),
        "padj": _series_to_clean_list(pvalues),
        "meta": {
            "source": "pydeseq",
            "entry": entry_name,
            "pval_key": pvalue_col,
            "logfoldchange_key": logfoldchange_col,
            "score_key": score_col if score_col in df.columns else None,
        },
    }
    entry.update(extras)

    payload = adata.uns.setdefault(key, {"entries": {}})
    payload.setdefault("entries", {})
    payload["source_type"] = "volcano"
    payload["entries"][entry_name] = entry
    if default_entry or not payload.get("default_entry"):
        payload["default_entry"] = entry_name
    return payload


def _coerce_results_dataframe(results: pd.DataFrame | str | Path) -> pd.DataFrame:
    if isinstance(results, pd.DataFrame):
        return results.copy()
    path = Path(results)
    return pd.read_csv(path)


def _require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise KeyError(f"Missing required PyDESeq2 result column(s): {missing}")


def _resolve_gene_ids(
    df: pd.DataFrame,
    gene_id_col: str | None,
    reserved_columns: set[str],
) -> tuple[list[str], str | None]:
    if gene_id_col and gene_id_col in df.columns:
        return _values_to_strings(df[gene_id_col]), gene_id_col
    if not isinstance(df.index, pd.RangeIndex):
        return _values_to_strings(df.index), None

    for candidate in ("Ensembl_gene_id_v", "gene_id", "Geneid", "geneid", "id"):
        if candidate in df.columns and candidate not in reserved_columns:
            return _values_to_strings(df[candidate]), candidate

    if len(df.columns) > 0:
        first_col = str(df.columns[0])
        if first_col not in reserved_columns and pd.api.types.is_object_dtype(df[first_col]):
            return _values_to_strings(df[first_col]), first_col

    return _values_to_strings(df.index), None


def _resolve_gene_labels(df: pd.DataFrame, gene_name_col: str, gene_ids: list[str]) -> list[str]:
    if gene_name_col not in df.columns:
        return gene_ids

    labels = []
    for label, gene_id in zip(df[gene_name_col], gene_ids, strict=True):
        if pd.isna(label) or str(label).strip() == "":
            labels.append(gene_id)
        else:
            labels.append(str(label))
    return labels


def _extra_columns(df: pd.DataFrame, gene_ids: list[str], reserved: set[str]) -> dict[str, list[Any]]:
    extras: dict[str, list[Any]] = {"gene_id": gene_ids}
    for column in df.columns:
        if column in reserved:
            continue
        extras[str(column)] = _series_to_clean_list(df[column])
    return extras


def _series_to_clean_list(series: pd.Series) -> list[Any]:
    cleaned = series.astype(object).where(pd.notna(series), None)
    return cleaned.tolist()


def _values_to_strings(values) -> list[str]:
    return [str(value) for value in values]
