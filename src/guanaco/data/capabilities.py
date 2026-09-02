"""Fast, metadata-only detection of visualizations supported by a dataset."""

from __future__ import annotations

import re
from collections.abc import Mapping

import numpy as np
import pandas as pd

from guanaco.data.capability_schema import PLOT_KEYS
from guanaco.data.ligand_receptor import discover_ligand_receptor_results

_PEAK_RE = re.compile(r"^\s*([^:\s]+)\s*:\s*([0-9,]+)\s*-\s*([0-9,]+)\s*$")
_NHOOD_SUFFIX = "_nhood_enrichment"
_CO_OCCURRENCE_SUFFIX = "_co_occurrence"


def _obs_col(obs, column):
    values = obs[column]
    return values.to_series() if hasattr(values, "to_series") else values


def _discrete_labels(adata, *, max_unique=50) -> list[str]:
    obs = getattr(adata, "obs", None)
    if obs is None or getattr(obs, "shape", (0,))[0] == 0:
        return []
    labels = []
    for column in obs.columns:
        dtype = obs.dtypes[column]
        if isinstance(dtype, pd.CategoricalDtype):
            cardinality = len(dtype.categories)
        elif (
            pd.api.types.is_object_dtype(dtype)
            or pd.api.types.is_string_dtype(dtype)
            or pd.api.types.is_bool_dtype(dtype)
        ):
            cardinality = int(_obs_col(obs, column).nunique(dropna=True))
        else:
            continue
        if cardinality < max_unique:
            labels.append(str(column))
    return labels


def _has_feature_matrix(adata) -> bool:
    return bool(
        adata is not None
        and getattr(adata, "n_obs", 0) > 0
        and getattr(adata, "n_vars", 0) > 0
    )


def _has_continuous_observation(adata) -> bool:
    if adata is None or getattr(adata, "obs", None) is None:
        return False
    supported_dtypes = {"float32", "float64", "int32", "int64"}
    return any(
        str(adata.obs.dtypes[column]) in supported_dtypes
        and _obs_col(adata.obs, column).nunique(dropna=True) > 50
        for column in adata.obs.columns
    )


def _has_genomic_peak_features(
    adata, *, sample_size: int = 2000, min_hits: int = 5
) -> bool:
    if adata is None or getattr(adata, "n_vars", 0) == 0:
        return False
    hits = sum(
        _PEAK_RE.match(str(name)) is not None
        for name in adata.var_names[: min(sample_size, adata.n_vars)]
    )
    var = getattr(adata, "var", None)
    has_coordinates = var is not None and {"chrom", "start", "end"}.issubset(
        set(var.columns)
    )
    return hits >= min_hits or has_coordinates


def _has_spatial_relationships(adata) -> bool:
    if adata is None:
        return False
    for uns_key in adata.uns:
        if not str(uns_key).endswith(_NHOOD_SUFFIX):
            continue
        cluster_key = str(uns_key)[: -len(_NHOOD_SUFFIX)]
        if cluster_key not in adata.obs.columns:
            continue
        values = _obs_col(adata.obs, cluster_key)
        if not isinstance(values.dtype, pd.CategoricalDtype):
            continue
        n_clusters = len(values.cat.categories)
        nhood = adata.uns.get(str(uns_key))
        occurrence = adata.uns.get(f"{cluster_key}{_CO_OCCURRENCE_SUFFIX}")
        if not isinstance(nhood, Mapping) or not isinstance(occurrence, Mapping):
            continue
        zscore = np.asarray(nhood.get("zscore"))
        count = np.asarray(nhood.get("count"))
        occ = np.asarray(occurrence.get("occ"))
        interval = np.asarray(occurrence.get("interval"))
        if (
            n_clusters > 0
            and zscore.shape == (n_clusters, n_clusters)
            and count.shape == zscore.shape
            and occ.ndim == 3
            and occ.shape[:2] == zscore.shape
            and interval.ndim == 1
            and occ.shape[2] == max(0, interval.size - 1)
        ):
            return True
    return False


def has_plot_capability(
    key: str,
    adata,
    *,
    has_igv: bool = False,
    modality_name: str | None = None,
    feature_data_available: bool | None = None,
    discrete_data_available: bool | None = None,
) -> bool:
    """Return whether one canonical plot can produce a useful data view."""
    if key == "igv":
        return has_igv
    if key in {"multiomics-composition", "cross-modal-concordance"}:
        return modality_name == "multiomics"
    if adata is None:
        return False

    has_features = (
        _has_feature_matrix(adata)
        if feature_data_available is None
        else feature_data_available
    )
    has_discrete = (
        bool(_discrete_labels(adata))
        if discrete_data_available is None
        else discrete_data_available
    )
    if key in {"dotplot", "heatmap", "violin", "ridge"}:
        return has_features and has_discrete
    if key == "pseudotime":
        return has_features and _has_continuous_observation(adata)
    if key in {"split-violin", "stacked-bar"}:
        return has_discrete
    if key == "ligand-receptor":
        return bool(discover_ligand_receptor_results(adata))
    if key == "peak-browser":
        return _has_genomic_peak_features(adata)
    if key == "paga":
        paga = adata.uns.get("paga")
        return isinstance(paga, Mapping) and "connectivities" in paga
    if key == "volcano":
        return "volcano" in adata.uns or "rank_genes_groups" in adata.uns
    if key == "spatial-relationships":
        return _has_spatial_relationships(adata)
    return key in PLOT_KEYS and has_features


def detect_plot_capabilities(
    adata,
    *,
    has_igv: bool = False,
    modality_name: str | None = None,
    feature_data_available: bool | None = None,
    discrete_data_available: bool | None = None,
) -> tuple[str, ...]:
    """Return every compatible plot key in stable display order."""
    return tuple(
        key
        for key in PLOT_KEYS
        if has_plot_capability(
            key,
            adata,
            has_igv=has_igv,
            modality_name=modality_name,
            feature_data_available=feature_data_available,
            discrete_data_available=discrete_data_available,
        )
    )
