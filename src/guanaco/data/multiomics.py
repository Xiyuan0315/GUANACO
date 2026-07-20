"""Lightweight, cell-aligned access to several MuData modalities.

The joint visualization never concatenates complete modality matrices.  It keeps
only a small shared ``obs``/``obsm`` frame in memory and materializes the feature
columns requested by the active plot.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from guanaco.data.loader import get_discrete_labels, obs_col
from guanaco.utils.embeddings import embedding_to_numpy, is_embedding_obsm
from guanaco.utils.gene_extraction_utils import (
    extract_gene_expression,
    prewarm_gene_cache,
)


JOINT_TAB_ID = "multiomics"


@dataclass(frozen=True)
class _FeatureRef:
    modality: str
    key: str


@dataclass(frozen=True)
class _ObsRef:
    storage: str
    key: str
    modality: str | None = None


def _modality_label(modality: str) -> str:
    return str(modality).upper()


def expose_global_obs_to_modalities(mdata) -> None:
    """Expose MuData-level metadata to modality views without persisting copies.

    Unqualified columns such as ``condition`` are shared by every modality.
    Qualified columns such as ``rna:cluster`` are exposed only to their owning
    modality as ``cluster``. Existing local columns always take precedence.
    """
    global_ids = pd.Index(mdata.obs_names.astype(str))
    if not global_ids.is_unique:
        raise ValueError("MuData global obs_names contain duplicate cell IDs.")

    modalities = tuple(str(name) for name in mdata.mod.keys())
    for modality in modalities:
        adata = mdata.mod[modality]
        modality_ids = pd.Index(adata.obs_names.astype(str))
        indexer = global_ids.get_indexer(modality_ids)
        if np.any(indexer < 0):
            continue

        for column in mdata.obs.columns:
            raw_column = str(column)
            owner, separator, local_column = raw_column.partition(":")
            if separator and owner in modalities:
                if owner != modality:
                    continue
                target = local_column
            else:
                target = raw_column

            if target in adata.obs.columns:
                continue
            values = obs_col(mdata.obs, raw_column).iloc[indexer].copy()
            values.index = adata.obs_names
            adata.obs[target] = values


def _same_values(left: pd.Series, right: pd.Series) -> bool:
    """Value comparison tolerant of equivalent categorical/string dtypes."""
    if len(left) != len(right):
        return False
    left_values = left.astype("string").fillna("<NA>").to_numpy()
    right_values = right.astype("string").fillna("<NA>").to_numpy()
    return bool(np.array_equal(left_values, right_values))


def _unique_label(candidate: str, occupied: set[str], suffix: str) -> str:
    if candidate not in occupied:
        return candidate
    labelled = f"{candidate} ({suffix})"
    counter = 2
    while labelled in occupied:
        labelled = f"{candidate} ({suffix} {counter})"
        counter += 1
    return labelled


class MultiOmicsSource:
    """A same-cell MuData view used by scatter and marker visualizations."""

    def __init__(self, mdata):
        self.mdata = mdata
        self.modalities = tuple(str(name) for name in mdata.mod.keys())
        if len(self.modalities) < 2:
            raise ValueError("A multi-omics view requires at least two modalities.")

        self.cell_ids = pd.Index(mdata.obs_names.astype(str), name=mdata.obs_names.name)
        if not self.cell_ids.is_unique:
            raise ValueError("MuData global obs_names contain duplicate cell IDs.")

        self._row_indexers: dict[str, np.ndarray] = {}
        for modality in self.modalities:
            mod_ids = pd.Index(mdata.mod[modality].obs_names.astype(str))
            if not mod_ids.is_unique:
                raise ValueError(f"Modality '{modality}' contains duplicate cell IDs.")
            indexer = mod_ids.get_indexer(self.cell_ids)
            if len(mod_ids) != len(self.cell_ids) or np.any(indexer < 0):
                raise ValueError(
                    f"Modality '{modality}' does not contain the complete MuData cell set."
                )
            self._row_indexers[modality] = indexer.astype(np.int64, copy=False)

        self._feature_refs: OrderedDict[str, _FeatureRef] = OrderedDict()
        self._obs_refs: OrderedDict[str, _ObsRef] = OrderedDict()
        self._embedding_refs: OrderedDict[str, tuple[str, str]] = OrderedDict()
        self._raw_feature_lookup: dict[tuple[str, str], str] = {}
        self._materialized_cache: OrderedDict[tuple[str, ...], ad.AnnData] = (
            OrderedDict()
        )
        self._materialized_cache_size = 8
        self._feature_score_cache: OrderedDict[tuple[str, ...], np.ndarray] = (
            OrderedDict()
        )
        self._feature_score_cache_size = 16

        self._build_feature_catalog()
        obs = self._build_obs_catalog()
        obsm = self._build_embedding_catalog()
        self.base_adata = ad.AnnData(
            X=np.empty((len(self.cell_ids), 0), dtype=np.float32),
            obs=obs,
            obsm=obsm,
        )
        self.base_adata.obs_names = self.cell_ids.copy()

    @property
    def label(self) -> str:
        return " + ".join(_modality_label(mod) for mod in self.modalities)

    @property
    def feature_names(self) -> list[str]:
        return list(self._feature_refs)

    @property
    def obs_names(self) -> list[str]:
        return list(self._obs_refs)

    @property
    def embedding_names(self) -> list[str]:
        return list(self._embedding_refs)

    def preferred_embedding(self, modality: str) -> str | None:
        candidates = [
            label
            for label, (owner, _raw_key) in self._embedding_refs.items()
            if owner == modality
        ]
        return next(
            (label for label in candidates if label.endswith("UMAP")),
            candidates[0] if candidates else None,
        )

    @property
    def discrete_obs_names(self) -> list[str]:
        return get_discrete_labels(self.base_adata)

    def _build_feature_catalog(self) -> None:
        occupied: set[str] = set()
        for modality in self.modalities:
            mod_label = _modality_label(modality)
            for raw_key in self.mdata.mod[modality].var_names.astype(str):
                label = _unique_label(f"{mod_label} · {raw_key}", occupied, "feature")
                occupied.add(label)
                self._feature_refs[label] = _FeatureRef(modality, str(raw_key))
                self._raw_feature_lookup[(modality.lower(), str(raw_key).lower())] = (
                    label
                )

    def _aligned_mod_obs(self, modality: str, column: str) -> pd.Series:
        series = obs_col(self.mdata.mod[modality].obs, column)
        aligned = series.iloc[self._row_indexers[modality]].copy()
        aligned.index = self.cell_ids
        return aligned

    def _add_obs(
        self,
        frame: pd.DataFrame,
        occupied: set[str],
        candidate: str,
        values: pd.Series,
        ref: _ObsRef,
    ) -> str:
        label = _unique_label(candidate, occupied, "metadata")
        occupied.add(label)
        values = values.copy()
        values.index = self.cell_ids
        frame[label] = values
        self._obs_refs[label] = ref
        return label

    def _build_obs_catalog(self) -> pd.DataFrame:
        frame = pd.DataFrame(index=self.cell_ids)
        occupied = set(self._feature_refs)
        modality_columns: dict[tuple[str, str], tuple[str, pd.Series]] = {}

        for modality in self.modalities:
            for column in self.mdata.mod[modality].obs.columns:
                raw_column = str(column)
                values = self._aligned_mod_obs(modality, raw_column)
                label = self._add_obs(
                    frame,
                    occupied,
                    f"{_modality_label(modality)} · {raw_column}",
                    values,
                    _ObsRef("modality", raw_column, modality),
                )
                modality_columns[(modality, raw_column)] = (label, values)

        global_obs = self.mdata.obs
        for column in global_obs.columns:
            raw_column = str(column)
            values = obs_col(global_obs, raw_column).copy()
            values.index = self.cell_ids

            prefix, separator, local_column = raw_column.partition(":")
            if separator and prefix in self.modalities:
                existing = modality_columns.get((prefix, local_column))
                if existing is not None and _same_values(existing[1], values):
                    continue
                if existing is None:
                    self._add_obs(
                        frame,
                        occupied,
                        f"{_modality_label(prefix)} · {local_column}",
                        values,
                        _ObsRef("global", raw_column, prefix),
                    )
                    continue

            global_label = raw_column if not separator else f"MuData · {raw_column}"
            self._add_obs(
                frame,
                occupied,
                global_label,
                values,
                _ObsRef("global", raw_column),
            )

        return frame

    def _build_embedding_catalog(self) -> dict[str, np.ndarray]:
        obsm: dict[str, np.ndarray] = {}
        occupied: set[str] = set()
        for modality in self.modalities:
            for raw_key in self.mdata.mod[modality].obsm.keys():
                stored_values = self.mdata.mod[modality].obsm[raw_key]
                if not is_embedding_obsm(raw_key, stored_values):
                    continue
                display_key = str(raw_key).removeprefix("X_").upper()
                label = _unique_label(
                    f"{_modality_label(modality)} · {display_key}",
                    occupied,
                    "embedding",
                )
                occupied.add(label)
                values = embedding_to_numpy(stored_values)
                obsm[label] = values[self._row_indexers[modality]]
                self._embedding_refs[label] = (modality, str(raw_key))
        if not obsm:
            raise ValueError("No modality contains an embedding in obsm.")
        return obsm

    def is_feature(self, value: str | None) -> bool:
        return bool(value) and value in self._feature_refs

    def feature_modality(self, value: str | None) -> str | None:
        ref = self._feature_refs.get(value or "")
        return ref.modality if ref is not None else None

    def first_feature(self, modality: str) -> str | None:
        return next(
            (
                label
                for label, ref in self._feature_refs.items()
                if ref.modality == modality
            ),
            None,
        )

    def search_features(
        self,
        modality: str,
        query: str | None,
        *,
        limit: int = 20,
    ) -> list[str]:
        """Return a bounded feature search scoped to one modality."""
        if modality not in self.modalities:
            return []
        needle = str(query or "").strip().lower()
        matches: list[str] = []
        for raw_key in self.mdata.mod[modality].var_names.astype(str):
            if needle and needle not in raw_key.lower():
                continue
            label = self._raw_feature_lookup[(modality.lower(), raw_key.lower())]
            matches.append(label)
            if len(matches) >= limit:
                break
        return matches

    def search_all_features(
        self,
        query: str | None,
        *,
        limit: int = 20,
    ) -> list[str]:
        """Return a bounded search across every modality."""
        needle = str(query or "").strip().lower()
        if not needle:
            return self.feature_preview(
                per_modality=max(1, limit // len(self.modalities))
            )[:limit]
        matches: list[str] = []
        for label in self._feature_refs:
            if needle in label.lower():
                matches.append(label)
                if len(matches) >= limit:
                    break
        return matches

    def is_obs(self, value: str | None) -> bool:
        return bool(value) and value in self._obs_refs

    def feature_preview(self, per_modality: int = 10) -> list[str]:
        selected: list[str] = []
        counts = {modality: 0 for modality in self.modalities}
        remaining = set(self.modalities)
        for label, ref in self._feature_refs.items():
            if counts[ref.modality] >= per_modality:
                continue
            selected.append(label)
            counts[ref.modality] += 1
            if counts[ref.modality] >= per_modality:
                remaining.discard(ref.modality)
                if not remaining:
                    break
        return selected

    def default_features(
        self,
        modality_markers: dict[str, Iterable[str] | None] | None = None,
        *,
        per_modality: int = 3,
    ) -> list[str]:
        modality_markers = modality_markers or {}
        selected: list[str] = []
        for modality in self.modalities:
            raw_markers = list(modality_markers.get(modality) or [])
            candidates = [raw_markers]
            candidates.append(
                self.mdata.mod[modality].var_names[:per_modality].astype(str).tolist()
            )
            modality_count = 0
            for candidate_group in candidates:
                for raw_key in candidate_group:
                    label = self._raw_feature_lookup.get(
                        (modality.lower(), str(raw_key).lower())
                    )
                    if label and label not in selected:
                        selected.append(label)
                        modality_count += 1
                    if modality_count >= per_modality:
                        break
                if modality_count >= per_modality:
                    break
        return selected

    def resolve_text_feature(self, value: str) -> str | None:
        text = str(value).strip()
        if text in self._feature_refs:
            return text
        modality, separator, raw_key = text.partition("::")
        if separator:
            return self._raw_feature_lookup.get(
                (modality.strip().lower(), raw_key.strip().lower())
            )
        matches = [
            label
            for (modality_key, raw_lower), label in self._raw_feature_lookup.items()
            if raw_lower == text.lower()
        ]
        return matches[0] if len(set(matches)) == 1 else None

    def score_features(self, features: Iterable[str]) -> np.ndarray:
        """Return one raw feature or a Scanpy set score aligned to global cells."""
        feature_names = tuple(
            sorted(
                name
                for name in dict.fromkeys(features)
                if name in self._feature_refs
            )
        )
        if not feature_names:
            raise ValueError("Select at least one valid feature in each set.")

        modalities = {self._feature_refs[name].modality for name in feature_names}
        if len(modalities) != 1:
            raise ValueError("Each feature set must come from one matrix.")

        cached = self._feature_score_cache.get(feature_names)
        if cached is not None:
            self._feature_score_cache.move_to_end(feature_names)
            return cached

        modality = next(iter(modalities))
        mod_adata = self.mdata.mod[modality]
        raw_keys = [self._feature_refs[name].key for name in feature_names]
        if mod_adata.X is None:
            raise ValueError(f"Matrix '{modality}' does not contain feature values.")

        if len(feature_names) == 1:
            prewarm_gene_cache(mod_adata, raw_keys)
            raw_values = extract_gene_expression(mod_adata, raw_keys[0])
            values = np.asarray(
                raw_values[self._row_indexers[modality]],
                dtype=np.float64,
            )
            values.setflags(write=False)
            self._cache_feature_score(feature_names, values)
            return values

        # score_genes only writes its result to obs. A lightweight wrapper lets it
        # operate on the source matrix without mutating shared MuData metadata or
        # copying the complete expression matrix.
        matrix = mod_adata.X
        if hasattr(matrix, "to_memory"):
            matrix = matrix.to_memory()
        elif not isinstance(matrix, np.ndarray) and not sparse.issparse(matrix):
            matrix = np.asarray(matrix)

        score_adata = ad.AnnData(
            X=matrix,
            obs=pd.DataFrame(index=mod_adata.obs_names.copy()),
            var=pd.DataFrame(index=mod_adata.var_names.copy()),
        )
        score_name = "_guanaco_feature_set_score"
        import scanpy as sc

        try:
            sc.tl.score_genes(
                score_adata,
                raw_keys,
                score_name=score_name,
                random_state=0,
                use_raw=False,
            )
        except Exception as exc:
            raise ValueError(
                f"Could not calculate a feature-set score for matrix '{modality}'."
            ) from exc
        values = score_adata.obs[score_name].to_numpy(dtype=np.float64, copy=True)
        values = values[self._row_indexers[modality]]
        values.setflags(write=False)

        self._cache_feature_score(feature_names, values)
        return values

    def _cache_feature_score(
        self,
        feature_names: tuple[str, ...],
        values: np.ndarray,
    ) -> None:
        self._feature_score_cache[feature_names] = values
        self._feature_score_cache.move_to_end(feature_names)
        while len(self._feature_score_cache) > self._feature_score_cache_size:
            self._feature_score_cache.popitem(last=False)

    def materialize(
        self,
        features: Iterable[str] | None = None,
        *,
        embedding: str | None = None,
    ) -> ad.AnnData:
        """Return a small AnnData containing only requested feature columns."""
        feature_names = tuple(
            name for name in dict.fromkeys(features or ()) if name in self._feature_refs
        )
        cache_key = (*feature_names, "__embedding__", embedding or "")
        cached = self._materialized_cache.get(cache_key)
        if cached is None:
            by_modality: dict[str, list[tuple[str, str]]] = {
                modality: [] for modality in self.modalities
            }
            for display_name in feature_names:
                ref = self._feature_refs[display_name]
                by_modality[ref.modality].append((display_name, ref.key))

            columns: dict[str, np.ndarray] = {}
            for modality, requested in by_modality.items():
                if not requested:
                    continue
                mod_adata = self.mdata.mod[modality]
                raw_keys = [raw_key for _, raw_key in requested]
                prewarm_gene_cache(mod_adata, raw_keys)
                indexer = self._row_indexers[modality]
                for display_name, raw_key in requested:
                    values = extract_gene_expression(mod_adata, raw_key)
                    columns[display_name] = np.asarray(
                        values[indexer], dtype=np.float32
                    )

            if feature_names:
                matrix = np.column_stack([columns[name] for name in feature_names])
            else:
                matrix = np.empty((len(self.cell_ids), 0), dtype=np.float32)

            selected_obsm = (
                {embedding: np.asarray(self.base_adata.obsm[embedding])}
                if embedding in self._embedding_refs
                else {}
            )
            cached = ad.AnnData(
                X=matrix,
                obs=self.base_adata.obs.copy(),
                obsm=selected_obsm,
            )
            cached.obs_names = self.cell_ids.copy()
            cached.var_names = pd.Index(feature_names)
            if embedding and embedding in self._embedding_refs:
                modality, raw_embedding = self._embedding_refs[embedding]
                source = self.mdata.mod[modality]
                if raw_embedding == "spatial" and "spatial" in source.uns:
                    cached.uns["spatial"] = source.uns["spatial"]
            self._materialized_cache[cache_key] = cached
            self._materialized_cache.move_to_end(cache_key)
            while len(self._materialized_cache) > self._materialized_cache_size:
                self._materialized_cache.popitem(last=False)
        else:
            self._materialized_cache.move_to_end(cache_key)

        return cached


def try_build_multiomics_source(mdata) -> tuple[MultiOmicsSource | None, str | None]:
    """Build a joint view, returning a user-readable reason when unavailable."""
    try:
        return MultiOmicsSource(mdata), None
    except ValueError as exc:
        return None, str(exc)
