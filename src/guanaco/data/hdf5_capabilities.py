"""Fast visualization-capability inspection for local H5AD and H5MU files.

The wizard only needs structural facts: dimensions, column encodings, a small
feature-name sample, and result-object schemas.  Reading those facts directly
from HDF5 avoids constructing AnnData/MuData objects and, critically, avoids
importing the Muon/Scanpy runtime just to populate a list of checkboxes.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import h5py


_PEAK_RE = re.compile(r"^\s*([^:\s]+)\s*:\s*([0-9,]+)\s*-\s*([0-9,]+)\s*$")
_NHOOD_SUFFIX = "_nhood_enrichment"
_CO_OCCURRENCE_SUFFIX = "_co_occurrence"
_IDENTITY_ALIASES = {
    "source": ("source", "sender", "source_group", "source_cell", "celltype_a"),
    "target": ("target", "receiver", "target_group", "target_cell", "celltype_b"),
    "ligand": ("ligand_complex", "ligand", "ligand_name", "gene_a", "partner_a"),
    "receptor": (
        "receptor_complex",
        "receptor",
        "receptor_name",
        "gene_b",
        "partner_b",
    ),
}


@dataclass(frozen=True)
class HDF5CapabilitySummary:
    """Capability facts needed by the configuration wizard."""

    modalities: tuple[str, ...] | None
    available_plots: frozenset[str]


def _text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _attr_texts(group: h5py.Group, key: str) -> tuple[str, ...]:
    values = group.attrs.get(key, ())
    if isinstance(values, (str, bytes)):
        values = (values,)
    return tuple(_text(value) for value in values)


def _encoding(obj) -> str:
    return _text(obj.attrs.get("encoding-type", "")).lower()


def _frame_columns(frame: h5py.Group | None) -> tuple[str, ...]:
    if not isinstance(frame, h5py.Group):
        return ()
    configured = _attr_texts(frame, "column-order")
    if configured:
        return configured
    index_name = _text(frame.attrs.get("_index", "_index"))
    return tuple(name for name in frame.keys() if name not in {"_index", index_name})


def _frame_length(frame: h5py.Group | None) -> int:
    if not isinstance(frame, h5py.Group):
        return 0
    index_name = _text(frame.attrs.get("_index", "_index"))
    index = frame.get(index_name) or frame.get("_index")
    if isinstance(index, h5py.Dataset) and index.shape:
        return int(index.shape[0])
    for column in frame.values():
        if isinstance(column, h5py.Dataset) and column.shape:
            return int(column.shape[0])
        if isinstance(column, h5py.Group):
            codes = column.get("codes")
            if isinstance(codes, h5py.Dataset) and codes.shape:
                return int(codes.shape[0])
    return 0


def _sample_values(dataset: h5py.Dataset, *, limit: int = 512):
    if dataset.ndim != 1 or not dataset.shape or dataset.shape[0] == 0:
        return ()
    size = int(dataset.shape[0])
    if size <= limit:
        return dataset[:]

    block = max(1, limit // 3)
    middle = max(block, (size - block) // 2)
    slices = (
        dataset[:block],
        dataset[middle : middle + block],
        dataset[size - block :],
    )
    return tuple(value for values in slices for value in values)


def _categorical_cardinality(column) -> int | None:
    if not isinstance(column, h5py.Group) or _encoding(column) != "categorical":
        return None
    categories = column.get("categories")
    if isinstance(categories, h5py.Dataset) and categories.shape:
        return int(categories.shape[0])
    return None


def _is_discrete_column(column, *, max_unique: int = 50) -> bool:
    cardinality = _categorical_cardinality(column)
    if cardinality is not None:
        return cardinality < max_unique
    if isinstance(column, h5py.Group):
        return "boolean" in _encoding(column)
    if not isinstance(column, h5py.Dataset) or column.ndim != 1:
        return False
    if column.dtype.kind == "b":
        return True
    if column.dtype.kind not in {"O", "S", "U"} and "string" not in _encoding(column):
        return False
    values = {_text(value) for value in _sample_values(column)}
    return 0 < len(values) < max_unique


def _has_many_numeric_values(column, *, min_unique: int = 51) -> bool:
    if not isinstance(column, h5py.Dataset) or column.ndim != 1:
        return False
    if column.dtype.kind not in {"i", "u", "f"}:
        return False
    values = set()
    for value in _sample_values(column):
        try:
            if value != value:  # NaN
                continue
        except Exception:
            continue
        values.add(value.item() if hasattr(value, "item") else value)
        if len(values) >= min_unique:
            return True
    return False


def _discrete_columns(obs: h5py.Group | None) -> set[str]:
    if not isinstance(obs, h5py.Group):
        return set()
    return {
        name
        for name in _frame_columns(obs)
        if name in obs and _is_discrete_column(obs[name])
    }


def _has_continuous_observation(obs: h5py.Group | None) -> bool:
    if not isinstance(obs, h5py.Group):
        return False
    return any(
        name in obs and _has_many_numeric_values(obs[name])
        for name in _frame_columns(obs)
    )


def _has_peak_features(var: h5py.Group | None, n_vars: int) -> bool:
    if not isinstance(var, h5py.Group) or n_vars == 0:
        return False
    if {"chrom", "start", "end"}.issubset(set(_frame_columns(var))):
        return True
    index_name = _text(var.attrs.get("_index", "_index"))
    index = var.get(index_name) or var.get("_index")
    if not isinstance(index, h5py.Dataset):
        return False
    sample_size = min(2000, n_vars)
    hits = sum(_PEAK_RE.match(_text(name)) is not None for name in index[:sample_size])
    return hits >= 5


def _normalized_column_names(frame: h5py.Group) -> set[str]:
    return {
        name.strip().lower().replace(" ", "_").replace(".", "_").replace("-", "_")
        for name in _frame_columns(frame)
    }


def _has_ligand_receptor_result(uns: h5py.Group | None) -> bool:
    if not isinstance(uns, h5py.Group):
        return False
    for value in uns.values():
        if not isinstance(value, h5py.Group):
            continue
        if _encoding(value) == "dataframe":
            columns = _normalized_column_names(value)
            if all(
                any(alias in columns for alias in aliases)
                for aliases in _IDENTITY_ALIASES.values()
            ):
                return True
        for table_name in ("means", "significant_means"):
            table = value.get(table_name)
            if isinstance(table, h5py.Group) and _encoding(table) == "dataframe":
                return True
    return False


def _shape(group: h5py.Group, key: str) -> tuple[int, ...] | None:
    value = group.get(key)
    return tuple(value.shape) if isinstance(value, h5py.Dataset) else None


def _has_spatial_relationships(
    obs: h5py.Group | None,
    uns: h5py.Group | None,
    discrete_columns: set[str],
) -> bool:
    if not isinstance(obs, h5py.Group) or not isinstance(uns, h5py.Group):
        return False
    for uns_key in uns.keys():
        if not str(uns_key).endswith(_NHOOD_SUFFIX):
            continue
        cluster_key = str(uns_key)[: -len(_NHOOD_SUFFIX)]
        if cluster_key not in discrete_columns or cluster_key not in obs:
            continue
        n_clusters = _categorical_cardinality(obs[cluster_key])
        nhood = uns.get(str(uns_key))
        occurrence = uns.get(f"{cluster_key}{_CO_OCCURRENCE_SUFFIX}")
        if (
            not n_clusters
            or not isinstance(nhood, h5py.Group)
            or not isinstance(occurrence, h5py.Group)
        ):
            continue
        zscore = _shape(nhood, "zscore")
        count = _shape(nhood, "count")
        occ = _shape(occurrence, "occ")
        interval = _shape(occurrence, "interval")
        if (
            zscore == (n_clusters, n_clusters)
            and count == zscore
            and occ is not None
            and len(occ) == 3
            and occ[:2] == zscore
            and interval is not None
            and len(interval) == 1
            and occ[2] == max(0, interval[0] - 1)
        ):
            return True
    return False


def _scan_anndata_group(
    group: h5py.Group,
    modality_name: str | None,
    *,
    shared_obs: h5py.Group | None = None,
) -> set[str]:
    obs = group.get("obs")
    var = group.get("var")
    uns = group.get("uns")
    n_obs = _frame_length(obs)
    n_vars = _frame_length(var)
    has_features = n_obs > 0 and n_vars > 0
    local_discrete = _discrete_columns(obs)
    shared_discrete = _discrete_columns(shared_obs)
    discrete = local_discrete | shared_discrete
    has_discrete = bool(discrete)
    has_peaks = _has_peak_features(var, n_vars)

    plots: set[str] = set()
    if has_features and has_discrete:
        plots.update(("dotplot", "heatmap", "violin", "ridge"))
    if has_features and (
        _has_continuous_observation(obs) or _has_continuous_observation(shared_obs)
    ):
        plots.add("pseudotime")
    if has_discrete:
        plots.update(("split-violin", "stacked-bar"))
    if _has_ligand_receptor_result(uns):
        plots.add("ligand-receptor")
    if isinstance(uns, h5py.Group):
        paga = uns.get("paga")
        if isinstance(paga, h5py.Group) and "connectivities" in paga:
            plots.add("paga")
        if "volcano" in uns or "rank_genes_groups" in uns:
            plots.add("volcano")
    if _has_spatial_relationships(obs, uns, local_discrete) or (
        shared_obs is not obs
        and _has_spatial_relationships(shared_obs, uns, shared_discrete)
    ):
        plots.add("spatial-relationships")
    if has_peaks:
        plots.add("peak-browser")
    return plots


def _modality_names(mod_group: h5py.Group) -> tuple[str, ...]:
    configured = _attr_texts(mod_group, "mod-order")
    names = [name for name in configured if name in mod_group]
    names.extend(name for name in mod_group.keys() if name not in names)
    return tuple(names)


def inspect_hdf5_capabilities(path: str | Path) -> HDF5CapabilitySummary:
    """Inspect one local H5AD/H5MU using HDF5 metadata only."""
    source = Path(path)
    with h5py.File(source, "r") as handle:
        if source.suffix.lower() == ".h5mu":
            mod_group = handle.get("mod")
            if not isinstance(mod_group, h5py.Group):
                raise ValueError("MuData file does not contain a 'mod' group.")
            modalities = _modality_names(mod_group)
            shared_obs = handle.get("obs")
            plots: set[str] = set()
            for name in modalities:
                plots.update(
                    _scan_anndata_group(
                        mod_group[name],
                        name,
                        shared_obs=(
                            shared_obs if isinstance(shared_obs, h5py.Group) else None
                        ),
                    )
                )
            if len(modalities) >= 2:
                plots.add("multiomics-composition")
                plots.add("cross-modal-concordance")
            return HDF5CapabilitySummary(modalities, frozenset(plots))

        plots = _scan_anndata_group(handle, "rna")
        return HDF5CapabilitySummary(None, frozenset(plots))
