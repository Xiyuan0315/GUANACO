import time
import gc
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.sparse import issparse


class GeneExpressionCache:
    """O(1) LRU + TTL cache for 1D gene vectors."""

    def __init__(self, max_size=64, max_age_seconds=1800):
        self.max_size = int(max_size)
        self.max_age_seconds = float(max_age_seconds)
        self._data = OrderedDict()  # key -> (timestamp, value)

    def _adata_id(self, adata):
        # backed: filename is stable; else use id
        if getattr(adata, "isbacked", False) and getattr(adata, "filename", None):
            return ("backed", adata.filename)
        return ("mem", id(adata))

    def _make_key(self, adata, gene, layer=None, use_raw=False):
        # include shape to distinguish filtered/unfiltered views
        return (self._adata_id(adata), adata.shape, gene, layer, bool(use_raw))

    def get_or_compute(self, adata, gene, layer, use_raw, compute_fn):
        key = self._make_key(adata, gene, layer, use_raw)
        now = time.time()

        if key in self._data:
            ts, val = self._data[key]
            if now - ts < self.max_age_seconds:
                self._data.move_to_end(key)  # mark as recently used
                return val
            # expired
            del self._data[key]

        val = compute_fn()
        self._data[key] = (now, val)
        self._data.move_to_end(key)

        # evict LRU
        while len(self._data) > self.max_size:
            self._data.popitem(last=False)

        # periodic garbage collection
        if len(self._data) % 20 == 0:
            gc.collect()

        return val

    def clear(self):
        self._data.clear()
        gc.collect()

    def info(self):
        return {
            "size": len(self._data),
            "max_size": self.max_size,
            "keys": list(self._data.keys()),
        }


# --------------------------
# Internal helper functions
# --------------------------

def _get_data_source(adata, use_raw: bool):
    """Return adata.raw (if requested & present) else adata."""
    return adata.raw if use_raw and getattr(adata, "raw", None) is not None else adata


def _get_matrix(data_source, layer):
    """Return matrix for layer or X."""
    if layer is None:
        return data_source.X
    if layer not in data_source.layers:
        raise KeyError(f"Layer '{layer}' not found.")
    return data_source.layers[layer]


def _gene_indexer(var_names: pd.Index, genes):
    """
    Vectorized gene lookup: returns integer positions; -1 for missing.
    genes can be a str or list-like.
    """
    if isinstance(genes, str):
        genes = [genes]
    idx = var_names.get_indexer(genes)
    return np.asarray(idx, dtype=np.int64), list(genes)


def _raise_if_missing(idx, genes):
    if np.any(idx == -1):
        missing = np.asarray(genes, dtype=object)[idx == -1]
        raise KeyError(f"Genes not found: {missing.tolist()}")


def _take_columns_backed_aware(X, col_idx):
    """
    Read X[:, col_idx] efficiently.

    For many backed / on-disk arrays, reading columns in increasing order
    can be faster. We read sorted columns then restore original order.
    """
    col_idx = np.asarray(col_idx, dtype=np.int64)
    if col_idx.size == 0:
        return X[:, :0]

    order = np.argsort(col_idx)
    sorted_idx = col_idx[order]

    sub = X[:, sorted_idx]

    inv = np.empty_like(order)
    inv[order] = np.arange(order.size)
    return sub[:, inv]


# --------------------------
# Public API
# --------------------------

_gene_cache = GeneExpressionCache(max_size=64, max_age_seconds=1800)


def extract_gene_expression(
    adata,
    gene,
    layer=None,
    use_raw=False,
    use_cache=True,
    dtype=None,
):
    """
    Fast 1D extraction for a single gene.

    Parameters
    ----------
    adata : AnnData
    gene : str
    layer : str | None
    use_raw : bool
    use_cache : bool
    dtype : numpy dtype | None
        If provided (e.g., np.float32), cast output without extra copy when possible.

    Returns
    -------
    np.ndarray (1D)
    """
    def _compute():
        data_source = _get_data_source(adata, use_raw)
        X = _get_matrix(data_source, layer)

        idx, genes = _gene_indexer(data_source.var_names, gene)
        _raise_if_missing(idx, genes)
        j = int(idx[0])

        col = X[:, j]
        if issparse(col):
            out = col.toarray().ravel()
        else:
            out = np.asarray(col).ravel()

        if dtype is not None:
            out = out.astype(dtype, copy=False)
        return out

    if use_cache:
        return _gene_cache.get_or_compute(adata, gene, layer, use_raw, _compute)
    return _compute()


def extract_multiple_genes(
    adata,
    genes,
    layer=None,
    use_raw=False,
    cell_indices=None,
    allow_missing=False,
    to_dense=True,
    dtype=np.float32,
):
    """
    Extract multiple genes efficiently.

    Parameters
    ----------
    adata : AnnData
    genes : list[str] | str
    layer : str | None
    use_raw : bool
    cell_indices : array-like | None
        Subset of cells to include (indices or boolean mask).
    allow_missing : bool
        If False, raise on any missing gene.
        If True, drop missing genes.
    to_dense : bool
        If True, convert to dense (DataFrame requires dense).
    dtype : numpy dtype | None
        Cast dense arrays to dtype (default float32).

    Returns
    -------
    pd.DataFrame
        rows = cells, columns = genes
    """
    data_source = _get_data_source(adata, use_raw)
    X = _get_matrix(data_source, layer)

    idx, gene_list = _gene_indexer(data_source.var_names, genes)

    if not allow_missing:
        _raise_if_missing(idx, gene_list)
        keep_mask = np.ones_like(idx, dtype=bool)
    else:
        keep_mask = idx != -1
        idx = idx[keep_mask]
        gene_list = list(np.asarray(gene_list, dtype=object)[keep_mask])
        if idx.size == 0:
            raise KeyError("No valid genes found.")

    # Prefer row slicing first (helps backed mode)
    if cell_indices is not None:
        cell_indices = np.asarray(cell_indices)
        X_sub = X[cell_indices, :]
    else:
        X_sub = X

    # Backed-aware column slicing
    gene_data = _take_columns_backed_aware(X_sub, idx)

    # Convert to dense if requested
    if to_dense:
        if issparse(gene_data):
            gene_data = gene_data.toarray()
        elif hasattr(gene_data, "compute"):  # dask arrays
            gene_data = gene_data.compute()
        else:
            gene_data = np.asarray(gene_data)

        if dtype is not None:
            gene_data = gene_data.astype(dtype, copy=False)

    obs_index = adata.obs_names if cell_indices is None else adata.obs_names[cell_indices]
    return pd.DataFrame(gene_data, columns=gene_list, index=obs_index)


def apply_transformation(expr, method="log1p", copy=True, clip_percentile=99):
    """
    Apply common transformations for visualization.

    method:
      - "log1p" / "log"
      - "zscore" / "z_score" (with clipping to percentiles)
    """
    if copy:
        expr = expr.copy()

    if method in ("zscore", "z_score"):
        if isinstance(expr, pd.DataFrame):
            mean = expr.mean(axis=0)
            std = expr.std(axis=0).replace(0, np.nan)
            expr = (expr - mean) / std
            expr = expr.fillna(0.0)

            arr = expr.to_numpy()
            upper = np.nanpercentile(arr, clip_percentile)
            lower = np.nanpercentile(arr, 100 - clip_percentile)
            expr = expr.clip(lower=lower, upper=upper)
        else:
            expr = np.asarray(expr)
            mean = expr.mean()
            std = expr.std()
            expr = (expr - mean) / std if std > 0 else (expr - mean)

            upper = np.percentile(expr, clip_percentile)
            lower = np.percentile(expr, 100 - clip_percentile)
            expr = np.clip(expr, lower, upper)

    elif method in ("log1p", "log"):
        expr = np.log1p(expr)

    return expr


def clear_gene_cache():
    """Clear the global gene expression cache."""
    _gene_cache.clear()


def get_cache_info():
    """Get information about the current cache state."""
    return _gene_cache.info()


def bin_cells_for_heatmap(df, gene_columns, groupby, n_bins, continuous_key=None):
    """
    Unified binning function for heatmap visualization to reduce memory usage.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with cells as rows, genes as columns
    gene_columns : list[str]
        Gene column names
    groupby : str
        Column name for cell type/group
    n_bins : int
        Number of bins total (distributed across groups)
    continuous_key : str | None
        If provided, sort by this continuous variable first (continuous binning)

    Returns
    -------
    pd.DataFrame
        Binned dataframe with reduced number of rows
    """
    if len(df) == 0:
        cols = [groupby, *gene_columns] if continuous_key is None else [groupby, continuous_key, *gene_columns]
        return df[cols].copy()

    n_bins = max(1, min(int(n_bins), len(df)))

    # Use float32 view for gene matrix to reduce memory
    gene_data = df[gene_columns].to_numpy(dtype=np.float32, copy=False)

    # --- Continuous ordering binning (single global sequence) ---
    if continuous_key is not None:
        df_sorted = df.sort_values(continuous_key)
        n_cells = len(df_sorted)

        edges = np.unique(np.linspace(0, n_cells, num=n_bins + 1, dtype=int))
        starts = edges[:-1]
        ends = edges[1:]
        valid = ends > starts
        starts = starts[valid]
        ends = ends[valid]

        sorted_gene_data = df_sorted[gene_columns].to_numpy(dtype=np.float32, copy=False)
        sorted_group_data = df_sorted[groupby].to_numpy()
        sorted_continuous = df_sorted[continuous_key].to_numpy(dtype=np.float32, copy=False)

        gene_sums = np.add.reduceat(sorted_gene_data, starts, axis=0)
        lengths = (ends - starts).reshape(-1, 1).astype(np.float32)
        gene_means = gene_sums / lengths

        continuous_sums = np.add.reduceat(sorted_continuous, starts)
        continuous_means = continuous_sums / (ends - starts)

        # Mode group label per bin (fast enough at bin-count scale)
        group_modes = []
        for s, e in zip(starts, ends):
            values, counts = np.unique(sorted_group_data[s:e], return_counts=True)
            group_modes.append(values[np.argmax(counts)])

        result_df = pd.DataFrame(gene_means, columns=gene_columns)
        result_df.insert(0, continuous_key, continuous_means.astype(np.float32))
        result_df.insert(0, groupby, group_modes)
        return result_df

    # --- Group-wise binning ---
    group_codes, unique_groups = pd.factorize(df[groupby], sort=False)

    binned_blocks = []
    binned_group_labels = []
    total_cells = len(df)

    for code, group in enumerate(unique_groups):
        group_idx = np.flatnonzero(group_codes == code)
        n_cells_in_group = group_idx.size
        if n_cells_in_group == 0:
            continue

        group_gene_data = gene_data[group_idx]

        # For tiny groups, keep as-is to avoid over-averaging
        if n_cells_in_group <= 10:
            binned_blocks.append(group_gene_data)
            binned_group_labels.extend([group] * n_cells_in_group)
            continue

        # Allocate bins proportional to group size
        group_bins = max(1, int(n_bins * n_cells_in_group / total_cells))
        group_bins = min(group_bins, n_cells_in_group)

        edges = np.unique(np.linspace(0, n_cells_in_group, num=group_bins + 1, dtype=int))
        starts = edges[:-1]
        ends = edges[1:]
        valid = ends > starts
        starts = starts[valid]
        ends = ends[valid]

        group_sums = np.add.reduceat(group_gene_data, starts, axis=0)
        lengths = (ends - starts).reshape(-1, 1).astype(np.float32)
        group_means = group_sums / lengths

        binned_blocks.append(group_means)
        binned_group_labels.extend([group] * group_means.shape[0])

    if not binned_blocks:
        return pd.DataFrame(columns=[groupby, *gene_columns])

    binned_gene_data = np.vstack(binned_blocks)
    result_df = pd.DataFrame(binned_gene_data, columns=gene_columns)
    result_df.insert(0, groupby, binned_group_labels)
    return result_df
