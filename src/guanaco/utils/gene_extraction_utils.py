import time
import gc
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.sparse import issparse


class GeneExpressionCache:
    """O(1) LRU + TTL cache for 1D gene vectors, bounded by item count *and* bytes.

    Each entry is an ``n_obs``-length vector, so its size scales with the dataset.
    The byte budgets make the cache self-scaling: large datasets keep fewer
    vectors resident, which matters in backed mode where all cells are served.
    """

    def __init__(
        self,
        max_size=64,
        max_age_seconds=1800,
        max_bytes=256 * 1024 * 1024,
        max_pinned_bytes=128 * 1024 * 1024,
    ):
        self.max_size = int(max_size)
        self.max_age_seconds = float(max_age_seconds)
        self.max_bytes = int(max_bytes)
        self.max_pinned_bytes = int(max_pinned_bytes)
        self._data = OrderedDict()  # key -> (timestamp, value)
        self._pinned = {}  # key -> value; never expires (still bounded by bytes)
        self._data_bytes = 0
        self._pinned_bytes = 0

    @staticmethod
    def _nbytes(value):
        """Best-effort byte size of a cached vector (numpy array)."""
        return int(getattr(value, "nbytes", 0) or 0)

    def _adata_id(self, adata):
        # backed: filename is stable; else use id
        if getattr(adata, "isbacked", False) and getattr(adata, "filename", None):
            return ("backed", adata.filename)
        return ("mem", id(adata))

    def _make_key(self, adata, gene, layer=None, use_raw=False):
        # include shape to distinguish filtered/unfiltered views
        return (self._adata_id(adata), adata.shape, gene, layer, bool(use_raw))

    def _pop_data(self, key):
        """Remove an LRU entry and decrement the byte counter."""
        item = self._data.pop(key, None)
        if item is not None:
            self._data_bytes -= self._nbytes(item[1])

    def _evict_data(self):
        """Evict LRU entries until within both the item and byte budgets."""
        while self._data and (
            len(self._data) > self.max_size or self._data_bytes > self.max_bytes
        ):
            old_key, (_, old_val) = self._data.popitem(last=False)
            self._data_bytes -= self._nbytes(old_val)

    def _store_data(self, key, val, now):
        self._pop_data(key)  # avoid double-counting if replacing
        self._data[key] = (now, val)
        self._data_bytes += self._nbytes(val)
        self._data.move_to_end(key)
        self._evict_data()

    def get_or_compute(self, adata, gene, layer, use_raw, compute_fn):
        key = self._make_key(adata, gene, layer, use_raw)

        # Pinned genes are always served from memory and never expire/evict.
        pinned = self._pinned.get(key)
        if pinned is not None:
            return pinned

        now = time.time()

        if key in self._data:
            ts, val = self._data[key]
            if now - ts < self.max_age_seconds:
                self._data.move_to_end(key)  # mark as recently used
                return val
            # expired
            self._pop_data(key)

        val = compute_fn()
        self._store_data(key, val, now)
        return val

    def peek(self, adata, gene, layer=None, use_raw=False):
        """Return a cached vector without computing it; None if absent or expired.

        Lets a batched reader serve cache hits and collect only the misses, so it
        can read those in a single column slice.
        """
        key = self._make_key(adata, gene, layer, use_raw)
        pinned = self._pinned.get(key)
        if pinned is not None:
            return pinned
        if key in self._data:
            ts, val = self._data[key]
            if time.time() - ts < self.max_age_seconds:
                self._data.move_to_end(key)
                return val
            self._pop_data(key)
        return None

    def store(self, adata, gene, layer, use_raw, value):
        """Insert a precomputed vector into the LRU cache (no-op if already pinned)."""
        key = self._make_key(adata, gene, layer, use_raw)
        if key in self._pinned:
            return
        self._store_data(key, value, time.time())

    def pin(self, adata, gene, layer, use_raw, compute_fn):
        """Compute (if needed) and keep a gene vector pinned in memory.

        Pinning is skipped once ``max_pinned_bytes`` is reached; the vector is then
        stored as an ordinary (evictable) LRU entry instead, so prewarming a huge
        dataset can't grow without bound.
        """
        key = self._make_key(adata, gene, layer, use_raw)
        if key in self._pinned:
            return self._pinned[key]
        val = compute_fn()
        nb = self._nbytes(val)
        if self._pinned_bytes + nb > self.max_pinned_bytes:
            # Over the pin budget: keep it, but as an evictable LRU entry.
            self._store_data(key, val, time.time())
            return val
        self._pinned[key] = val
        self._pinned_bytes += nb
        # Drop any LRU copy so we don't keep two references to the same vector.
        self._pop_data(key)
        return val

    def clear(self):
        self._data.clear()
        self._pinned.clear()
        self._data_bytes = 0
        self._pinned_bytes = 0
        gc.collect()

    def info(self):
        return {
            "size": len(self._data),
            "pinned": len(self._pinned),
            "max_size": self.max_size,
            "data_mb": round(self._data_bytes / 1024 / 1024, 1),
            "pinned_mb": round(self._pinned_bytes / 1024 / 1024, 1),
            "max_mb": round(self.max_bytes / 1024 / 1024, 1),
            "keys": list(self._data.keys()),
        }


# --------------------------
# Internal helper functions
# --------------------------

def _densify(block):
    """Return a dense ``numpy`` array from a numpy / scipy-sparse / dask block.

    Cloud-backed (lazy) ``X``/``layers`` are dask arrays whose chunks may be
    scipy-sparse. ``.compute()`` reads only the requested slice off cloud/disk but
    yields a sparse matrix, which still has to be densified. Already-dense or
    in-memory blocks pass straight through.
    """
    if hasattr(block, "compute"):  # dask array -> pull just this slice into memory
        # Use the synchronous scheduler: under a gunicorn sync worker, dask's default
        # threaded scheduler calls fsspec's async HTTP from worker threads while its
        # event loop lives in another thread, which deadlocks (the worker hangs in
        # `compute()` until gunicorn times out and SIGKILLs it -- looks like an OOM
        # but isn't). Running inline in the calling thread reads the same single
        # column with no thread hand-off, so cloud-backed reads work under gunicorn.
        block = block.compute(scheduler="synchronous")
    if issparse(block):
        block = block.toarray()
    return np.asarray(block)


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
    return np.asarray(idx, dtype=np.int64)


def _raise_if_missing(idx, genes):
    if np.any(idx == -1):
        gene_list = [genes] if isinstance(genes, str) else list(genes)
        missing = np.asarray(gene_list, dtype=object)[idx == -1]
        raise KeyError(f"Genes not found: {missing.tolist()}")


def _compute_gene_vector(adata, gene, layer=None, use_raw=False, dtype=None):
    """Read a single gene's 1D expression vector from X/layer (disk or memory)."""
    data_source = _get_data_source(adata, use_raw)
    X = _get_matrix(data_source, layer)

    idx = _gene_indexer(data_source.var_names, gene)
    _raise_if_missing(idx, gene)
    j = int(idx[0])

    # X may be a lazy/cloud-backed dask array: this slice reads only column j.
    out = _densify(X[:, j]).ravel()

    if dtype is not None:
        out = out.astype(dtype, copy=False)
    return out


def _compute_gene_block(adata, genes, layer=None, use_raw=False, dtype=None):
    """Read several genes' columns in ONE slice: returns {gene: 1D vector}.

    ``X[:, idxs]`` is a single column slice -- one ``indptr`` scan for a sparse
    matrix, one read for a cloud-backed dask array -- instead of one slice per
    gene. Genes must already be known to exist in ``var_names`` (callers filter).
    """
    data_source = _get_data_source(adata, use_raw)
    X = _get_matrix(data_source, layer)

    idxs = np.asarray(data_source.var_names.get_indexer(genes), dtype=np.int64)
    block = _densify(X[:, idxs])
    if block.ndim == 1:  # single gene -> keep 2D so the column split below works
        block = block.reshape(-1, 1)

    out = {}
    for col, gene in enumerate(genes):
        vec = np.ascontiguousarray(block[:, col]).ravel()
        if dtype is not None:
            vec = vec.astype(dtype, copy=False)
        out[gene] = vec
    return out


# --------------------------
# Public API
# --------------------------

# Byte ceiling lowered from the 256 MB default: the interactive workflow shows
# < ~20 genes at a time (added/removed one by one), so a smaller resident working
# set is enough and caps memory on very large datasets. The heatmap additionally
# uses a binned-per-gene cache (see plots/heatmap.py) so it no longer depends on
# full-length columns staying resident.
_gene_cache = GeneExpressionCache(
    max_size=64, max_age_seconds=1800, max_bytes=64 * 1024 * 1024
)

# Default cap on how many config genes we pre-load per dataset.
PIN_GENE_LIMIT = 50


def pin_genes(adata, genes, layer=None, use_raw=False, max_genes=PIN_GENE_LIMIT, dtype=np.float32):
    """Pre-load (pin) gene vectors into memory so the first access is instant.

    Intended for backed datasets, where the first read of a gene otherwise hits
    disk. All genes are read in ONE batched column slice rather than one slice per
    gene -- the difference between one and N ``indptr`` scans / cloud reads at
    startup. Pinned vectors are never evicted by the LRU cache and never expire.
    Genes missing from ``adata.var_names`` are skipped. Returns the count pinned.
    """
    if isinstance(genes, str):
        genes = [genes]
    if not genes:
        return 0

    var_names = adata.var_names
    to_pin = []
    for gene in genes:
        if max_genes is not None and len(to_pin) >= max_genes:
            break
        if gene in var_names and gene not in to_pin:
            to_pin.append(gene)
    if not to_pin:
        return 0

    # One slice for all genes; then pin each precomputed vector (the pin's compute_fn
    # just returns it, so no further reads happen).
    vectors = _compute_gene_block(adata, to_pin, layer, use_raw, dtype)
    pinned = 0
    for gene in to_pin:
        vec = vectors.get(gene)
        if vec is None:
            continue
        _gene_cache.pin(adata, gene, layer, use_raw, lambda v=vec: v)
        pinned += 1
    return pinned


def prewarm_gene_cache(adata, genes, layer=None, use_raw=False, dtype=np.float32):
    """Batch-read any not-yet-cached genes in one slice and populate the LRU cache.

    Call this before a per-gene extraction loop: the genes already cached (or
    pinned) are left untouched, the rest are read in a single ``X[:, idxs]`` slice
    and stored, so the following ``extract_gene_expression`` calls all hit the
    cache instead of issuing one column slice each. Missing genes are skipped.
    """
    if isinstance(genes, str):
        genes = [genes]
    if not genes:
        return

    data_source = _get_data_source(adata, use_raw)
    var_names = data_source.var_names

    misses = []
    seen = set()
    for gene in genes:
        if gene in seen or gene not in var_names:
            continue
        seen.add(gene)
        if _gene_cache.peek(adata, gene, layer, use_raw) is None:
            misses.append(gene)
    if not misses:
        return

    vectors = _compute_gene_block(adata, misses, layer, use_raw, dtype)
    for gene, vec in vectors.items():
        _gene_cache.store(adata, gene, layer, use_raw, vec)


def extract_gene_expression(
    adata,
    gene,
    layer=None,
    use_raw=False,
    use_cache=True,
    dtype=np.float32,
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
        Cast output to this dtype without an extra copy when possible. Defaults to
        ``np.float32`` so cached vectors are half the size of float64 (display
        plots don't need float64 precision); pass ``None`` to keep the source
        dtype.

    Returns
    -------
    np.ndarray (1D)
    """
    def _compute():
        return _compute_gene_vector(adata, gene, layer, use_raw, dtype)

    if use_cache:
        return _gene_cache.get_or_compute(adata, gene, layer, use_raw, _compute)
    return _compute()


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
        # Number of original cells aggregated into each bin (for true-count widths).
        result_df['__bin_size__'] = (ends - starts).astype(np.int64)
        return result_df

    # --- Group-wise binning ---
    gene_data = df[gene_columns].to_numpy(dtype=np.float32, copy=False)
    group_codes, unique_groups = pd.factorize(df[groupby], sort=False)

    binned_blocks = []
    binned_group_labels = []
    binned_sizes = []
    total_cells = len(df)

    for code, group in enumerate(unique_groups):
        group_idx = np.flatnonzero(group_codes == code)
        n_cells_in_group = group_idx.size
        if n_cells_in_group == 0:
            continue

        group_gene_data = gene_data[group_idx]

        # For tiny groups, keep as-is to avoid over-averaging (each row = 1 cell)
        if n_cells_in_group <= 10:
            binned_blocks.append(group_gene_data)
            binned_group_labels.extend([group] * n_cells_in_group)
            binned_sizes.extend([1] * n_cells_in_group)
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
        binned_sizes.extend((ends - starts).astype(np.int64).tolist())

    if not binned_blocks:
        return pd.DataFrame(columns=[groupby, *gene_columns, '__bin_size__'])

    binned_gene_data = np.vstack(binned_blocks)
    result_df = pd.DataFrame(binned_gene_data, columns=gene_columns)
    result_df.insert(0, groupby, binned_group_labels)
    # Number of original cells aggregated into each bin (for true-count widths).
    result_df['__bin_size__'] = np.asarray(binned_sizes, dtype=np.int64)
    return result_df
