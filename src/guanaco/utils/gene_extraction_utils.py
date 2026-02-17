import numpy as np
import pandas as pd
from scipy.sparse import issparse
import time
import gc


class GeneExpressionCache:
    """Cache for gene expression data to avoid repeated extractions."""
    
    def __init__(self, max_size=64, max_age_seconds=1800):  # Reduced to 64 items, 30-min expiration
        self.cache = {}
        self.max_size = max_size
        self.max_age_seconds = max_age_seconds
        self.timestamps = {}
        self.access_order = []
        
    def _make_key(self, adata, gene, layer=None, use_raw=False):
        """Create a unique key for caching."""
        # For backed mode, use the filename as a stable identifier
        if hasattr(adata, 'isbacked') and adata.isbacked and hasattr(adata, 'filename'):
            adata_id = adata.filename
        else:
            adata_id = id(adata)
        
        # Include data shape to differentiate between filtered and unfiltered data
        data_shape = adata.shape
        
        return (adata_id, gene, layer, use_raw, data_shape)
    
    def get(self, adata, gene, layer=None, use_raw=False):
        """Get gene expression from cache or compute if not cached."""
        key = self._make_key(adata, gene, layer, use_raw)
        current_time = time.time()
        
        # Check if cached and not expired
        if key in self.cache:
            if current_time - self.timestamps.get(key, 0) < self.max_age_seconds:
                # Move to end of access order
                if key in self.access_order:
                    self.access_order.remove(key)
                self.access_order.append(key)
                return self.cache[key]  # Return reference, not copy (saves memory)
            else:
                # Expired - remove from cache
                del self.cache[key]
                del self.timestamps[key]
                if key in self.access_order:
                    self.access_order.remove(key)
        
        # Compute if not in cache - use extract without cache to avoid recursion
        expr = extract_gene_expression(adata, gene, layer, use_raw, use_cache=False)
        
        # Add to cache with size limit
        if len(self.cache) >= self.max_size:
            # Remove least recently used
            if self.access_order:
                lru_key = self.access_order.pop(0)
                del self.cache[lru_key]
                del self.timestamps[lru_key]
        
        self.cache[key] = expr
        self.timestamps[key] = current_time
        self.access_order.append(key)
        
        # Periodic garbage collection
        if len(self.cache) % 20 == 0:
            gc.collect()
        
        return expr
    
    def clear(self):
        """Clear the cache and free memory."""
        self.cache.clear()
        self.timestamps.clear()
        self.access_order.clear()
        gc.collect()


def extract_gene_expression(adata, gene, layer=None, use_raw=False, use_cache=True):
    """
    Optimized gene expression extraction with caching support.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    gene : str
        Gene name to extract
    layer : str, optional
        Layer to extract from (default: None uses .X)
    use_raw : bool
        Whether to use raw data (default: False)
    use_cache : bool
        Whether to use caching (default: True)
    
    Returns:
    --------
    np.ndarray
        Gene expression values as 1D array
    """
    # Use cache if enabled
    if use_cache:
        return _gene_cache.get(adata, gene, layer, use_raw)
    
    # Direct extraction without cache
    # Use raw data if requested
    data_source = adata.raw if use_raw and adata.raw is not None else adata
    
    # Get gene index
    if gene not in data_source.var_names:
        raise ValueError(f"Gene '{gene}' not found in data")
    
    gene_idx = data_source.var_names.get_loc(gene)
    
    # Get expression matrix
    if layer is not None:
        if layer not in data_source.layers:
            raise ValueError(f"Layer '{layer}' not found")
        X = data_source.layers[layer]
    else:
        X = data_source.X
    
    # Extract gene expression
    gene_data = X[:, gene_idx]
    
    # Handle sparse matrices efficiently
    if issparse(gene_data):
        # Only convert the specific column to dense
        gene_expr = gene_data.toarray().flatten()
    else:
        gene_expr = gene_data.flatten() if hasattr(gene_data, 'flatten') else np.asarray(gene_data).flatten()
    
    return gene_expr


def extract_multiple_genes(adata, genes, layer=None, use_raw=False, cell_indices=None):
    """
    Extract multiple genes efficiently in a single operation.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    genes : list of str
        List of gene names to extract
    layer : str, optional
        Layer to extract from (default: None uses .X)
    use_raw : bool
        Whether to use raw data (default: False)
    cell_indices : array-like, optional
        Optional cell indices to subset before extracting genes.
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with genes as columns and cells as rows
    """
    # Use raw data if requested
    data_source = adata.raw if use_raw and adata.raw is not None else adata
    
    # Get gene indices
    gene_indices = []
    valid_genes = []
    for gene in genes:
        if gene in data_source.var_names:
            gene_indices.append(data_source.var_names.get_loc(gene))
            valid_genes.append(gene)
        else:
            print(f"Warning: Gene '{gene}' not found in data")
    
    if not gene_indices:
        raise ValueError("No valid genes found")
    
    # Get expression matrix
    if layer is not None:
        if layer not in data_source.layers:
            raise ValueError(f"Layer '{layer}' not found")
        X = data_source.layers[layer]
    else:
        X = data_source.X
    
    if cell_indices is not None:
        try:
            # Prefer row-first slicing to avoid loading all cells when backed/on-disk.
            gene_data = X[cell_indices, :][:, gene_indices]
        except Exception:
            gene_data = X[:, gene_indices]
            gene_data = gene_data[cell_indices, :]
    else:
        # Extract all genes at once
        gene_data = X[:, gene_indices]
    
    # Convert to dense if sparse
    if issparse(gene_data):
        gene_data = gene_data.toarray()
    elif hasattr(gene_data, 'compute'):  # Handle dask arrays
        gene_data = gene_data.compute()
    
    obs_index = adata.obs_names if cell_indices is None else adata.obs_names[cell_indices]
    # Create DataFrame
    df = pd.DataFrame(gene_data, columns=valid_genes, index=obs_index)
    
    return df


def apply_transformation(expr, method='log1p', copy=True, clip_percentile=99):
    if copy:
        expr = expr.copy()
    
    if method in ['zscore', 'z_score']:
        # Z-score normalization
        if isinstance(expr, pd.DataFrame):
            expr = expr.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
        else:
            mean = np.mean(expr)
            std = np.std(expr)
            expr = (expr - mean) / std if std > 0 else expr - mean
        
        # Clip to ±99% percentile
        upper = np.percentile(expr, clip_percentile)
        lower = np.percentile(expr, 100 - clip_percentile)
        expr = np.clip(expr, lower, upper)

    elif method in ['log1p', 'log']:
        expr = np.log1p(expr)

    return expr

# Global cache instance
_gene_cache = GeneExpressionCache()


def clear_gene_cache():
    """Clear the global gene expression cache."""
    _gene_cache.clear()


def get_cache_info():
    """Get information about the current cache state."""
    return {
        'size': len(_gene_cache.cache),
        'max_size': _gene_cache.max_size,
        'keys': list(_gene_cache.cache.keys())
    }


def bin_cells_for_heatmap(df, gene_columns, groupby, n_bins, continuous_key=None):
    """
    Unified binning function for heatmap visualization to reduce memory usage.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with cells as rows, genes as columns
    gene_columns : list
        List of gene column names
    groupby : str
        Column name for cell type/group
    n_bins : int
        Number of bins to create
    continuous_key : str, optional
        If provided, sort by this continuous variable first (for heatmap2)
    
    Returns:
    --------
    pd.DataFrame
        Binned dataframe with reduced number of rows
    """
    if len(df) == 0:
        return df[[groupby, *gene_columns]].copy() if continuous_key is None else df[[groupby, continuous_key, *gene_columns]].copy()

    n_bins = max(1, min(int(n_bins), len(df)))
    gene_data = df[gene_columns].to_numpy(dtype=np.float32, copy=False)

    # If continuous ordering is requested, use continuous binning approach
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

        group_modes = []
        for s, e in zip(starts, ends):
            values, counts = np.unique(sorted_group_data[s:e], return_counts=True)
            group_modes.append(values[np.argmax(counts)])

        result_df = pd.DataFrame(gene_means, columns=gene_columns)
        result_df.insert(0, continuous_key, continuous_means.astype(np.float32))
        result_df.insert(0, groupby, group_modes)
        return result_df

    # Factorize once to avoid fragile object equality (e.g. tuple group keys).
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
        if n_cells_in_group <= 10:
            binned_blocks.append(group_gene_data)
            binned_group_labels.extend([group] * n_cells_in_group)
            continue

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
