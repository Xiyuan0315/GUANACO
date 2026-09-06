"""Legacy memory helpers; new data access uses gene_extraction_utils / data.loader."""

import gc
import numpy as np
from scipy import sparse
from functools import wraps
import psutil
import os

from guanaco.utils.gene_extraction_utils import extract_gene_expression

def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def sparse_safe_slice(matrix, indices, axis=0):
    """Efficiently slice sparse matrix without converting to dense."""
    return matrix[indices] if axis == 0 else matrix[:, indices]

def sparse_to_dense_batched(sparse_matrix, batch_size=1000):
    """Convert sparse matrix to dense in batches to reduce peak memory."""
    if not sparse.issparse(sparse_matrix):
        return sparse_matrix
    
    n_rows = sparse_matrix.shape[0]
    result = np.empty(sparse_matrix.shape, dtype=sparse_matrix.dtype)
    
    for i in range(0, n_rows, batch_size):
        end_idx = min(i + batch_size, n_rows)
        result[i:end_idx] = sparse_matrix[i:end_idx].toarray()
    
    return result

def memory_efficient_gene_expression(adata, gene_name, transformation=None):
    """Compatibility wrapper preserving source dtype and legacy transformations."""
    gene_expr = extract_gene_expression(adata, gene_name, use_cache=False, dtype=None)
    # Legacy z-scores use population std, no clipping, and keep constants unchanged.
    
    if transformation == 'log':
        gene_expr = np.log1p(gene_expr)
    elif transformation == 'z_score':
        mean = np.mean(gene_expr)
        std = np.std(gene_expr)
        if std > 0:
            gene_expr = (gene_expr - mean) / std
    
    return gene_expr

def cleanup_memory():
    """Force garbage collection to free memory."""
    gc.collect()

def memory_monitor(func):
    """Decorator to monitor memory usage of functions."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_mem = get_memory_usage()
        result = func(*args, **kwargs)
        end_mem = get_memory_usage()
        
        if end_mem - start_mem > 100:  # Log if > 100MB increase
            print(f"Memory increased by {end_mem - start_mem:.1f}MB in {func.__name__}")
        
        return result
    return wrapper

class LazyAnnData:
    """Lazy wrapper for AnnData that loads data on demand."""
    
    def __init__(self, file_path, max_cells=None, seed=None):
        self.file_path = file_path
        self.max_cells = max_cells
        self.seed = seed
        self._adata = None
    
    @property
    def adata(self):
        """Load AnnData on first access."""
        if self._adata is None:
            from guanaco.data.loader import load_adata
            self._adata = load_adata(
                self.file_path,
                max_cells=self.max_cells,
                seed=self.seed
            )
        return self._adata
    
    def __getattr__(self, name):
        """Delegate attribute access to the loaded AnnData."""
        return getattr(self.adata, name)
