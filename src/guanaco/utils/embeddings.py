"""Shared recognition and materialization helpers for embedding coordinates."""

from __future__ import annotations

import numpy as np
from scipy.sparse import issparse


_EXACT_EMBEDDING_NAMES = {
    "diffmap",
    "fa",
    "harmony",
    "latent",
    "lsi",
    "mde",
    "pca",
    "phate",
    "scvi",
    "spatial",
    "tsne",
    "umap",
}
_EMBEDDING_NAME_PARTS = {
    "diffmap",
    "draw_graph",
    "harmony",
    "latent",
    "lsi",
    "mde",
    "pca",
    "phate",
    "scvi",
    "tsne",
    "umap",
}


def is_embedding_obsm(key: str, values=None) -> bool:
    """Return whether an ``obsm`` entry is recognizably coordinate data."""
    normalized = str(key).strip().lower().replace("-", "_")
    descriptor = normalized.removeprefix("x_")
    recognized = (
        normalized == "spatial"
        or descriptor in _EXACT_EMBEDDING_NAMES
        or any(part in descriptor for part in _EMBEDDING_NAME_PARTS)
    )
    if not recognized:
        return False
    if values is None:
        return True
    shape = getattr(values, "shape", ())
    return len(shape) == 2 and shape[1] >= 2


def embedding_to_numpy(values) -> np.ndarray:
    """Materialize dense or sparse embedding coordinates as a 2D array."""
    if hasattr(values, "to_memory"):
        values = values.to_memory()
    if issparse(values):
        values = values.toarray()
    elif hasattr(values, "compute"):
        values = values.compute()
    array = np.asarray(values)
    if array.ndim != 2 or array.shape[1] < 2:
        raise ValueError("Embedding coordinates must be a two-dimensional matrix.")
    return array
