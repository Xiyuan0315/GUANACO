"""Bounded column reads and grouped summaries shared by dot plots and ATAC."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse

from guanaco.utils.gene_extraction_utils import iter_feature_blocks


@dataclass
class FeatureSummary:
    means: np.ndarray
    fractions: np.ndarray
    mean: np.ndarray
    std: np.ndarray


def grouped_feature_stats(
    adata, indices, codes, n_groups, *, rows=None, layer=None,
    transformation=None, skip_nan=True, moments=False, read_block=None, metric="both",
):
    """Summarize selected cells without a cells-by-all-features dense matrix.

    ``codes`` describes selected rows; -1 excludes a row from group summaries,
    but not from overall moments. Detection treats NaN as not positive, matching
    (values > 0).mean(). ATAC uses skip_nan=False, dot plots use pandas semantics.
    Only elementwise transforms belong here; global percentile clipping stays
    with the caller because clipping each block would change the result.
    """
    indices = np.asarray(indices, dtype=np.int64)
    codes = np.asarray(codes, dtype=np.int64)
    counts = np.bincount(codes[codes >= 0], minlength=n_groups)
    n_cells = len(codes)
    result = FeatureSummary(
        np.full((n_groups, len(indices)), np.nan),
        np.full((n_groups, len(indices)), np.nan),
        np.full(len(indices), np.nan), np.full(len(indices), np.nan),
    )
    if not n_cells:
        return result
    group_cat = pd.Categorical.from_codes(codes, categories=range(n_groups))
    members = None
    for start, block in iter_feature_blocks(adata, indices, layer, read_block=read_block):
        if sparse.issparse(block):
            block = block.tocsc(copy=True)
            block.sum_duplicates()
        if rows is not None:
            block = block[rows, :]
        block = block.astype(np.float32, copy=False)
        if transformation in ("log", "log1p"):
            with np.errstate(invalid="ignore", divide="ignore"):
                if sparse.issparse(block):
                    block.data = np.log1p(block.data)
                else:
                    block = np.log1p(block)

        if sparse.issparse(block) and np.isfinite(block.data).all():
            # Work only on stored values. Include implicit zeros in denominators
            # and variance, and use centered moments to avoid cancellation.
            for j in range(block.shape[1]):
                lo, hi = block.indptr[j:j + 2]
                values = block.data[lo:hi].astype(np.float64)
                groups = codes[block.indices[lo:hi]]
                included = groups >= 0
                col = start + j
                if metric != "detection":
                    sums = np.bincount(groups[included], weights=values[included], minlength=n_groups)
                    np.divide(sums, counts, out=result.means[:, col], where=counts > 0)
                if metric != "mean":
                    positive = np.bincount(groups[included & (values > 0)], minlength=n_groups)
                    np.divide(positive, counts, out=result.fractions[:, col], where=counts > 0)
                if moments and n_cells:
                    mean = values.sum() / n_cells
                    result.mean[col] = mean
                    if n_cells > 1:
                        variance = ((values - mean) ** 2).sum() + (n_cells - len(values)) * mean ** 2
                        result.std[col] = np.sqrt(variance / (n_cells - 1))
        else:
            # Dense/backed and non-finite data use the original NaN semantics,
            # but only for the current bounded block.
            values = block.toarray() if sparse.issparse(block) else np.asarray(block)
            frame = pd.DataFrame(values)
            stop = start + block.shape[1]
            if metric != "detection":
                if skip_nan:
                    means = frame.groupby(group_cat, observed=True, sort=False).mean().reindex(range(n_groups))
                    result.means[:, start:stop] = means.to_numpy()
                else:
                    if members is None:
                        members = [np.flatnonzero(codes == group) for group in range(n_groups)]
                    for group, positions in enumerate(members):
                        if len(positions):
                            result.means[group, start:stop] = values[positions].mean(axis=0)
            if metric != "mean":
                fractions = (frame > 0).groupby(group_cat, observed=True, sort=False).mean().reindex(range(n_groups))
                result.fractions[:, start:stop] = fractions.to_numpy()
            if moments:
                result.mean[start:stop] = frame.mean().to_numpy()
                result.std[start:stop] = frame.std().to_numpy()
    return result
