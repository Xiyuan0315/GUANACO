"""Compatibility tests for legacy entry points backed by shared readers."""

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from guanaco.utils.memory_utils import (
    LazyAnnData,
    memory_efficient_gene_expression,
    sparse_safe_slice,
    sparse_to_dense_batched,
)


@pytest.mark.parametrize("storage", [np.asarray, sparse.csr_matrix, sparse.csc_matrix])
@pytest.mark.parametrize("transformation", [None, "log", "z_score", "unknown"])
@pytest.mark.parametrize("constant", [False, True])
def test_legacy_expression_preserves_values_dtype_and_transformation(storage, transformation, constant):
    values = np.full(4, 3.0) if constant else np.array([0.0, 1.0, 2.0, 100.0])
    data = ad.AnnData(storage(values[:, None]))
    expected = values.copy()
    if transformation == "log":
        expected = np.log1p(expected)
    elif transformation == "z_score" and not constant:
        expected = (expected - expected.mean()) / expected.std()

    actual = memory_efficient_gene_expression(data, "0", transformation)

    np.testing.assert_allclose(actual, expected)
    assert actual.dtype == values.dtype
    original = data.X.toarray() if sparse.issparse(data.X) else data.X
    np.testing.assert_array_equal(original[:, 0], values)


def test_legacy_expression_missing_gene_raises():
    with pytest.raises(KeyError):
        memory_efficient_gene_expression(ad.AnnData(np.ones((2, 1))), "missing")


@pytest.mark.parametrize("storage", [np.asarray, sparse.csr_matrix, sparse.csc_matrix])
@pytest.mark.parametrize("mode", ["backed", "lazy"])
def test_legacy_expression_uses_shared_reader_for_filtered_disk_and_lazy_data(tmp_path, storage, mode):
    values = np.arange(12, dtype=np.float64).reshape(4, 3)
    data = ad.AnnData(storage(values))
    if mode == "backed":
        path = tmp_path / "expression.h5ad"
        data.write_h5ad(path)
        data = ad.read_h5ad(path, backed="r")
    else:
        da = pytest.importorskip("dask.array")
        data.X = da.from_array(data.X, chunks=(2, 3))
    try:
        actual = memory_efficient_gene_expression(data[[3, 1]], "2", "log")
        np.testing.assert_allclose(actual, np.log1p(values[[3, 1], 2]))
        assert actual.dtype == values.dtype
    finally:
        if mode == "backed":
            data.file.close()


@pytest.mark.parametrize("storage", [np.asarray, sparse.csr_matrix, sparse.csc_matrix])
@pytest.mark.parametrize("axis", [0, 1])
def test_legacy_slice_preserves_order_and_sparse_storage(storage, axis):
    values = np.arange(12).reshape(3, 4)
    actual = sparse_safe_slice(storage(values), [2, 0, 2], axis=axis)
    expected = values[[2, 0, 2]] if axis == 0 else values[:, [2, 0, 2]]
    assert sparse.issparse(actual) == (storage is not np.asarray)
    np.testing.assert_array_equal(actual.toarray() if sparse.issparse(actual) else actual, expected)


@pytest.mark.parametrize("storage", [np.asarray, sparse.csr_matrix, sparse.csc_matrix])
def test_legacy_batched_conversion_preserves_values(storage):
    expected = np.arange(12, dtype=np.float32).reshape(4, 3)
    np.testing.assert_array_equal(sparse_to_dense_batched(storage(expected), batch_size=3), expected)


def test_lazy_wrapper_defers_loading_and_reuses_loaded_data(tmp_path):
    path = tmp_path / "later.h5ad"
    wrapper = LazyAnnData(path, max_cells=2, seed=0)
    # Construction must not attempt to read a file that does not yet exist.
    ad.AnnData(np.arange(12, dtype=np.float32).reshape(4, 3)).write_h5ad(path)
    loaded = wrapper.adata
    assert wrapper.adata is loaded
    assert wrapper.n_obs == 2
    np.testing.assert_array_equal(wrapper.X.toarray() if sparse.issparse(wrapper.X) else wrapper.X, loaded.X)
