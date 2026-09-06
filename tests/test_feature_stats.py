"""Numerical equivalence and bounded reads, rather than timing-sensitive tests."""

import weakref

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from guanaco.utils import feature_stats as stats
from guanaco.utils import gene_extraction_utils as genes


@pytest.mark.parametrize("consumer", ["gene", "heatmap", "violin"])
def test_caches_distinguish_equal_sized_backed_views(tmp_path, monkeypatch, consumer):
    from guanaco.pages.matrix.plots.heatmap import plot_unified_heatmap
    from guanaco.pages.matrix.plots.violin1 import plot_violin1

    x = np.arange(12, dtype=np.float32).reshape(6, 2)
    path = tmp_path / "views.h5ad"
    ad.AnnData(x, obs=pd.DataFrame({"group": ["A"] * 6}, index=list("abcdef"))).write_h5ad(path)
    data = ad.read_h5ad(path, backed="r")
    monkeypatch.setattr(genes, "_gene_cache", genes.GeneExpressionCache())

    def read(view):
        if consumer == "heatmap":
            fig = plot_unified_heatmap(view, ["0"], "group", labels=["A"], color_config=["red"])
            return np.asarray(fig.data[0].z).ravel()
        if consumer == "violin":
            fig = plot_violin1(view, ["0"], "group", labels=["A"])
            return np.asarray(fig.data[0].y)
        return genes.extract_gene_expression(view, "0")

    try:
        first, second = data[:3], data[3:]
        np.testing.assert_array_equal(read(first), x[:3, 0])
        np.testing.assert_array_equal(read(second), x[3:, 0])
        np.testing.assert_array_equal(read(first), x[:3, 0])
    finally:
        data.file.close()


def test_dataset_cache_identity_does_not_retain_data():
    import gc

    data = ad.AnnData(np.ones((2, 2)))
    ref = weakref.ref(data)
    token = genes.dataset_cache_token(data)
    assert genes.dataset_cache_token(data) == token
    del data
    gc.collect()
    assert ref() is None
    assert genes.dataset_cache_token(ad.AnnData(np.ones((2, 2)))) != token


@pytest.mark.parametrize("storage", [np.asarray, sparse.csr_matrix, sparse.csc_matrix])
@pytest.mark.parametrize("transform", [None, "log1p"])
@pytest.mark.parametrize("nonfinite", [False, True])
@pytest.mark.parametrize("skip_nan", [False, True])
def test_grouped_stats_match_dense_reference(monkeypatch, storage, transform, nonfinite, skip_nan):
    x = np.array([[0, 2, 0], [1, 0, 3], [2, 4, 0], [0, 0, 0], [5, 2, 1]], dtype=np.float32)
    if nonfinite:
        x[0, 1], x[2, 0] = np.nan, np.inf
    data = ad.AnnData(storage(x))
    rows = np.array([4, 0, 2, 3])
    codes = np.array([1, 0, 1, -1])
    cols = [2, 0, 1]
    reads = []
    monkeypatch.setattr(genes, "GENE_READ_BUDGET", data.n_obs * 8)

    def read(indices):
        reads.append(list(indices))
        return data.X[:, indices]

    actual = stats.grouped_feature_stats(
        data, cols, codes, 3, rows=rows, transformation=transform,
        skip_nan=skip_nan, moments=True, read_block=read,
    )
    assert reads == [[2], [0], [1]]
    values = x[rows][:, cols]
    if transform:
        values = np.log1p(values)
    frame = pd.DataFrame(values)
    cat = pd.Categorical.from_codes(codes, categories=range(3))
    expected_mean = frame.groupby(cat, observed=True).mean().reindex(range(3)).to_numpy()
    if not skip_nan:
        for group in [0, 1]:
            expected_mean[group] = values[codes == group].mean(axis=0)
    expected_fraction = (frame > 0).groupby(cat, observed=True).mean().reindex(range(3))
    np.testing.assert_allclose(actual.means, expected_mean, rtol=1e-6, equal_nan=True)
    np.testing.assert_allclose(actual.fractions, expected_fraction, equal_nan=True)
    np.testing.assert_allclose(actual.mean, frame.mean(), rtol=1e-6, equal_nan=True)
    np.testing.assert_allclose(actual.std, frame.std(), rtol=1e-6, equal_nan=True)


def test_finite_sparse_summary_does_not_densify(monkeypatch):
    # Duplicate entries and explicit zeros must count as one cell, not two hits.
    x = sparse.csc_matrix((np.array([1., 2., 0., -2.]), [0, 0, 1, 2], [0, 4]), shape=(4, 1))
    data = ad.AnnData(x)

    def forbid(*args, **kwargs):
        raise AssertionError("A sparse summary must not allocate a dense cell matrix")

    monkeypatch.setattr(sparse.csc_matrix, "toarray", forbid)
    result = stats.grouped_feature_stats(data, [0], [0, 0, 1, 1], 2, moments=True)
    np.testing.assert_allclose(result.means[:, 0], [1.5, -1])
    np.testing.assert_allclose(result.fractions[:, 0], [0.5, 0])
    np.testing.assert_allclose(result.std, np.std([3, 0, -2, 0], ddof=1))


@pytest.mark.parametrize("storage", [np.asarray, sparse.csc_matrix])
@pytest.mark.parametrize("metric", ["mean", "detection"])
def test_summary_only_computes_requested_metric(storage, metric):
    data = ad.AnnData(storage(np.array([[1, 0], [3, 1]], dtype=np.float32)))
    result = stats.grouped_feature_stats(data, [0, 1], [0, 0], 1, metric=metric, skip_nan=False)
    if metric == "mean":
        np.testing.assert_allclose(result.means, [[2, 0.5]])
        assert np.isnan(result.fractions).all()
    else:
        np.testing.assert_allclose(result.fractions, [[1, 0.5]])
        assert np.isnan(result.means).all()


@pytest.mark.parametrize("reader", ["genes", "blocks"])
def test_batched_reads_protect_future_cache_hits(monkeypatch, reader):
    data = ad.AnnData(np.arange(24, dtype=np.float32).reshape(4, 6))
    expected = data.X.copy()
    if reader == "blocks":
        da = pytest.importorskip("dask.array")
        data.X = da.from_array(data.X, chunks=(4, 1))
    cache = genes.GeneExpressionCache(max_size=1)
    monkeypatch.setattr(genes, "_gene_cache", cache)
    monkeypatch.setattr(genes, "GENE_READ_BUDGET", data.n_obs * 8 * 2)
    cache.store(data, "5", None, False, expected[:, 5].copy())
    original = genes._compute_gene_block
    reads = []

    def read(adata, names, *args):
        reads.append(list(names))
        return original(adata, names, *args)

    monkeypatch.setattr(genes, "_compute_gene_block", read)
    if reader == "genes":
        result = np.column_stack([v for _, v in genes.iter_gene_expression(data, data.var_names)])
    else:
        result = np.column_stack([b for _, b in genes.iter_feature_blocks(data, range(data.n_vars))])
    assert reads == [["0", "1"], ["2", "3"], ["4"]]
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("storage", [sparse.csr_matrix, sparse.csc_matrix])
def test_feature_blocks_preserve_sparse_storage_and_request_order(monkeypatch, storage):
    x = np.arange(20, dtype=np.float32).reshape(4, 5)
    data = ad.AnnData(storage(x))
    monkeypatch.setattr(genes, "GENE_READ_BUDGET", data.n_obs * 8 * 2)
    blocks = list(genes.iter_feature_blocks(data, [4, 1, 4, 0, 2]))
    assert [start for start, _ in blocks] == [0, 2, 4]
    assert all(sparse.issparse(block) for _, block in blocks)
    np.testing.assert_array_equal(sparse.hstack([b for _, b in blocks]).toarray(), x[:, [4, 1, 4, 0, 2]])


def test_custom_block_reader_keeps_its_own_cache_policy(monkeypatch):
    da = pytest.importorskip("dask.array")
    x = np.arange(12, dtype=np.float32).reshape(4, 3)
    data = ad.AnnData(da.from_array(x, chunks=(4, 1)))

    def forbid(*args, **kwargs):
        raise AssertionError("A supplied reader must not use the default gene cache")

    monkeypatch.setattr(genes, "iter_gene_expression", forbid)
    blocks = list(genes.iter_feature_blocks(data, [2, 0], read_block=lambda cols: sparse.csc_matrix(x[:, cols])))
    np.testing.assert_array_equal(blocks[0][1].toarray(), x[:, [2, 0]])


def test_feature_iterator_does_not_retain_consumed_dense_block():
    data = ad.AnnData(np.ones((4, 3), dtype=np.float32))
    blocks = genes.iter_feature_blocks(data, [0, 1])
    _, block = next(blocks)
    reference = weakref.ref(block)
    del block
    assert reference() is None, "The iterator must not keep an extra dense block alive"


@pytest.mark.parametrize("backed", [False, True])
@pytest.mark.parametrize("source", ["X", "layer", "raw"])
@pytest.mark.parametrize("storage", [np.asarray, sparse.csr_matrix])
def test_feature_reads_support_views_layers_raw_and_backed(tmp_path, backed, source, storage):
    x = np.arange(30, dtype=np.float32).reshape(5, 6)
    data = ad.AnnData(storage(x))
    data.layers["counts"] = sparse.csc_matrix(x * 2)
    data.raw = data.copy()
    if backed:
        path = tmp_path / "features.h5ad"
        data.write_h5ad(path)
        data = ad.read_h5ad(path, backed="r")
    view = data[[4, 1, 3], :]
    try:
        result = list(genes.iter_gene_expression(
            view, ["4", "1", "4"], layer="counts" if source == "layer" else None,
            use_raw=source == "raw",
        ))
        expected = x[[4, 1, 3]][:, [4, 1, 4]] * (2 if source == "layer" else 1)
        np.testing.assert_array_equal(np.column_stack([v for _, v in result]), expected)
        if source != "raw":
            summary = stats.grouped_feature_stats(
                view, [4, 1, 4], [1, 0, 1], 2, layer="counts" if source == "layer" else None,
            )
            np.testing.assert_allclose(summary.means, [expected[1], expected[[0, 2]].mean(axis=0)])
    finally:
        genes.clear_gene_cache()
        if backed:
            data.file.close()


def test_prewarm_does_not_load_more_genes_than_cache_can_hold(monkeypatch):
    data = ad.AnnData(np.ones((3, 6), dtype=np.float32))
    monkeypatch.setattr(genes, "_gene_cache", genes.GeneExpressionCache(max_size=2))
    genes.prewarm_gene_cache(data, data.var_names.tolist())
    assert [key[2] for key in genes._gene_cache._data] == ["0", "1"]


@pytest.mark.parametrize("storage", [np.asarray, sparse.csr_matrix, sparse.csc_matrix])
def test_lazy_feature_read_and_summary(storage):
    da = pytest.importorskip("dask.array")
    x = np.array([[0, 2, 0], [1, 0, 3], [2, 4, 0]], dtype=np.float32)
    data = ad.AnnData(da.from_array(storage(x), chunks=(2, 1)))
    result = stats.grouped_feature_stats(data, [2, 0], [1, 0, 1], 2)
    np.testing.assert_allclose(result.means, [[3, 1], [0, 1]])
