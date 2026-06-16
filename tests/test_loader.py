"""Tests for sc_data loading across .h5ad, .h5mu, local .zarr and cloud .zarr."""

from pathlib import Path
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from guanaco.data import loader


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _toy_adata(n_obs: int = 40, n_vars: int = 6) -> ad.AnnData:
    rng = np.random.default_rng(0)
    X = rng.random((n_obs, n_vars)).astype("float32")
    obs = {
        "cell_type": np.array(["A", "B"] * (n_obs // 2)),
    }
    return ad.AnnData(X=X, obs=obs)


@pytest.fixture
def h5ad_file(tmp_path: Path) -> Path:
    path = tmp_path / "data.h5ad"
    _toy_adata().write_h5ad(path)
    return path


@pytest.fixture
def h5mu_file(tmp_path: Path) -> Path:
    import muon as mu

    mdata = mu.MuData({"rna": _toy_adata()})
    path = tmp_path / "data.h5mu"
    mdata.write_h5mu(path)
    return path


@pytest.fixture
def zarr_store(tmp_path: Path) -> Path:
    path = tmp_path / "data.zarr"
    _toy_adata().write_zarr(path)
    return path


def _sparse_adata(n_obs: int = 120, n_vars: int = 10, fmt: str = "csc") -> ad.AnnData:
    from scipy.sparse import csc_matrix, csr_matrix

    rng = np.random.default_rng(1)
    dense = (rng.random((n_obs, n_vars)) > 0.6) * rng.random((n_obs, n_vars))
    dense = dense.astype("float32")
    mat = csc_matrix(dense) if fmt == "csc" else csr_matrix(dense)
    adata = ad.AnnData(
        X=mat,
        obs={"cell_type": np.array(["A", "B"] * (n_obs // 2))},
        var={"gene_ids": np.arange(n_vars)},
    )
    adata.uns["_dense"] = dense  # ground truth for assertions
    return adata


@pytest.fixture
def csc_zarr_store(tmp_path: Path):
    path = tmp_path / "csc.zarr"
    adata = _sparse_adata(fmt="csc")
    dense = adata.uns.pop("_dense")
    adata.write_zarr(path)
    return path, dense


@pytest.fixture
def csr_zarr_store(tmp_path: Path) -> Path:
    path = tmp_path / "csr.zarr"
    adata = _sparse_adata(fmt="csr")
    adata.uns.pop("_dense")
    adata.write_zarr(path)
    return path


@pytest.fixture
def csr_with_csc_layer_store(tmp_path: Path):
    """A store like the public demo: CSR X plus a gene-major CSC copy as a layer."""
    from scipy.sparse import csc_matrix

    path = tmp_path / "csr_plus_csc.zarr"
    adata = _sparse_adata(fmt="csr")
    dense = adata.uns.pop("_dense")
    adata.layers["X_csc"] = csc_matrix(dense)
    adata.write_zarr(path)
    return path, dense


# ---------------------------------------------------------------------------
# Remote-URI detection helpers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "uri",
    [
        "s3://bucket/path/dataset.zarr",
        "gs://bucket/path/dataset.zarr",
        "https://example.com/path/dataset.zarr",
        "http://example.com/path/dataset.zarr",
    ],
)
def test_is_remote_uri_true(uri):
    assert loader._is_remote_uri(uri) is True


@pytest.mark.parametrize(
    "value",
    ["/abs/local/data.zarr", "/abs/local/data.h5ad", "relative/data.zarr", Path("/abs/x.zarr")],
)
def test_is_remote_uri_false(value):
    assert loader._is_remote_uri(value) is False


@pytest.mark.parametrize(
    "value,expected",
    [
        ("/abs/data.h5ad", ".h5ad"),
        ("/abs/data.h5mu", ".h5mu"),
        ("/abs/data.zarr", ".zarr"),
        ("/abs/data.zarr/", ".zarr"),
        ("s3://bucket/path/dataset.zarr", ".zarr"),
        ("https://example.com/x/dataset.zarr?versionId=abc", ".zarr"),
    ],
)
def test_source_suffix(value, expected):
    assert loader._source_suffix(value) == expected


# ---------------------------------------------------------------------------
# Local loading still works
# ---------------------------------------------------------------------------

def test_load_h5ad(h5ad_file):
    adata = loader.load_adata(str(h5ad_file), max_cells=None)
    assert adata.n_obs == 40


def test_load_h5mu(h5mu_file):
    mdata = loader.load_adata(str(h5mu_file), max_cells=None)
    assert hasattr(mdata, "mod")
    assert mdata.mod["rna"].n_obs == 40


def test_load_local_zarr(zarr_store):
    adata = loader.load_adata(str(zarr_store), max_cells=None)
    assert adata.n_obs == 40


def test_local_zarr_downsamples(zarr_store):
    adata = loader.load_adata(str(zarr_store), max_cells=10, seed=0)
    assert adata.n_obs == 10


def test_backed_csc_zarr_keeps_x_lazy(csc_zarr_store):
    """Cloud-backed mode keeps X lazy (dask) but materializes obs/var into memory."""
    path, _ = csc_zarr_store
    adata = loader.load_adata(str(path), backed=True)
    # X stays lazy; metadata is in-memory pandas.
    assert hasattr(adata.X, "compute"), f"expected lazy dask X, got {type(adata.X)}"
    assert isinstance(adata.obs, pd.DataFrame)
    assert isinstance(adata.var, pd.DataFrame)
    assert adata.n_obs == 120


def test_backed_csc_zarr_reads_only_requested_genes(csc_zarr_store):
    """Gene columns are pulled correctly on demand from the lazy backed store."""
    from guanaco.utils.gene_extraction_utils import (
        extract_gene_expression,
        clear_gene_cache,
    )

    path, dense = csc_zarr_store
    adata = loader.load_adata(str(path), backed=True)
    clear_gene_cache()

    # var names are g-less ints by default; use the actual names.
    names = adata.var_names.tolist()
    j = 3
    vec = extract_gene_expression(adata, names[j])
    assert np.allclose(vec, dense[:, j])

    # Multiple genes are read one cached column at a time.
    cols = [names[1], names[4], names[7]]
    got = np.column_stack([extract_gene_expression(adata, c) for c in cols])
    assert np.allclose(got, dense[:, [1, 4, 7]])


def test_backed_csc_zarr_caches_genes(csc_zarr_store):
    """A re-requested gene is served from cache, not re-read."""
    from guanaco.utils.gene_extraction_utils import (
        extract_gene_expression,
        clear_gene_cache,
        get_cache_info,
    )

    path, _ = csc_zarr_store
    adata = loader.load_adata(str(path), backed=True)
    clear_gene_cache()
    name = adata.var_names[2]
    extract_gene_expression(adata, name)
    extract_gene_expression(adata, name)  # second hit -> cache
    info = get_cache_info()
    # exactly one resident vector for the single distinct gene requested
    assert info["size"] == 1


def test_backed_csr_zarr_loads(csr_zarr_store):
    """Any encoding loads in backed mode: X stays lazy, metadata is in memory."""
    adata = loader.load_adata(str(csr_zarr_store), backed=True)
    assert hasattr(adata.X, "compute")
    assert isinstance(adata.obs, pd.DataFrame)
    assert isinstance(adata.var, pd.DataFrame)


def test_backed_csr_zarr_reads_genes(csr_zarr_store):
    """A backed gene read returns correct values regardless of encoding."""
    from guanaco.utils.gene_extraction_utils import extract_gene_expression, clear_gene_cache

    adata = loader.load_adata(str(csr_zarr_store), backed=True)
    clear_gene_cache()
    name = adata.var_names[2]
    vec = extract_gene_expression(adata, name)
    assert vec.shape[0] == adata.n_obs


def test_nonbacked_csr_zarr_allowed(csr_zarr_store):
    """CSR is fine for the in-memory (non-backed) path -- it gets materialized."""
    adata = loader.load_adata(str(csr_zarr_store), backed=False, max_cells=None)
    assert adata.n_obs == 120


def _x_chunktype(adata):
    return type(adata.X._meta).__name__


def test_backed_auto_swaps_csc_layer(csr_with_csc_layer_store):
    """Backed mode auto-serves expression from the CSC layer when X is CSR."""
    path, dense = csr_with_csc_layer_store
    adata = loader.load_adata(str(path), backed=True)
    # The CSC layer is swapped into X (now gene-major) and consumed.
    assert _x_chunktype(adata) == "csc_matrix"
    assert "X_csc" not in adata.layers
    # X is rechunked to one gene per column so a single-gene read fetches one column
    # (not read_lazy's default ~1000-gene chunk).
    assert adata.X.chunksize == (adata.n_obs, 1)
    # Reads are still correct.
    j = 4
    col = adata.X[:, j].compute()
    out = col.toarray().ravel()
    assert np.allclose(out, dense[:, j])


def test_backed_explicit_expression_layer(csr_with_csc_layer_store):
    """An explicit expression_layer is honoured."""
    path, dense = csr_with_csc_layer_store
    adata = loader.load_adata(str(path), backed=True, expression_layer="X_csc")
    assert _x_chunktype(adata) == "csc_matrix"
    assert "X_csc" not in adata.layers


def test_backed_unknown_expression_layer_falls_back(csr_with_csc_layer_store, capsys):
    """A bad expression_layer name warns and falls back to X (no crash)."""
    path, _ = csr_with_csc_layer_store
    adata = loader.load_adata(str(path), backed=True, expression_layer="nope")
    assert "not found" in capsys.readouterr().out
    # X unchanged (still CSR), layer still present.
    assert _x_chunktype(adata) == "csr_matrix"


def test_local_zarr_uses_read_lazy(zarr_store):
    """The AnnData zarr path goes through anndata.experimental.read_lazy."""
    import anndata.experimental as aexp

    real = aexp.read_lazy
    calls = {"n": 0}

    def spy(store, *a, **k):
        calls["n"] += 1
        return real(store, *a, **k)

    with mock.patch.object(aexp, "read_lazy", side_effect=spy):
        adata = loader.load_adata(str(zarr_store), max_cells=None)

    assert calls["n"] >= 1
    assert adata.n_obs == 40


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

def test_unsupported_extension_raises(tmp_path):
    bad = tmp_path / "data.txt"
    bad.write_text("nope")
    with pytest.raises(ValueError, match="Unsupported file extension"):
        loader.load_adata(str(bad))


def test_relative_path_rejected():
    with pytest.raises(ValueError, match="must be absolute"):
        loader.load_adata("relative/data.h5ad")


def test_remote_h5ad_downloaded_then_read_locally(tmp_path):
    """A remote .h5ad/.h5mu is no longer rejected: it's downloaded to a local cache
    and read from disk (HDF5 can't stream over the network like .zarr)."""
    fake_local = tmp_path / "cached.h5ad"
    sentinel = object()
    with mock.patch.object(
        loader, "_download_remote_to_cache", return_value=fake_local
    ) as dl, mock.patch.object(loader.ad, "read_h5ad", return_value=sentinel):
        result = loader.load_adata("s3://bucket/path/dataset.h5ad", backed=True)

    dl.assert_called_once_with("s3://bucket/path/dataset.h5ad")
    assert result is sentinel


def test_remote_unsupported_suffix_rejected():
    with pytest.raises(ValueError, match="Remote sc_data must be a .zarr store or an"):
        loader.load_adata("s3://bucket/path/dataset.csv")


# ---------------------------------------------------------------------------
# Cloud .zarr (mocked -- no real network/fsspec required)
# ---------------------------------------------------------------------------

def test_remote_zarr_not_path_mangled():
    """A cloud URI reaches the zarr opener verbatim, never collapsed to a Path."""
    import anndata.experimental as aexp

    captured = {}

    def fake_open(store):
        captured["store"] = store
        return object()  # sentinel; read_lazy is mocked so it is never consumed

    with mock.patch.object(loader, "_zarr_is_mudata", return_value=False), mock.patch.object(
        loader, "_open_zarr_group", side_effect=fake_open
    ), mock.patch.object(aexp, "read_lazy", return_value=_toy_adata()):
        adata = loader.load_adata("s3://my-bucket/path/my_dataset.zarr", max_cells=None)

    assert captured["store"] == "s3://my-bucket/path/my_dataset.zarr"
    assert adata.n_obs == 40


def test_remote_backed_zarr_uses_uri_verbatim():
    """Backed mode opens the cloud URI unchanged (no Path mangling), keeping X remote."""
    import anndata.experimental as aexp

    captured = {}

    def fake_open(store):
        captured["store"] = store
        return object()

    with mock.patch.object(loader, "_zarr_is_mudata", return_value=False), mock.patch.object(
        loader, "_open_zarr_group", side_effect=fake_open
    ), mock.patch.object(aexp, "read_lazy", return_value=_toy_adata()):
        adata = loader.load_adata("s3://my-bucket/path/my_dataset.zarr", backed=True)

    assert captured["store"] == "s3://my-bucket/path/my_dataset.zarr"
    assert adata.n_obs == 40


def test_remote_zarr_mudata_uses_muon():
    import muon as mu

    captured = {}

    def fake_read_zarr(store):
        captured["store"] = store
        return mu.MuData({"rna": _toy_adata()})

    with mock.patch.object(loader, "_zarr_is_mudata", return_value=True), mock.patch.object(
        loader.mu, "read_zarr", side_effect=fake_read_zarr
    ):
        mdata = loader.load_adata("gs://my-bucket/path/my_dataset.zarr", max_cells=None)

    assert captured["store"] == "gs://my-bucket/path/my_dataset.zarr"
    assert hasattr(mdata, "mod")


# ---------------------------------------------------------------------------
# initialize_data + lazy loading with cloud URIs
# ---------------------------------------------------------------------------

def test_initialize_data_accepts_cloud_uri_lazily():
    cfg = {
        "Demo": {"sc_data": "s3://bucket/path/dataset.zarr", "description": "d"},
    }
    datasets = loader.initialize_data(cfg=cfg, lazy_load=True)
    bundle = datasets["Demo"]
    # Original cloud URI stored verbatim for later lazy load.
    assert bundle.adata_path == "s3://bucket/path/dataset.zarr"
    assert bundle._adata is None


def test_lazy_load_calls_load_adata_with_cloud_uri():
    cfg = {
        "Demo": {"sc_data": "s3://bucket/path/dataset.zarr", "description": "d"},
    }
    datasets = loader.initialize_data(cfg=cfg, lazy_load=True)
    bundle = datasets["Demo"]

    with mock.patch.object(loader, "load_adata", return_value=_toy_adata()) as m:
        _ = bundle.adata

    m.assert_called_once()
    assert m.call_args.args[0] == "s3://bucket/path/dataset.zarr"


def test_initialize_data_rejects_relative_sc_data():
    cfg = {"Demo": {"sc_data": "relative/data.zarr"}}
    with pytest.raises(ValueError, match="must be absolute"):
        loader.initialize_data(cfg=cfg, lazy_load=True)


# ---------------------------------------------------------------------------
# Discrete label detection
# ---------------------------------------------------------------------------

def _labels_adata(n_obs: int = 200) -> ad.AnnData:
    rng = np.random.default_rng(3)
    obs = pd.DataFrame(
        {
            # categorical with few levels -> a discrete label
            "cell_type": pd.Categorical(np.array(["A", "B", "C"])[rng.integers(0, 3, n_obs)]),
            # plain object with few levels -> a discrete label
            "condition": np.array(["ctrl", "treat"])[rng.integers(0, 2, n_obs)],
            # continuous -> not a label
            "n_genes": rng.random(n_obs).astype("float32"),
            # high-cardinality categorical -> not a label
            "barcode": pd.Categorical([f"bc{i}" for i in range(n_obs)]),
        }
    )
    return ad.AnnData(X=rng.random((n_obs, 4)).astype("float32"), obs=obs)


def test_get_discrete_labels_selects_low_cardinality():
    labels = loader.get_discrete_labels(_labels_adata(), max_unique=50)
    assert labels == ["condition", "cell_type"]  # sorted by ascending cardinality (2, 3)


def test_get_discrete_labels_uses_declared_categories_for_categoricals():
    # An unused category still counts toward cardinality (O(1) metadata read),
    # so a categorical declaring >= max_unique levels is excluded.
    obs = pd.DataFrame(
        {"many": pd.Categorical(["A", "B"], categories=[f"c{i}" for i in range(60)])}
    )
    adata = ad.AnnData(X=np.zeros((2, 2), dtype="float32"), obs=obs)
    assert loader.get_discrete_labels(adata, max_unique=50) == []


def test_get_discrete_labels_empty_obs():
    adata = ad.AnnData(X=np.zeros((3, 2), dtype="float32"))
    assert loader.get_discrete_labels(adata) == []
