"""Tests for sc_data loading across .h5ad, .h5mu, local .zarr and cloud .zarr."""

from pathlib import Path
from unittest import mock

import anndata as ad
import numpy as np
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


def test_zarr_backed_warns_and_loads_all(zarr_store, capsys):
    # backed mode has no zarr equivalent: warn, and serve *all* cells (ignore max_cells).
    adata = loader.load_adata(str(zarr_store), max_cells=10, backed=True)
    assert adata.n_obs == 40
    assert "backed mode is not supported for .zarr" in capsys.readouterr().out


def test_local_zarr_uses_read_lazy(zarr_store):
    """The AnnData zarr path goes through anndata.experimental.read_lazy."""
    import anndata.experimental as aexp

    real = aexp.read_lazy
    calls = {}

    def spy(store, *a, **k):
        calls["store"] = str(store)
        return real(store, *a, **k)

    with mock.patch.object(aexp, "read_lazy", side_effect=spy):
        adata = loader.load_adata(str(zarr_store), max_cells=None)

    assert calls["store"] == str(zarr_store)
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


def test_remote_non_zarr_rejected():
    with pytest.raises(ValueError, match="Remote sc_data must be a .zarr store"):
        loader.load_adata("s3://bucket/path/dataset.h5ad")


# ---------------------------------------------------------------------------
# Cloud .zarr (mocked -- no real network/fsspec required)
# ---------------------------------------------------------------------------

def test_remote_zarr_not_path_mangled():
    """A cloud URI is passed to read_lazy verbatim, never collapsed to a Path."""
    import anndata.experimental as aexp

    captured = {}

    def fake_read_lazy(store, *a, **k):
        captured["store"] = store
        return _toy_adata()

    with mock.patch.object(loader, "_zarr_is_mudata", return_value=False), mock.patch.object(
        aexp, "read_lazy", side_effect=fake_read_lazy
    ):
        adata = loader.load_adata("s3://my-bucket/path/my_dataset.zarr", max_cells=None)

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
