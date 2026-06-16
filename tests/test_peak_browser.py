import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from guanaco import widget as pl
from guanaco.pages.matrix.plots import atac_browser


@pytest.fixture(autouse=True)
def _clear_atac_caches():
    # atac_browser caches the peak index by id(adata); across tests Python can
    # reuse an object id, so clear the module caches to keep each test isolated.
    atac_browser._peak_index_cache.clear()
    atac_browser._signal_cache._store.clear()
    yield


def _peak_adata(n_obs=40, n_peaks=30):
    rng = np.random.default_rng(0)
    x = rng.random((n_obs, n_peaks), dtype=np.float64).astype(np.float32)
    var_names = [f"chr1:{1000 + i * 1000}-{1500 + i * 1000}" for i in range(n_peaks)]
    obs = pd.DataFrame(
        {"CellType": np.array(["A", "B"] * (n_obs // 2))[:n_obs]},
        index=[f"c{i}" for i in range(n_obs)],
    )
    return AnnData(X=x, obs=obs, var=pd.DataFrame(index=var_names))


def _bar_traces(fig):
    return [t for t in fig.data if t.type == "bar"]


def test_peak_browser_draws_one_track_per_group():
    fig = pl.peak_browser(_peak_adata(), groupby="CellType", return_fig=True)
    assert len(_bar_traces(fig)) == 2  # one accessibility track per cell type


def test_peak_browser_accepts_a_locus_string():
    fig = pl.peak_browser(
        _peak_adata(), region="chr1:1,000-10,000", groupby="CellType", return_fig=True
    )
    assert len(_bar_traces(fig)) >= 1


def test_peak_browser_detection_metric_builds():
    fig = pl.peak_browser(_peak_adata(), groupby="CellType", metric="detection", return_fig=True)
    assert len(_bar_traces(fig)) >= 1


def test_default_y_mode_is_shared():
    import inspect

    assert inspect.signature(pl.peak_browser).parameters["y_mode"].default == "shared"


def test_peak_browser_y_mode_auto_builds():
    fig = pl.peak_browser(_peak_adata(), groupby="CellType", y_mode="auto", return_fig=True)
    assert len(_bar_traces(fig)) == 2


def test_peak_browser_labels_subset_groups():
    fig = pl.peak_browser(
        _peak_adata(), groupby="CellType", labels=["A"], return_fig=True
    )
    assert len(_bar_traces(fig)) == 1


def test_peak_browser_without_peaks_raises():
    adata = AnnData(
        X=np.zeros((4, 3), dtype="float32"),
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC"]),
    )
    with pytest.raises(ValueError):
        pl.peak_browser(adata)


def test_peak_browser_rejects_unparseable_region_without_annotation():
    with pytest.raises(ValueError, match="gene name"):
        pl.peak_browser(_peak_adata(), region="not-a-locus")


def _mini_gtf(tmp_path):
    # Place the gene inside the fixture's peak coordinate range (1 kb – ~30 kb); the
    # browser clamps navigation to where peaks actually exist, so a gene far outside
    # the data's peaks would be pulled back and not appear.
    p = tmp_path / "mini.gtf"
    p.write_text(
        'chr1\tHAVANA\tgene\t5000\t25000\t.\t+\t.\tgene_id "G1"; gene_name "MYGENE"; gene_type "protein_coding";\n'
        'chr1\tHAVANA\ttranscript\t5000\t25000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "MYGENE"; transcript_type "protein_coding";\n'
        'chr1\tHAVANA\texon\t5000\t8000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    )
    return str(p)


def test_peak_browser_searches_by_gene_name(tmp_path):
    fig = pl.peak_browser(
        _peak_adata(), region="MYGENE", groupby="CellType",
        gene_annotation=_mini_gtf(tmp_path), return_fig=True,
    )
    # Navigated to the gene and drew a gene-model track.
    assert any(t.name == "Transcripts" for t in fig.data)


def test_gene_name_without_annotation_raises():
    with pytest.raises(ValueError, match="gene name"):
        pl.peak_browser(_peak_adata(), region="MYGENE")


def test_unknown_gene_name_raises(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        pl.peak_browser(_peak_adata(), region="NOSUCHGENE", gene_annotation=_mini_gtf(tmp_path))


def test_peak_index_cache_evicts_when_adata_is_collected():
    import gc

    adata = _peak_adata()
    atac_browser.build_peak_index(adata)
    assert len(atac_browser._peak_index_cache) == 1

    del adata
    gc.collect()
    # weakref.finalize must drop the entry so a recycled id() can't read it back.
    assert len(atac_browser._peak_index_cache) == 0
