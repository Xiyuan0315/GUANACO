import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from guanaco.pages.matrix.plots.atac_browser import (
    build_peak_index,
    compute_atac_signal,
    has_genomic_peak_features,
    parse_locus,
    plot_atac_browser,
)
from guanaco.pages.matrix.plots import gene_annotation as gene_annotation_module
from guanaco.pages.matrix.plots.gene_annotation import (
    find_gene_region,
    load_gene_annotation,
    query_gene_models,
    resolve_annotation_source,
)


def _small_atac():
    x = sparse.csr_matrix(
        np.array(
            [
                [1, 0, 2, 0],
                [0, 1, 0, 3],
                [2, 0, 0, 1],
            ],
            dtype=np.float32,
        )
    )
    return AnnData(
        x,
        obs=pd.DataFrame({"cell_type": ["T", "B", "T"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["chr1:100-200", "chr1:300-400", "chr2:100-200", "not_a_peak"]),
    )


def test_parse_locus():
    assert parse_locus("chr1:1,000-2,500") == ("chr1", 1000, 2500)
    assert parse_locus("chr1:2500-1000") == ("chr1", 1000, 2500)
    assert parse_locus("TP53") is None


def test_peak_index_from_var_names():
    adata = _small_atac()
    assert has_genomic_peak_features(adata, min_hits=2)
    index = build_peak_index(adata)
    assert index.chroms == ("chr1", "chr2")
    assert index.var_indices["chr1"].tolist() == [0, 1]


def test_compute_signal_respects_selected_cells_and_groupby():
    adata = _small_atac()
    payload = compute_atac_signal(
        adata,
        {"chrom": "chr1", "start": 0, "end": 500},
        selected_cells=["c1", "c3"],
        groupby="cell_type",
        metric="mean",
    )
    assert payload["n_cells"] == 2
    assert len(payload["signals"]) == 1
    assert payload["signals"][0]["name"] == "T"
    assert np.allclose(payload["signals"][0]["values"], [1.5, 0.0])


def test_signal_keeps_tracks_in_peak_free_window():
    # Peaks clustered at the two ends of chr1 leave 400k-500k as an interior gap.
    x = sparse.csr_matrix(
        np.array([[1, 2, 0, 0], [0, 0, 1, 3], [2, 1, 0, 1]], dtype=np.float32)
    )
    adata = AnnData(
        x,
        obs=pd.DataFrame({"cell_type": ["T", "B", "T"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame(
            index=["chr1:1000-2000", "chr1:2000-3000", "chr1:900000-901000", "chr1:901000-902000"]
        ),
    )
    payload = compute_atac_signal(
        adata,
        {"chrom": "chr1", "start": 400000, "end": 500000},
        groupby="cell_type",
        labels=["T", "B"],
        group_order=["B", "T"],
    )
    # No peaks in the window, but every selected track must still be emitted (empty)
    # so panning into a peak-free stretch doesn't collapse the figure to gene-only.
    assert payload["peaks"]["total"] == 0
    assert [s["name"] for s in payload["signals"]] == ["B", "T"]
    assert all(len(s["values"]) == 0 for s in payload["signals"])


def test_plot_atac_browser_returns_figure():
    adata = _small_atac()
    payload = compute_atac_signal(
        adata,
        {"chrom": "chr1", "start": 0, "end": 500},
        metric="detection",
    )
    fig = plot_atac_browser(payload)
    assert fig.layout.dragmode == "pan"
    assert len(fig.data) >= 1
    assert fig.layout.showlegend is False
    assert fig.layout.title.text is None


def test_resolve_annotation_source(monkeypatch, tmp_path):
    # None and plain local paths pass through untouched (no download attempted).
    assert resolve_annotation_source(None) is None
    assert resolve_annotation_source("/data/foo.gtf.gz") == "/data/foo.gtf.gz"

    # Genome ids (and their assembly aliases) and URLs route through the cached
    # downloader -- stubbed here so the test never touches the network.
    fetched = []

    def fake_download(url):
        fetched.append(url)
        return tmp_path / url.split("/")[-1]

    monkeypatch.setattr(gene_annotation_module, "_download_to_cache", fake_download)

    assert str(resolve_annotation_source("hg38")).endswith("gencode.v50.basic.annotation.gtf.gz")
    assert str(resolve_annotation_source("GRCh38")).endswith("gencode.v50.basic.annotation.gtf.gz")
    assert str(resolve_annotation_source("mm10")).endswith("gencode.vM25.basic.annotation.gtf.gz")
    assert str(resolve_annotation_source("https://example.org/custom.gtf.gz")).endswith("custom.gtf.gz")
    assert len(fetched) == 4


def test_gene_annotation_models_from_gtf(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\ttest\tgene\t80\t420\t.\t+\t.\tgene_id "GENE1"; gene_name "ABC1"; gene_type "protein_coding";',
                'chr1\ttest\ttranscript\t90\t400\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; transcript_name "ABC1-201"; transcript_type "protein_coding"; tag "basic";',
                'chr1\ttest\texon\t90\t150\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";',
                'chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";',
            ]
        )
    )
    index = load_gene_annotation(str(gtf))
    models = query_gene_models(index, "chr1", 0, 500)
    assert len(models) == 1
    assert models[0]["gene_name"] == "ABC1"
    assert models[0]["transcripts"][0]["transcript_id"] == "TX1"
    assert models[0]["transcripts"][0]["exons"] == [(90, 150), (300, 400)]
    # Framed by the transcript extent (TX1: 90-400), not the wider gene record
    # (80-420) -- so the search isn't thrown off by isoforms outside the basic set.
    assert find_gene_region(index, "ABC1", flank=10) == {"chrom": "chr1", "start": 80, "end": 410}


def test_plot_atac_browser_can_draw_gene_models(tmp_path):
    adata = _small_atac()
    payload = compute_atac_signal(
        adata,
        {"chrom": "chr1", "start": 0, "end": 500},
        metric="mean",
    )
    gene_models = [
        {
            "gene_id": "GENE1",
            "gene_name": "ABC1",
            "gene_type": "protein_coding",
            "chrom": "chr1",
            "start": 80,
            "end": 420,
            "strand": "+",
            "transcripts": [
                {
                    "transcript_id": "TX1",
                    "transcript_name": "ABC1-201",
                    "transcript_type": "protein_coding",
                    "start": 90,
                    "end": 400,
                    "strand": "+",
                    "exons": [(90, 150), (300, 400)],
                }
            ],
        }
    ]
    fig = plot_atac_browser(payload, gene_models=gene_models)
    assert len(fig.data) >= 4
    assert any(annotation.text == "Genes" for annotation in fig.layout.annotations)
    assert any(annotation.text == "Genes" and annotation.bgcolor for annotation in fig.layout.annotations)
    assert fig.layout.yaxis.title.text is None
