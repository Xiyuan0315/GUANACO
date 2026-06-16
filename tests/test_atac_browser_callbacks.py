from guanaco.pages.matrix.callbacks.atac_browser_callbacks import (
    _gene_track_label,
    _gene_model_unavailable,
    _range_from_relayout,
)


def test_gene_track_label_defaults_to_genes_when_no_path():
    assert _gene_track_label(None) == "Genes"
    assert _gene_track_label("") == "Genes"


def test_gene_track_label_strips_directory_and_known_extensions():
    assert _gene_track_label("/data/refs/gencode.v44.gtf") == "gencode.v44"
    assert _gene_track_label("/data/refs/gencode.v44.gff3.gz") == "gencode.v44"
    assert _gene_track_label("ann.gff") == "ann"


def test_range_from_relayout_reads_separate_range_keys():
    relayout = {"xaxis.range[0]": "100.0", "xaxis.range[1]": "250.9"}
    # Coordinates are truncated to int basepairs.
    assert _range_from_relayout(relayout) == (100, 250)


def test_range_from_relayout_reads_list_range_key():
    assert _range_from_relayout({"xaxis.range": [10.0, 80.0]}) == (10, 80)


def test_range_from_relayout_none_for_irrelevant_or_empty_events():
    assert _range_from_relayout(None) is None
    assert _range_from_relayout({}) is None
    assert _range_from_relayout({"yaxis.range[0]": 1}) is None


def test_gene_model_unavailable_formats_error():
    assert _gene_model_unavailable("bad file") == "Gene model unavailable: bad file"
