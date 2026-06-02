import guanaco
import pandas as pd
from guanaco import cli
from guanaco.pages.matrix.plots.volcano import (
    deg_csv,
    deg_rows,
    load_volcano_payload,
    plot_volcano,
    volcano_degs_filename,
    volcano_entry_options,
)
from guanaco.pages.matrix.plots.grn_demo import (
    default_grn_edge_threshold,
    grn_context_options,
    grn_dataframe,
    grn_has_weight,
    grn_weight_range,
    grn_weight_step,
    has_grn_data,
)
from guanaco.utils.volcano_utils import save_pydeseq_results_to_adata_uns
from guanaco.utils.colors import (
    continuous_colormap_options,
    continuous_colormap_overview,
    discrete_palette_names,
    plot_continuous_colormap_overview,
    discrete_palette_overview,
    plot_discrete_palette_overview,
    resolve_discrete_palette,
)


def test_version_attribute_exists():
    assert isinstance(guanaco.__version__, str)
    assert guanaco.__version__


def test_cli_entrypoint_is_callable():
    assert callable(cli.main)


def test_continuous_colormap_overview_shows_original_colormaps_only():
    overview = continuous_colormap_overview()

    assert set(overview) == {"linear", "diverging"}
    assert any(entry["label"] == "plotly/linear_viridis" for entry in overview["linear"])
    assert any(entry["label"] == "plotly/diverging_rdbu" for entry in overview["diverging"])
    assert not any(
        entry["name"].endswith("_r")
        for entries in overview.values()
        for entry in entries
    )
    assert all(
        entry["label"].split("/", 1)[1].startswith(f"{group_name}_")
        for group_name in overview
        for entry in overview[group_name]
    )


def test_continuous_colormap_options_include_reversed_variants():
    options = continuous_colormap_options()
    labels = {option["label"] for option in options}
    values = {option["value"] for option in options}

    assert "plotly/linear_viridis" in labels
    assert "plotly/linear_viridis_r" in labels
    assert "cmc/diverging_vik" in labels
    assert "cmc/diverging_vik_r" in labels
    assert "viridis_r" in values
    assert "cmc:vik_r" in values
    assert all("/" in option["label"] for option in options)


def test_discrete_palette_names_use_source_prefix_for_generated_palettes():
    names = set(discrete_palette_names())

    assert "okabe_ito" in names
    assert "tol_vibrant" in names
    assert "polychrome_colorsafe" in names
    assert "krzywinski_24" in names
    assert "krzywinski_24_main" not in names
    assert "ibm" in names
    assert "Okabe_Ito" not in names
    assert "glasbey/default" in names
    assert "glasbey/safe" in names
    assert "Plotly" not in names
    assert "plotly/Plotly" in names
    assert not any(name.startswith("plotly/") and name.endswith("_r") for name in names)
    assert "cc/glasbey" in names
    assert "cmc/batlowS" in names
    assert "cmc/grayCS" in names
    assert not any(name.startswith("colorcet/") for name in names)


def test_crameri_categorical_palette_resolves_as_discrete_colors():
    palette = resolve_discrete_palette("cmc/batlowS")

    assert len(palette) == 100
    assert all(color.startswith("rgb(") for color in palette)


def test_dynamic_glasbey_palette_generates_requested_size():
    default_palette = resolve_discrete_palette("glasbey/default", 37)
    safe_palette = resolve_discrete_palette("glasbey/safe", 37)

    assert len(default_palette) == 37
    assert len(safe_palette) == 37
    assert default_palette != safe_palette
    assert all(color.startswith("#") for color in default_palette + safe_palette)


def test_discrete_palette_overview_groups_palettes_by_source():
    overview = discrete_palette_overview(dynamic_n_colors=12)

    assert {source: len(entries) for source, entries in overview.items()} == {
        "custom": 9,
        "plotly": 19,
        "cc": 13,
        "cmc": 16,
        "glasbey": 2,
    }
    assert any(entry["label"] == "okabe_ito" for entry in overview["custom"])
    assert any(entry["label"] == "plotly/Plotly" for entry in overview["plotly"])
    assert any(entry["label"] == "glasbey/safe" and entry["n_colors"] == 12 for entry in overview["glasbey"])


def test_plot_discrete_palette_overview_returns_figure():
    fig = plot_discrete_palette_overview(source="custom")

    try:
        assert len(fig.axes) == 1
        assert fig.axes[0].get_title() == "custom (9)"
    finally:
        fig.clear()


class _FakeAnnData:
    def __init__(self, uns):
        self.uns = uns


def _fake_volcano_adata():
    return _FakeAnnData(
        {
            "volcano": {
                "entries": {
                    "treated_vs_control": {
                        "gene": ["GeneA", "GeneB", "GeneC"],
                        "logfoldchange": [2.0, -1.5, 0.2],
                        "score": [5.0, -4.0, 0.5],
                        "padj": [0.01, 0.04, 0.8],
                        "symbol": ["A", "B", "C"],
                    }
                },
                "default_entry": "treated_vs_control",
            }
        }
    )


def test_volcano_payload_options_and_deg_export():
    adata = _fake_volcano_adata()
    options, default_entry = volcano_entry_options(adata)
    payload = load_volcano_payload(adata)
    entry = payload["entries"]["treated_vs_control"]

    assert options == [{"label": "treated_vs_control", "value": "treated_vs_control"}]
    assert default_entry == "treated_vs_control"
    assert payload["source_type"] == "volcano"
    assert len(deg_rows(entry, "logfoldchange", 0.05, 1.0)) == 2
    assert "symbol" in deg_csv(entry, "logfoldchange", 0.05, 1.0).splitlines()[0]
    assert (
        volcano_degs_filename("treated_vs_control", entry, "logfoldchange", 0.05, 0.5)
        == "treated_vs_control_padj005_log2fc05_deg.csv"
    )


def test_plot_volcano_returns_figure():
    payload = load_volcano_payload(_fake_volcano_adata())
    entry = payload["entries"]["treated_vs_control"]
    fig = plot_volcano("treated_vs_control", entry, top_n=2)

    assert len(fig.data) == 4
    assert fig.layout.title.text == "Volcano Plot: treated_vs_control"
    assert fig.layout.xaxis.title.text == "log fold change"


def test_save_pydeseq_dataframe_to_adata_uns_for_volcano():
    adata = _FakeAnnData({})
    results = pd.DataFrame(
        {
            "baseMean": [3083.4, 131.8, 6307318.0],
            "log2FoldChange": [-0.05, -0.88, 0.045],
            "lfcSE": [0.16, 0.46, 0.05],
            "stat": [-0.32, -1.92, 0.86],
            "pvalue": [0.75, 0.054, 0.38],
            "padj": [0.91, 0.34, 0.72],
            "gene_name": ["TSPAN6", "TNMD", None],
        },
        index=pd.Index(["ENSG00000000003.15", "ENSG00000000005.6", "N_ambiguous"], name="Ensembl_gene_id_v"),
    )

    save_pydeseq_results_to_adata_uns(adata, results, "treated_vs_control")
    payload = load_volcano_payload(adata)
    entry = payload["entries"]["treated_vs_control"]

    assert payload["default_entry"] == "treated_vs_control"
    assert entry["gene"] == ["TSPAN6", "TNMD", "N_ambiguous"]
    assert entry["logfoldchange"] == [-0.05, -0.88, 0.045]
    assert entry["padj"] == [0.91, 0.34, 0.72]
    assert entry["extra_columns"]["gene_id"] == ["ENSG00000000003.15", "ENSG00000000005.6", "N_ambiguous"]
    assert entry["extra_columns"]["baseMean"] == [3083.4, 131.8, 6307318.0]
    assert entry["meta"]["pval_key"] == "padj"


def test_save_pydeseq_csv_to_adata_uns_infers_gene_id_column(tmp_path):
    adata = _FakeAnnData({})
    csv_path = tmp_path / "deseq_results.csv"
    pd.DataFrame(
        {
            "Ensembl_gene_id_v": ["ENSG1", "ENSG2"],
            "baseMean": [10.0, 20.0],
            "log2FoldChange": [1.2, -0.7],
            "stat": [3.1, -2.2],
            "pvalue": [0.001, 0.02],
            "padj": [0.01, 0.04],
            "gene_name": ["GeneA", "GeneB"],
        }
    ).to_csv(csv_path, index=False)

    save_pydeseq_results_to_adata_uns(adata, csv_path, "csv_contrast")
    payload = load_volcano_payload(adata)
    entry = payload["entries"]["csv_contrast"]

    assert entry["gene"] == ["GeneA", "GeneB"]
    assert entry["extra_columns"]["gene_id"] == ["ENSG1", "ENSG2"]
    assert "Ensembl_gene_id_v" not in entry["extra_columns"]


def test_grn_dataframe_contract_supports_context_alias_and_weight():
    adata = _FakeAnnData(
        {
            "grn": pd.DataFrame(
                {
                    "source": ["STAT1", "NFKB1"],
                    "target": ["IRF1", "TNF"],
                    "regulation": ["+", "-"],
                    "weight": [0.92, 0.84],
                    "cell_type": ["Monocyte", "Macrophage"],
                }
            )
        }
    )

    options, context_column, has_weight = grn_context_options(adata)

    assert has_grn_data(adata)
    assert grn_dataframe(adata).shape == (2, 5)
    assert context_column == "cell_type"
    assert has_weight
    assert grn_has_weight(grn_dataframe(adata))
    assert grn_weight_range(grn_dataframe(adata)) == (0.84, 0.92)
    assert default_grn_edge_threshold(grn_dataframe(adata)) == 0.84
    assert grn_weight_step(grn_dataframe(adata)) == 0.05
    assert options == [
        {"label": "All", "value": "All"},
        {"label": "Macrophage", "value": "Macrophage"},
        {"label": "Monocyte", "value": "Monocyte"},
    ]


def test_grn_weight_threshold_supports_importance_scale():
    adata = _FakeAnnData(
        {
            "grn": pd.DataFrame(
                {
                    "source": ["P4hb", "Klf1"],
                    "target": ["Mpo", "Car2"],
                    "regulation": ["+", "+"],
                    "weight": [14.179311, 3.230082],
                }
            )
        }
    )
    grn_df = grn_dataframe(adata)

    assert grn_weight_range(grn_df) == (3.230082, 14.179311)
    assert default_grn_edge_threshold(grn_df) == 3.230082
    assert grn_weight_step(grn_df) == 1.0


def test_grn_dataframe_requires_source_target_regulation():
    adata = _FakeAnnData({"grn": pd.DataFrame({"source": ["STAT1"], "target": ["IRF1"]})})

    assert not has_grn_data(adata)

    try:
        grn_dataframe(adata)
    except ValueError as exc:
        assert "regulation" in str(exc)
    else:
        raise AssertionError("Expected missing regulation column to fail")


def test_plot_continuous_colormap_overview_returns_figure():
    fig = plot_continuous_colormap_overview(colormap_type="linear", source="plotly")

    try:
        assert len(fig.axes) == 1
        assert fig.axes[0].get_title() == "Linear colormaps (15)"
    finally:
        fig.clear()
