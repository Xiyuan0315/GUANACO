import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from dash import Dash, no_update

from guanaco.data.network_sources import (
    NetworkEdge,
    fetch_network_edges,
    normalize_organism,
)
from guanaco.pages.matrix.analysis.network import (
    build_network_payload,
    parse_gene_list,
    trrust_key_regulator_enrichment,
)
from guanaco.pages.matrix.callbacks.network_callbacks import (
    _build_view,
    _network_status,
    register_network_callbacks,
)
from guanaco.pages.matrix.plots.network import (
    _edge_width,
    _elements,
    network_layout,
    network_stylesheet,
)
from guanaco.pages.visualizations.layout import generate_visualization_sections
from guanaco.pages.visualizations.registry import resolve_plot_components


def _adata():
    return AnnData(
        X=np.asarray(
            [
                [0.0, 1.0, 2.0],
                [2.0, 3.0, 4.0],
                [4.0, 5.0, 6.0],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["TP53", "EGFR", "CDK2"]),
    )


def _walk(component):
    if component is None:
        return
    if isinstance(component, (list, tuple)):
        for child in component:
            yield from _walk(child)
        return
    yield component
    children = getattr(component, "children", None)
    if children is not None:
        yield from _walk(children)


def _by_id(component, component_id):
    return next(
        item for item in _walk(component) if getattr(item, "id", None) == component_id
    )


def test_parse_gene_list_accepts_pasted_signatures_and_deduplicates():
    assert parse_gene_list("TP53, EGFR\nCDK2;tp53") == ["TP53", "EGFR", "CDK2"]


def test_network_is_available_only_for_rna_modalities():
    source = _adata()
    assert resolve_plot_components(
        source,
        ["network"],
        modality_name="rna",
    ) == ("network",)
    assert (
        resolve_plot_components(
            source,
            ["network"],
            modality_name="atac",
        )
        == ()
    )
    assert resolve_plot_components(source, ["grn"], modality_name="rna") == ()


def test_network_layout_is_a_self_contained_exploratory_tab():
    sections = generate_visualization_sections(
        _adata(),
        ["TP53"],
        [],
        "p",
        optional_plot_components=["network"],
        modality_name="rna",
    )
    tabs = _by_id(sections, "p-exploratory-tabs")
    assert [tab.value for tab in tabs.children] == ["network-tab"]
    assert [tab.label for tab in tabs.children] == ["Network"]
    assert _by_id(sections, "p-network-genes") is not None
    assert _by_id(sections, "p-network-type").value == "ppi"
    assert _by_id(sections, "p-network-view").value == "first-order"
    string_mode = _by_id(sections, "p-network-string-physical")
    assert string_mode.value == []
    assert string_mode.options == [{"label": "Physical", "value": "physical"}]
    assert (
        getattr(_by_id(sections, "p-network-string-type-wrap"), "style", None) is None
    )
    tf_direction = _by_id(sections, "p-network-tf-direction")
    assert tf_direction.value == "targets"
    assert [option["value"] for option in tf_direction.options] == [
        "targets",
        "regulators",
        "both",
    ]
    assert _by_id(sections, "p-network-tf-direction-wrap").style == {"display": "none"}
    assert not any(
        getattr(item, "id", None) == "p-network-max-nodes" for item in _walk(sections)
    )
    confidence = _by_id(sections, "p-network-string-confidence")
    assert confidence.value == 0.7
    assert confidence.min == 0.4
    assert confidence.max == 0.9
    assert confidence.updatemode == "mouseup"
    assert confidence.className == "dbc-slider"
    assert (
        getattr(_by_id(sections, "p-network-string-confidence-wrap"), "style", None)
        is None
    )
    assert not any(
        getattr(item, "id", None) == "p-network-mirna-regulators"
        for item in _walk(sections)
    )
    layout = _by_id(sections, "p-network-layout")
    assert layout.value == "cose"
    assert [option["value"] for option in layout.options] == [
        "cose",
        "breadthfirst",
        "circle",
        "concentric",
    ]
    assert _by_id(sections, "p-network-build") is not None
    primary_controls = _by_id(sections, "p-network-primary-controls")
    primary_ids = {
        getattr(component, "id", None) for component in _walk(primary_controls)
    }
    assert {
        "p-network-genes",
        "p-network-type",
        "p-network-string-physical",
        "p-network-tf-direction",
        "p-network-build",
    } <= primary_ids
    assert "p-network-view" not in primary_ids
    secondary_controls = _by_id(sections, "p-network-secondary-controls")
    secondary_ids = {
        getattr(component, "id", None) for component in _walk(secondary_controls)
    }
    assert {
        "p-network-view",
        "p-network-layout",
        "p-network-string-confidence",
    } <= secondary_ids
    assert "p-network-string-physical" not in secondary_ids
    assert "p-network-tf-direction" not in secondary_ids
    assert "p-network-download-svg" not in secondary_ids
    assert _by_id(sections, "p-network-options-toggle").children == "▸ More options"
    assert _by_id(sections, "p-network-options-collapse").is_open is False
    graph_shell = _by_id(sections, "p-network-graph-shell")
    assert _by_id(graph_shell, "p-network-download-svg") is not None
    assert _by_id(sections, "p-network-loading").target_components == {
        "p-network-graph": "children"
    }
    regulator_section = _by_id(sections, "p-network-regulator-section")
    assert regulator_section.style == {"display": "none"}
    assert _by_id(sections, "p-network-regulator-results") is not None
    assert not any(
        getattr(item, "id", None) == "p-single-cell-genes-selection"
        for item in _walk(tabs)
    )


def test_string_response_is_normalized_and_supports_functional_or_physical_networks():
    captured = {}

    def request_json(url, params, **kwargs):
        captured.update(url=url, params=params, kwargs=kwargs)
        return [
            {
                "preferredName_A": "TP53",
                "preferredName_B": "MDM2",
                "score": 0.98,
            }
        ]

    edges = fetch_network_edges(
        "ppi",
        ["TP53"],
        "human",
        request_json=request_json,
        use_cache=False,
    )
    assert edges == [
        NetworkEdge(
            source="TP53",
            target="MDM2",
            score=0.98,
            database="STRING v12 functional",
        )
    ]
    assert captured["params"]["network_type"] == "functional"
    assert captured["params"]["required_score"] == 400
    assert captured["params"]["show_query_node_labels"] == 1
    assert captured["params"]["species"] == 9606
    assert captured["kwargs"]["method"] == "POST"

    physical = fetch_network_edges(
        "ppi",
        ["TP53"],
        "human",
        request_json=request_json,
        ppi_network_type="physical",
        use_cache=False,
    )
    assert physical[0].database == "STRING v12 physical"
    assert captured["params"]["network_type"] == "physical"


def test_trrust_response_is_filtered_in_both_directions_and_normalized():
    captured = []

    def request_text(url):
        captured.append(url)
        return "\n".join(
            [
                "STAT1\tIRF1\tActivation\t12345678",
                "JAK1\tSTAT1\tRepression\t23456789",
                "TP53\tCDKN1A\tUnknown\t34567890",
            ]
        )

    edges = fetch_network_edges(
        "tf-gene",
        ["STAT1"],
        "human",
        request_text=request_text,
        use_cache=False,
    )
    assert captured == ["https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"]
    assert edges == [
        NetworkEdge(
            "STAT1",
            "IRF1",
            directed=True,
            effect="activation",
            database="TRRUST v2",
            references="12345678",
        ),
        NetworkEdge(
            "JAK1",
            "STAT1",
            directed=True,
            effect="inhibition",
            database="TRRUST v2",
            references="23456789",
        ),
    ]


def test_trrust_key_regulator_enrichment_matches_published_example():
    query = ["CD74", "CD52", "HLA-DRA", "HLA-DRB1", "CIITA"]
    target_sets = {
        "ILF3": (["HLA-DRA", "HLA-DRB1"], 7),
        "RFX1": (["HLA-DRA", "CIITA"], 15),
        "RFXANK": (["HLA-DRA", "HLA-DRB1"], 15),
        "RFXAP": (["HLA-DRA", "HLA-DRB1"], 15),
        "RFX5": (["HLA-DRA", "HLA-DRB1"], 18),
        "CIITA": (["HLA-DRA", "HLA-DRB1"], 31),
    }
    used = set(query)
    edges = []
    for tf, (overlap, target_count) in target_sets.items():
        targets = list(overlap)
        while len(targets) < target_count:
            target = f"BACKGROUND_{tf}_{len(targets)}"
            targets.append(target)
            used.add(target)
        edges.extend(NetworkEdge(tf, target, directed=True) for target in targets)

    background = list(used)
    background.extend(f"UNIVERSE_{index}" for index in range(18_862 - len(background)))
    adata = AnnData(
        X=np.zeros((1, len(background)), dtype=np.float32),
        var=pd.DataFrame(index=background),
    )
    enrichment = trrust_key_regulator_enrichment(adata, query, edges)
    results = {row["tf"]: row for row in enrichment["results"]}

    assert enrichment["background_size"] == 18_862
    assert enrichment["valid_query_count"] == 5
    assert len(results) == 6
    assert results["ILF3"]["p_value"] == pytest.approx(1.1799579e-6)
    assert results["ILF3"]["fdr"] == pytest.approx(7.0797474e-6)


def test_kegg_metabolite_network_uses_organism_reaction_chain():
    responses = {
        "https://rest.kegg.jp/list/hsa": (
            "hsa:3098\tCDS\t10:69270000..69401882\t"
            "HK1, CNSHA5, HK; hexokinase-1 isoform HKI"
        ),
        "https://rest.kegg.jp/link/ec/hsa:3098": "hsa:3098\tec:2.7.1.1",
        "https://rest.kegg.jp/link/reaction/ec:2.7.1.1": ("ec:2.7.1.1\trn:R01786"),
        "https://rest.kegg.jp/link/compound/rn:R01786": (
            "rn:R01786\tcpd:C00002\nrn:R01786\tcpd:C00267"
        ),
        "https://rest.kegg.jp/list/cpd:C00002+cpd:C00267": (
            "C00002\tATP; Adenosine 5'-triphosphate\nC00267\talpha-D-Glucose"
        ),
    }
    requested = []

    def request_text(url):
        requested.append(url)
        return responses[url]

    edges = fetch_network_edges(
        "metabolite",
        ["HK1"],
        "human",
        request_text=request_text,
        request_json=lambda *_a, **_k: pytest.fail("OmniPath should not be called"),
        use_cache=False,
    )

    assert requested == list(responses)
    assert edges == [
        NetworkEdge(
            "ATP (C00002)",
            "HK1",
            source_type="metabolite",
            directed=False,
            database="KEGG reactions",
            references="R01786",
        ),
        NetworkEdge(
            "alpha-D-Glucose (C00267)",
            "HK1",
            source_type="metabolite",
            directed=False,
            database="KEGG reactions",
            references="R01786",
        ),
    ]


def test_mirtarbase_query_merges_evidence_and_filters_exact_target():
    rows = "\n".join(
        [
            "MTI ID,Species (miRNA),Species (Target),miRNA,Target Gene,Reporter Assay,"
            "Western Blot,qPCR,Microarray,NGS,pSILAC,CLIP-Seq,CLASH,Other,Sum,Papers,Score",
            "MIRT1,Homo sapiens,Homo sapiens,hsa-miR-21-5p,PTEN,0,0,0,0,1,0,0,1,0,2,1,Tier D: Low",
            "MIRT2,Homo sapiens,Homo sapiens,hsa-miR-21-5p,PTEN,1,1,1,0,0,0,0,0,0,3,2,Tier B: High",
            "MIRT3,Homo sapiens,Homo sapiens,hsa-miR-1-3p,SPTEN,0,0,0,0,1,0,0,0,0,1,1,Tier D: Low",
            "MIRT4,Homo sapiens,Mus musculus,hsa-miR-1-3p,Pten,1,0,0,0,0,0,0,0,0,1,1,Tier A: Very High",
        ]
    )
    requested = []

    def request_text(url):
        requested.append(url)
        return rows

    edges = fetch_network_edges(
        "mirna",
        ["PTEN"],
        "human",
        mirtarbase_request=request_text,
        request_json=lambda *_a, **_k: pytest.fail("OmniPath should not be called"),
        use_cache=False,
    )
    assert requested == [
        "https://awi.cuhk.edu.cn/miRTarBase/search/results/download/"
        "?mode=target&keyword=PTEN&species=hsa"
    ]
    assert edges == [
        NetworkEdge(
            "hsa-miR-21-5p",
            "PTEN",
            source_type="mirna",
            directed=True,
            score=3.0,
            database="miRTarBase",
            references="MIRT1, MIRT2",
        )
    ]


def test_mirtarbase_uses_species_specific_mouse_file():
    requested = []
    header = (
        "MTI ID,Species (miRNA),Species (Target),miRNA,Target Gene,Reporter Assay,"
        "Western Blot,qPCR,Microarray,NGS,pSILAC,CLIP-Seq,CLASH,Other,Sum,Papers,Score\n"
    )

    def request_text(url):
        requested.append(url)
        return (
            header + "MIRT1,Mus musculus,Mus musculus,mmu-miR-21a-5p,Pten,1,0,0,0,"
            "0,0,0,0,0,1,1,Tier C: Middle"
        )

    edges = fetch_network_edges(
        "mirna",
        ["Pten"],
        "mouse",
        mirtarbase_request=request_text,
        use_cache=False,
    )
    assert requested[0].endswith("?mode=target&keyword=Pten&species=mmu")
    assert edges[0].source == "mmu-miR-21a-5p"


def test_trrust_network_rejects_rat():
    with pytest.raises(ValueError, match="human and mouse"):
        fetch_network_edges(
            "tf-gene", ["Tp53"], "rat", request_text=lambda *_a, **_k: ""
        )


def test_network_payload_colors_only_measured_gene_nodes():
    graph = build_network_payload(
        _adata(),
        ["TP53", "EGFR"],
        "tf-gene",
        [
            NetworkEdge("TP53", "EGFR", directed=True, effect="activation", score=3),
            NetworkEdge("MIR21", "TP53", source_type="mirna", directed=True, score=2),
            NetworkEdge("MDM2", "TP53", directed=True, score=1),
        ],
    )
    nodes = {node["id"]: node for node in graph["nodes"]}
    assert nodes["TP53"]["measured"] is True
    assert nodes["TP53"]["expression"] == pytest.approx(2.0)
    assert nodes["EGFR"]["expression"] == pytest.approx(3.0)
    assert nodes["MIR21"]["entity_type"] == "mirna"
    assert nodes["MIR21"]["measured"] is False
    assert nodes["MDM2"]["measured"] is False
    assert "MDM2" in graph["missing_expression"]
    assert nodes["TP53"]["entity_type"] == "tf"


def test_database_only_components_are_removed_but_isolated_queries_remain():
    graph = build_network_payload(
        _adata(),
        ["TP53", "CDK2"],
        "ppi",
        [
            NetworkEdge("TP53", "EGFR", score=0.9),
            NetworkEdge("EXTERNAL_A", "EXTERNAL_B", score=0.95),
        ],
    )

    nodes = {node["id"]: node for node in graph["nodes"]}
    assert set(nodes) == {"TP53", "EGFR", "CDK2"}
    assert nodes["CDK2"]["is_query"] is True
    assert nodes["CDK2"]["degree"] == 0
    assert {(edge["source"], edge["target"]) for edge in graph["edges"]} == {
        ("TP53", "EGFR")
    }


def test_ppi_node_size_uses_network_degree():
    graph = build_network_payload(
        _adata(),
        ["TP53"],
        "ppi",
        [
            NetworkEdge("TP53", "EGFR", score=0.9),
            NetworkEdge("TP53", "CDK2", score=0.8),
        ],
    )
    nodes = {node["id"]: node for node in graph["nodes"]}
    assert nodes["TP53"]["degree"] == 2
    assert nodes["TP53"]["node_size"] == pytest.approx(44)
    assert nodes["EGFR"]["node_size"] == pytest.approx(
        10 + 34 * np.log1p(1) / np.log1p(2)
    )
    assert nodes["TP53"]["node_size"] > nodes["EGFR"]["node_size"]

    rules = {
        rule["selector"]: rule["style"]
        for rule in network_stylesheet()
    }
    assert rules["node[degree >= 0]"]["width"] == "data(node_size)"
    assert rules["node[degree >= 0]"]["height"] == "data(node_size)"


def test_tf_gene_node_sizes_use_degree():
    adata = AnnData(
        X=np.asarray(
            [
                [0, 1, 1, 1],
                [0, 1, 1, 1],
                [0, 0, 0, 1],
                [1, 0, 0, 1],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["TF_LOW", "TF_HIGH", "TARGET_HALF", "TARGET_HIGH"]),
    )
    graph = build_network_payload(
        adata,
        ["TF_HIGH"],
        "tf-gene",
        [
            NetworkEdge("TF_LOW", "TARGET_HALF", directed=True),
            NetworkEdge("TF_HIGH", "TARGET_HALF", directed=True),
            NetworkEdge("TF_HIGH", "TARGET_HIGH", directed=True),
        ],
    )
    nodes = {node["id"]: node for node in graph["nodes"]}

    assert nodes["TF_HIGH"]["degree"] == 2
    assert nodes["TF_HIGH"]["node_size"] == pytest.approx(44)
    assert nodes["TARGET_HALF"]["node_size"] == pytest.approx(44)
    assert nodes["TF_HIGH"]["node_size"] > nodes["TARGET_HIGH"]["node_size"]

    rules = {
        rule["selector"]: rule["style"]
        for rule in network_stylesheet()
    }
    assert rules["node[degree >= 0]"]["width"] == "data(node_size)"


def test_metabolite_node_sizes_use_network_degree():
    adata = AnnData(
        X=np.ones((1, 1), dtype=np.float32),
        var=pd.DataFrame(index=["HK1"]),
    )
    graph = build_network_payload(
        adata,
        ["HK1"],
        "metabolite",
        [
            NetworkEdge("ATP", "HK1", source_type="metabolite"),
            NetworkEdge("Glucose", "HK1", source_type="metabolite"),
        ],
    )
    nodes = {node["id"]: node for node in graph["nodes"]}

    assert nodes["HK1"]["degree"] == 2
    assert nodes["HK1"]["node_size"] == pytest.approx(44)
    assert nodes["ATP"]["node_size"] < nodes["HK1"]["node_size"]


def test_network_expression_colors_follow_continuous_colormap():
    graph = build_network_payload(
        _adata(),
        ["TP53"],
        "ppi",
        [NetworkEdge("TP53", "EGFR", score=0.9)],
    )

    def node_colors(colormap):
        return {
            element["data"]["id"]: element["data"]["expression_color"]
            for element in _elements(graph, colormap)
            if "source" not in element["data"]
        }

    assert node_colors("viridis")["TP53"] != node_colors("plasma")["TP53"]


def test_only_input_genes_receive_the_query_outline_class():
    graph = build_network_payload(
        _adata(),
        ["TP53"],
        "ppi",
        [NetworkEdge("TP53", "EGFR", score=0.9)],
    )
    elements = {
        element["data"]["id"]: element
        for element in _elements(graph)
        if "source" not in element["data"]
    }

    assert "network-query-node" in elements["TP53"]["classes"].split()
    assert "network-query-node" not in elements["EGFR"]["classes"].split()


def test_zoom_labels_keep_queries_and_limit_early_hub_labels():
    graph = build_network_payload(
        _adata(),
        ["TP53"],
        "ppi",
        [NetworkEdge("TP53", f"GENE_{index}", score=0.9) for index in range(15)],
    )
    elements = {
        element["data"]["id"]: element
        for element in _elements(graph)
        if "source" not in element["data"]
    }
    query_classes = elements["TP53"]["classes"].split()
    hub_count = sum(
        "network-hub-node" in element.get("classes", "").split()
        for node_id, element in elements.items()
        if node_id != "TP53"
    )

    assert "network-query-node" in query_classes
    assert hub_count == 4

    rules = {rule["selector"]: rule["style"] for rule in network_stylesheet()}
    assert rules["node"]["min-zoomed-font-size"] == 24
    assert rules["node.network-hub-node"]["min-zoomed-font-size"] == 14
    assert rules["node.network-query-node"]["min-zoomed-font-size"] == 0


def test_tf_force_layout_separates_hubs_more_than_other_networks():
    tf_layout = network_layout("cose", "tf-gene")
    ppi_layout = network_layout("cose", "ppi")

    assert tf_layout["nodeRepulsion"] == 18000
    assert tf_layout["idealEdgeLength"] == 130
    assert tf_layout["componentSpacing"] == 120
    assert ppi_layout["nodeRepulsion"] == 6500
    assert ppi_layout["idealEdgeLength"] == 90


def test_tf_direction_filters_cached_edges_without_refetching():
    source = {
        "network_type": "tf-gene",
        "query_genes": ["STAT1"],
        "edges": [
            NetworkEdge("STAT1", "IRF1", directed=True).to_dict(),
            NetworkEdge("JAK1", "STAT1", directed=True).to_dict(),
        ],
    }

    targets = _build_view(_adata(), source, "first-order", tf_direction="targets")
    regulators = _build_view(_adata(), source, "first-order", tf_direction="regulators")
    both = _build_view(_adata(), source, "first-order", tf_direction="both")

    assert {(edge["source"], edge["target"]) for edge in targets["edges"]} == {
        ("STAT1", "IRF1")
    }
    assert {(edge["source"], edge["target"]) for edge in regulators["edges"]} == {
        ("JAK1", "STAT1")
    }
    assert len(both["edges"]) == 2


def test_mirna_network_prioritizes_shared_regulators_with_a_uniform_color():
    adata = AnnData(
        X=np.zeros((1, 3), dtype=np.float32),
        var=pd.DataFrame(index=["A", "B", "C"]),
    )
    edges = (
        [
            NetworkEdge("miR-1", target, source_type="mirna", directed=True, score=4)
            for target in ["A", "B", "C"]
        ]
        + [
            NetworkEdge("miR-2", target, source_type="mirna", directed=True, score=3)
            for target in ["A", "B"]
        ]
        + [NetworkEdge("miR-3", "A", source_type="mirna", directed=True, score=4)]
    )

    shared = build_network_payload(
        adata,
        ["A", "B", "C"],
        "mirna",
        edges,
        mirna_regulator_filter="shared-20",
    )
    all_validated = build_network_payload(
        adata,
        ["A", "B", "C"],
        "mirna",
        edges,
        mirna_regulator_filter="all",
    )

    assert {edge["source"] for edge in shared["edges"]} == {"miR-1", "miR-2"}
    assert shared["retained_regulator_count"] == 2
    assert shared["mirna_min_targets"] == 2
    nodes = {node["id"]: node for node in shared["nodes"]}
    assert nodes["miR-1"]["degree"] == 3
    assert nodes["miR-1"]["node_size"] == pytest.approx(44)
    assert nodes["miR-1"]["node_size"] > nodes["miR-2"]["node_size"]
    mirna_rule = {
        rule["selector"]: rule["style"]
        for rule in network_stylesheet()
    }['node[entity_type = "mirna"], node[entity_type = "microrna"]']
    assert mirna_rule["background-color"] == "#7b4f9d"
    assert "width" not in mirna_rule
    all_node_rule = {
        rule["selector"]: rule["style"]
        for rule in network_stylesheet()
    }["node[degree >= 0]"]
    assert all_node_rule["width"] == "data(node_size)"
    assert {edge["source"] for edge in all_validated["edges"]} == {
        "miR-1",
        "miR-2",
        "miR-3",
    }


def test_string_ppi_stylesheet_removes_arrows():
    rules = {rule["selector"]: rule["style"] for rule in network_stylesheet()}
    assert rules["node"]["border-width"] == 0
    assert rules["node.network-query-node"]["border-width"] == 3
    assert rules["edge.network-undirected-edge"]["target-arrow-shape"] == "none"

    graph = build_network_payload(
        _adata(),
        ["TP53"],
        "ppi",
        [NetworkEdge("TP53", "EGFR", database="STRING v12 functional")],
    )
    edge = next(element for element in _elements(graph) if "source" in element["data"])
    assert edge["classes"] == "network-undirected-edge"


def test_string_ppi_edge_width_uses_absolute_confidence_scale():
    assert _edge_width(0.4, [0.4, 0.7], network_type="ppi") == pytest.approx(1.2)
    assert _edge_width(0.7, [0.4, 0.7], network_type="ppi") == pytest.approx(3.2)
    assert _edge_width(0.7, [0.7, 0.99], network_type="ppi") == pytest.approx(3.2)
    assert _edge_width(1.0, [0.7, 1.0], network_type="ppi") == pytest.approx(5.2)


def test_network_views_keep_full_source_and_extract_input_or_minimum_edges():
    edges = [
        NetworkEdge("TP53", "CDK2", directed=True),
        NetworkEdge("CDK2", "EGFR", directed=True),
        NetworkEdge("TP53", "EGFR", directed=True),
    ]
    first_order = build_network_payload(
        _adata(), ["TP53", "EGFR"], "tf-gene", edges, view_mode="first-order"
    )
    input_only = build_network_payload(
        _adata(), ["TP53", "EGFR"], "tf-gene", edges, view_mode="input-only"
    )
    minimum = build_network_payload(
        _adata(), ["TP53", "EGFR"], "tf-gene", edges, view_mode="minimum"
    )

    assert first_order["available_edge_count"] == 3
    assert {(edge["source"], edge["target"]) for edge in input_only["edges"]} == {
        ("TP53", "EGFR")
    }
    assert {(edge["source"], edge["target"]) for edge in minimum["edges"]} == {
        ("TP53", "EGFR")
    }


def test_network_status_reports_exact_limits_instead_of_confidence():
    graph = build_network_payload(
        _adata(),
        ["TP53"],
        "tf-gene",
        [NetworkEdge("TP53", "EGFR"), NetworkEdge("TP53", "CDK2")],
        max_nodes=2,
        max_edges=8,
    )
    status = _network_status(graph)
    assert "Showing 2 of 3 nodes and 1 of 2 interactions" in status
    assert "highest-confidence" not in status


def test_network_callback_queries_only_on_build_and_keeps_previous_graph_on_error():
    app = Dash(__name__)
    calls = []
    adata = _adata()

    def resolve_filtered(filtered):
        indices = (filtered or {}).get("cell_indices")
        return adata[indices] if indices is not None else adata

    def edge_fetcher(network_type, genes, organism, *, ppi_network_type="functional"):
        calls.append((network_type, genes, organism, ppi_network_type))
        return [
            NetworkEdge("TP53", "EGFR", score=0.9),
            NetworkEdge("TP53", "CDK2", score=0.8),
        ]

    register_network_callbacks(
        app,
        adata,
        "p",
        organism="human",
        resolve_plot_adata_from_filter=resolve_filtered,
        hash_list_signature=lambda values: str(values),
        edge_fetcher=edge_fetcher,
    )
    callback_entry = next(
        value
        for key, value in app.callback_map.items()
        if "p-network-graph.children" in key
    )
    assert {item["id"] for item in callback_entry["inputs"]} == {"p-network-build"}
    callback = callback_entry["callback"].__wrapped__
    (
        component,
        source,
        graph,
        rendered_key,
        highlight,
        status,
        regulator_style,
        regulator_component,
    ) = callback(
        1,
        "TP53, EGFR",
        "ppi",
        [],
        "targets",
        "first-order",
        0.85,
        "cose",
        "viridis",
        {"cell_indices": ["c1", "c2"]},
        ["c1"],
        {"len": 1, "hash": "c1"},
        "dataset-exploration",
        "network-tab",
        None,
    )
    assert calls == [("ppi", ["TP53", "EGFR"], "human", "functional")]
    assert graph["edges"][0]["source"] == "TP53"
    assert graph["ppi_network_type"] == "functional"
    nodes = {node["id"]: node for node in graph["nodes"]}
    assert nodes["TP53"]["expression"] == pytest.approx(0.0)
    assert nodes["EGFR"]["expression"] == pytest.approx(1.0)
    assert len(graph["edges"]) == 1
    assert len(source["edges"]) == 2
    assert rendered_key
    assert highlight is None
    assert status == ""
    assert regulator_style == {"display": "none"}
    assert regulator_component is None
    assert component.className == "network-graph-content"

    view_entry = next(
        value
        for value in app.callback_map.values()
        if {item["id"] for item in value["inputs"]}
        == {
            "p-network-view",
            "p-network-string-confidence",
            "p-network-tf-direction",
        }
    )
    view_callback = view_entry["callback"].__wrapped__
    _, input_graph, input_highlight, input_status = view_callback(
        "input-only",
        0.7,
        "targets",
        source,
        "cose",
        "viridis",
        {"cell_indices": ["c1", "c2"]},
        ["c1"],
    )
    assert len(input_graph["edges"]) == 1
    assert input_highlight is None
    assert input_status == ""
    assert calls == [("ppi", ["TP53", "EGFR"], "human", "functional")]

    appearance_entry = next(
        value
        for value in app.callback_map.values()
        if {item["id"] for item in value["inputs"]}
        == {"p-network-layout", "p-scatter-color-map-dropdown"}
    )
    appearance_callback = appearance_entry["callback"].__wrapped__
    relaid_component, relaid_highlight = appearance_callback(
        "circle", "viridis", graph
    )
    assert (
        _by_id(relaid_component, "p-network-cytoscape-view").layout["name"] == "circle"
    )
    assert relaid_highlight is None
    assert calls == [("ppi", ["TP53", "EGFR"], "human", "functional")]

    def failing_fetcher(*_args, **_kwargs):
        raise ValueError("database unavailable")

    failing = Dash(__name__ + "-failing")
    failing_adata = _adata()
    register_network_callbacks(
        failing,
        failing_adata,
        "q",
        resolve_plot_adata_from_filter=lambda _filtered: failing_adata,
        hash_list_signature=lambda values: str(values),
        edge_fetcher=failing_fetcher,
    )
    failing_callback = next(
        value
        for key, value in failing.callback_map.items()
        if "q-network-graph.children" in key
    )["callback"].__wrapped__
    result = failing_callback(
        1,
        "TP53",
        "ppi",
        [],
        "targets",
        "first-order",
        0.7,
        "cose",
        "viridis",
        {"cell_indices": None},
        None,
        None,
        "dataset-exploration",
        "network-tab",
        None,
    )
    assert result[:5] == (no_update, no_update, no_update, no_update, no_update)
    assert result[5] == "database unavailable"
    assert result[6:] == (no_update, no_update)


def test_trrust_build_adds_key_regulator_results():
    app = Dash(__name__ + "-trrust")
    query = ["TP53", "EGFR", "CDK2", "STAT1", "IRF1"]
    trrust_edges = [
        NetworkEdge("TF1", gene, directed=True, database="TRRUST v2") for gene in query
    ]
    calls = []

    def edge_fetcher(_network_type, _genes, _organism):
        return trrust_edges

    def regulator_fetcher(organism):
        calls.append(organism)
        return trrust_edges

    adata = AnnData(
        X=np.zeros((1, len(query)), dtype=np.float32),
        var=pd.DataFrame(index=query),
    )
    register_network_callbacks(
        app,
        adata,
        "r",
        organism="human",
        resolve_plot_adata_from_filter=lambda _filtered: adata,
        hash_list_signature=lambda values: str(values),
        edge_fetcher=edge_fetcher,
        regulator_fetcher=regulator_fetcher,
    )
    callback = next(
        value
        for key, value in app.callback_map.items()
        if "r-network-graph.children" in key
    )["callback"].__wrapped__
    result = callback(
        1,
        ", ".join(query),
        "tf-gene",
        [],
        "targets",
        "first-order",
        0.7,
        "cose",
        "viridis",
        {"cell_indices": None},
        None,
        None,
        "dataset-exploration",
        "network-tab",
        None,
    )

    assert calls == ["human"]
    assert result[6] == {"display": "block"}
    assert result[7].className == "network-enrichment-results"


def test_normalize_organism_accepts_names_genomes_and_taxon_ids():
    assert normalize_organism("hg38") == ("human", 9606)
    assert normalize_organism("mouse") == ("mouse", 10090)
    assert normalize_organism(10116) == ("rat", 10116)
