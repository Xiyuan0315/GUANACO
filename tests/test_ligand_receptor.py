import base64

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from dash import Dash
from dash.exceptions import PreventUpdate

from guanaco.data.ligand_receptor import (
    discover_ligand_receptor_results,
    load_ligand_receptor_result,
    parse_uploaded_ligand_receptor_results,
)
from guanaco.pages.matrix.callbacks import ligand_receptor_callbacks
from guanaco.pages.matrix.callbacks.ligand_receptor_callbacks import (
    _selection_candidate,
    register_ligand_receptor_callbacks,
)
from guanaco.pages.matrix.plots.ligand_receptor import (
    build_ligand_receptor_view,
    ligand_receptor_detail,
    ligand_receptor_highlight_stylesheet,
    plot_ligand_receptor_dotplot,
)
from guanaco.pages.visualizations.layout import generate_visualization_sections
from guanaco.pages.visualizations.registry import resolve_plot_components


def _liana_table():
    return pd.DataFrame(
        {
            "source": ["Monocyte", "Monocyte", "B cell"],
            "target": ["T cell", "T cell", "T cell"],
            "ligand_complex": ["IL1B", "CXCL10", "CD40"],
            "receptor_complex": ["IL1R1", "CXCR3", "CD40LG"],
            "magnitude_rank": [0.02, 0.001, 0.05],
            "specificity_rank": [0.03, 0.004, 0.2],
        }
    )


def _adata():
    adata = AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(
                    ["Monocyte", "Monocyte", "B cell", "B cell", "T cell", "T cell"]
                )
            },
            index=[f"cell-{index}" for index in range(6)],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    adata.uns["liana_res"] = _liana_table()
    return adata


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


def test_liana_and_cellchat_long_tables_use_the_canonical_adapter():
    adata = _adata()
    options = discover_ligand_receptor_results(adata)
    assert options == [
        {"label": "liana_res · interaction table", "value": "liana_res"}
    ]

    payload = load_ligand_receptor_result(adata, "liana_res")
    assert payload["format"] == "interaction table"
    assert payload["default_magnitude"] == "magnitude_rank"
    assert payload["default_specificity"] == "specificity_rank"
    assert payload["metrics"][0]["direction"] == "lower"
    assert payload["records"][0]["ligand"] == "IL1B"

    adata.uns["cellchat"] = pd.DataFrame(
        {
            "source": ["A"],
            "target": ["B"],
            "ligand": ["L"],
            "receptor": ["R"],
            "prob": [0.7],
            "pval": [0.01],
        }
    )
    cellchat = load_ligand_receptor_result(adata, "cellchat")
    assert cellchat["default_magnitude"] == "prob"
    assert cellchat["default_specificity"] == "pval"


def test_squidpy_mapping_is_converted_from_pair_matrices():
    adata = _adata()
    interaction_index = pd.MultiIndex.from_tuples(
        [("L1", "R1"), ("L2", "R2")],
        names=["source", "target"],
    )
    group_columns = pd.MultiIndex.from_tuples(
        [("A", "B"), ("B", "A")],
        names=["source", "target"],
    )
    adata.uns["cluster_ligrec"] = {
        "means": pd.DataFrame(
            [[1.2, 0.0], [2.3, 1.1]],
            index=interaction_index,
            columns=group_columns,
        ),
        "pvalues": pd.DataFrame(
            [[0.01, np.nan], [0.02, 0.2]],
            index=interaction_index,
            columns=group_columns,
        ),
        "metadata": pd.DataFrame(),
    }
    payload = load_ligand_receptor_result(adata, "cluster_ligrec")
    assert payload["format"] == "Squidpy"
    assert payload["default_magnitude"] == "mean"
    assert payload["default_specificity"] == "pvalue"
    assert {
        (row["source"], row["target"], row["ligand"], row["receptor"])
        for row in payload["records"]
    } >= {("A", "B", "L1", "R1"), ("B", "A", "L2", "R2")}


def test_cellphonedb_wide_files_and_canonical_upload_are_supported():
    means = pd.DataFrame(
        {
            "interacting_pair": ["L|R"],
            "gene_a": ["L"],
            "gene_b": ["R"],
            "receptor_a": [False],
            "receptor_b": [True],
            "A|B": [2.5],
        }
    )
    pvalues = means.copy()
    pvalues["A|B"] = 0.01
    adata = _adata()
    adata.uns["cpdb"] = {"means": means, "pvalues": pvalues}
    payload = load_ligand_receptor_result(adata, "cpdb")
    assert payload["format"] == "CellPhoneDB"
    assert payload["records"][0] == {
        "source": "A",
        "target": "B",
        "ligand": "L",
        "receptor": "R",
        "mean": 2.5,
        "pvalue": 0.01,
    }

    csv = _liana_table().to_csv(index=False).encode()
    encoded = "data:text/csv;base64," + base64.b64encode(csv).decode()
    uploaded = parse_uploaded_ligand_receptor_results(encoded, "results.csv")
    assert uploaded["format"] == "uploaded interaction table"
    assert len(uploaded["records"]) == 3


def test_circle_edges_and_dotplot_share_one_filtered_view():
    payload = load_ligand_receptor_result(_adata(), "liana_res")
    view = build_ligand_receptor_view(
        payload,
        magnitude="magnitude_rank",
        specificity="specificity_rank",
        senders=["Monocyte"],
    )
    assert len(view["edges"]) == 1
    assert view["edges"][0]["pair_count"] == 2
    assert view["interactions"][0]["ligand"] == "CXCL10"

    figure = plot_ligand_receptor_dotplot(
        view,
        selected_edge=view["edges"][0]["id"],
    )
    assert len(figure.data[0].x) == 2
    assert len(figure.data) == 4
    assert figure.layout.legend.title.text == "specificity_rank"
    assert figure.data[0].marker.colorbar.y == 0.7
    assert figure.layout.legend.y == 0.36
    assert figure.data[0].marker.reversescale is True
    assert figure.data[0].customdata[0][1] == view["edges"][0]["id"]

    stylesheet = ligand_receptor_highlight_stylesheet(
        view,
        edge_id=view["edges"][0]["id"],
    )
    selectors = {rule["selector"] for rule in stylesheet}
    assert f'edge[id = "{view["edges"][0]["id"]}"]' in selectors

    thresholded = build_ligand_receptor_view(
        payload,
        magnitude="magnitude_rank",
        specificity="specificity_rank",
        magnitude_range=[0, 0.01],
        specificity_range=[0, 0.01],
        senders=["Monocyte"],
    )
    assert [row["ligand"] for row in thresholded["interactions"]] == ["CXCL10"]


def test_layout_and_callbacks_are_result_driven_not_expression_driven(monkeypatch):
    adata = _adata()
    assert resolve_plot_components(
        adata,
        ["ligand-receptor"],
        modality_name="atac",
    ) == ("ligand-receptor",)

    sections = generate_visualization_sections(
        adata,
        ["G1"],
        ["cell_type"],
        "p",
        optional_plot_components=["ligand-receptor"],
        modality_name="rna",
    )
    assert _by_id(sections, "p-lr-magnitude").value == "magnitude_rank"
    assert _by_id(sections, "p-lr-senders").options
    assert _by_id(sections, "p-lr-options-collapse").is_open is False
    assert _by_id(sections, "p-lr-options-toggle").children == "▸ More options"
    assert _by_id(sections, "p-lr-magnitude-range").value == [0.001, 0.05]
    assert _by_id(sections, "p-lr-specificity-range").value == [0.004, 0.2]
    assert _by_id(sections, "p-lr-selection-store") is not None
    assert _by_id(sections, "p-lr-dotplot") is not None
    component_ids = {
        getattr(component, "id", None) for component in _walk(sections)
    }
    assert "p-lr-result-source" not in component_ids
    assert "p-lr-result-store" not in component_ids
    assert "p-lr-upload" not in component_ids
    assert "p-lr-run" not in component_ids
    assert "p-lr-groupby" not in component_ids

    app = Dash(__name__)
    register_ligand_receptor_callbacks(
        app,
        adata,
        "p",
        color_config=["#111111", "#777777", "#bbbbbb"],
    )
    payload = load_ligand_receptor_result(adata, "liana_res")
    assert payload["default_magnitude"] == "magnitude_rank"

    render_entry = next(
        value
        for key, value in app.callback_map.items()
        if key.startswith("..p-lr-network.children")
    )
    render_inputs = {item["id"] for item in render_entry["inputs"]}
    assert "p-lr-result-store" not in render_inputs
    assert "p-global-filtered-data" not in render_inputs
    assert "p-selected-cells-store" not in render_inputs
    component, view = render_entry["callback"].__wrapped__(
        "magnitude_rank",
        "specificity_rank",
        [0.001, 0.05],
        [0.004, 0.2],
        [],
        [],
        "plotly/Plotly",
    )
    assert component.className == "lr-network-content"
    assert len(view["interactions"]) == 3

    interaction = view["interactions"][0]
    selection = _selection_candidate(
        "p-lr-dotplot.clickData",
        None,
        None,
        {
            "points": [
                {
                    "customdata": [
                        interaction["interaction_id"],
                        interaction["edge_id"],
                    ]
                }
            ]
        },
    )
    assert selection == {
        "kind": "interaction",
        "interaction_id": interaction["interaction_id"],
        "edge_id": interaction["edge_id"],
    }
    selection_entry = app.callback_map["p-lr-selection-store.data"]
    monkeypatch.setattr(
        ligand_receptor_callbacks,
        "_triggered_property",
        lambda: "p-lr-dotplot.clickData",
    )
    selected = selection_entry["callback"].__wrapped__(
        view,
        None,
        None,
        {
            "points": [
                {
                    "customdata": [
                        interaction["interaction_id"],
                        interaction["edge_id"],
                    ]
                }
            ]
        },
        None,
    )
    assert selected == selection
    dotplot_entry = app.callback_map["p-lr-dotplot.figure"]
    monkeypatch.setattr(
        ligand_receptor_callbacks,
        "_triggered_property",
        lambda: "p-lr-selection-store.data",
    )
    with pytest.raises(PreventUpdate):
        dotplot_entry["callback"].__wrapped__(
            view,
            selected,
            "plotly/Plotly",
        )
    monkeypatch.setattr(
        ligand_receptor_callbacks,
        "_triggered_property",
        lambda: "p-lr-dotplot.clickData",
    )
    assert (
        selection_entry["callback"].__wrapped__(
            view,
            None,
            None,
            {
                "points": [
                    {
                        "customdata": [
                            interaction["interaction_id"],
                            interaction["edge_id"],
                        ]
                    }
                ]
            },
            selected,
        )
        is None
    )
    detail_entry = app.callback_map["p-lr-detail.children"]
    detail = detail_entry["callback"].__wrapped__(selection, view)
    assert detail.className == "lr-detail-content"
    assert detail.children[0].children == "1 selected interaction"
    all_details = detail_entry["callback"].__wrapped__(None, view)
    assert all_details.children[0].children == "3 interactions"
    node_details = detail_entry["callback"].__wrapped__(
        {"kind": "node", "node_id": "Monocyte"},
        view,
    )
    assert node_details.children[0].children == "2 selected interactions"
    assert len(node_details.children[1].children.data) == 2

    highlight_entry = app.callback_map["p-lr-cytoscape-view.stylesheet"]
    reset_stylesheet = highlight_entry["callback"].__wrapped__(None, view)
    assert reset_stylesheet == ligand_receptor_highlight_stylesheet(view)

    many_view = dict(view)
    many_view["interactions"] = [
        {
            **view["interactions"][index % len(view["interactions"])],
            "interaction_id": f"many-{index}",
        }
        for index in range(25)
    ]
    top_figure = plot_ligand_receptor_dotplot(many_view)
    assert len(top_figure.data[0].x) == 20
    assert top_figure.layout.title.text.endswith("Top 20 of 25")
    top_details = ligand_receptor_detail(many_view)
    assert top_details.children[0].children == "25 interactions"
    table = top_details.children[1].children
    assert len(table.data) == 25
    assert table.page_size == 25
    assert table.page_action == "native"
