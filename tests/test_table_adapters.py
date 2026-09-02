import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from guanaco.linking.data import DataStore
from guanaco.linking.model import MarkMembers, ViewState, view
from guanaco.linking.table_adapters import (
    NetworkAdapter,
    TableHeatmapAdapter,
    TablePlotAdapter,
)


def _table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0],
            "y": [2.0, 1.0, 3.0],
            "group": ["a", "a", "b"],
        },
        index=["row-1", "row-2", "row-3"],
    )


def _row_ids(figure: go.Figure) -> list[str]:
    return [str(row[0]) for trace in figure.data for row in trace.customdata]


@pytest.mark.parametrize(
    ("plot", "trace_type"),
    [
        ("plotly.scatter", "scatter"),
        ("plotly.bar", "bar"),
        ("plotly.line", "scatter"),
        ("plotly.area", "scatter"),
    ],
)
def test_table_plot_marks_keep_dataframe_row_ids(plot, trace_type) -> None:
    frame = _table()
    figure = TablePlotAdapter().render(
        view(plot, id="table", x="x", y="y", group="group"),
        DataStore.from_data(frame),
    )

    assert {trace.type for trace in figure.data} == {trace_type}
    assert set(_row_ids(figure)) == set(frame.index)


def test_table_plot_filters_highlights_and_decodes_the_same_index() -> None:
    store = DataStore.from_data(_table())
    spec = view("plotly.scatter", id="table", x="x", y="y")
    adapter = TablePlotAdapter()

    filtered = adapter.render(spec, store, ViewState(rows=("row-2",)))
    highlighted = adapter.render(
        spec,
        store,
        ViewState(cells=("row-2",), highlight_axes=("cell",)),
    )
    selected = adapter.decode_event(
        "select",
        {
            "points": [
                {"customdata": ["row-3"]},
                {"customdata": ["row-1"]},
                {"customdata": ["row-3"]},
            ]
        },
        spec,
        store,
    )

    assert _row_ids(filtered) == ["row-2"]
    assert _row_ids(highlighted) == ["row-1", "row-2", "row-3"]
    assert list(highlighted.data[0].selectedpoints) == [1]
    assert selected == MarkMembers(
        rows=("row-3", "row-1"),
    )
    assert adapter.emits == {"row", "cell", "feature"}
    assert adapter.accepts == {"row", "cell", "feature"}


def test_filtered_scatter_keeps_full_numeric_axes_by_default() -> None:
    store = DataStore.from_data(_table())
    spec = view("plotly.scatter", id="table", x="x", y="y")
    adapter = TablePlotAdapter()

    full = adapter.render(spec, store)
    filtered = adapter.render(spec, store, ViewState(rows=("row-2",)))
    automatic = adapter.render(
        view("plotly.scatter", id="automatic", x="x", y="y", fixed_axes=False),
        store,
        ViewState(rows=("row-2",)),
    )

    assert tuple(filtered.layout.xaxis.range) == tuple(full.layout.xaxis.range)
    assert tuple(filtered.layout.yaxis.range) == tuple(full.layout.yaxis.range)
    assert full.layout.xaxis.range[0] < _table()["x"].min()
    assert full.layout.xaxis.range[1] > _table()["x"].max()
    assert automatic.layout.xaxis.range is None
    assert automatic.layout.yaxis.range is None


def test_filtered_scatter_drops_unused_categories_in_declared_order() -> None:
    frame = _table()
    frame["category"] = pd.Categorical(
        ["alpha", "beta", "gamma"],
        categories=["gamma", "beta", "alpha"],
        ordered=True,
    )
    frame["other"] = pd.Categorical(
        ["first", "second", "third"],
        categories=["third", "second", "first"],
        ordered=True,
    )
    store = DataStore.from_data(frame)
    spec = view("plotly.scatter", id="table", x="category", y="other")

    full = TablePlotAdapter().render(spec, store)
    filtered = TablePlotAdapter().render(spec, store, ViewState(rows=("row-1",)))

    assert tuple(full.layout.xaxis.categoryarray) == ("gamma", "beta", "alpha")
    assert tuple(filtered.layout.xaxis.categoryarray) == ("alpha",)
    assert tuple(filtered.layout.yaxis.categoryarray) == ("first",)
    assert tuple(filtered.layout.xaxis.range) == (-0.5, 0.5)
    assert tuple(filtered.layout.yaxis.range) == (-0.5, 0.5)


def test_scatter_supports_explicit_equal_axes_and_one_global_color_scale() -> None:
    frame = _table().assign(score=[-0.7, 0.1, 0.9])
    spec = view(
        "plotly.scatter",
        id="relationship",
        x="x",
        y="y",
        group="group",
        color="score",
        x_range=(-10, 10),
        y_range=(-5, 15),
        equal_axes=True,
    )

    figure = TablePlotAdapter().render(spec, DataStore.from_data(frame))

    assert tuple(figure.layout.xaxis.range) == (-10, 10)
    assert tuple(figure.layout.yaxis.range) == (-5, 15)
    assert figure.layout.yaxis.scaleanchor == "x"
    assert {trace.marker.coloraxis for trace in figure.data} == {"coloraxis"}
    assert figure.layout.coloraxis.cmin == pytest.approx(-0.7)
    assert figure.layout.coloraxis.cmax == pytest.approx(0.9)
    assert figure.layout.coloraxis.colorbar.title.text == "score"


def test_scatter_adds_a_generic_numpy_layout_image() -> None:
    pixels = np.zeros((3, 4, 3), dtype=np.uint8)
    pixels[:, :, 1] = 180
    spec = view(
        "plotly.scatter",
        id="spatial",
        x="x",
        y="y",
        layout_image={
            "source": pixels,
            "x": -0.1,
            "y": 3.1,
            "sizex": 2.2,
            "sizey": 2.2,
            "opacity": 0.65,
        },
    )

    image = (
        TablePlotAdapter().render(spec, DataStore.from_data(_table())).layout.images[0]
    )

    assert str(image.source).startswith("data:image/png;base64,")
    assert image.x == pytest.approx(-0.1)
    assert image.y == pytest.approx(3.1)
    assert image.sizex == pytest.approx(2.2)
    assert image.sizey == pytest.approx(2.2)
    assert image.opacity == pytest.approx(0.65)
    assert image.layer == "below"


def _heatmap_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "column": ["A", "B", "A", "A"],
            "row": ["R1", "R1", "R2", "R1"],
            "value": [0.2, 0.4, 0.8, 1.0],
        },
        index=["tile-1", "tile-2", "tile-3", "tile-4"],
    )


def _heatmap_spec(**options):
    return view(
        "plotly.heatmap",
        id="heatmap",
        x="column",
        y="row",
        value="value",
        **options,
    )


@pytest.mark.parametrize(
    ("aggregate", "expected"),
    [
        ("count", 2.0),
        ("sum", 1.2),
        ("mean", 0.6),
        ("median", 0.6),
        ("min", 0.2),
        ("max", 1.0),
    ],
)
def test_long_heatmap_aggregates_duplicate_atomic_rows(aggregate, expected) -> None:
    figure = TableHeatmapAdapter().render(
        _heatmap_spec(aggregate=aggregate),
        DataStore.from_data(_heatmap_table()),
    )

    assert figure.data[0].z[0][0] == pytest.approx(expected)


def test_heatmap_tile_decodes_all_rows_and_filter_reaggregates() -> None:
    store = DataStore.from_data(_heatmap_table())
    spec = _heatmap_spec(aggregate="sum")
    adapter = TableHeatmapAdapter()

    members = adapter.decode_event(
        "click", {"points": [{"x": "A", "y": "R1"}]}, spec, store
    )
    filtered = adapter.render(spec, store, ViewState(rows=("tile-4",)))

    expected = ("tile-1", "tile-4")
    assert members == MarkMembers(rows=expected)
    assert filtered.data[0].z[0][0] == pytest.approx(1.0)


def test_heatmap_highlight_marks_tiles_without_filtering_them() -> None:
    store = DataStore.from_data(_heatmap_table())
    spec = _heatmap_spec()

    figure = TableHeatmapAdapter().render(
        spec,
        store,
        ViewState(rows=("tile-2",), highlight_axes=("row",)),
    )

    assert len(figure.data) == 2
    assert list(figure.data[0].x) == ["A", "B"]
    assert list(figure.data[1].x) == ["B"]
    assert list(figure.data[1].y) == ["R1"]


def _network_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sender": ["NK", "NK", "NK", "Mono"],
            "receiver": ["B", "B", "B", "T"],
            "score": [1.0, 3.0, 5.0, 8.0],
        },
        index=["interaction-3", "interaction-1", "interaction-2", "interaction-4"],
    )


def _network_edges(component):
    cytoscape = component.children[1]
    return [element for element in cytoscape.elements if "source" in element["data"]]


def test_network_groups_long_rows_and_decodes_an_edge_without_number_labels() -> None:
    store = DataStore.from_data(_network_table())
    spec = view(
        "network",
        id="network",
        source="sender",
        target="receiver",
        edge_weight="score",
        aggregate="mean",
        height=480,
    )
    adapter = NetworkAdapter()

    component = adapter.render(spec, store, component_id="network-component")
    edges = _network_edges(component)
    nk_to_b = next(edge for edge in edges if edge["data"]["weight"] == 3)
    members = adapter.decode_event("click", nk_to_b["data"], spec, store)

    expected = ("interaction-3", "interaction-1", "interaction-2")
    assert members == MarkMembers(rows=expected)
    assert component.style["height"] == "480px"
    assert component.children[1].id == "network-component"
    assert "label" not in nk_to_b["data"]
    edge_style = next(
        rule["style"]
        for rule in component.children[1].stylesheet
        if rule["selector"] == "edge"
    )
    assert "label" not in edge_style


def test_network_filter_reaggregates_and_highlight_keeps_all_edges() -> None:
    store = DataStore.from_data(_network_table())
    spec = view("network", id="network", source="sender", target="receiver")
    adapter = NetworkAdapter()

    filtered = adapter.render(
        spec,
        store,
        ViewState(rows=("interaction-1", "interaction-4")),
    )
    highlighted = adapter.render(
        spec,
        store,
        ViewState(cells=("interaction-2",), highlight_axes=("cell",)),
    )

    filtered_edges = _network_edges(filtered)
    highlighted_edges = _network_edges(highlighted)
    assert len(filtered_edges) == 2
    assert {edge["data"]["weight"] for edge in filtered_edges} == {1}
    assert len(highlighted_edges) == 2
    assert (
        sum(
            edge.get("classes") == "guanaco-linked-selected"
            for edge in highlighted_edges
        )
        == 1
    )


def test_table_target_can_filter_one_axis_and_highlight_another() -> None:
    store = DataStore.from_data(_table())
    spec = view("plotly.scatter", id="detail", x="x", y="y")
    state = ViewState(
        rows=("row-1", "row-2"),
        cells=("row-2", "row-3"),
        highlight_axes=("cell",),
    )

    figure = TablePlotAdapter().render(spec, store, state)
    selected = {
        str(trace.customdata[position][0])
        for trace in figure.data
        for position in (trace.selectedpoints or ())
    }

    assert _row_ids(figure) == ["row-1", "row-2"]
    assert selected == {"row-2"}
