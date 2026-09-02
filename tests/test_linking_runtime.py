from __future__ import annotations

from collections.abc import Mapping

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from guanaco.linking.base import PlotAdapter
from guanaco.linking.data import DataStore
from guanaco.linking.model import MarkMembers, ViewSpec, link, view
from guanaco.linking.registry import PlotRegistry
from guanaco.linking.engine import state_for_target
from guanaco.linking.runtime import LinkedView, _coalesce_events


class IndexAdapter(PlotAdapter):
    events = {"click": "clickData", "select": "selectedData"}

    def __init__(
        self,
        *,
        emits=("row", "cell", "feature"),
        accepts=("row", "cell", "feature"),
    ):
        self.emits = frozenset(emits)
        self.accepts = frozenset(accepts)

    def validate(self, spec: ViewSpec, store: DataStore) -> None:
        store.source(spec.data)

    def render(self, spec, store, state=None, *, component_id=None):
        del store, state, component_id
        return go.Figure(
            go.Scatter(x=[0, 1], y=[0, 1], customdata=[["a"], ["b"]]),
            layout={"title": spec.title},
        )

    def decode_event(
        self,
        event: str,
        payload: Mapping | None,
        spec: ViewSpec,
        store: DataStore,
    ) -> MarkMembers | None:
        del event, spec, store
        ids = tuple(
            str(point["customdata"][0])
            for point in ((payload or {}).get("points") or ())
        )
        return MarkMembers(rows=ids, cells=ids, features=ids) if ids else None


def _registry(**adapters) -> PlotRegistry:
    registry = PlotRegistry()
    for name, adapter in adapters.items():
        registry.register(name, adapter)
    return registry


def _adata() -> ad.AnnData:
    value = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["CD4", "IL7R"]),
    )
    value.obsm["X_umap"] = np.asarray([[0.0, 0.0], [1.0, 1.0]])
    return value


def test_link_defaults_follow_the_two_data_domains() -> None:
    anndata_view = LinkedView(
        _adata(),
        views=[view("plot", "umap"), view("plot", "violin")],
        links=[link("umap", "violin")],
        registry=_registry(plot=IndexAdapter(emits=("cell",))),
    )
    table = pd.DataFrame({"x": [1, 2]}, index=["row-1", "row-2"])
    table_view = LinkedView(
        table,
        views=[view("plot", "overview"), view("plot", "detail")],
        links=[link("overview", "detail")],
        registry=_registry(plot=IndexAdapter(emits=("row",))),
    )

    assert anndata_view.links[0].selection_by == "cell"
    assert anndata_view.links[0].resolved_action == "highlight"
    assert table_view.links[0].selection_by == "row"
    assert table_view.links[0].resolved_action == "filter"


def test_table_to_anndata_requires_the_index_meaning() -> None:
    data = {
        "scores": pd.DataFrame({"x": [1, 2]}, index=["cell-1", "cell-2"]),
        "adata": _adata(),
    }
    specs = [
        view("plot", "scores", data="scores"),
        view("plot", "umap", data="adata"),
    ]
    registry = _registry(plot=IndexAdapter())

    with pytest.raises(ValueError, match="by='cell'.*by='feature'"):
        LinkedView(data, views=specs, links=[link("scores", "umap")], registry=registry)

    workspace = LinkedView(
        data,
        views=specs,
        links=[link("scores", "umap", by="cell")],
        registry=registry,
    )
    assert workspace.links[0].selection_by == "cell"


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("duplicate", "Duplicate links"),
        ("no_events", "no selectable browser event"),
        ("no_emits", "emits no IDs"),
        ("no_accepts", "accepts no identity axes"),
        ("no_overlap", "shares no row IDs"),
    ],
)
def test_runtime_rejects_ambiguous_or_incomplete_link_contracts(case, message) -> None:
    data = pd.DataFrame({"x": [1]}, index=["row-1"])
    specs = [view("source", "source"), view("target", "target")]
    links = [link("source", "target")]
    source = IndexAdapter(emits=("row",))
    target = IndexAdapter(accepts=("row",))

    if case == "duplicate":
        links *= 2
    elif case == "no_events":
        source.events = {}
    elif case == "no_emits":
        source.emits = frozenset()
    elif case == "no_accepts":
        target.accepts = frozenset()
    else:
        data = {
            "left": pd.DataFrame({"x": [1]}, index=["row-left"]),
            "right": pd.DataFrame({"x": [2]}, index=["row-right"]),
        }
        specs = [
            view("source", "source", data="left"),
            view("target", "target", data="right"),
        ]

    with pytest.raises(ValueError, match=message):
        LinkedView(
            data,
            views=specs,
            links=links,
            registry=_registry(source=source, target=target),
        )


def test_details_are_terminal_not_new_link_sources() -> None:
    table = pd.DataFrame({"x": [1]}, index=["row-1"])

    with pytest.raises(ValueError, match="overview.*terminal detail"):
        LinkedView(
            table,
            views=[view("plot", name) for name in ("a", "b", "c")],
            links=[link("a", "b"), link("b", "c")],
            registry=_registry(plot=IndexAdapter(emits=("row",))),
        )


def test_callback_batch_prefers_a_non_empty_sibling_event() -> None:
    selected = MarkMembers(cells=("cell-1",))

    assert _coalesce_events([("umap", selected), ("umap", None)]) == [
        ("umap", selected)
    ]


def test_notebook_can_retrieve_and_clear_the_latest_selection() -> None:
    workspace = LinkedView(
        _adata(),
        views=[view("plot", "umap"), view("plot", "violin")],
        links=[link("umap", "violin")],
        registry=_registry(plot=IndexAdapter(emits=("cell",))),
    )

    workspace._record_members("umap", MarkMembers(cells=("cell-2",)))
    selected = workspace.get_selection("umap")

    assert selected is not None
    assert selected.by == "cell"
    assert selected.ids == ("cell-2",)

    workspace._record_members("umap", None)
    assert workspace.get_selection("umap") is None


def test_notebook_retrieves_typed_ids_from_a_custom_table_adapter() -> None:
    frame = pd.DataFrame({"x": [1, 2]}, index=["cell-1", "cell-2"])
    workspace = LinkedView(
        frame,
        views=[view("plot", "source"), view("plot", "detail")],
        links=[link("source", "detail", by="cell")],
        registry=_registry(plot=IndexAdapter(emits=("cell",))),
    )

    workspace._record_members("source", MarkMembers(cells=("cell-2",)))

    assert workspace.get_selection("source").ids == ("cell-2",)


def test_logical_table_key_expands_one_parent_to_many_detail_rows() -> None:
    overview = pd.DataFrame(
        {
            "pair_id": ["pair-1"],
            "source": ["A"],
            "target": ["B"],
            "value": [1.0],
        },
        index=["pair-1"],
    )
    details = pd.DataFrame(
        {"pair_id": ["pair-1", "pair-1", "other"], "x": [0, 1, 2], "y": [0, 1, 2]},
        index=["spot-1", "spot-2", "spot-3"],
    )
    workspace = LinkedView(
        {"overview": overview, "details": details},
        views=[
            view(
                "plotly.heatmap",
                "overview",
                data="overview",
                x="source",
                y="target",
                value="value",
            ),
            view("plotly.scatter", "details", data="details", x="x", y="y"),
        ],
        links=[link("overview", "details", key="pair_id")],
    )

    workspace._record_members("overview", MarkMembers(rows=("pair-1",)))

    target = state_for_target(
        workspace._state,
        "details",
        workspace.links,
        workspace._table_views,
        workspace._row_key_maps,
    )
    assert target.rows == ("spot-1", "spot-2")
    assert workspace.get_selection("overview").ids == ("pair-1",)
    assert workspace.describe()["links"][0]["key"] == "pair_id"


def test_source_state_keeps_only_id_axes_used_by_outgoing_links() -> None:
    workspace = LinkedView(
        _adata(),
        views=[view("plot", "source"), view("plot", "detail")],
        links=[link("source", "detail", by="feature")],
        registry=_registry(plot=IndexAdapter()),
    )

    workspace._record_members(
        "source", MarkMembers(cells=("cell-1",), features=("CD4",))
    )

    assert workspace._state["sources"]["source"]["members"] == {"features": ["CD4"]}


def test_source_and_detail_have_distinct_interaction_configuration() -> None:
    table = pd.DataFrame({"x": [1, 2]}, index=["a", "b"])
    workspace = LinkedView(
        table,
        views=[
            view("plot", "source", title="Source"),
            view("plot", "detail", title="A long centered detail title"),
        ],
        links=[link("source", "detail")],
        registry=_registry(plot=IndexAdapter(emits=("row",))),
    )

    source = workspace._component(workspace._view_by_id["source"])
    detail = workspace._component(workspace._view_by_id["detail"])

    assert source.figure.layout.dragmode == "lasso"
    assert source.config.get("displayModeBar", True) is True
    assert detail.figure.layout.dragmode == "zoom"
    assert detail.config["displayModeBar"] is False
    assert detail.config["scrollZoom"] is True
    assert detail.figure.layout.title.x == 0.5
