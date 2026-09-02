from guanaco.linking.engine import reduce_members, state_for_target
from guanaco.linking.model import MarkMembers, link


def _route(state, source, members):
    return reduce_members(state, source_view=source, members=members)


def test_one_to_one_routes_rows_as_a_filter() -> None:
    links = [link("network", "pairs")]

    state = _route(
        None,
        "network",
        MarkMembers(rows=("id1", "id2")),
    )
    target = state_for_target(state, "pairs", links)

    assert target.rows == ("id1", "id2")
    assert not target.is_highlighted("row")


def test_one_to_many_projects_cell_and_feature_from_the_same_mark() -> None:
    links = [
        link("matrix", "cells", by="cell"),
        link("matrix", "expression", by="feature"),
    ]

    state = _route(
        None,
        "matrix",
        MarkMembers(cells=("cell-1", "cell-2"), features=("IL7R",)),
    )

    cells = state_for_target(state, "cells", links)
    expression = state_for_target(state, "expression", links)
    assert cells.cells == ("cell-1", "cell-2")
    assert cells.is_highlighted("cell")
    assert expression.features == ("IL7R",)


def test_many_to_one_uses_latest_source_then_falls_back_when_cleared() -> None:
    links = [link("source_a", "detail"), link("source_b", "detail")]
    state = _route(None, "source_a", MarkMembers(rows=("row-a",)))
    state = _route(state, "source_b", MarkMembers(rows=("row-b",)))

    assert state_for_target(state, "detail", links).rows == ("row-b",)

    state = _route(state, "source_b", None)

    assert state_for_target(state, "detail", links).rows == ("row-a",)


def test_many_to_many_routes_each_source_to_each_terminal_detail() -> None:
    links = [
        link("source_a", "detail_1"),
        link("source_a", "detail_2"),
        link("source_b", "detail_1"),
        link("source_b", "detail_2"),
    ]
    state = _route(None, "source_a", MarkMembers(rows=("row-a",)))
    state = _route(state, "source_b", MarkMembers(rows=("row-b",)))

    assert state_for_target(state, "detail_1", links).rows == ("row-b",)
    assert state_for_target(state, "detail_2", links).rows == ("row-b",)


def test_filter_and_highlight_can_coexist_on_different_identity_axes() -> None:
    links = [
        link("rows", "detail"),
        link("cells", "detail", by="cell"),
    ]
    state = _route(None, "rows", MarkMembers(rows=("row-1", "row-2")))
    state = _route(
        state,
        "cells",
        MarkMembers(cells=("cell-2", "cell-3")),
    )

    target = state_for_target(state, "detail", links)
    assert target.rows == ("row-1", "row-2")
    assert target.cells == ("cell-2", "cell-3")
    assert not target.is_highlighted("row")
    assert target.is_highlighted("cell")


def test_reducer_records_source_members_and_genuine_clear() -> None:
    links = [link("umap", "detail", by="cell")]
    state = _route(
        None,
        "umap",
        MarkMembers(cells=("cell-2", "cell-4")),
    )

    assert state["sources"]["umap"]["members"]["cells"] == ["cell-2", "cell-4"]
    assert state_for_target(state, "detail", links).cells == ("cell-2", "cell-4")

    state = _route(state, "umap", None)

    assert "umap" not in state["sources"]
    assert state_for_target(state, "detail", links).is_empty


def test_reducer_never_guesses_that_rows_are_cells() -> None:
    links = [link("scores", "umap", by="cell")]
    state = _route(
        None,
        "scores",
        MarkMembers(rows=("cell-1",)),
    )

    assert state_for_target(state, "umap", links).cells is None
    assert state["sources"]["scores"]["members"]["rows"] == ["cell-1"]
