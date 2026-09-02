from __future__ import annotations

import pytest

from guanaco.linking.model import (
    LinkSpec,
    MarkMembers,
    Selection,
    ViewSpec,
    ViewState,
    link,
    view,
)


def test_view_keeps_data_and_plot_options() -> None:
    spec = view(
        "plotly.scatter",
        "points",
        "external",
        "External points",
        x="score",
        y="effect",
    )

    assert spec == ViewSpec(
        plot="plotly.scatter",
        id="points",
        data="external",
        title="External points",
        options={"x": "score", "y": "effect"},
    )


@pytest.mark.parametrize("field", ["plot", "id"])
def test_view_requires_non_empty_names(field: str) -> None:
    values = {"plot": "scatter", "id": "points", field: "  "}

    with pytest.raises(ValueError, match=field):
        ViewSpec(**values)


def test_view_data_is_a_source_name_not_a_data_object() -> None:
    with pytest.raises(TypeError, match="named data source"):
        view("plotly.scatter", "points", data=object())


@pytest.mark.parametrize(
    ("by", "action", "selection_by", "resolved_action"),
    [
        (None, None, "row", "filter"),
        ("cell", None, "cell", "highlight"),
        ("cell", "filter", "cell", "filter"),
        ("feature", None, "feature", None),
    ],
)
def test_link_resolves_identity_and_display_action(
    by: str | None,
    action: str | None,
    selection_by: str,
    resolved_action: str | None,
) -> None:
    spec = link("source", "target", by=by, action=action)

    assert spec == LinkSpec("source", "target", by, action)
    assert spec.selection_by == selection_by
    assert spec.resolved_action == resolved_action


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"by": "row"}, "`by` must be None"),
        ({"by": "gene"}, "`by` must be None"),
        ({"action": "select"}, "`action` must be None"),
        ({"by": "feature", "action": "filter"}, "Feature links"),
    ],
)
def test_link_rejects_invalid_vocabulary(kwargs: dict[str, str], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        link("source", "target", **kwargs)


def test_selection_and_mark_members_normalize_identity_ids() -> None:
    selection = Selection(by="cell", ids=(7, "cell-8"), source_view="umap")
    members = MarkMembers(
        rows=(1, 2),
        cells=("cell-1", "cell-2"),
        features=("CD3D",),
    )

    assert selection.ids == ("7", "cell-8")
    assert members.project("row") == ("1", "2")
    assert members.project("cell") == ("cell-1", "cell-2")
    assert members.project("feature") == ("CD3D",)
    assert MarkMembers.from_dict(members.to_dict()) == members


def test_view_state_distinguishes_unset_empty_and_highlighted_axes() -> None:
    empty = ViewState.empty()
    active = ViewState(rows=(), cells=("cell-1",), highlight_axes=("cell",))

    assert empty.is_empty
    assert empty.rows is None
    assert not active.is_empty
    assert active.rows == ()
    assert active.features is None
    assert active.is_highlighted("cell")
    assert not active.is_highlighted("row")


@pytest.mark.parametrize(
    "factory",
    [
        lambda: Selection(by="rows", ids=(), source_view="source"),
        lambda: Selection(by="row", ids=(None,), source_view="source"),
        lambda: MarkMembers(rows=("",)),
    ],
)
def test_identity_objects_reject_invalid_axes_and_ids(factory) -> None:
    with pytest.raises(ValueError):
        factory()
