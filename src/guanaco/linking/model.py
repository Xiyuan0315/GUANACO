"""Public, plot-agnostic types for index-linked views."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Literal, Mapping


SelectionBy = Literal["row", "cell", "feature"]
LinkBy = Literal["cell", "feature"]
ActionName = Literal["highlight", "filter"]

_AXES = frozenset({"row", "cell", "feature"})
_LINK_AXES = frozenset({"cell", "feature"})
_ACTIONS = frozenset({"highlight", "filter"})


def _name(value: object, label: str) -> str:
    result = str(value)
    if not result.strip():
        raise ValueError(f"`{label}` must be a non-empty string.")
    return result


def _ids(values: Iterable[object], label: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"`{label}` must be an iterable of IDs, not a string.")
    result: list[str] = []
    for value in values:
        if value is None:
            raise ValueError(f"`{label}` cannot contain None.")
        item = str(value)
        if not item.strip():
            raise ValueError(f"`{label}` cannot contain an empty ID.")
        result.append(item)
    return tuple(result)


def _axis(value: object) -> SelectionBy:
    if value not in _AXES:
        raise ValueError(
            f"`by` must be one of 'row', 'cell', or 'feature'; received {value!r}."
        )
    return value  # type: ignore[return-value]


@dataclass(frozen=True)
class ViewSpec:
    """One plot in a linked workspace."""

    plot: str
    id: str
    data: str | None = None
    title: str | None = None
    options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "plot", _name(self.plot, "plot"))
        object.__setattr__(self, "id", _name(self.id, "id"))
        if self.data is not None:
            if not isinstance(self.data, str):
                raise TypeError("`data` must be a named data source, not an object.")
            object.__setattr__(self, "data", _name(self.data, "data"))
        object.__setattr__(self, "options", dict(self.options))


@dataclass(frozen=True)
class LinkSpec:
    """A directed link that carries IDs, optionally through a table key."""

    source: str
    target: str
    by: LinkBy | None = None
    action: ActionName | None = None
    key: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "source", _name(self.source, "source"))
        object.__setattr__(self, "target", _name(self.target, "target"))
        if self.by is not None and self.by not in _LINK_AXES:
            raise ValueError(
                f"`by` must be None, 'cell', or 'feature'; received {self.by!r}."
            )
        if self.action is not None and self.action not in _ACTIONS:
            raise ValueError(
                "`action` must be None, 'highlight', or 'filter'; "
                f"received {self.action!r}."
            )
        if self.key is not None and (
            not isinstance(self.key, str) or not self.key.strip()
        ):
            raise ValueError("`key` must be a non-empty table column name.")
        if self.by == "feature" and self.action is not None:
            raise ValueError(
                "Feature links do not accept `action`; the selected feature "
                "becomes the target's feature context."
            )

    @property
    def selection_by(self) -> SelectionBy:
        return self.by or "row"

    @property
    def resolved_action(self) -> ActionName | None:
        if self.action is not None:
            return self.action
        return {"row": "filter", "cell": "highlight", "feature": None}[
            self.selection_by
        ]


@dataclass(frozen=True)
class Selection:
    """Stable IDs emitted by a source view."""

    by: SelectionBy
    ids: tuple[str, ...]
    source_view: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "by", _axis(self.by))
        object.__setattr__(self, "ids", _ids(self.ids, "ids"))
        object.__setattr__(self, "source_view", _name(self.source_view, "source_view"))


@dataclass(frozen=True)
class MarkMembers:
    """The row, cell, and feature IDs represented by one or more marks."""

    rows: tuple[str, ...] = ()
    cells: tuple[str, ...] = ()
    features: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("rows", "cells", "features"):
            object.__setattr__(self, name, _ids(getattr(self, name), name))

    def project(self, by: SelectionBy) -> tuple[str, ...]:
        return getattr(self, f"{_axis(by)}s")

    def to_dict(self) -> dict[str, list[str]]:
        return {
            name: list(values)
            for name in ("rows", "cells", "features")
            if (values := getattr(self, name))
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> MarkMembers:
        return cls(
            **{
                name: tuple(value.get(name) or ())
                for name in ("rows", "cells", "features")
            }
        )


@dataclass(frozen=True)
class ViewState:
    """Selections applied to a target; ``None`` means unconstrained."""

    rows: tuple[str, ...] | None = None
    cells: tuple[str, ...] | None = None
    features: tuple[str, ...] | None = None
    highlight_axes: tuple[SelectionBy, ...] = ()

    def __post_init__(self) -> None:
        for name in ("rows", "cells", "features"):
            values = getattr(self, name)
            if values is not None:
                object.__setattr__(self, name, _ids(values, name))
        object.__setattr__(
            self,
            "highlight_axes",
            tuple(dict.fromkeys(_axis(value) for value in self.highlight_axes)),
        )

    @classmethod
    def empty(cls) -> ViewState:
        return cls()

    @property
    def is_empty(self) -> bool:
        return all(
            getattr(self, name) is None for name in ("rows", "cells", "features")
        )

    def is_highlighted(self, by: SelectionBy) -> bool:
        return _axis(by) in self.highlight_axes


def view(
    plot: str,
    id: str,
    data: str | None = None,
    title: str | None = None,
    **options: Any,
) -> ViewSpec:
    return ViewSpec(plot=plot, id=id, data=data, title=title, options=options)


def link(
    source: str,
    target: str,
    *,
    by: LinkBy | None = None,
    action: ActionName | None = None,
    key: str | None = None,
) -> LinkSpec:
    return LinkSpec(source=source, target=target, by=by, action=action, key=key)
