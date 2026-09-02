"""Normalize linked-view inputs to named AnnData or table sources."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, Literal

import pandas as pd


def _is_mudata(value: Any) -> bool:
    mod = getattr(value, "mod", None)
    return hasattr(mod, "keys") and hasattr(mod, "__getitem__")


def _is_anndata(value: Any) -> bool:
    return not _is_mudata(value) and all(
        hasattr(value, name) for name in ("obs", "var", "obs_names", "var_names")
    )


def _name(value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Data source names must be non-empty strings.")
    return value


def _stable_ids(values: Any, label: str) -> pd.Index:
    index = pd.Index(values)
    has_null = (
        index.to_frame(index=False).isna().to_numpy().any()
        if isinstance(index, pd.MultiIndex)
        else index.isna().any()
    )
    if has_null:
        raise ValueError(f"{label} must not contain null IDs.")
    if not index.is_unique:
        raise ValueError(f"{label} must contain unique IDs.")
    result = pd.Index((str(value) for value in index), dtype="object", name=index.name)
    if any(not value.strip() for value in result):
        raise ValueError(f"{label} must not contain empty or whitespace-only IDs.")
    if not result.is_unique:
        raise ValueError(f"{label} must remain unique after conversion to strings.")
    return result


@dataclass(frozen=True, slots=True)
class DataSource:
    """One named AnnData or DataFrame payload and its stable IDs."""

    name: str
    kind: Literal["anndata", "table"]
    data: Any
    row_ids: pd.Index | None = field(init=False, default=None, repr=False)
    cell_ids: pd.Index | None = field(init=False, default=None, repr=False)
    feature_ids: pd.Index | None = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        _name(self.name)
        if self.kind == "table":
            if not isinstance(self.data, pd.DataFrame):
                raise TypeError(
                    f"Table source {self.name!r} must be a pandas DataFrame."
                )
            object.__setattr__(
                self,
                "row_ids",
                _stable_ids(self.data.index, f"Table source {self.name!r} index"),
            )
            return
        if self.kind != "anndata":
            raise ValueError(f"Unsupported data source kind: {self.kind!r}.")
        if not _is_anndata(self.data):
            raise TypeError(f"Source {self.name!r} is not AnnData-like.")
        object.__setattr__(
            self,
            "cell_ids",
            _stable_ids(self.data.obs_names, f"AnnData source {self.name!r} obs_names"),
        )
        object.__setattr__(
            self,
            "feature_ids",
            _stable_ids(self.data.var_names, f"AnnData source {self.name!r} var_names"),
        )


@dataclass(frozen=True, slots=True)
class DataStore:
    """Named sources used by one linked workspace."""

    sources: Mapping[str, DataSource]
    default_source: str | None = None

    def __post_init__(self) -> None:
        sources = dict(self.sources)
        if not sources:
            raise ValueError("A data store requires at least one data source.")
        for name, source in sources.items():
            _name(name)
            if not isinstance(source, DataSource):
                raise TypeError(f"Data source {name!r} must be a DataSource.")
            if source.name != name:
                raise ValueError(
                    f"Data source key {name!r} does not match source name {source.name!r}."
                )
        if self.default_source is not None and self.default_source not in sources:
            raise ValueError(f"Unknown default data source: {self.default_source!r}.")
        object.__setattr__(self, "sources", sources)

    @classmethod
    def from_data(cls, data: Any) -> DataStore:
        if isinstance(data, cls):
            return data
        if _is_anndata(data) or isinstance(data, pd.DataFrame):
            kind = "anndata" if _is_anndata(data) else "table"
            source = DataSource("main", kind, data)
            return cls({"main": source}, "main")
        if _is_mudata(data):
            data = data.mod
        if not isinstance(data, Mapping):
            raise TypeError(
                "Linked-view data must be AnnData, MuData, a pandas DataFrame, or "
                "a mapping of named AnnData/DataFrame sources."
            )
        sources: dict[str, DataSource] = {}
        for raw_name, value in data.items():
            name = _name(raw_name)
            if name in sources:
                raise ValueError(f"Duplicate data source name: {name!r}.")
            if _is_anndata(value):
                kind = "anndata"
            elif isinstance(value, pd.DataFrame):
                kind = "table"
            else:
                raise TypeError(
                    f"Named source {name!r} must be AnnData-like or a pandas DataFrame."
                )
            sources[name] = DataSource(name, kind, value)
        if not sources:
            raise ValueError("A data mapping requires at least one source.")
        default = next(iter(sources)) if len(sources) == 1 else None
        return cls(sources, default)

    def source(self, name: str | None = None) -> DataSource:
        if name is None:
            name = self.default_source
            if name is None:
                available = ", ".join(repr(item) for item in self.sources)
                raise ValueError(
                    "Linked-view data has multiple sources; specify `data=` on the "
                    f"view. Available sources: {available}."
                )
        try:
            return self.sources[name]
        except KeyError as error:
            raise ValueError(f"Unknown linked-view data source: {name!r}.") from error
