"""Small adapter protocol shared by the index-based linking runtime."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping
from types import MappingProxyType
from typing import Any, TypeAlias

import plotly.graph_objects as go
from dash.development.base_component import Component

from .data import DataStore
from .model import MarkMembers, SelectionBy, ViewSpec, ViewState


Rendered: TypeAlias = go.Figure | Component


class PlotAdapter(ABC):
    """Translate browser marks to indices and indices back to one plot."""

    events: Mapping[str, str] = MappingProxyType({})
    emits: frozenset[str] = frozenset()
    accepts: frozenset[str] = frozenset()

    def identity_ids(
        self,
        by: SelectionBy,
        spec: ViewSpec,
        store: DataStore,
    ) -> tuple[str, ...] | None:
        """Optionally refine the IDs this plot can consume on one axis.

        Most adapters use the container's native obs/var/index domain. A plot
        may widen or narrow that domain when its visual grammar supports extra
        stable names, such as embedding color from an AnnData ``obs`` column.
        """

        del by, spec, store
        return None

    @abstractmethod
    def validate(self, spec: ViewSpec, store: DataStore) -> None:
        """Validate one view before the Dash application is compiled."""

    @abstractmethod
    def render(
        self,
        spec: ViewSpec,
        store: DataStore,
        state: ViewState | None = None,
        *,
        component_id: str | None = None,
    ) -> Rendered:
        """Render the view for an optional index selection state."""

    @abstractmethod
    def decode_event(
        self,
        event: str,
        payload: Mapping[str, Any] | None,
        spec: ViewSpec,
        store: DataStore,
    ) -> MarkMembers | None:
        """Decode one browser payload into atomic mark memberships."""


__all__ = ["PlotAdapter", "Rendered"]
