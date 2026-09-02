"""Dash runtime for directed index-linked views.

Plots own rendering and mark decoding. This module only validates shared
indices, routes source events, and redraws terminal detail views. A table key
may map atomic table rows to native AnnData cell/feature IDs for mixed links.
"""

from __future__ import annotations

import re
import threading
import uuid
from collections.abc import Sequence
from dataclasses import replace
from html import escape as escape_html
from typing import Any

import pandas as pd
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, ctx, dcc, html, no_update
from dash.development.base_component import Component

from guanaco.utils.plot_style import apply_guanaco_figure_style

from .data import DataSource, DataStore
from .engine import reduce_members, state_for_target
from .model import LinkSpec, MarkMembers, Selection, SelectionBy, ViewSpec, ViewState
from .registry import PlotRegistry, default_plot_registry


_SOURCE_CONFIG = {
    "responsive": True,
    "displaylogo": False,
    "scrollZoom": False,
    "doubleClick": "reset",
    "modeBarButtonsToRemove": ["autoScale2d", "zoomIn2d", "zoomOut2d"],
    "toImageButtonOptions": {"format": "svg", "filename": "linked_view"},
}
_DETAIL_CONFIG = {
    "responsive": True,
    "displaylogo": False,
    "displayModeBar": False,
    "editable": False,
    "scrollZoom": True,
    "doubleClick": "reset",
}


def _safe_prefix(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "-", str(value)).strip("-") or "linked"


def _merge_members(values: Sequence[MarkMembers]) -> MarkMembers:
    def union(name: str) -> tuple[str, ...]:
        return tuple(dict.fromkeys(x for value in values for x in getattr(value, name)))

    return MarkMembers(
        rows=union("rows"), cells=union("cells"), features=union("features")
    )


def _coalesce_events(
    events: Sequence[tuple[str, MarkMembers | None]],
) -> list[tuple[str, MarkMembers | None]]:
    """Merge simultaneous Plotly click/select properties for each source."""

    grouped: dict[str, list[MarkMembers | None]] = {}
    for source, members in events:
        grouped.setdefault(source, []).append(members)
    return [
        (source, _merge_members(nonempty) if nonempty else None)
        for source, values in grouped.items()
        for nonempty in ([value for value in values if value is not None],)
    ]


class LinkedView:
    """Plots connected by stable DataFrame, ``obs_names``, or ``var_names`` IDs."""

    def __init__(
        self,
        data: Any,
        *,
        views: Sequence[ViewSpec],
        links: Sequence[LinkSpec],
        layout: str = "overview-detail",
        registry: PlotRegistry | None = None,
        prefix: str | None = None,
        title: str = "GUANACO Linked View",
    ) -> None:
        self.registry = registry or default_plot_registry()
        self.store = DataStore.from_data(data)
        self.views = tuple(views)
        self._view_by_id = {spec.id: spec for spec in self.views}
        self._table_views = frozenset(
            spec.id for spec in self.views if self._source(spec).kind == "table"
        )
        self.links = tuple(links)
        self.layout_mode = str(layout)
        self.prefix = _safe_prefix(prefix or f"linked-{uuid.uuid4().hex[:8]}")
        self.title = str(title)
        self._registered_apps: set[int] = set()
        self._component_kinds: dict[str, str] = {}
        self._marimo_server: dict[str, Any] | None = None
        self._state_lock = threading.RLock()
        self._state: dict[str, Any] = {"sources": {}}
        self._validate()
        self._row_key_maps = self._build_row_key_maps()

    def _source(self, spec: ViewSpec) -> DataSource:
        return self.store.source(spec.data)

    def _build_row_key_map(
        self, spec: ViewSpec, key: str | None
    ) -> dict[str, str]:
        """Map atomic table row IDs to one link's optional logical key."""

        source = self._source(spec)
        if source.kind != "table":
            return {}
        raw_ids = tuple(str(value) for value in source.row_ids)
        if key is None:
            return dict(zip(raw_ids, raw_ids, strict=True))
        if key not in source.data.columns:
            raise ValueError(
                f"View {spec.id!r} is missing link-key column {key!r}."
            )
        values = source.data[key]
        if values.isna().any():
            raise ValueError(
                f"View {spec.id!r} link-key column {key!r} cannot contain nulls."
            )
        keys = tuple(str(value) for value in values)
        if any(not value.strip() for value in keys):
            raise ValueError(
                f"View {spec.id!r} link-key column {key!r} cannot contain empty IDs."
            )
        return dict(zip(raw_ids, keys, strict=True))

    def _build_row_key_maps(
        self,
    ) -> dict[tuple[str, str | None], dict[str, str]]:
        maps: dict[tuple[str, str | None], dict[str, str]] = {}
        for link in self.links:
            for view_id in (link.source, link.target):
                spec = self._view_by_id[view_id]
                maps.setdefault(
                    (view_id, link.key), self._build_row_key_map(spec, link.key)
                )
        return maps

    def _validate(self) -> None:
        if not self.views:
            raise ValueError("A linked view requires at least one view.")
        if len(self._view_by_id) != len(self.views):
            raise ValueError("Linked view IDs must be unique.")
        if self.layout_mode not in {"overview-detail", "grid", "column"}:
            raise ValueError("`layout` must be 'overview-detail', 'grid', or 'column'.")
        for spec in self.views:
            self.registry.get(spec.plot).validate(spec, self.store)

        self.links = tuple(self._resolve_link(link) for link in self.links)
        link_keys = [
            (link.source, link.target, link.selection_by) for link in self.links
        ]
        if len(link_keys) != len(set(link_keys)):
            raise ValueError("Duplicate links must be removed from a linked view.")
        overlap = {link.source for link in self.links} & {
            link.target for link in self.links
        }
        if overlap:
            name = sorted(overlap)[0]
            raise ValueError(
                f"View {name!r} cannot be both overview and detail; links are "
                "one-way and acyclic: overview → terminal detail."
            )
        for link in self.links:
            self._validate_link(link)

    def _resolve_link(self, link: LinkSpec) -> LinkSpec:
        if link.source not in self._view_by_id:
            raise ValueError(f"Unknown link source view: {link.source!r}.")
        if link.target not in self._view_by_id:
            raise ValueError(f"Unknown link target view: {link.target!r}.")
        if link.source == link.target:
            raise ValueError("A linked view cannot link a plot to itself.")
        if link.by is not None:
            return link
        kinds = {
            self._source(self._view_by_id[link.source]).kind,
            self._source(self._view_by_id[link.target]).kind,
        }
        if kinds == {"anndata"}:
            return replace(link, by="cell")
        if kinds == {"table"}:
            return link
        raise ValueError(
            f"Link {link.source!r} → {link.target!r} crosses table and AnnData; "
            "declare `by='cell'` or `by='feature'`."
        )

    def _identity_ids(
        self, spec: ViewSpec, by: SelectionBy, key: str | None = None
    ) -> pd.Index | None:
        adapter = self.registry.get(spec.plot)
        custom = adapter.identity_ids(by, spec, self.store)
        if custom is not None:
            return pd.Index(custom, dtype="object")
        source = self._source(spec)
        if source.kind == "table":
            return pd.Index(
                tuple(
                    dict.fromkeys(self._build_row_key_map(spec, key).values())
                ),
                dtype="object",
            )
        return source.cell_ids if by == "cell" else source.feature_ids

    def _validate_link(self, link: LinkSpec) -> None:
        source_spec = self._view_by_id[link.source]
        target_spec = self._view_by_id[link.target]
        source = self._source(source_spec)
        target = self._source(target_spec)
        source_adapter = self.registry.get(source_spec.plot)
        target_adapter = self.registry.get(target_spec.plot)
        by = link.selection_by

        if by == "row" and (source.kind != "table" or target.kind != "table"):
            raise ValueError("Row links require two DataFrame-backed views.")
        if link.key is not None:
            table_endpoints = [source, target]
            if not any(item.kind == "table" for item in table_endpoints):
                raise ValueError(
                    "`key` requires at least one table-backed view; use native "
                    "AnnData obs_names or var_names for an AnnData-only link."
                )
            for endpoint in (source_spec, target_spec):
                if self._source(endpoint).kind == "table":
                    self._build_row_key_map(endpoint, link.key)
        if (
            source.kind == "anndata"
            and by == "cell"
            and source_spec.options.get("render_backend") == "datashader"
        ):
            raise ValueError(
                "A selectable cell source requires scattergl, not datashader."
            )
        if not source_adapter.events:
            raise ValueError(f"View {link.source!r} has no selectable browser event.")
        if by not in source_adapter.emits:
            choices = ", ".join(sorted(source_adapter.emits))
            raise ValueError(
                f"View {link.source!r} emits {choices or 'no'} IDs, not {by!r} IDs."
            )
        if by not in target_adapter.accepts:
            choices = ", ".join(sorted(target_adapter.accepts)) or "no identity axes"
            raise ValueError(
                f"View {link.target!r} cannot consume {by!r} IDs; it accepts {choices}."
            )
        source_ids = self._identity_ids(source_spec, by, link.key)
        target_ids = self._identity_ids(target_spec, by, link.key)
        if source_ids is None or target_ids is None:
            raise ValueError(
                f"Link {link.source!r} → {link.target!r} has no {by} index."
            )
        if (
            len(source_ids)
            and len(target_ids)
            and not len(source_ids.intersection(target_ids, sort=False))
        ):
            raise ValueError(
                f"Link {link.source!r} → {link.target!r} shares no {by} IDs."
            )

    def _selection_axis(
        self, view_id: str, by: SelectionBy | None = None
    ) -> SelectionBy:
        axes = tuple(
            dict.fromkeys(
                link.selection_by for link in self.links if link.source == view_id
            )
        )
        if by is not None:
            if by not in {"row", "cell", "feature"}:
                raise ValueError("`by` must be 'row', 'cell', or 'feature'.")
            if by not in axes:
                raise ValueError(f"View {view_id!r} has no outgoing {by!r} link.")
            return by
        if len(axes) == 1:
            return axes[0]
        if not axes:
            raise ValueError(f"View {view_id!r} is not a link source.")
        raise ValueError(f"View {view_id!r} emits several ID axes; pass `by=`.")

    def _record_members(self, view_id: str, members: MarkMembers | None) -> None:
        """Update the same reducer state used by Dash (also useful in tests)."""

        with self._state_lock:
            self._state = reduce_members(
                self._state,
                source_view=view_id,
                members=self._linked_members(view_id, members),
            )

    def _linked_members(
        self, view_id: str, members: MarkMembers | None
    ) -> MarkMembers | None:
        """Keep only identities consumed by this source's outgoing links."""

        if members is None:
            return None
        axes = {link.selection_by for link in self.links if link.source == view_id}
        return MarkMembers(
            rows=members.rows if view_id in self._table_views else (),
            cells=members.cells if "cell" in axes else (),
            features=members.features if "feature" in axes else (),
        )

    def get_selection(
        self, view_id: str, *, by: SelectionBy | None = None
    ) -> Selection | None:
        if view_id not in self._view_by_id:
            raise ValueError(f"Unknown linked view: {view_id!r}.")
        axis = self._selection_axis(view_id, by)
        with self._state_lock:
            snapshot = (self._state.get("sources") or {}).get(view_id)
        if not snapshot:
            return None
        members = MarkMembers.from_dict(snapshot.get("members") or {})
        ids = (
            members.rows or members.project(axis)
            if view_id in self._table_views
            else members.project(axis)
        )
        if view_id in self._table_views:
            keys = {
                link.key
                for link in self.links
                if link.source == view_id and link.selection_by == axis
            }
            if len(keys) == 1:
                row_keys = self._row_key_maps.get((view_id, next(iter(keys))), {})
                ids = tuple(dict.fromkeys(row_keys.get(value, value) for value in ids))
        return Selection(axis, ids, view_id) if ids else None

    def describe(self) -> dict[str, Any]:
        return {
            "views": [
                {
                    "id": spec.id,
                    "plot": spec.plot,
                    "data": self._source(spec).name,
                    "data_kind": self._source(spec).kind,
                }
                for spec in self.views
            ],
            "links": [self._describe_link(link) for link in self.links],
            "layout": self.layout_mode,
        }

    @staticmethod
    def _describe_link(link: LinkSpec) -> dict[str, Any]:
        result = {
            "source": link.source,
            "target": link.target,
            "by": link.selection_by,
            "action": link.resolved_action,
        }
        if link.key is not None:
            result["key"] = link.key
        return result

    def _interactive_id(self, view_id: str) -> str:
        return f"{self.prefix}-view-{view_id}"

    def _container_id(self, view_id: str) -> str:
        return f"{self.prefix}-container-{view_id}"

    @property
    def _state_id(self) -> str:
        return f"{self.prefix}-state"

    def _is_source(self, view_id: str) -> bool:
        return any(link.source == view_id for link in self.links)

    def _render(self, spec: ViewSpec, state: ViewState | None = None):
        adapter = self.registry.get(spec.plot)
        render_spec = (
            replace(spec, options={**spec.options, "_link_source": True})
            if self._is_source(spec.id)
            else spec
        )
        rendered = adapter.render(
            render_spec,
            self.store,
            state or ViewState.empty(),
            component_id=self._interactive_id(spec.id),
        )
        if isinstance(rendered, go.Figure):
            figure = apply_guanaco_figure_style(rendered)
            if self._is_source(spec.id):
                selectable = "select" in adapter.events
                figure.update_layout(
                    dragmode="lasso" if selectable else "zoom",
                    clickmode="event+select" if selectable else "event",
                    uirevision=f"{spec.id}:axes",
                )
            else:
                figure.update_layout(
                    dragmode="zoom", clickmode="event", uirevision=f"{spec.id}:axes"
                )
            return figure
        if not isinstance(rendered, Component):
            raise TypeError(f"Unsupported render result: {type(rendered).__name__}.")
        return rendered

    def _component(self, spec: ViewSpec, state: ViewState | None = None):
        rendered = self._render(spec, state)
        if isinstance(rendered, go.Figure):
            self._component_kinds[spec.id] = "figure"
            height = spec.options.get("height", "600px")
            height = f"{height}px" if isinstance(height, (int, float)) else str(height)
            return dcc.Graph(
                id=self._interactive_id(spec.id),
                figure=rendered,
                config=dict(
                    _SOURCE_CONFIG if self._is_source(spec.id) else _DETAIL_CONFIG
                ),
                responsive=True,
                style={"height": height, "width": "100%"},
            )
        self._component_kinds[spec.id] = "component"
        return rendered

    def component(self):
        cards = [
            html.Section(
                html.Div(self._component(spec), id=self._container_id(spec.id)),
                style={
                    "background": "white",
                    "border": "1px solid #E5E7EB",
                    "borderRadius": "12px",
                    "minWidth": 0,
                    "overflow": "hidden",
                },
            )
            for spec in self.views
        ]
        columns = "minmax(0, 1fr)"
        if self.layout_mode == "overview-detail" and len(cards) == 2:
            columns = "minmax(0, 1.15fr) minmax(0, 0.85fr)"
        elif self.layout_mode == "grid":
            columns = "repeat(auto-fit, minmax(min(100%, 460px), 1fr))"
        return html.Div(
            [
                dcc.Store(id=self._state_id, data=self._state),
                html.H2(
                    self.title,
                    style={"fontSize": "1.35rem", "margin": "0 0 16px"},
                ),
                html.Div(
                    cards,
                    style={
                        "display": "grid",
                        "gap": "16px",
                        "gridTemplateColumns": columns,
                    },
                ),
            ],
            style={
                "background": "#F8FAFC",
                "boxSizing": "border-box",
                "fontFamily": "Inter, ui-sans-serif, system-ui, sans-serif",
                "minHeight": "100vh",
                "padding": "22px",
            },
        )

    def register(self, app: Dash) -> None:
        """Register one event router plus one redraw callback per detail view."""

        if id(app) in self._registered_apps:
            return
        self._registered_apps.add(id(app))
        inputs: list[Input] = []
        lookup: dict[str, tuple[str, str]] = {}
        for source_id in dict.fromkeys(link.source for link in self.links):
            spec = self._view_by_id[source_id]
            adapter = self.registry.get(spec.plot)
            for event, prop in adapter.events.items():
                component_id = self._interactive_id(source_id)
                inputs.append(Input(component_id, prop))
                lookup[f"{component_id}.{prop}"] = (source_id, event)

        if inputs:

            @app.callback(
                Output(self._state_id, "data"),
                inputs,
                State(self._state_id, "data"),
                prevent_initial_call=True,
            )
            def route_events(*values):
                current = values[-1] if values else None
                decoded: list[tuple[str, MarkMembers | None]] = []
                for triggered in ctx.triggered or ():
                    binding = lookup.get(str(triggered.get("prop_id", "")))
                    if binding is None:
                        continue
                    source_id, event = binding
                    spec = self._view_by_id[source_id]
                    members = self.registry.get(spec.plot).decode_event(
                        event, triggered.get("value"), spec, self.store
                    )
                    decoded.append((source_id, members))
                if not decoded:
                    return no_update
                next_state = dict(current or self._state)
                for source_id, members in _coalesce_events(decoded):
                    next_state = reduce_members(
                        next_state,
                        source_view=source_id,
                        members=self._linked_members(source_id, members),
                    )
                with self._state_lock:
                    self._state = next_state
                return next_state

        for target_id in dict.fromkeys(link.target for link in self.links):
            spec = self._view_by_id[target_id]
            kind = self._component_kinds.get(target_id)
            if kind is None:
                kind = (
                    "figure"
                    if isinstance(self._render(spec), go.Figure)
                    else "component"
                )
                self._component_kinds[target_id] = kind

            def register_target(spec=spec, kind=kind):
                output = (
                    Output(self._interactive_id(spec.id), "figure")
                    if kind == "figure"
                    else Output(self._container_id(spec.id), "children")
                )

                @app.callback(
                    output,
                    Input(self._state_id, "data"),
                    prevent_initial_call=True,
                )
                def redraw(state):
                    target_state = state_for_target(
                        state,
                        spec.id,
                        self.links,
                        self._table_views,
                        self._row_key_maps,
                    )
                    return (
                        self._render(spec, target_state)
                        if kind == "figure"
                        else self._component(spec, target_state)
                    )

            register_target()

    def create_app(self, *, name: str = __name__) -> Dash:
        app = Dash(name, suppress_callback_exceptions=True)
        app.title = self.title
        app.layout = self.component()
        self.register(app)
        return app

    def show(self, **run_kwargs: Any) -> LinkedView:
        self.create_app().run(**run_kwargs)
        return self

    def show_jupyter(
        self,
        *,
        mode: str = "inline",
        width: str = "100%",
        height: int = 700,
        port: int = 8050,
        **run_kwargs: Any,
    ) -> LinkedView:
        run_kwargs.setdefault("debug", False)
        self.create_app().run(
            port=port,
            jupyter_mode=mode,
            jupyter_width=width,
            jupyter_height=height,
            **run_kwargs,
        )
        return self

    def show_marimo(
        self,
        *,
        host: str = "127.0.0.1",
        port: int | None = None,
        width: str = "100%",
        height: int | str = 720,
        server_url: str | None = None,
    ):
        try:
            import marimo as mo
        except ImportError as error:
            raise ImportError(
                "Marimo display requires `pip install guanaco-viz[notebook]`."
            ) from error
        requested_port = int(port or 0)
        active = self._marimo_server
        if active is None or not active["thread"].is_alive():
            from werkzeug.serving import make_server

            server = make_server(
                host, requested_port, self.create_app().server, threaded=True
            )
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            active = self._marimo_server = {
                "server": server,
                "thread": thread,
                "port": int(server.server_port),
            }
        elif requested_port and requested_port != active["port"]:
            raise RuntimeError(
                f"Already served on port {active['port']}; call `stop_marimo()` first."
            )
        display_host = "127.0.0.1" if host in {"0.0.0.0", "::"} else host
        url = server_url or f"http://{display_host}:{active['port']}/"
        height_css = f"{height}px" if isinstance(height, int) else str(height)
        return mo.Html(
            f'<iframe src="{escape_html(url, quote=True)}" '
            f'title="{escape_html(self.title, quote=True)}" '
            f'style="width:{escape_html(str(width), quote=True)};'
            f'height:{escape_html(height_css, quote=True)};border:0"></iframe>'
        )

    def stop_marimo(self) -> None:
        active = self._marimo_server
        if active is not None:
            active["server"].shutdown()
            active["server"].server_close()
            active["thread"].join(timeout=2)
            self._marimo_server = None


def linked_view(
    data: Any,
    *,
    views: Sequence[ViewSpec],
    links: Sequence[LinkSpec],
    layout: str = "overview-detail",
    registry: PlotRegistry | None = None,
    prefix: str | None = None,
    title: str = "GUANACO Linked View",
) -> LinkedView:
    return LinkedView(
        data,
        views=views,
        links=links,
        layout=layout,
        registry=registry,
        prefix=prefix,
        title=title,
    )


__all__ = ["LinkedView", "linked_view"]
