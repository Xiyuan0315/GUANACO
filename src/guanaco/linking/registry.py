"""Registry of plot renderers that speak the index-link adapter protocol."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from .base import PlotAdapter


class PlotRegistry:
    """Map public plot names to reusable rendering/event adapters."""

    def __init__(self) -> None:
        self._adapters: dict[str, PlotAdapter] = {}

    def register(
        self,
        names: str | Iterable[str],
        adapter: PlotAdapter,
        *,
        replace: bool = False,
    ) -> PlotRegistry:
        raw_names = (names,) if isinstance(names, str) else tuple(names)
        if not raw_names:
            raise ValueError("Register at least one plot name.")
        keys = tuple(str(name).strip() for name in raw_names)
        if any(not key for key in keys):
            raise ValueError("Plot names must be non-empty strings.")
        if len(set(keys)) != len(keys):
            raise ValueError("Plot names must be unique within one registration.")
        conflicts = [name for name in keys if name in self._adapters]
        if conflicts and not replace:
            raise ValueError(f"A plot adapter named {conflicts[0]!r} already exists.")
        for name in keys:
            self._adapters[name] = adapter
        return self

    def get(self, name: str) -> PlotAdapter:
        try:
            return self._adapters[str(name)]
        except KeyError as error:
            available = ", ".join(sorted(self._adapters))
            raise ValueError(
                f"Unknown linked plot type {name!r}. Available: {available}."
            ) from error

    def contract(self, name: str) -> dict[str, Any]:
        adapter = self.get(name)
        return {
            "emits": tuple(sorted(getattr(adapter, "emits", ()))),
            "accepts": tuple(sorted(getattr(adapter, "accepts", ()))),
        }

    def contracts(self) -> dict[str, dict[str, Any]]:
        return {name: self.contract(name) for name in self._adapters}


def default_plot_registry() -> PlotRegistry:
    """Build GUANACO's native and generic Plotly adapter registry."""

    from .native_adapters import (
        CoexpressionAdapter,
        CompositionAdapter,
        EmbeddingAdapter,
        FeatureDistributionAdapter,
        FeatureGroupMatrixAdapter,
        HeatmapAdapter,
        PseudotimeAdapter,
        ViolinAdapter,
        VolcanoAdapter,
    )
    from .table_adapters import NetworkAdapter, TableHeatmapAdapter, TablePlotAdapter

    registry = PlotRegistry()
    registry.register(
        ("embedding", "umap", "pca", "tsne", "spatial"),
        EmbeddingAdapter(),
    )
    registry.register("coexpression", CoexpressionAdapter())
    registry.register("violin", ViolinAdapter())
    registry.register("ridge", FeatureDistributionAdapter("ridge"))
    registry.register("violin_grouped", FeatureDistributionAdapter("split_violin"))
    registry.register("dotplot", FeatureGroupMatrixAdapter("dotplot"))
    registry.register("matrixplot", FeatureGroupMatrixAdapter("matrixplot"))
    registry.register("heatmap", HeatmapAdapter())
    registry.register("stacked_bar", CompositionAdapter())
    registry.register("pseudotime", PseudotimeAdapter())
    registry.register("volcano", VolcanoAdapter())

    table_plot = TablePlotAdapter()
    registry.register(
        (
            "plotly.scatter",
            "plotly.bar",
            "plotly.line",
            "plotly.area",
        ),
        table_plot,
    )
    registry.register("plotly.heatmap", TableHeatmapAdapter())
    registry.register("network", NetworkAdapter())
    return registry


__all__ = ["PlotRegistry", "default_plot_registry"]
