"""Data models, loaders, and runtime dataset registry."""

from .loader import (
    DEFAULT_COLORS,
    DatasetBundle,
    get_discrete_labels,
    get_modality_variables,
    get_ref_track,
    initialize_data,
    load_adata,
    load_config,
    load_tracks_from_s3,
)
__all__ = [
    "DEFAULT_COLORS",
    "DatasetBundle",
    "color_config",
    "datasets",
    "get_discrete_labels",
    "get_modality_variables",
    "get_ref_track",
    "initialize_data",
    "load_adata",
    "load_config",
    "load_tracks_from_s3",
]


def __getattr__(name):
    if name in {"color_config", "datasets"}:
        from . import registry

        return getattr(registry, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
