"""Data models, loaders, and runtime dataset registry.

Keep this package initializer lightweight.  Importing a metadata helper such as
``guanaco.data.hdf5_capabilities`` must not initialize Muon and Scanpy.  Public
loader attributes remain available through lazy attribute resolution.
"""

from importlib import import_module


_LOADER_EXPORTS = {
    "DatasetBundle",
    "get_discrete_labels",
    "get_modality_variables",
    "get_ref_track",
    "initialize_data",
    "load_adata",
    "load_config",
    "load_tracks_from_s3",
}

__all__ = [
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
    if name in _LOADER_EXPORTS:
        return getattr(import_module(".loader", __name__), name)
    if name in {"color_config", "datasets"}:
        from . import registry

        return getattr(registry, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
