"""GUANACO package."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("guanaco-viz")
except PackageNotFoundError:
    __version__ = "0+unknown"

# Scanpy-style plotting namespace for notebooks: ``import guanaco as gc; gc.pl.umap(...)``.
from guanaco import widget as pl  # noqa: E402

__all__ = ["pl", "__version__"]
