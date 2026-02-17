"""GUANACO package."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("guanaco-viz")
except PackageNotFoundError:
    __version__ = "0+unknown"
