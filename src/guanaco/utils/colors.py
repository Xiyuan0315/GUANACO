"""Color palette registry for matrix plots.

Continuous colormaps:
- Plotly named color scales are restricted to the approved perceptually uniform
  set in `PLOTLY_PERCEPTUALLY_UNIFORM`.
- Colorcet continuous palettes use `cc:<name>` and are limited to linear and
  diverging original names.
- Crameri scientific colormaps use `cmc:<name>`, for example `cmc:batlow`.
- Dropdown labels use `source/type_name`, for example `plotly/linear_viridis`,
  `cc/diverging_bwr_40_95_c42`, and `cmc/diverging_vik`.
- Dropdown options include reversed variants with `_r` suffixes. Overview
  helpers show only original colormaps.

Discrete palettes:
- Custom colorblind-friendly palettes are defined in `cvd_color.json` and keep
  their configured names.
- Plotly qualitative palettes are added as `plotly/<name>`.
- Colorcet Glasbey palettes are added as `cc/<name>`.
- Crameri categorical palettes are added as `cmc/<name>`.
- Dynamic Glasbey palettes are available as `glasbey/default` and
  `glasbey/safe`.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Literal

import plotly.express as px

COLORBLIND_PALETTE_PATH = Path(__file__).with_name("cvd_color.json")
DEFAULT_CONTINUOUS_COLORMAP = "Viridis"
COLORCET_PREFIX = "cc:"
CMCRAMERI_PREFIX = "cmc:"
DYNAMIC_DISCRETE_PALETTES = {
    "glasbey/default": {"colorblind_safe": False},
    "glasbey/safe": {"colorblind_safe": True},
}
CRAMERI_SEQUENTIAL = [
    "batlow",
    "batlowK",
    "batlowW",
    "acton",
    "bamako",
    "bilbao",
    "buda",
    "davos",
    "devon",
    "grayC",
    "hawaii",
    "imola",
    "lajolla",
    "lapaz",
    "lipari",
    "navia",
    "nuuk",
    "oslo",
    "tokyo",
    "turku",
    "glasgow",
]
CRAMERI_DIVERGING = [
    "vik",
    "roma",
    "broc",
    "cork",
    "bam",
    "berlin",
    "lisbon",
    "tofino",
    "vanimo",
    "managua",
]
CRAMERI_CATEGORICAL = [
    "batlowS",
    "devonS",
    "davosS",
    "osloS",
    "lapazS",
    "actonS",
    "lajollaS",
    "bilbaoS",
    "grayCS",
    "tokyoS",
    "turkuS",
    "bamakoS",
    "nuukS",
    "hawaiiS",
    "budaS",
    "imolaS",
]
PLOTLY_PERCEPTUALLY_UNIFORM = {
    "sequential": {
        "viridis": "viridis",
        "plasma": "plasma",
        "inferno": "inferno",
        "magma": "magma",
        "cividis": "cividis",
        "thermal": "thermal",
        "haline": "haline",
        "solar": "solar",
        "ice": "ice",
        "deep": "deep",
        "dense": "dense",
        "matter": "matter",
        "speed": "speed",
        "amp": "amp",
        "tempo": "tempo",
    },
    "diverging": {
        "balance": "balance",
        "curl": "curl",
        "delta": "delta",
        "rdbu": "rdbu",
        "brbg": "brbg",
        "puor": "puor",
        "piyg": "piyg",
        "prgn": "prgn",
    },
}


@lru_cache(maxsize=1)
def discrete_palette_registry() -> dict[str, list[str]]:
    with open(COLORBLIND_PALETTE_PATH, "r") as f:
        palette_json = json.load(f)

    palettes = dict(palette_json.get("color_palettes", {}))
    palettes.update(_plotly_qualitative_palettes(palettes))
    palettes.update(_colorcet_discrete_palettes(palettes))
    palettes.update(_cmcrameri_discrete_palettes(palettes))
    return palettes


def discrete_palette_config() -> dict[str, dict[str, list[str]]]:
    return {"color_palettes": discrete_palette_registry()}


def discrete_palette_names() -> list[str]:
    return list(discrete_palette_registry().keys()) + list(DYNAMIC_DISCRETE_PALETTES)


def discrete_palette_options() -> list[dict[str, str]]:
    return [{"label": name, "value": name} for name in discrete_palette_names()]


def discrete_palette_overview(dynamic_n_colors: int = 24) -> dict[str, list[dict[str, object]]]:
    """Return available discrete palettes grouped by source."""
    overview = {"custom": [], "plotly": [], "cc": [], "cmc": [], "glasbey": []}
    for name, palette in discrete_palette_registry().items():
        source = _discrete_palette_source(name)
        overview[source].append(_discrete_palette_entry(name, palette, dynamic=False))

    for name in DYNAMIC_DISCRETE_PALETTES:
        palette = resolve_discrete_palette(name, dynamic_n_colors)
        overview["glasbey"].append(_discrete_palette_entry(name, palette, dynamic=True))

    return overview


def plot_discrete_palette_overview(
    source: Literal["custom", "plotly", "cc", "cmc", "glasbey"] | None = None,
    *,
    dynamic_n_colors: int = 24,
    max_colors_per_palette: int = 32,
    row_height: float = 0.24,
    max_label_length: int = 48,
    figsize: tuple[float, float] | None = None,
):
    """Plot available discrete palettes and return a Matplotlib figure."""
    import matplotlib.pyplot as plt

    grouped_entries = _filtered_discrete_palette_overview(source, dynamic_n_colors)
    group_names = [group for group, entries in grouped_entries.items() if entries]
    if not group_names:
        raise ValueError("No discrete palettes match the requested filters.")

    max_rows = max(len(grouped_entries[group]) for group in group_names)
    if figsize is None:
        width = max(8, min(16, 3.8 * len(group_names)))
        height = max(2.5, max_rows * row_height + 0.9)
        figsize = (width, height)

    fig, axes = plt.subplots(
        ncols=len(group_names),
        figsize=figsize,
        squeeze=False,
        constrained_layout=True,
    )

    for ax, group_name in zip(axes[0], group_names, strict=True):
        entries = grouped_entries[group_name]
        for row_index, entry in enumerate(entries):
            palette = entry["palette"][:max_colors_per_palette]
            for color_index, color in enumerate(palette):
                ax.barh(
                    row_index,
                    1,
                    left=color_index,
                    color=_normalize_matplotlib_color(color),
                    edgecolor="none",
                    height=0.8,
                )

        ax.set_title(f"{group_name} ({len(entries)})")
        ax.set_xlim(0, max_colors_per_palette)
        ax.set_ylim(len(entries) - 0.5, -0.5)
        ax.set_xticks([])
        ax.set_yticks(range(len(entries)))
        ax.set_yticklabels(
            [_shorten_label(entry["label"], max_label_length) for entry in entries],
            fontsize=7,
        )
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

    return fig


def _filtered_discrete_palette_overview(
    source: str | None,
    dynamic_n_colors: int,
) -> dict[str, list[dict[str, object]]]:
    overview = discrete_palette_overview(dynamic_n_colors=dynamic_n_colors)
    if source is not None:
        if source not in overview:
            raise ValueError("source must be 'custom', 'plotly', 'cc', 'cmc', 'glasbey', or None.")
        overview = {source: overview[source]}
    return overview


def _discrete_palette_entry(name: str, palette: list[str], dynamic: bool) -> dict[str, object]:
    return {
        "source": _discrete_palette_source(name),
        "name": name.split("/", 1)[1] if "/" in name else name,
        "label": name,
        "value": name,
        "n_colors": len(palette),
        "dynamic": dynamic,
        "palette": palette,
    }


def _discrete_palette_source(name: str) -> str:
    if "/" not in name:
        return "custom"
    return name.split("/", 1)[0]


def get_discrete_palette(name: str | None, default: list[str] | None = None) -> list[str] | None:
    if name is None:
        return default
    return discrete_palette_registry().get(name, default)


def resolve_discrete_palette(
    name: str | None,
    n_colors: int | None = None,
    default: list[str] | None = None,
) -> list[str] | None:
    """Return a fixed or dynamic discrete palette.

    Dynamic Glasbey palettes are generated for ``n_colors`` and cached for the
    current Python session.
    """
    if name is None:
        return default
    if name in DYNAMIC_DISCRETE_PALETTES:
        return _dynamic_glasbey_palette(name, n_colors or 0)
    return discrete_palette_registry().get(name, default)


@lru_cache(maxsize=128)
def _dynamic_glasbey_palette(name: str, n_colors: int) -> list[str]:
    n_colors = max(1, int(n_colors))
    colorblind_safe = bool(DYNAMIC_DISCRETE_PALETTES[name]["colorblind_safe"])
    try:
        import glasbey
    except Exception:
        return _fallback_glasbey_palette(n_colors)

    return list(
        glasbey.create_palette(
            n_colors,
            as_hex=True,
            colorblind_safe=colorblind_safe,
            optimize_palette=False,
        )
    )


def _fallback_glasbey_palette(n_colors: int) -> list[str]:
    try:
        import colorcet as cc
    except Exception:
        return px.colors.qualitative.Plotly

    palette = list(getattr(cc, "glasbey", []))
    if not palette:
        return px.colors.qualitative.Plotly
    repeats = (n_colors // len(palette)) + 1
    return (palette * repeats)[:n_colors]


def continuous_colormap_options() -> list[dict[str, str]]:
    return [
        {"label": entry["label"], "value": entry["value"]}
        for entry in _continuous_colormap_entries(include_reversed=True)
    ]


def continuous_colormap_overview() -> dict[str, list[dict[str, str]]]:
    """Return available continuous colormaps grouped by semantic type.

    Each entry uses the same label shown in Dash dropdowns and keeps the stored
    value used internally by Plotly figure builders.
    """
    overview = {"linear": [], "diverging": []}

    for group_name, group in PLOTLY_PERCEPTUALLY_UNIFORM.items():
        colormap_type = _plotly_colormap_label_group(group_name)
        overview[colormap_type].extend(
            _continuous_colormap_entry("plotly", colormap_type, name, value)
            for name, value in group.items()
        )

    for name in _colorcet_continuous_original_names():
        colormap_type = _colorcet_colormap_type(name)
        overview[colormap_type].append(
            _continuous_colormap_entry(
                "cc",
                colormap_type,
                name,
                f"{COLORCET_PREFIX}{name}",
            )
        )

    for colormap_type, name, value in _cmcrameri_continuous_entries():
        overview[colormap_type].append(
            _continuous_colormap_entry("cmc", colormap_type, name, value)
        )

    return overview


def _continuous_colormap_entries(include_reversed: bool = False) -> list[dict[str, str]]:
    entries = [
        entry
        for group_entries in continuous_colormap_overview().values()
        for entry in group_entries
    ]
    if not include_reversed:
        return entries

    with_reversed = []
    for entry in entries:
        with_reversed.append(entry)
        with_reversed.append(_reversed_continuous_colormap_entry(entry))
    return with_reversed


def _reversed_continuous_colormap_entry(entry: dict[str, str]) -> dict[str, str]:
    return {
        **entry,
        "name": f"{entry['name']}_r",
        "label": f"{entry['label']}_r",
        "value": f"{entry['value']}_r",
    }


def plot_continuous_colormap_overview(
    colormap_type: Literal["linear", "diverging"] | None = None,
    source: Literal["plotly", "cc", "cmc"] | None = None,
    *,
    samples: int = 256,
    row_height: float = 0.24,
    max_label_length: int = 52,
    figsize: tuple[float, float] | None = None,
):
    """Plot available continuous colormaps and return a Matplotlib figure.

    Parameters
    ----------
    colormap_type
        Plot only ``"linear"`` or ``"diverging"`` maps. By default both groups
        are shown as separate panels.
    source
        Optional source filter: ``"plotly"``, ``"cc"``, or ``"cmc"``.
    samples
        Number of color samples per gradient.
    row_height
        Figure height used per colormap row.
    max_label_length
        Long labels are shortened to this length in the plot margin.
    figsize
        Optional Matplotlib figure size.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import ListedColormap

    grouped_entries = _filtered_continuous_colormap_overview(colormap_type, source)
    group_names = [group for group, entries in grouped_entries.items() if entries]
    if not group_names:
        raise ValueError("No continuous colormaps match the requested filters.")

    max_rows = max(len(grouped_entries[group]) for group in group_names)
    if figsize is None:
        width = 14 if len(group_names) == 2 else 7
        height = max(2.5, max_rows * row_height + 0.8)
        figsize = (width, height)

    fig, axes = plt.subplots(
        ncols=len(group_names),
        figsize=figsize,
        squeeze=False,
        constrained_layout=True,
    )
    gradient = np.linspace(0, 1, max(2, samples)).reshape(1, -1)

    for ax, group_name in zip(axes[0], group_names, strict=True):
        entries = grouped_entries[group_name]
        for row_index, entry in enumerate(entries):
            cmap = ListedColormap(_sample_continuous_colormap(entry["value"], samples))
            ax.imshow(
                gradient,
                aspect="auto",
                cmap=cmap,
                extent=(0, 1, row_index - 0.4, row_index + 0.4),
            )

        ax.set_title(f"{group_name.capitalize()} colormaps ({len(entries)})")
        ax.set_xlim(0, 1)
        ax.set_ylim(len(entries) - 0.5, -0.5)
        ax.set_xticks([])
        ax.set_yticks(range(len(entries)))
        ax.set_yticklabels(
            [_shorten_label(entry["label"], max_label_length) for entry in entries],
            fontsize=7,
        )
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

    return fig


def _filtered_continuous_colormap_overview(
    colormap_type: str | None,
    source: str | None,
) -> dict[str, list[dict[str, str]]]:
    overview = continuous_colormap_overview()
    if colormap_type is not None:
        if colormap_type not in overview:
            raise ValueError("colormap_type must be 'linear', 'diverging', or None.")
        overview = {colormap_type: overview[colormap_type]}

    if source is not None:
        if source not in {"plotly", "cc", "cmc"}:
            raise ValueError("source must be 'plotly', 'cc', 'cmc', or None.")
        overview = {
            group_name: [entry for entry in entries if entry["source"] == source]
            for group_name, entries in overview.items()
        }

    return overview


def _sample_continuous_colormap(color_map: str, samples: int) -> list[tuple[float, float, float]]:
    from matplotlib.colors import to_rgb
    from plotly.colors import sample_colorscale

    positions = [i / (max(2, samples) - 1) for i in range(max(2, samples))]
    colorscale = resolve_continuous_colorscale(color_map)
    sampled_colors = sample_colorscale(colorscale, positions)
    return [_plotly_color_to_rgb(color, to_rgb) for color in sampled_colors]


def _plotly_color_to_rgb(color: str, to_rgb) -> tuple[float, float, float]:
    if color.startswith("rgb"):
        channels = color[color.find("(") + 1 : color.rfind(")")].split(",")[:3]
        return tuple(float(channel.strip()) / 255 for channel in channels)
    return to_rgb(color)


def _normalize_matplotlib_color(color: str):
    if isinstance(color, str) and color.startswith("rgb"):
        channels = color[color.find("(") + 1 : color.rfind(")")].split(",")[:3]
        return tuple(float(channel.strip()) / 255 for channel in channels)
    return color


def _shorten_label(label: str, max_length: int) -> str:
    if len(label) <= max_length:
        return label
    return f"{label[: max_length - 1]}..."


def _continuous_colormap_entry(
    source: str, colormap_type: str, name: str, value: str
) -> dict[str, str]:
    label_name = name
    type_prefix = f"{colormap_type}_"
    if source != "cc" or not name.startswith(type_prefix):
        label_name = f"{type_prefix}{name}"

    return {
        "source": source,
        "type": colormap_type,
        "name": name,
        "label": f"{source}/{label_name}",
        "value": value,
    }


def _plotly_colormap_label_group(group_name: str) -> str:
    return "linear" if group_name == "sequential" else group_name


def resolve_continuous_colorscale(color_map):
    if not isinstance(color_map, str):
        return color_map
    if color_map.startswith(COLORCET_PREFIX):
        palette = _colorcet_palette(color_map.removeprefix(COLORCET_PREFIX).strip())
        return _palette_to_colorscale(palette)

    if color_map.startswith(CMCRAMERI_PREFIX):
        palette = _cmcrameri_palette(color_map.removeprefix(CMCRAMERI_PREFIX).strip())
        return _palette_to_colorscale(palette)

    return color_map


def _palette_to_colorscale(palette: list[str]):
    if len(palette) < 2:
        return DEFAULT_CONTINUOUS_COLORMAP

    denom = len(palette) - 1
    return [[i / denom, color] for i, color in enumerate(palette)]


def _plotly_qualitative_palettes(existing: dict[str, list[str]]) -> dict[str, list[str]]:
    palettes = {}
    for name in dir(px.colors.qualitative):
        key = f"plotly/{name}"
        if name.startswith("_") or name.endswith("_r") or key in existing:
            continue
        palette = getattr(px.colors.qualitative, name, None)
        if isinstance(palette, (list, tuple)) and palette and all(isinstance(c, str) for c in palette):
            palettes[key] = list(palette)
    return palettes


def _colorcet_discrete_palettes(existing: dict[str, list[str]]) -> dict[str, list[str]]:
    palettes = {}
    try:
        import colorcet as cc
    except Exception:
        return palettes

    for name, palette in cc.palette.items():
        if "glasbey" not in name.lower():
            continue
        key = f"cc/{name}"
        if key not in existing:
            palettes[key] = list(palette)
    return palettes


def _cmcrameri_discrete_palettes(existing: dict[str, list[str]]) -> dict[str, list[str]]:
    palettes = {}
    try:
        import cmcrameri.cm as cmc
    except Exception:
        return palettes

    cmaps = getattr(cmc, "cmaps", None)
    available = set(cmaps.keys()) if isinstance(cmaps, dict) else set(dir(cmc))
    for name in CRAMERI_CATEGORICAL:
        key = f"cmc/{name}"
        if key in existing or name not in available:
            continue
        cmap = cmaps.get(name) if isinstance(cmaps, dict) else getattr(cmc, name, None)
        if cmap is not None:
            palettes[key] = [_rgba_to_rgb_string(cmap(i / (cmap.N - 1))) for i in range(cmap.N)]
    return palettes


def _colorcet_continuous_original_names() -> list[str]:
    try:
        import colorcet as cc
    except Exception:
        return []

    linear = cc.all_original_names(group="linear")
    diverging = cc.all_original_names(group="diverging")
    return sorted(name for name in set(linear + diverging) if "rainbow" not in name.lower())


def _colorcet_colormap_type(name: str) -> str:
    base_name = name.removesuffix("_r")
    if base_name.startswith("diverging_"):
        return "diverging"
    return "linear"


def _colorcet_palette(name: str) -> list[str]:
    if not name:
        return []

    try:
        import colorcet as cc
    except Exception:
        return []

    cmap = cc.cm.get(name)
    if cmap is None:
        return []

    return [_rgba_to_rgb_string(cmap(i / 255)) for i in range(256)]


def _cmcrameri_continuous_entries() -> list[tuple[str, str, str]]:
    try:
        import cmcrameri.cm as cmc
    except Exception:
        return []

    cmaps = getattr(cmc, "cmaps", None)
    available = set(cmaps.keys()) if isinstance(cmaps, dict) else set(dir(cmc))
    options = []
    for group_name, names in (
        ("linear", CRAMERI_SEQUENTIAL),
        ("diverging", CRAMERI_DIVERGING),
    ):
        options.extend(
            (group_name, name, f"{CMCRAMERI_PREFIX}{name}")
            for name in names
            if name in available
        )
    return options


def _cmcrameri_palette(name: str) -> list[str]:
    if not name:
        return []

    try:
        import cmcrameri.cm as cmc
    except Exception:
        return []

    cmaps = getattr(cmc, "cmaps", None)
    cmap = cmaps.get(name) if isinstance(cmaps, dict) else getattr(cmc, name, None)
    if cmap is None and name.endswith("_r"):
        base_name = name.removesuffix("_r")
        cmap = cmaps.get(base_name) if isinstance(cmaps, dict) else getattr(cmc, base_name, None)
        if cmap is not None:
            return [_rgba_to_rgb_string(cmap(i / 255)) for i in range(255, -1, -1)]
    if cmap is None:
        return []

    return [_rgba_to_rgb_string(cmap(i / 255)) for i in range(256)]


def _rgba_to_rgb_string(rgba) -> str:
    r, g, b = [round(float(channel) * 255) for channel in rgba[:3]]
    return f"rgb({r},{g},{b})"
