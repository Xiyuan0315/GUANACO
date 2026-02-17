import json
from pathlib import Path

import dash_draggable
import plotly.express as px
from dash import html

from guanaco.utils.ui_helpers import (
    labeled_dropdown,
    labeled_radioitems,
    graph_flex_container,
)


def _palette_names():
    cvd_color_path = Path(__file__).resolve().parents[1] / "cvd_color.json"
    with open(cvd_color_path, "r") as f:
        palette_json = json.load(f)

    plotly_qualitative_palettes = {}
    for name in dir(px.colors.qualitative):
        if name.startswith("_"):
            continue
        palette = getattr(px.colors.qualitative, name, None)
        if isinstance(palette, (list, tuple)) and palette and all(isinstance(c, str) for c in palette):
            plotly_qualitative_palettes[name] = list(palette)
    for k, v in plotly_qualitative_palettes.items():
        if k not in palette_json.get("color_palettes", {}):
            palette_json["color_palettes"][k] = v
    try:
        import colorcet as cc

        for name, palette in cc.palette.items():
            if "glasbey" not in name.lower():
                continue
            key = f"colorcet/{name}"
            if key not in palette_json.get("color_palettes", {}):
                palette_json["color_palettes"][key] = list(palette)
    except Exception:
        # colorcet is optional; keep existing palettes if unavailable.
        pass
    return list(palette_json["color_palettes"].keys())


def _continuous_colormap_options():
    options = [{"label": scale, "value": scale} for scale in px.colors.named_colorscales()]
    try:
        import colorcet as cc

        cc_names = sorted(
            name for name in cc.palette.keys()
            if "glasbey" not in name.lower()
        )
        options.extend(
            {"label": f"colorcet/{name}", "value": f"cc:{name}"}
            for name in cc_names
        )
    except Exception:
        # colorcet is optional; keep Plotly maps if unavailable.
        pass
    return options


def generate_heatmap_layout(adata, prefix):
    label_list = adata.obs.columns.to_list()
    colorscales = _continuous_colormap_options()
    palette_names = _palette_names()

    heatmap_transformation_selection = labeled_radioitems(
        "Transformation:",
        f"{prefix}-heatmap-transformation",
        [
            {"label": "None", "value": "None"},
            {"label": "Log", "value": "log"},
            {"label": "Z-score (across cell)", "value": "z_score"},
        ],
        value="None",
        inline=True,
        radio_style={"marginBottom": "10px"},
    )

    heatmap_secondary_dropdown = labeled_dropdown(
        "Secondary Annotation:",
        f"{prefix}-heatmap-label-dropdown",
        [{"label": "None", "value": "None"}] + [{"label": label, "value": label} for label in label_list],
        value="None",
        clearable=True,
        label_style={"fontWeight": "bold", "marginTop": "10px"},
        dropdown_style={"width": "200px"},
        wrapper_style={"marginBottom": "10px"},
    )

    heatmap_color_map_dropdown = labeled_dropdown(
        "Continuous ColorMap:",
        f"{prefix}-heatmap-colorscale-dropdown",
        colorscales,
        value="viridis",
        clearable=False,
        dropdown_style={"width": "200px", "marginBottom": "10px"},
    )

    secondary_annotation_colormap_dropdown = html.Div(
        labeled_dropdown(
            "Secondary Annotation ColorMap:",
            f"{prefix}-heatmap-secondary-colormap-dropdown",
            [{"label": name, "value": name} for name in palette_names],
            value="Pastel",
            placeholder="Select colormap for secondary annotation",
            clearable=False,
            dropdown_style={"width": "200px", "marginBottom": "10px"},
        ),
        id=f"{prefix}-heatmap-secondary-colormap-wrapper",
        style={"display": "none"},
    )

    draggable_container = dash_draggable.GridLayout(
        id=f"{prefix}-draggable-heatmap",
        className="grid-layout-no-border",
        children=[graph_flex_container(f"{prefix}-heatmap")],
    )

    return html.Div(
        [
            heatmap_transformation_selection,
            heatmap_secondary_dropdown,
            heatmap_color_map_dropdown,
            secondary_annotation_colormap_dropdown,
            draggable_container,
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )
