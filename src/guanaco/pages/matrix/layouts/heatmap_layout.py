from dash import dcc, html

from guanaco.utils.render_guard import rendered_key_store
from guanaco.utils.colors import discrete_palette_options
from guanaco.utils.plot_config import heatmap_config
from guanaco.utils.ui_helpers import (
    labeled_dropdown,
    labeled_radioitems,
    graph_flex_container,
    responsive_graph_grid,
)


def generate_heatmap_layout(adata, prefix):
    label_list = adata.obs.columns.to_list()
    palette_options = discrete_palette_options()

    heatmap_transformation_selection = labeled_radioitems(
        "Standardization:",
        f"{prefix}-heatmap-standardization",
        [
            {"label": "None", "value": "None"},
            {"label": "Min–max", "value": "minmax"},
            {"label": "Z-score", "value": "zscore"},
        ],
        value="minmax",
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

    secondary_annotation_colormap_dropdown = html.Div(
        labeled_dropdown(
            "Secondary Annotation ColorMap:",
            f"{prefix}-heatmap-secondary-colormap-dropdown",
            palette_options,
            value="Pastel",
            placeholder="Select colormap for secondary annotation",
            clearable=False,
            dropdown_style={"width": "200px", "marginBottom": "10px"},
        ),
        id=f"{prefix}-heatmap-secondary-colormap-wrapper",
        style={"display": "none"},
    )

    # Width-responsive grid (scales with the screen instead of capping at ~1200px),
    # matching the dotplot: opens at 3/4 width x ~510px, drag-resize bounded by the
    # helper's min/max. The container id lets the per-item min/max reach
    # react-grid-layout. Fresh grid id (was "-draggable-heatmap").
    grid_item_id = f"{prefix}-heatmap-grid-item"
    draggable_container = responsive_graph_grid(
        f"{prefix}-heatmap-grid",
        grid_item_id,
        graph_flex_container(f"{prefix}-heatmap", config=heatmap_config, container_id=grid_item_id),
    )

    return html.Div(
        [
            rendered_key_store(prefix, "heatmap"),
            dcc.Store(id=f"{prefix}-heatmap-reset-link"),
            heatmap_transformation_selection,
            heatmap_secondary_dropdown,
            secondary_annotation_colormap_dropdown,
            draggable_container,
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )
