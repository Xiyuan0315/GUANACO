from pathlib import Path
import re

from dash import Dash, dcc, html
from dash_draggable import ResponsiveGridLayout

from guanaco.utils.ui_helpers import responsive_graph_grid


def test_dash4_control_theme_is_served_as_a_global_asset():
    assets = Path(__file__).parents[1] / "src/guanaco/assets"
    app = Dash(__name__, assets_folder=str(assets))
    app.layout = html.Div([
        dcc.Dropdown(id="features", options=["Cortex_1", "Cortex_2"], multi=True),
        dcc.Slider(id="size", min=0, max=10, value=3),
        dcc.RangeSlider(id="range", min=0, max=10, value=[2, 8]),
    ])
    client = app.server.test_client()
    assert "scientific_style.css" in client.get("/").get_data(as_text=True)
    response = client.get("/assets/scientific_style.css")
    assert response.status_code == 200
    css = response.get_data(as_text=True)
    # Root scope also reaches portaled menus. html:root outranks the async
    # component's :root defaults without an order-dependent !important patch.
    rule = re.search(r"html:root\s*\{([^}]+)\}", css)
    assert rule is not None
    assert "--Dash-Fill-Interactive-Strong: #343a40;" in rule[1]
    assert "--Dash-Fill-Primary-Hover: #f1f3f5;" in rule[1]
    assert "--Dash-Fill-Disabled: #dee2e6;" in rule[1]


def test_control_theme_does_not_keep_obsolete_dash3_selectors():
    css = (Path(__file__).parents[1] / "src/guanaco/assets/scientific_style.css").read_text()
    assert ".Select-" not in css
    assert ".rc-slider-" not in css


def test_responsive_graph_grid_uses_dash_draggable_with_size_constraints():
    child = html.Div(id="plot-item")

    wrapper = responsive_graph_grid(
        "plot-grid",
        "plot-item",
        child,
        w=8,
        h=11,
        min_w=5,
        min_h=6,
        max_w=12,
        max_h=18,
    )

    grid = wrapper
    assert isinstance(grid, ResponsiveGridLayout)
    assert grid.className == "grid-layout-no-border"
    assert grid.id == "plot-grid"
    assert grid.height == 30
    assert grid.isDraggable is True
    assert grid.isResizable is True
    assert grid.resizeHandles == ["se"]
    assert grid.save is True
    assert grid.clearSavedLayout is False
    assert grid.gridCols == {"lg": 12, "md": 12, "sm": 12, "xs": 12, "xxs": 12}
    assert grid.layouts["lg"] == [
        {
            "i": "plot-item",
            "x": 0,
            "y": 0,
            "w": 8,
            "h": 11,
            "minW": 5,
            "minH": 6,
            "maxW": 12,
            "maxH": 18,
        }
    ]
    assert grid.layouts["sm"][0]["w"] == 12
    grid_item = grid.children[0]
    assert grid_item.id == "plot-item"
    assert child.id == "plot-item"


def test_vendored_dash_draggable_uses_dash4_public_child_layout_api():
    source = (
        Path(__file__).parents[1]
        / "vendor/dash_draggable/src/lib/components/ResponsiveGridLayout.react.js"
    ).read_text()

    assert "window.dash_component_api.getLayout(child.props.componentPath)" in source
    assert "_dashprivate_layout" not in source
    assert "this.layouts = all_layouts" in source
    assert "saveToLs(`${id}-layouts`, all_layouts)" in source
