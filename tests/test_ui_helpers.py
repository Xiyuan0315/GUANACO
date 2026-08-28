from pathlib import Path

from dash import html
from dash_draggable import ResponsiveGridLayout

from guanaco.utils.ui_helpers import responsive_graph_grid


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
