from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_draggable

from guanaco.utils.plot_config import common_config


LOADING_OVERLAY_STYLE = {"visibility": "visible", "filter": "blur(2px)"}

# Standard react-grid-layout responsive breakpoints (px) and a uniform 12-col
# grid at each. ResponsiveGridLayout (unlike the fixed-1200px GridLayout) is
# width-aware via WidthProvider, so an item's width is a fraction of the *actual*
# container and scales to fill large screens instead of capping at ~1200px.
_GRID_BREAKPOINTS = {"lg": 1200, "md": 996, "sm": 768, "xs": 480, "xxs": 0}
_GRID_COLS = {bp: 12 for bp in _GRID_BREAKPOINTS}


def responsive_graph_grid(grid_id, item_id, child, *, w=9, h=13,
                          min_w=6, min_h=7, max_w=12, max_h=20):
    """A width-responsive dash_draggable grid holding one resizable graph item.

    w/h are the default open size and min_*/max_* bound drag-resize, all in grid
    units (12 cols; rows are 30px). Because the grid is responsive, w is a
    fraction of the container width, so the plot scales with the screen. Small
    breakpoints open full-width (max_w) for readability. ``child`` must carry
    ``id=item_id`` so dash_draggable keys it (key = child id) and the per-item
    minW/minH/maxW/maxH pass through to react-grid-layout.
    """
    def item(width):
        return {"i": item_id, "x": 0, "y": 0, "w": width, "h": h,
                "minW": min_w, "minH": min_h, "maxW": max_w, "maxH": max_h}

    layouts = {
        "lg": [item(w)],
        "md": [item(w)],
        "sm": [item(max_w)],
        "xs": [item(max_w)],
        "xxs": [item(max_w)],
    }
    return dash_draggable.ResponsiveGridLayout(
        id=grid_id,
        className="grid-layout-no-border",
        layouts=layouts,
        breakpoints=dict(_GRID_BREAKPOINTS),
        gridCols=dict(_GRID_COLS),
        height=30,
        isResizable=True,
        isDraggable=True,
        style={"backgroundColor": "transparent", "padding": "0px", "border": "none", "boxShadow": "none"},
        children=[child],
    )


def labeled_dropdown(
    label,
    dropdown_id,
    options,
    *,
    value=None,
    placeholder=None,
    clearable=True,
    multi=False,
    label_style=None,
    label_id=None,
    dropdown_style=None,
    wrapper_style=None,
):
    label_kwargs = {'style': label_style or {'fontWeight': 'bold', 'marginBottom': '5px'}}
    if label_id is not None:
        label_kwargs['id'] = label_id
    return html.Div(
        [
            html.Label(label, **label_kwargs),
            dcc.Dropdown(
                id=dropdown_id,
                options=options,
                value=value,
                placeholder=placeholder,
                clearable=clearable,
                multi=multi,
                style=dropdown_style or {},
            ),
        ],
        style=wrapper_style or {},
    )


def labeled_radioitems(
    label,
    radio_id,
    options,
    *,
    value=None,
    inline=True,
    label_style=None,
    radio_style=None,
    wrapper_style=None,
):
    return html.Div(
        [
            html.Label(label, style=label_style or {'fontWeight': 'bold', 'marginBottom': '5px'}),
            dbc.RadioItems(
                id=radio_id,
                options=options,
                value=value,
                inline=inline,
                style=radio_style or {},
            ),
        ],
        style=wrapper_style or {},
    )


def switch_checklist(check_id, label):
    return html.Div(
        [
            dbc.Checklist(
                id=check_id,
                options=[{'label': label, 'value': 'show'}],
                value=[],
                switch=True,
            )
        ]
    )


def graph_flex_container(graph_id, *, graph_style=None, config=None, container_id=None):
    # container_id lets a dash_draggable GridLayout key this wrapper (key = child
    # id) so a `layout` entry can target it by `i` and pass minW/minH/maxW/maxH
    # through to react-grid-layout.
    div_kwargs = {"id": container_id} if container_id else {}
    return html.Div(
        children=[
            dcc.Graph(
                id=graph_id,
                config=config or common_config,
                responsive=True,
                style=graph_style or {"min-height": "0", "flex-grow": "1"},
            )
        ],
        style={
            "height": '100%',
            "width": '100%',
            "display": "flex",
            "flex-direction": "column",
            "flex-grow": "0",
        },
        **div_kwargs,
    )


def component_flex_container(component_id, *, style=None):
    return html.Div(
        id=component_id,
        style=style or {
            "height": "100%",
            "width": "100%",
            "display": "flex",
            "flex-direction": "column",
            "flex-grow": "0",
        },
    )
