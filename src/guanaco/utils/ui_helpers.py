from dash import dcc, html
import dash_bootstrap_components as dbc

from guanaco.utils.plot_config import common_config


LOADING_OVERLAY_STYLE = {"visibility": "visible", "filter": "blur(2px)"}


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


def graph_flex_container(graph_id, *, graph_style=None, config=None):
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
