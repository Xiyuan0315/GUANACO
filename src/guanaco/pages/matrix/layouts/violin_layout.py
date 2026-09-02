import dash_bootstrap_components as dbc
from dash import dcc, html

from guanaco.utils.plot_config import common_config
from guanaco.utils.ui_helpers import (
    LOADING_OVERLAY_STYLE,
    graph_flex_container,
    labeled_dropdown,
    responsive_graph_grid,
    switch_checklist,
)


def generate_violin_layout(default_gene_markers, discrete_label_list, prefix):
    """Stacked violin tab: one violin per selected gene, stacked as subplot rows.

    Driven by the shared left-panel controls (genes / annotation / labels); the only
    local option is "Show Box Plot".
    """
    violin_show_box1 = switch_checklist(f"{prefix}-show-box1", "Show Box Plot")

    return html.Div(
        [
            html.Div(
                [violin_show_box1],
                style={"marginBottom": "15px"},
            ),
            dcc.Store(id=f"{prefix}-violin-plot-cache-store"),
            dcc.Store(id=f"{prefix}-violin1-rendered-key"),
            # No dcc.Loading wrapper here: the stacked-violin figure is served
            # from cache on tab switches (the gene/label data is cached server
            # side), so the spinner only ever *flashed* on re-entry without
            # signalling real work. Rendering the graph directly removes that
            # distracting flash.
            dcc.Graph(
                id=f"{prefix}-violin-plot1",
                config=common_config,
                responsive=True,
                style={"width": "100%", "minHeight": "400px"},
            ),
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )


def generate_ridge_layout(default_gene_markers, prefix):
    """Single-gene ridgeline plot driven by the shared annotation and label controls."""
    genes = list(default_gene_markers or [])
    gene_selection = labeled_dropdown(
        "Gene:",
        f"{prefix}-ridge-gene-selection",
        [{"label": gene, "value": gene} for gene in genes],
        value=genes[0] if genes else None,
        clearable=False,
        placeholder="Select one of the chosen genes",
        wrapper_style={"flex": "1", "minWidth": "240px"},
    )
    show_box = switch_checklist(f"{prefix}-ridge-show-box", "Show Box Plot")

    return html.Div(
        [
            html.Div(
                [gene_selection, show_box],
                style={
                    "display": "flex",
                    "alignItems": "end",
                    "gap": "20px",
                    "marginBottom": "15px",
                },
            ),
            dcc.Store(id=f"{prefix}-ridge-rendered-key"),
            dcc.Loading(
                id=f"{prefix}-ridge-loading",
                type="circle",
                overlay_style=LOADING_OVERLAY_STYLE,
                children=[
                    responsive_graph_grid(
                        f"{prefix}-ridge-grid",
                        f"{prefix}-ridge-grid-item",
                        graph_flex_container(
                            f"{prefix}-ridge-plot",
                            container_id=f"{prefix}-ridge-grid-item",
                            config=common_config,
                        ),
                    )
                ],
            ),
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )


def generate_split_violin_layout(
    adata, default_gene_markers, discrete_label_list, prefix
):
    """Comparative violin tab: one gene compared across obs1 (and optional obs2).

    Only the essentials -- gene + primary grouping -- are shown up front. The optional
    secondary grouping, analysis mode, statistical test and box-plot toggle live in a
    collapsible "Grouping & statistics" panel so the default view stays uncluttered.
    """
    violin_show_box2 = switch_checklist(f"{prefix}-show-box2", "Show Box Plot")
    default_meta = discrete_label_list[0] if discrete_label_list else None

    meta1_selection = labeled_dropdown(
        "Group by (Obs1):",
        f"{prefix}-meta1-selection",
        [{"label": meta, "value": meta} for meta in discrete_label_list],
        value=default_meta,
        clearable=False,
        placeholder="Select primary grouping",
        wrapper_style={"flex": "1"},
    )

    meta2_selection = labeled_dropdown(
        "Compare by (Obs2):",
        f"{prefix}-meta2-selection",
        [{"label": "None", "value": "none"}]
        + [{"label": meta, "value": meta} for meta in discrete_label_list],
        value="none",
        clearable=False,
        placeholder="Optional secondary grouping",
        wrapper_style={"flex": "1"},
    )

    mode_selection = labeled_dropdown(
        "Analysis mode:",
        f"{prefix}-mode-selection",
        [
            {"label": "Mode 1: One variable", "value": "mode1"},
            {"label": "Mode 2: Split/grouped by Obs2", "value": "mode2"},
            {"label": "Mode 3: Linear model (obs1 + obs2)", "value": "mode3"},
            {"label": "Mode 4: Mixed model (obs2 as random effect)", "value": "mode4"},
        ],
        value="mode1",
        clearable=False,
        wrapper_style={"flex": "1"},
    )

    test_method_selection = labeled_dropdown(
        "Statistical test:",
        f"{prefix}-test-method-selection",
        [
            {"label": "Auto", "value": "auto"},
            {"label": "None", "value": "none"},
            {"label": "Mann-Whitney U", "value": "mwu-test"},
            {"label": "T-test", "value": "ttest"},
            {"label": "Kruskal-Wallis", "value": "kw-test"},
            {"label": "ANOVA", "value": "anova"},
            {"label": "Linear Model", "value": "linear-model"},
            {
                "label": "Linear Model with Interaction",
                "value": "linear-model-interaction",
            },
            {"label": "Mixed Model", "value": "mixed-model"},
        ],
        value="auto",
        clearable=False,
        wrapper_style={"flex": "1"},
    )

    violin2_gene_selection = labeled_dropdown(
        "Select Gene",
        f"{prefix}-violin2-gene-selection",
        [{"label": gene, "value": gene} for gene in default_gene_markers],
        value=default_gene_markers[0] if default_gene_markers else None,
        wrapper_style={"flex": "1"},
    )

    data_layer_selection = dcc.Store(
        id=f"{prefix}-violin2-data-layer",
        data="X",
    )

    advanced_toggle = dbc.Button(
        "▸ More options",
        id=f"{prefix}-split-violin-options-toggle",
        color="link",
        size="sm",
        style={
            "padding": "2px 0",
            "textDecoration": "none",
            "fontWeight": "normal",
            "marginBottom": "10px",
        },
    )

    advanced_panel = dbc.Collapse(
        html.Div(
            [
                html.Div(
                    [
                        mode_selection,
                        meta2_selection,
                        test_method_selection,
                    ],
                    style={"display": "flex", "marginBottom": "10px", "gap": "10px"},
                ),
                html.Div(
                    id=f"{prefix}-mode-explanation",
                    style={
                        "fontSize": "12px",
                        "color": "gray",
                        "marginBottom": "10px",
                        "fontStyle": "italic",
                    },
                ),
            ],
            style={
                "padding": "12px",
                "border": "1px solid #e9ecef",
                "borderRadius": "6px",
                "background": "#fbfbfc",
                "marginBottom": "12px",
            },
        ),
        id=f"{prefix}-split-violin-options-collapse",
        is_open=False,
    )

    return html.Div(
        [
            data_layer_selection,
            html.Div(
                [violin2_gene_selection, meta1_selection],
                style={"display": "flex", "gap": "10px", "marginBottom": "10px"},
            ),
            html.Div(
                [violin_show_box2],
                style={"marginBottom": "15px"},
            ),
            advanced_toggle,
            advanced_panel,
            dcc.Loading(
                id=f"{prefix}-loading-violin2",
                type="circle",
                overlay_style=LOADING_OVERLAY_STYLE,
                children=[
                    responsive_graph_grid(
                        f"{prefix}-violin2-grid",
                        f"{prefix}-violin2-grid-item",
                        html.Div(
                            id=f"{prefix}-violin2-grid-item",
                            children=dcc.Graph(
                                id=f"{prefix}-violin-plot2",
                                config=common_config,
                                responsive=True,
                                style={
                                    "min-height": "0",
                                    "flex-grow": "1",
                                    "border": "none",
                                    "box-shadow": "none",
                                },
                            ),
                            style={
                                "height": "100%",
                                "width": "100%",
                                "display": "flex",
                                "flex-direction": "column",
                                "flex-grow": "0",
                                "border": "none",
                                "box-shadow": "none",
                            },
                        ),
                    )
                ],
            ),
        ],
        style={"padding": "20px", "marginBottom": "15px"},
    )
