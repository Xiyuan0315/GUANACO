import dash_draggable
from dash import dcc, html

from guanaco.config import common_config
from guanaco.utils.ui_helpers import labeled_dropdown, labeled_radioitems, switch_checklist


def generate_violin_layout(default_gene_markers, discrete_label_list, prefix):
    violin1_transformation_selection = labeled_radioitems(
        "Transformation:",
        f"{prefix}-violin-log-or-zscore",
        [
            {"label": "None", "value": "False"},
            {"label": "Log", "value": "log"},
        ],
        value="False",
        inline=True,
        radio_style={"marginLeft": "10px"},
        wrapper_style={"marginBottom": "15px"},
    )

    violin2_transformation_selection = labeled_radioitems(
        "Transformation:",
        f"{prefix}-violin2-log-or-zscore",
        [
            {"label": "None", "value": "False"},
            {"label": "Log", "value": "log"},
        ],
        value="False",
        inline=False,
        wrapper_style={"marginBottom": "10px"},
    )

    violin_show_box1 = switch_checklist(f"{prefix}-show-box1", "Show Box Plot")
    violin_show_scatter1 = switch_checklist(f"{prefix}-show-scatter1", "Show Scatter Points")

    violin1_more_options = html.Div(
        [
            html.Label("More Options:", style={"fontWeight": "bold"}),
            html.Div([violin_show_box1, violin_show_scatter1], style={"display": "flex", "gap": "20px"}),
        ],
        style={"marginBottom": "15px"},
    )

    violin_show_box2 = switch_checklist(f"{prefix}-show-box2", "Show Box Plot")
    violin_show_scatter2 = switch_checklist(f"{prefix}-show-scatter2", "Show Scatter Points")

    meta1_selection = labeled_dropdown(
        "Obs1 (Primary):",
        f"{prefix}-meta1-selection",
        [{"label": meta, "value": meta} for meta in discrete_label_list],
        value=discrete_label_list[0],
        clearable=False,
        placeholder="Select primary metadata",
        wrapper_style={"flex": "1"},
    )

    meta2_selection = labeled_dropdown(
        "Obs2 (Secondary):",
        f"{prefix}-meta2-selection",
        [{"label": "None", "value": "none"}] + [{"label": meta, "value": meta} for meta in discrete_label_list],
        value="none",
        clearable=False,
        placeholder="Select secondary metadata (optional)",
        wrapper_style={"flex": "1"},
    )

    mode_selection = labeled_dropdown(
        "Analysis Mode:",
        f"{prefix}-mode-selection",
        [
            {"label": "Mode 1: One metadata only", "value": "mode1"},
            {"label": "Mode 2: Facet by meta1, compare meta2", "value": "mode2"},
            {"label": "Mode 3: Linear model (meta1 + meta2)", "value": "mode3"},
            {"label": "Mode 4: Mixed model (meta1 + (1|meta2))", "value": "mode4"},
        ],
        value="mode1",
        clearable=False,
        wrapper_style={"flex": "1"},
    )

    test_method_selection = labeled_dropdown(
        "Statistical Test:",
        f"{prefix}-test-method-selection",
        [
            {"label": "None", "value": "none"},
            {"label": "Mann-Whitney U", "value": "mwu-test"},
            {"label": "T-test", "value": "ttest"},
            {"label": "Kruskal-Wallis", "value": "kw-test"},
            {"label": "ANOVA", "value": "anova"},
            {"label": "Linear Model", "value": "linear-model"},
            {"label": "Linear Model with Interaction", "value": "linear-model-interaction"},
            {"label": "Mixed Model", "value": "mixed-model"},
        ],
        value="none",
        clearable=False,
        wrapper_style={"flex": "1"},
    )

    violin2_gene_selection = labeled_dropdown(
        "Select Gene",
        f"{prefix}-violin2-gene-selection",
        [{"label": gene, "value": gene} for gene in default_gene_markers],
        value=default_gene_markers[0] if default_gene_markers else None,
    )

    return html.Div(
        [
            html.Div(
                [
                    violin1_transformation_selection,
                    violin1_more_options,
                    dcc.Store(id=f"{prefix}-violin-plot-cache-store"),
                    dcc.Loading(
                        id="loading-violin1",
                        type="circle",
                        children=[
                            dcc.Graph(
                                id=f"{prefix}-violin-plot1",
                                config=common_config,
                                style={"width": "100%", "minHeight": "400px"},
                            )
                        ],
                    ),
                ],
                style={"marginBottom": "30px", "padding": "20px"},
            ),
            html.Div(
                [
                    html.H4(
                        "Split Violin/Grouped Violin",
                        style={"textAlign": "center", "margin": "10px 0", "fontWeight": "bold"},
                    ),
                    violin2_gene_selection,
                    violin2_transformation_selection,
                    html.Div([meta1_selection, meta2_selection], style={"display": "flex", "marginBottom": "10px", "gap": "10px"}),
                    html.Div([mode_selection, test_method_selection], style={"display": "flex", "marginBottom": "10px", "gap": "10px"}),
                    html.Div(
                        id=f"{prefix}-mode-explanation",
                        style={"fontSize": "12px", "color": "gray", "marginBottom": "10px", "fontStyle": "italic"},
                    ),
                    html.Div(
                        [
                            html.Label("More Options:", style={"fontWeight": "bold"}),
                            violin_show_box2,
                            violin_show_scatter2,
                        ],
                        style={"marginBottom": "10px"},
                    ),
                    dcc.Loading(
                        id="loading-violin2",
                        type="circle",
                        children=[
                            html.Div(
                                [
                                    dash_draggable.GridLayout(
                                        id=f"{prefix}-draggable-violin2",
                                        className="grid-layout-no-border",
                                        children=[
                                            html.Div(
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
                                            )
                                        ],
                                        isResizable=True,
                                        isDraggable=True,
                                        height=30,
                                        gridCols=12,
                                        style={"background": "transparent", "padding": "0px"},
                                    )
                                ]
                            )
                        ],
                        style={"padding": "10px", "display": "flex", "flex-direction": "column", "height": "100%"},
                    ),
                ],
                style={"padding": "10px"},
            ),
        ],
        style={"width": "100%"},
    )
