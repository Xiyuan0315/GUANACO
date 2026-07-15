from guanaco.config_wizard import (
    DEFAULT_PLOTS,
    EXPLORATORY_VISUALIZATION_PLOTS,
    MARKER_VISUALIZATION_PLOTS,
    OPTIONAL_PLOTS,
)


def test_wizard_plot_groups_match_dashboard_labels_and_keys():
    assert MARKER_VISUALIZATION_PLOTS == (
        ("Dotplot", "dotplot"),
        ("Heatmap", "heatmap"),
        ("Violin Plot", "violin"),
        ("Expression Trend", "expression-trend"),
    )
    assert EXPLORATORY_VISUALIZATION_PLOTS == (
        ("Comparative Violin", "split-violin"),
        ("Composition", "stacked-bar"),
        ("PAGA", "paga"),
        ("Volcano Plot", "volcano"),
        ("GRN", "grn"),
        ("Peak Browser", "peak-browser"),
    )


def test_wizard_serialization_order_contains_each_plot_key_once():
    values = [value for _label, value in OPTIONAL_PLOTS]

    assert OPTIONAL_PLOTS == (
        MARKER_VISUALIZATION_PLOTS + EXPLORATORY_VISUALIZATION_PLOTS
    )
    assert len(values) == len(set(values))
    assert DEFAULT_PLOTS <= set(values)
