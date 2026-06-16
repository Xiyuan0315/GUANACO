from guanaco.pages.matrix.callbacks.scatter_callbacks import _is_reset_relayout


def test_is_reset_relayout_true_for_xaxis_autorange():
    # Double-click "reset axes" emits an autorange relayout event.
    assert _is_reset_relayout({"xaxis.autorange": True}) is True


def test_is_reset_relayout_true_for_yaxis_autorange():
    assert _is_reset_relayout({"yaxis.autorange": True}) is True


def test_is_reset_relayout_false_for_zoom_range_event():
    # Zoom/pan emits explicit axis ranges, not autorange -> not a reset.
    zoom = {"xaxis.range[0]": 0.0, "xaxis.range[1]": 5.0}
    assert _is_reset_relayout(zoom) is False


def test_is_reset_relayout_false_for_empty_or_none():
    assert _is_reset_relayout(None) is False
    assert _is_reset_relayout({}) is False
