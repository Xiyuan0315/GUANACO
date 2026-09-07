"""Offline exports contain precomputed results, not implicit analyses or links."""

import json
import re

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

import guanaco as gc
from guanaco.static_export import StaticPanel, freeze_violin, write_dashboard


def payload(path):
    document = path.read_text()
    match = re.search(
        r'<script type="application/json" id="guanaco-static-data">(.*?)</script>',
        document,
        re.S,
    )
    return json.loads(match.group(1))


@pytest.fixture
def adata():
    return ad.AnnData(
        np.array(
            [[0.0, 1.0], [1.0, 2.0], [2.0, 3.0], [3.0, 0.0], [0.0, 1.0], [1.0, 0.0]]
        ),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(["A", "A", "A", "B", "B", "B"]),
                "condition": pd.Categorical(["x", "y", "x", "y", "x", "y"]),
            },
            index=[f"c{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )


def test_existing_plots_export_with_two_independent_presets(tmp_path, adata):
    figures = {}
    for groupby in ("cell_type", "condition"):
        figures[groupby] = {
            "dotplot": gc.pl.dotplot(adata, ["G1", "G2"], groupby, return_fig=True),
            "heatmap": gc.pl.heatmap(adata, ["G1", "G2"], groupby, return_fig=True),
            "violin": freeze_violin(
                gc.pl.violin(adata, ["G1", "G2"], groupby, return_fig=True)
            ),
        }
    panels = [
        StaticPanel(name, name, {g: figures[g][name] for g in figures})
        for name in ("dotplot", "heatmap", "violin")
    ]
    path = write_dashboard(tmp_path / "demo.html", panels)
    saved = payload(path)
    assert [p["id"] for p in saved] == ["dotplot", "heatmap", "violin"]
    for panel in saved:
        assert [p["label"] for p in panel["presets"]] == ["cell_type", "condition"]
        assert "links" not in panel
        assert all(
            t["type"] != "violin" for p in panel["presets"] for t in p["figure"]["data"]
        )
    document = path.read_text()
    assert "plotly.js v" in document
    assert not re.search(r"<script[^>]+src=", document)
    assert "connect-src 'none'" in document
    assert "window.dash_clientside" not in document


def test_freezing_does_not_modify_original_figure(adata):
    original = gc.pl.violin(adata, ["G1"], "cell_type", return_fig=True)
    before = original.to_json()
    frozen = freeze_violin(original)
    assert original.to_json() == before
    assert all(t.type != "violin" for t in frozen.data)
    assert frozen.layout.xaxis.type == "linear"
    assert list(frozen.layout.xaxis.ticktext) == ["A", "B"]
    summaries = [t.meta for t in frozen.data if t.meta]
    assert summaries[0]["sample_mean"] == 1.0
    assert summaries[0]["sample_count"] == 3


def violin(values, **kwargs):
    return go.Figure(
        go.Violin(
            y=values,
            x=["A"] * len(values),
            width=0.8,
            scalemode="width",
            points=False,
            spanmode="hard",
            **kwargs,
        )
    )


def test_gaussian_density_widths_are_precomputed_numerically():
    values = np.array([0.0, 1.0, 2.0])
    frozen = freeze_violin(violin(values, bandwidth=0.4), resolution=32)
    trace = frozen.data[0]
    grid = np.linspace(0, 2, 32)
    density = np.exp(-0.5 * ((grid[:, None] - values) / 0.4) ** 2).sum(axis=1)
    np.testing.assert_allclose(np.asarray(trace.x)[:32], -0.4 * density / density.max())
    np.testing.assert_allclose(np.asarray(trace.y)[:32], grid)


@pytest.mark.parametrize("values", [[0.0, 0.0, 0.0], [3.0]])
def test_constant_and_single_value_groups_are_finite(values):
    frozen = freeze_violin(violin(values))
    assert np.isfinite(frozen.data[0].x).all()
    assert np.isfinite(frozen.data[0].y).all()
    assert frozen.data[0].fill is None


def test_constant_group_retains_shared_bandwidth():
    fig = violin([0.0, 0.0], bandwidth=0.3)
    fig.update_traces(spanmode="soft")
    frozen = freeze_violin(fig)
    assert frozen.data[0].fill == "toself"
    assert min(frozen.data[0].y) == pytest.approx(-0.6)
    assert max(frozen.data[0].y) == pytest.approx(0.6)


def test_empty_or_nonfinite_groups_do_not_emit_nan():
    fig = violin([np.nan, np.inf])
    assert not freeze_violin(fig).data


@pytest.mark.parametrize(
    "change", [{"orientation": "h"}, {"box_visible": True}, {"side": "positive"}]
)
def test_unsupported_violin_modes_fail_explicitly(change):
    fig = violin([0.0, 1.0], bandwidth=0.4)
    fig.update_traces(**change)
    with pytest.raises(ValueError, match="Only vertical"):
        freeze_violin(fig)


def test_missing_bandwidth_fails_instead_of_silently_selecting_new_estimator():
    with pytest.raises(ValueError, match="explicit precomputed bandwidth"):
        freeze_violin(violin([0.0, 1.0]))


@pytest.mark.parametrize(
    "figure",
    [
        go.Figure(go.Violin(y=[1, 2])),
        go.Figure(go.Histogram(x=[1, 2])),
        go.Figure(go.Box(y=[1, 2])),
    ],
)
def test_statistical_traces_rejected_before_output(tmp_path, figure):
    target = tmp_path / "result.html"
    with pytest.raises(ValueError, match="Precompute"):
        write_dashboard(target, [StaticPanel("plot", "Plot", {"Default": figure})])
    assert not target.exists()


def test_script_delimiters_and_titles_are_escaped(tmp_path):
    attack = '</script><script>alert("not executed")</script>'
    fig = go.Figure(go.Scatter(x=[1], y=[2], text=[attack]))
    path = write_dashboard(
        tmp_path / "safe.html",
        [StaticPanel("plot", attack, {attack: fig})],
        title=attack,
    )
    assert attack not in path.read_text()
    assert payload(path)[0]["title"] == attack
    assert payload(path)[0]["presets"][0]["figure"]["data"][0]["text"] == [attack]


def test_existing_file_preserved_by_default(tmp_path):
    path = tmp_path / "existing.html"
    path.write_text("user content")
    with pytest.raises(FileExistsError):
        write_dashboard(path, [StaticPanel("plot", "Plot", {"Default": go.Figure()})])
    assert path.read_text() == "user content"


@pytest.mark.parametrize("ids", [["duplicate", "duplicate"], ["bad id"], ["</script>"]])
def test_invalid_ids_rejected(tmp_path, ids):
    with pytest.raises(ValueError, match="Panel IDs"):
        write_dashboard(
            tmp_path / "bad.html",
            [StaticPanel(id, id, {"Default": go.Figure()}) for id in ids],
        )


def test_empty_panels_and_presets_rejected(tmp_path):
    with pytest.raises(ValueError, match="At least one"):
        write_dashboard(tmp_path / "bad.html", [])
    with pytest.raises(ValueError, match="at least one precomputed preset"):
        write_dashboard(tmp_path / "bad.html", [StaticPanel("plot", "Plot", {})])


def test_missing_values_become_json_null(tmp_path):
    fig = go.Figure(go.Heatmap(z=[[1, None], [np.nan, 2]]))
    path = write_dashboard(
        tmp_path / "missing.html", [StaticPanel("heatmap", "Heatmap", {"Default": fig})]
    )
    assert payload(path)[0]["presets"][0]["figure"]["data"][0]["z"] == [
        [1, None],
        [None, 2],
    ]


def test_external_image_rejected_before_export(tmp_path):
    fig = go.Figure(go.Scatter(x=[1], y=[2]))
    fig.add_layout_image(source="https://example.org/private-image.png", x=0, y=1)
    with pytest.raises(ValueError, match="embedded data:image"):
        write_dashboard(
            tmp_path / "external.html", [StaticPanel("plot", "Plot", {"Default": fig})]
        )


def test_viewer_resources_are_not_auto_loaded_dash_assets():
    from pathlib import Path
    from guanaco import static_export

    root = Path(static_export.__file__).parent
    assert (root / "_static_viewer/viewer.js").is_file()
    assert not (root / "assets/static_viewer/viewer.js").exists()
    assert not (root / "assets/static_viewer/viewer.css").exists()
