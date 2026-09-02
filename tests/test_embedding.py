import base64
from io import BytesIO

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from anndata import AnnData
from PIL import Image

from guanaco.pages.matrix.plots import embedding as embedding_module
from guanaco.pages.matrix.plots.embedding import plot_coexpression_embedding, plot_embedding


def _embedding_adata():
    x = np.array(
        [
            [0.0, 0.0, 0.1],
            [1.0, 1.0, 0.8],
            [2.0, 0.5, 1.5],
        ],
        dtype=np.float32,
    )
    obs = pd.DataFrame(
        {
            "cell_type": ["T", "B", "T"],
            "score": [2.0, 1.0, 3.0],
        },
        index=["c1", "c2", "c3"],
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])
    adata = AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = np.array(
        [
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 1.0],
        ],
        dtype=np.float32,
    )
    adata.obsm["spatial"] = np.array(
        [
            [1.0, 1.0],
            [2.0, 3.0],
            [3.0, 2.0],
        ],
        dtype=np.float32,
    )
    adata.uns["spatial"] = {
        "lib1": {
            "images": {"hires": np.zeros((12, 20, 3), dtype=np.uint8)},
            "scalefactors": {"tissue_hires_scalef": 2.0},
        }
    }
    return adata


def test_spatial_background_downscales_bitmap_without_changing_coordinate_space(monkeypatch):
    monkeypatch.setattr(embedding_module, "SPATIAL_IMAGE_MAX_DIM", 10)

    fig = go.Figure()
    embedding_module._apply_spatial_background(
        fig,
        np.zeros((12, 20, 3), dtype=np.uint8),
        img_alpha=0.4,
    )

    image = fig.layout.images[0]
    encoded_image = image.source.split(",", 1)[1]
    with Image.open(BytesIO(base64.b64decode(encoded_image))) as downscaled:
        assert downscaled.size == (10, 6)
    assert image.sizex == 20
    assert image.sizey == 12
    assert image.opacity == 0.4
    assert tuple(fig.layout.xaxis.range) == (0, 20)
    assert tuple(fig.layout.yaxis.range) == (-12, 0)
    assert fig.layout.xaxis.autorange is False
    assert fig.layout.yaxis.autorange is False
    assert fig.layout.xaxis.constrain == "domain"
    assert fig.layout.yaxis.constrain == "domain"
    assert fig.layout.yaxis.scaleanchor == "x"


def test_continuous_embedding_order_uses_original_row_positions():
    fig = plot_embedding(_embedding_adata(), "X_umap", "score", mode="continuous", order="max")

    assert isinstance(fig.data[0], go.Scattergl)
    assert fig.data[0].customdata.tolist() == [1, 0, 2]


def test_categorical_embedding_customdata_exposes_source_row_position():
    fig = plot_embedding(
        _embedding_adata(),
        "X_umap",
        "cell_type",
        mode="categorical",
    )

    positions = set()
    for trace in fig.data:
        assert isinstance(trace, go.Scattergl)
        positions.update(np.asarray(trace.customdata).reshape(-1).tolist())
    assert positions == {0, 1, 2}


def test_spatial_categorical_embedding_keeps_original_image_extent():
    fig = plot_embedding(_embedding_adata(), "spatial", "cell_type", mode="categorical")

    assert all(isinstance(trace, go.Scattergl) for trace in fig.data)
    assert len(fig.layout.images) == 1
    assert fig.layout.images[0].sizex == 20
    assert fig.layout.images[0].sizey == 12
    assert tuple(fig.layout.xaxis.range) == (0, 20)
    assert tuple(fig.layout.yaxis.range) == (-12, 0)
    assert fig.layout.xaxis.constrain == "domain"
    assert fig.layout.yaxis.constrain == "domain"


def test_categorical_embedding_accepts_mixed_labels_and_drops_missing_values():
    adata = _embedding_adata()
    adata.obs["mixed_label"] = pd.Series(
        ["T", 2, np.nan],
        index=adata.obs_names,
        dtype=object,
    )

    fig = plot_embedding(
        adata,
        "X_umap",
        "mixed_label",
        mode="categorical",
    )

    assert {trace.name for trace in fig.data} == {"2", "T"}
    assert sum(len(trace.x) for trace in fig.data) == 2


def test_coexpression_embedding_customdata_uses_row_positions():
    fig = plot_coexpression_embedding(_embedding_adata(), "X_umap", "GeneA", "GeneB")

    assert fig.layout.title.text == "<b>Co-expression: GeneA & GeneB</b>"
    assert fig.data
    for trace in fig.data:
        assert isinstance(trace, go.Scattergl)
        assert set(trace.customdata.tolist()).issubset({0, 1, 2})
