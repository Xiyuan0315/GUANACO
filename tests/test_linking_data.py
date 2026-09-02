from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from guanaco.linking.data import DataSource, DataStore


def _adata(
    *,
    obs_names: list[object] | None = None,
    var_names: list[object] | None = None,
):
    obs_names = obs_names or ["cell-1", "cell-2"]
    var_names = var_names or ["CD3D", "MS4A1"]
    return SimpleNamespace(
        X=np.zeros((len(obs_names), len(var_names))),
        obs=pd.DataFrame(index=pd.Index(obs_names)),
        var=pd.DataFrame(index=pd.Index(var_names)),
        obs_names=pd.Index(obs_names),
        var_names=pd.Index(var_names),
    )


def test_anndata_and_table_expose_their_stable_identity_domains() -> None:
    cells = DataStore.from_data(_adata()).source()
    rows = DataStore.from_data(
        pd.DataFrame({"score": [0.8, 0.6]}, index=[1, 2])
    ).source()

    assert cells.kind == "anndata"
    assert cells.cell_ids.tolist() == ["cell-1", "cell-2"]
    assert cells.feature_ids.tolist() == ["CD3D", "MS4A1"]
    assert cells.row_ids is None
    assert rows.kind == "table"
    assert rows.row_ids.tolist() == ["1", "2"]
    assert rows.cell_ids is None
    assert rows.feature_ids is None


def test_mapping_and_mudata_create_named_sources_without_container_sources() -> None:
    mapping = DataStore.from_data(
        {
            "rna": _adata(),
            "scores": pd.DataFrame({"score": [1.0]}, index=["row-1"]),
        }
    )
    mudata = DataStore.from_data(
        SimpleNamespace(mod={"rna": _adata(), "protein": _adata()})
    )

    assert mapping.default_source is None
    assert mapping.source("rna").kind == "anndata"
    assert mapping.source("scores").kind == "table"
    assert tuple(mudata.sources) == ("rna", "protein")
    assert all(source.kind == "anndata" for source in mudata.sources.values())
    with pytest.raises(ValueError, match="multiple sources"):
        mapping.source()
    with pytest.raises(ValueError, match="Unknown linked-view data source"):
        mapping.source("missing")


@pytest.mark.parametrize(
    ("kind", "axis", "ids", "message"),
    [
        ("table", None, ["row-1", None], "must not contain null IDs"),
        ("table", None, ["row-1", "row-1"], "must contain unique IDs"),
        ("table", None, [1, "1"], "remain unique"),
        ("table", None, ["row-1", "  "], "empty or whitespace-only IDs"),
        ("anndata", "obs", ["cell-1", None], "must not contain null IDs"),
        ("anndata", "var", ["CD4", "CD4"], "must contain unique IDs"),
    ],
)
def test_browser_identity_indices_must_be_stable(kind, axis, ids, message) -> None:
    if kind == "table":
        value = pd.DataFrame({"value": range(len(ids))}, index=ids)
    else:
        value = _adata(**{f"{axis}_names": ids})

    with pytest.raises(ValueError, match=message):
        DataStore.from_data(value)


@pytest.mark.parametrize(
    ("factory", "error", "message"),
    [
        (lambda: DataStore.from_data({}), ValueError, "at least one source"),
        (
            lambda: DataStore.from_data({"bad": object()}),
            TypeError,
            "must be AnnData-like or a pandas DataFrame",
        ),
        (
            lambda: DataStore.from_data({"": pd.DataFrame(index=["row-1"])}),
            ValueError,
            "non-empty strings",
        ),
        (
            lambda: DataSource("payload", "raw", {"x": [1, 2]}),
            ValueError,
            "Unsupported data source kind",
        ),
    ],
)
def test_invalid_inputs_fail_at_the_data_boundary(factory, error, message) -> None:
    with pytest.raises(error, match=message):
        factory()
