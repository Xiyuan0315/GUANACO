import json
import subprocess
import sys

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

from guanaco.config_wizard import (
    DEFAULT_PLOTS,
    EXPLORATORY_VISUALIZATION_PLOTS,
    MARKER_VISUALIZATION_PLOTS,
    OPTIONAL_PLOTS,
    inspect_data_file,
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
        ("Network", "network"),
        ("Ligand–receptor", "ligand-receptor"),
        ("Spatial relationships", "spatial-relationships"),
        ("Peak Browser", "peak-browser"),
        ("Multi-omics coverage", "multiomics-composition"),
        ("Omics comparison", "cross-modal-concordance"),
    )


def test_wizard_serialization_order_contains_each_plot_key_once():
    values = [value for _label, value in OPTIONAL_PLOTS]

    assert OPTIONAL_PLOTS == (
        MARKER_VISUALIZATION_PLOTS + EXPLORATORY_VISUALIZATION_PLOTS
    )
    assert len(values) == len(set(values))
    assert DEFAULT_PLOTS <= set(values)


def test_wizard_inspects_metadata_and_offers_only_compatible_plots(tmp_path):
    n_obs = 60
    adata = AnnData(
        X=np.ones((n_obs, 3), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(["B", "T"] * (n_obs // 2)),
                "pseudotime": np.linspace(0, 1, n_obs, dtype=np.float32),
            },
            index=[f"cell-{index}" for index in range(n_obs)],
        ),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    adata.uns["liana_res"] = pd.DataFrame(
        {
            "source": ["B"],
            "target": ["T"],
            "ligand": ["G1"],
            "receptor": ["G2"],
            "score": [0.8],
        }
    )
    path = tmp_path / "capabilities.h5ad"
    adata.write_h5ad(path)

    inspection = inspect_data_file(str(path))

    assert inspection.modalities is None
    assert {
        "dotplot",
        "heatmap",
        "violin",
        "expression-trend",
        "split-violin",
        "stacked-bar",
        "network",
        "ligand-receptor",
    } <= inspection.available_plots
    assert {
        "paga",
        "volcano",
        "spatial-relationships",
        "peak-browser",
        "multiomics-composition",
        "cross-modal-concordance",
    }.isdisjoint(inspection.available_plots)


def test_wizard_detects_peak_features_without_atac_filename(tmp_path):
    adata = AnnData(
        X=np.ones((4, 5), dtype=np.float32),
        var=pd.DataFrame(
            index=[f"chr1:{start}-{start + 100}" for start in range(0, 500, 100)]
        ),
    )
    path = tmp_path / "generic.h5ad"
    adata.write_h5ad(path)

    inspection = inspect_data_file(str(path))

    assert "peak-browser" in inspection.available_plots
    assert "ligand-receptor" not in inspection.available_plots


def test_wizard_h5mu_modalities_inherit_shared_observation_metadata(tmp_path):
    patient_ids = pd.Index([f"patient-{index}" for index in range(60)])
    modalities = {
        name: AnnData(
            X=np.ones((60, 3), dtype=np.float32),
            obs=pd.DataFrame(index=patient_ids.copy()),
            var=pd.DataFrame(index=[f"{name}-feature-{index}" for index in range(3)]),
        )
        for name in ("rna", "drug")
    }
    mdata = MuData(modalities)
    mdata.obs["response"] = pd.Categorical(["Yes", "No"] * 30)
    mdata.obs["age"] = np.linspace(30, 90, 60, dtype=np.float32)
    path = tmp_path / "shared-clinical-metadata.h5mu"
    mdata.write_h5mu(path)

    inspection = inspect_data_file(str(path))

    assert inspection.modalities == ("rna", "drug")
    assert {
        "dotplot",
        "heatmap",
        "violin",
        "expression-trend",
        "split-violin",
        "stacked-bar",
        "network",
        "multiomics-composition",
        "cross-modal-concordance",
    } <= inspection.available_plots


def test_wizard_inspection_does_not_import_runtime_data_stack(tmp_path):
    path = tmp_path / "lightweight.h5ad"
    AnnData(
        X=np.ones((4, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {"group": pd.Categorical(["A", "A", "B", "B"])},
            index=[f"cell-{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    ).write_h5ad(path)
    script = """
import sys
from guanaco.config_wizard import inspect_data_file

inspection = inspect_data_file(sys.argv[1])
assert "dotplot" in inspection.available_plots
assert "anndata" not in sys.modules
assert "muon" not in sys.modules
assert "scanpy" not in sys.modules
"""

    subprocess.run(
        [sys.executable, "-c", script, str(path)],
        check=True,
        capture_output=True,
        text=True,
    )


def test_wizard_inspection_cache_invalidates_when_file_changes(tmp_path):
    path = tmp_path / "changing.h5ad"
    AnnData(
        X=np.ones((4, 5), dtype=np.float32),
        var=pd.DataFrame(index=[f"G{index}" for index in range(5)]),
    ).write_h5ad(path)
    initial = inspect_data_file(str(path))
    assert "peak-browser" not in initial.available_plots

    AnnData(
        X=np.ones((4, 5), dtype=np.float32),
        var=pd.DataFrame(
            index=[f"chr1:{start}-{start + 100}" for start in range(0, 500, 100)]
        ),
    ).write_h5ad(path)
    updated = inspect_data_file(str(path))

    assert "peak-browser" in updated.available_plots


def test_runtime_uses_wizard_config_even_if_loader_was_imported_early(tmp_path):
    config_path = tmp_path / "chosen.json"
    config_path.write_text(
        json.dumps({"title": "Chosen dataset", "settings": {}})
    )
    script = """
import os
import sys

os.environ.pop("GUANACO_CONFIG", None)
import guanaco.data.loader as loader
assert str(loader.JSON_PATH) == "guanaco.json"

os.environ["GUANACO_CONFIG"] = sys.argv[1]
import guanaco.data.registry as registry
assert str(registry.JSON_PATH) == sys.argv[1]
assert registry.cfg["title"] == "Chosen dataset"
"""

    subprocess.run(
        [sys.executable, "-c", script, str(config_path)],
        check=True,
        capture_output=True,
        text=True,
    )
