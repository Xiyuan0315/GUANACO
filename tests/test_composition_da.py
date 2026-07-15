import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from guanaco.pages.matrix.analysis import calculate_alr_welch
from guanaco.pages.matrix.plots.stacked_bar import (
    plot_composition_differential_abundance,
)


def _composition_adata():
    rows = []
    designs = {
        "A": [(5, 25, 70), (6, 24, 70), (4, 26, 70), (5, 25, 70)],
        "B": [(30, 25, 45), (32, 24, 44), (28, 26, 46), (31, 25, 44)],
        "C": [(60, 25, 15), (62, 24, 14), (58, 26, 16), (61, 25, 14)],
    }
    for group, samples in designs.items():
        for sample_number, population_counts in enumerate(samples, start=1):
            sample = f"{group}{sample_number}"
            for population, count in zip(
                ["Target", "Stable", "Reference"],
                population_counts,
                strict=True,
            ):
                rows.extend(
                    {
                        "condition": group,
                        "sample": sample,
                        "cell_type": population,
                    }
                    for _ in range(count)
                )

    obs = pd.DataFrame(rows)
    obs.index = [f"cell-{index}" for index in range(len(obs))]
    return AnnData(
        X=np.ones((len(obs), 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["G1"]),
    )


def test_alr_welch_compares_every_group_pair_and_population_globally():
    results = calculate_alr_welch(
        _composition_adata(),
        group_key="condition",
        population_key="cell_type",
        sample_key="sample",
        reference_population="Reference",
        group_order=["A", "B", "C"],
    )

    assert set(results["population"]) == {"Target", "Stable"}
    assert results["comparison"].drop_duplicates().tolist() == [
        "B − A",
        "C − A",
        "C − B",
    ]
    assert len(results) == 6
    assert (results[["n_a", "n_b"]] == 4).all().all()
    assert results["q_value"].between(0, 1).all()
    assert (results["q_value"] >= results["p_value"] - 1e-12).all()
    assert (
        results.loc[results["population"] == "Target", "effect"] > 0
    ).all()


def test_alr_welch_rejects_samples_that_cross_x_axis_groups():
    adata = _composition_adata()
    adata.obs.loc[adata.obs.index[0], "sample"] = "B1"

    with pytest.raises(ValueError, match="exactly one x-axis group"):
        calculate_alr_welch(
            adata,
            group_key="condition",
            population_key="cell_type",
            sample_key="sample",
            reference_population="Reference",
        )


def test_alr_welch_requires_replicated_samples_in_each_comparison():
    adata = _composition_adata()
    keep = (adata.obs["condition"] == "A") | (
        (adata.obs["condition"] == "B") & (adata.obs["sample"] == "B1")
    )

    with pytest.raises(ValueError, match="at least two samples"):
        calculate_alr_welch(
            adata[keep].copy(),
            group_key="condition",
            population_key="cell_type",
            sample_key="sample",
            reference_population="Reference",
        )


def test_differential_abundance_plot_places_pairwise_results_in_a_heatmap():
    results = calculate_alr_welch(
        _composition_adata(),
        group_key="condition",
        population_key="cell_type",
        sample_key="sample",
        reference_population="Reference",
        group_order=["A", "B", "C"],
    )
    fig = plot_composition_differential_abundance(results)

    trace = fig.data[0]
    assert trace.type == "heatmap"
    assert list(trace.x) == ["B − A", "C − A", "C − B"]
    assert set(trace.y) == {"Target", "Stable"}
    assert any(marker == "*" for row in trace.text for marker in row)
    assert fig.layout.title.text.startswith("ALR + pairwise Welch t-tests")
    assert "Reference population:" not in fig.layout.title.text
    assert "Differential abundance" not in fig.layout.title.text
