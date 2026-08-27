"""Sample-level differential abundance for categorical cell populations."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests

from guanaco.utils.obs_utils import obs_values


def _welch_result(group_a, group_b):
    """Return B-A effect, confidence interval, and two-sided Welch p-value."""
    group_a = np.asarray(group_a, dtype=float)
    group_b = np.asarray(group_b, dtype=float)
    n_a = len(group_a)
    n_b = len(group_b)
    mean_a = float(group_a.mean())
    mean_b = float(group_b.mean())
    effect = mean_b - mean_a
    var_a = float(group_a.var(ddof=1))
    var_b = float(group_b.var(ddof=1))
    term_a = var_a / n_a
    term_b = var_b / n_b
    variance = term_a + term_b

    if variance == 0:
        p_value = 1.0 if np.isclose(effect, 0) else 0.0
        return effect, effect, effect, p_value

    denominator = (term_a**2 / (n_a - 1)) + (term_b**2 / (n_b - 1))
    degrees_freedom = variance**2 / denominator if denominator else np.inf
    standard_error = np.sqrt(variance)
    statistic = effect / standard_error
    p_value = float(2 * student_t.sf(abs(statistic), degrees_freedom))
    margin = float(student_t.ppf(0.975, degrees_freedom) * standard_error)
    return effect, effect - margin, effect + margin, p_value


def calculate_alr_welch(
    adata,
    *,
    group_key,
    population_key,
    sample_key,
    reference_population,
    group_order=None,
    pseudocount=0.5,
    group_values=None,
):
    """Run all pairwise sample-level ALR Welch tests with global BH correction."""
    if not sample_key:
        raise ValueError(
            "Select a biological sample unit (for example donor or replicate)."
        )
    if group_key == population_key:
        raise ValueError(
            "Group bars by and Stack bars by must be different for differential abundance."
        )
    if sample_key in {group_key, population_key}:
        raise ValueError(
            "Sample unit must be different from the x-axis and population annotations."
        )
    if pseudocount <= 0:
        raise ValueError("Pseudocount must be greater than zero.")

    required = (group_key, population_key, sample_key)
    missing = [
        key
        for key in required
        if key not in adata.obs.columns
        and not (group_values and key in group_values)
    ]
    if missing:
        raise ValueError(f"Missing observation metadata: {', '.join(missing)}")

    frame = pd.DataFrame(
        {
            "group": obs_values(
                adata, group_key, group_values
            ).astype(object),
            "population": obs_values(
                adata, population_key, group_values
            ).astype(object),
            "sample": obs_values(
                adata, sample_key, group_values
            ).astype(object),
        }
    ).dropna()
    if frame.empty:
        raise ValueError("No complete cells are available for differential abundance.")
    frame = frame.astype(str)

    groups_per_sample = frame.groupby("sample", observed=True)["group"].nunique()
    invalid_samples = groups_per_sample[groups_per_sample > 1].index.tolist()
    if invalid_samples:
        preview = ", ".join(invalid_samples[:3])
        raise ValueError(
            "Each biological sample must belong to exactly one x-axis group. "
            f"Check: {preview}."
        )

    counts = pd.crosstab(frame["sample"], frame["population"])
    sample_groups = frame.groupby("sample", observed=True, sort=False)["group"].first()
    sample_groups = sample_groups.reindex(counts.index)
    if reference_population not in counts.columns:
        raise ValueError("Select a valid ALR reference population.")

    observed_groups = sample_groups.drop_duplicates().tolist()
    requested_order = [str(value) for value in group_order or []]
    ordered_groups = [
        value for value in requested_order if value in observed_groups
    ] + [value for value in observed_groups if value not in requested_order]
    if len(ordered_groups) < 2:
        raise ValueError("At least two x-axis groups are required.")

    proportions = counts.div(counts.sum(axis=1), axis=0)
    reference_counts = counts[reference_population].astype(float)
    records = []
    for population in counts.columns:
        if population == reference_population:
            continue
        alr = np.log(
            (counts[population].astype(float) + pseudocount)
            / (reference_counts + pseudocount)
        )
        for group_a, group_b in combinations(ordered_groups, 2):
            mask_a = sample_groups == group_a
            mask_b = sample_groups == group_b
            n_a = int(mask_a.sum())
            n_b = int(mask_b.sum())
            if n_a < 2 or n_b < 2:
                continue
            effect, ci_low, ci_high, p_value = _welch_result(
                alr[mask_a], alr[mask_b]
            )
            records.append(
                {
                    "population": str(population),
                    "group_a": group_a,
                    "group_b": group_b,
                    "comparison": f"{group_b} − {group_a}",
                    "effect": effect,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "p_value": p_value,
                    "n_a": n_a,
                    "n_b": n_b,
                    "mean_proportion_a": float(proportions.loc[mask_a, population].mean()),
                    "mean_proportion_b": float(proportions.loc[mask_b, population].mean()),
                }
            )

    if not records:
        raise ValueError(
            "No valid pairwise tests: every compared group needs at least two samples."
        )

    results = pd.DataFrame.from_records(records)
    results["q_value"] = multipletests(
        results["p_value"].to_numpy(), method="fdr_bh"
    )[1]
    return results
