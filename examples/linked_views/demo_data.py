"""Deterministic data preparation for the index-linked notebook demos.

The plotting examples deliberately receive only ordinary AnnData, MuData, and
pandas objects. Every external table has a unique, non-null string index; links
may additionally use a shared table column when one overview row expands to
many detail records.
"""

from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
CELL_TYPES = ("CD4 T", "CD8 T", "B", "NK", "Mono", "DC")
MARKERS = (
    "CD3D",
    "CD4",
    "IL7R",
    "CCL5",
    "MS4A1",
    "NKG7",
    "LST1",
    "FCGR3A",
    "CLEC10A",
    "CCR5",
    "CXCL8",
    "CXCR2",
    "IL7",
    "IL2RG",
    "TNF",
    "TNFRSF1A",
)


def _dense(values) -> np.ndarray:
    if hasattr(values, "toarray"):
        values = values.toarray()
    return np.asarray(values)


def _principal_coordinates(values) -> np.ndarray:
    """Return two deterministic principal coordinates without another dependency."""

    matrix = _dense(values).astype(float)
    matrix -= matrix.mean(axis=0, keepdims=True)
    left, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    coordinates = left[:, :2] * singular[:2]
    if coordinates.shape[1] == 1:
        coordinates = np.column_stack([coordinates[:, 0], np.zeros(len(matrix))])
    return coordinates.astype(np.float32)


def _first_existing(
    explicit: str | Path | None,
    environment_name: str,
    candidates: list[Path],
) -> Path | None:
    choices = [
        Path(explicit).expanduser() if explicit is not None else None,
        (
            Path(os.environ[environment_name]).expanduser()
            if os.environ.get(environment_name)
            else None
        ),
        *candidates,
    ]
    return next(
        (
            candidate.resolve()
            for candidate in choices
            if candidate and candidate.is_file()
        ),
        None,
    )


def make_single_cell(n_cells: int = 540, *, seed: int = 315) -> ad.AnnData:
    """Create a small zero-rich single-cell dataset with marker programs."""

    rng = np.random.default_rng(seed)
    labels = np.resize(np.asarray(CELL_TYPES, dtype=object), n_cells)
    rng.shuffle(labels)
    condition = np.where(rng.random(n_cells) < 0.48, "Control", "Stimulated")

    centers = {
        "CD4 T": (-3.2, 1.8),
        "CD8 T": (-2.3, -1.2),
        "B": (0.2, 3.0),
        "NK": (0.4, -2.9),
        "Mono": (3.2, -0.4),
        "DC": (3.4, 2.5),
    }
    coordinates = np.vstack(
        [np.asarray(centers[str(label)]) + rng.normal(0, 0.62, 2) for label in labels]
    ).astype(np.float32)
    coordinates[:, 0] += np.where(condition == "Stimulated", 0.34, -0.12)

    genes = pd.Index(MARKERS, name="feature")
    matrix = rng.gamma(1.25, 0.42, size=(n_cells, len(genes)))
    programs = {
        "CD4 T": {"CD3D": 3.8, "CD4": 3.4, "IL7R": 4.4, "IL2RG": 1.8},
        "CD8 T": {"CD3D": 3.2, "CCL5": 4.1, "NKG7": 2.5},
        "B": {"MS4A1": 5.2, "IL2RG": 1.0},
        "NK": {"NKG7": 5.3, "CCL5": 4.0, "CCR5": 1.7},
        "Mono": {"LST1": 4.8, "FCGR3A": 3.2, "CXCL8": 3.0, "TNF": 2.0},
        "DC": {"CLEC10A": 4.7, "IL7": 2.8, "TNFRSF1A": 1.4},
    }
    positions = {gene: index for index, gene in enumerate(genes)}
    for cell_type, program in programs.items():
        mask = labels == cell_type
        for gene, strength in program.items():
            matrix[mask, positions[gene]] += rng.gamma(
                2.2, strength / 2.2, size=int(mask.sum())
            )
    stimulated = condition == "Stimulated"
    matrix[stimulated, positions["TNF"]] += rng.gamma(
        2.0, 1.15, size=int(stimulated.sum())
    )
    matrix = np.log1p(rng.poisson(matrix)).astype(np.float32)

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(labels, categories=CELL_TYPES),
            "condition": pd.Categorical(
                condition, categories=["Control", "Stimulated"]
            ),
            "sample": pd.Categorical(
                [f"donor-{index % 6 + 1}" for index in range(n_cells)]
            ),
        },
        index=pd.Index(
            [f"cell-{index:04d}" for index in range(n_cells)], name="cell_id"
        ),
    )
    result = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=genes))
    result.obsm["X_umap"] = coordinates
    return result


def make_external_cell_table(adata: ad.AnnData, *, seed: int = 177) -> pd.DataFrame:
    """Make independent per-cell attributes indexed exactly like ``obs_names``."""

    rng = np.random.default_rng(seed)
    matrix = _dense(adata.X).astype(float)
    feature_position = {
        str(value): index for index, value in enumerate(adata.var_names)
    }
    library_complexity = np.log1p((matrix > 0).sum(axis=1))
    activation = (
        matrix[:, feature_position["TNF"]] + matrix[:, feature_position["CCL5"]]
    )
    activation += rng.normal(0, 0.18, adata.n_obs)
    result = pd.DataFrame(
        {
            "library_complexity": library_complexity,
            "activation_score": activation,
            "cell_type": adata.obs["cell_type"].astype(str).to_numpy(),
        },
        index=pd.Index(adata.obs_names.astype(str), name="cell_id"),
    )
    return result


def make_volcano_data(*, seed: int = 1208) -> tuple[ad.AnnData, pd.DataFrame]:
    """Return cell data plus a feature-indexed differential-result table."""

    adata = make_single_cell(560, seed=seed)
    matrix = _dense(adata.X).astype(float)
    stimulated = adata.obs["condition"].astype(str).eq("Stimulated").to_numpy()
    positions = {str(gene): index for index, gene in enumerate(adata.var_names)}
    rng = np.random.default_rng(seed + 1)
    for gene, strength in {"TNF": 3.8, "CCL5": 2.9, "CXCL8": 2.5}.items():
        matrix[stimulated, positions[gene]] += rng.gamma(
            2.2, strength / 2.2, int(stimulated.sum())
        )
    for gene, strength in {"IL7R": 2.8, "MS4A1": 2.3, "CLEC10A": 1.8}.items():
        matrix[~stimulated, positions[gene]] += rng.gamma(
            2.2, strength / 2.2, int((~stimulated).sum())
        )
    adata.X = matrix.astype(np.float32)

    treated = matrix[stimulated]
    control = matrix[~stimulated]
    fold_change = np.log2((treated.mean(axis=0) + 0.1) / (control.mean(axis=0) + 0.1))
    test = stats.ttest_ind(treated, control, axis=0, equal_var=False, nan_policy="omit")
    pvalue = np.nan_to_num(test.pvalue, nan=1.0, posinf=1.0, neginf=1.0)
    order = np.argsort(pvalue)
    ranked = pvalue[order] * len(pvalue) / np.arange(1, len(pvalue) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 1e-300, 1.0)
    significant = (adjusted <= 0.05) & (np.abs(fold_change) >= 0.5)
    labels = np.where(significant & (fold_change > 0), "Up", "Not significant")
    labels = np.where(significant & (fold_change < 0), "Down", labels)
    table = pd.DataFrame(
        {
            "feature": adata.var_names.astype(str),
            "log2_fold_change": fold_change,
            "adjusted_p_value": adjusted,
            "negative_log10_padj": -np.log10(adjusted),
            "result": labels,
        },
        index=pd.Index(adata.var_names.astype(str), name="feature_id"),
    )
    return adata, table


def make_pathway_demo(*, seed: int = 1018) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return pathway summary and long pathway-to-gene detail tables."""

    adata = make_single_cell(seed=seed)
    pathways = {
        "T-cell activation": ["CD3D", "CD4", "IL7R", "CCL5"],
        "Cytotoxicity": ["CCL5", "NKG7", "TNF"],
        "B-cell receptor": ["MS4A1", "CD3D", "IL2RG"],
        "Myeloid activation": ["LST1", "FCGR3A", "CXCL8"],
    }
    expression = pd.DataFrame(
        np.asarray(adata.X, dtype=float),
        index=adata.obs_names,
        columns=adata.var_names,
    )
    cell_types = adata.obs["cell_type"].astype(str)
    group_means = expression.assign(cell_type=cell_types).groupby(
        "cell_type", sort=False, observed=True
    )
    group_means = group_means.mean(numeric_only=True)
    means = expression.mean(axis=0)
    records: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    for pathway, genes in pathways.items():
        pathway_score = float(means.loc[genes].mean())
        summary.append(
            {
                "pathway": pathway,
                "pathway_score": pathway_score,
                "n_genes": len(genes),
            }
        )
        for cell_type, values in group_means.iterrows():
            for gene in genes:
                records.append(
                    {
                        "record_id": f"{pathway}::{cell_type}::{gene}",
                        "pathway": pathway,
                        "cell_type": cell_type,
                        "gene": gene,
                        "pathway_score": pathway_score,
                        "gene_mean": float(values[gene]),
                    }
                )
    pathway_table = pd.DataFrame.from_records(summary).set_index("pathway", drop=False)
    pathway_table.index.name = "pathway_id"
    gene_table = pd.DataFrame.from_records(records).set_index("record_id")
    return pathway_table, gene_table


def make_mudata_demo(*, seed: int = 930):
    """Small RNA/protein fallback with exactly matching cell indices."""

    import mudata as md

    base = make_single_cell(420, seed=seed)
    rng = np.random.default_rng(seed + 1)
    rna_features = ["CD4", "IL7R", "CCL5", "MS4A1", "NKG7", "LST1"]
    rna = base[:, rna_features].copy()
    rna.obsm["X_umap"] = np.asarray(base.obsm["X_umap"], dtype=np.float32)
    protein_features = pd.Index(["CD4 protein", "CD127", "CD20", "CD56", "CD14"])
    source_columns = [0, 1, 3, 4, 5]
    protein_matrix = np.maximum(
        _dense(rna.X)[:, source_columns] * np.asarray([0.95, 0.9, 0.85, 1.1, 1.0])
        + rng.normal(0, 0.55, size=(rna.n_obs, len(protein_features))),
        0,
    )
    protein = ad.AnnData(
        X=protein_matrix.astype(np.float32),
        obs=base.obs.copy(),
        var=pd.DataFrame(index=protein_features),
    )
    protein.obsm["X_pca"] = _principal_coordinates(protein.X)
    return md.MuData({"rna": rna, "protein": protein}), "CD4", "CD4 protein"


def load_pbmc_mudata_demo(
    path: str | Path | None = None,
    *,
    max_cells: int = 5_000,
    seed: int = 930,
):
    """Use the real PBMC RNA/ADT data when available, else a small fallback."""

    repository = Path(__file__).resolve().parents[2]
    selected = _first_existing(
        path,
        "GUANACO_PBMC_H5MU",
        [
            repository / "data" / "mdata_pbmc.h5mu",
            Path("/Users/xiyuanzhang/Documents/GUANACO_v2/data/mdata_pbmc.h5mu"),
        ],
    )
    if selected is None:
        mdata, rna_feature, protein_feature = make_mudata_demo(seed=seed)
        return mdata, rna_feature, protein_feature, None

    import mudata as md

    backed = md.read_h5mu(selected, backed="r")
    try:
        rna_source = backed.mod["rna"]
        protein_source = backed.mod["adt"]
        if not rna_source.obs_names.equals(protein_source.obs_names):
            raise ValueError("PBMC RNA and ADT cell indices do not match.")
        rng = np.random.default_rng(seed)
        positions = np.arange(rna_source.n_obs)
        if len(positions) > max_cells:
            positions = np.sort(rng.choice(positions, max_cells, replace=False))

        rna_features = [
            feature
            for feature in ("CD4", "IL7R", "LST1", "NKG7", "MS4A1")
            if feature in rna_source.var_names
        ]
        protein_features = [
            feature
            for feature in ("adt_CD4-1", "adt_CD4-2", "adt_CD14", "adt_CD19")
            if feature in protein_source.var_names
        ]
        if not rna_features or not protein_features:
            raise ValueError(
                "PBMC modalities do not contain the expected CD4 features."
            )
        obs = rna_source.obs.iloc[positions].copy()
        rna_matrix = _dense(rna_source[positions, rna_features].X).astype(np.float32)
        protein_matrix = _dense(protein_source[positions, protein_features].X).astype(
            np.float32
        )
        rna = ad.AnnData(
            X=rna_matrix,
            obs=obs.copy(),
            var=pd.DataFrame(index=pd.Index(rna_features)),
        )
        protein = ad.AnnData(
            X=protein_matrix,
            obs=obs.copy(),
            var=pd.DataFrame(index=pd.Index(protein_features)),
        )
        rna.obsm["X_umap"] = np.asarray(
            rna_source.obsm["X_umap"], dtype=np.float32
        )[positions, :2]
        protein.obsm["X_pca"] = (
            np.asarray(protein_source.obsm["X_pca"], dtype=np.float32)[positions, :2]
            if "X_pca" in protein_source.obsm
            else _principal_coordinates(protein_matrix)
        )
    finally:
        backed.file.close()

    result = md.MuData({"rna": rna, "protein": protein})
    return result, rna_features[0], protein_features[0], selected


def _fallback_liana(seed: int = 811) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    groups = ("NK", "B", "Mono", "CD4 T", "CD8 T", "DC")
    pairs = (
        ("CCL5", "CCR5"),
        ("IL7", "IL7R"),
        ("TNF", "TNFRSF1A"),
        ("CXCL8", "CXCR2"),
        ("IL7", "IL2RG"),
    )
    records = []
    for index in range(90):
        source = groups[index % len(groups)]
        target = groups[(index * 5 + 2) % len(groups)]
        ligand, receptor = pairs[index % len(pairs)]
        magnitude = float(rng.uniform(0.001, 0.09))
        specificity = float(rng.uniform(0.001, 0.09))
        records.append(
            {
                "source": source,
                "target": target,
                "ligand": ligand,
                "receptor": receptor,
                "cellphone_pvals": float(rng.uniform(0.0, 0.06)),
                "specificity_rank": specificity,
                "magnitude_rank": magnitude,
            }
        )
    result = pd.DataFrame.from_records(records)
    result.index = pd.Index(
        [f"lr-interaction-{index:06d}" for index in range(len(result))],
        name="interaction_id",
    )
    return result


def load_liana_long(path: str | Path | None = None) -> tuple[pd.DataFrame, Path | None]:
    """Load PBMC LIANA as one atomic long table indexed by interaction ID."""

    repository = Path(__file__).resolve().parents[2]
    selected = _first_existing(
        path,
        "GUANACO_LIANA_H5AD",
        [
            repository / "data" / "pbmc_liana_demo.h5ad",
            Path("/Users/xiyuanzhang/Documents/GUANACO_v2/data/pbmc_liana_demo.h5ad"),
        ],
    )
    if selected is None:
        result = _fallback_liana()
    else:
        from guanaco.data.ligand_receptor import (
            ligand_receptor_frame,
            load_default_ligand_receptor_result,
        )

        backed = ad.read_h5ad(selected, backed="r")
        try:
            payload = load_default_ligand_receptor_result(backed)
            if payload is None:
                raise ValueError(
                    f"{selected} does not contain a supported LIANA result."
                )
            result = ligand_receptor_frame(payload).copy()
        finally:
            backed.file.close()
        result.index = pd.Index(
            [f"lr-interaction-{index:06d}" for index in range(len(result))],
            name="interaction_id",
        )

    required = ("source", "target", "ligand", "receptor")
    missing = [column for column in required if column not in result]
    if missing:
        raise ValueError(f"LIANA table is missing {missing[0]!r}.")
    for column in ("cellphone_pvals", "specificity_rank", "magnitude_rank"):
        if column not in result:
            result[column] = 1.0
        result[column] = pd.to_numeric(result[column], errors="coerce").fillna(1.0)
    result["cell_pair"] = (
        result["source"].astype(str) + " → " + result["target"].astype(str)
    )
    result["ligand_receptor"] = (
        result["ligand"].astype(str) + " → " + result["receptor"].astype(str)
    )
    result["confidence"] = 1.0 - result["magnitude_rank"].clip(0.0, 1.0)
    return result, selected


def significant_liana(
    interactions: pd.DataFrame,
    *,
    pvalue: float = 0.05,
    rank: float = 0.02,
) -> pd.DataFrame:
    """One reusable significance filter for the LIANA network demo."""

    keep = (
        interactions["cellphone_pvals"].le(pvalue)
        & interactions["specificity_rank"].le(rank)
        & interactions["magnitude_rank"].le(rank)
    )
    result = interactions.loc[keep].sort_values(
        ["magnitude_rank", "specificity_rank"], kind="stable"
    )
    if result.empty:
        raise ValueError(
            "No LIANA interactions pass "
            f"cellphone_pvals <= {pvalue:g}, specificity_rank <= {rank:g}, "
            f"and magnitude_rank <= {rank:g}."
        )
    return result.copy()


def _fallback_spatial(seed: int = 824) -> tuple:
    rng = np.random.default_rng(seed)
    groups = np.asarray(["A", "B", "C", "D", "E"], dtype=object)
    centers = np.asarray([(0, 0), (4, 0), (1, 4), (5, 4), (8, 2)], dtype=float)
    labels = np.repeat(groups, 80)
    coordinates = np.vstack(
        [centers[index] + rng.normal(0, 0.72, size=(80, 2)) for index in range(5)]
    )
    coordinates = coordinates * np.asarray([72.0, 68.0]) + np.asarray([95.0, 80.0])
    spot_ids = np.asarray([f"spot-{index:04d}" for index in range(len(labels))])
    image = np.broadcast_to(
        np.asarray([0.92, 0.86, 0.90], dtype=np.float32), (600, 800, 3)
    ).copy()
    count = np.zeros((len(groups), len(groups)), dtype=int)
    zscore = np.full(count.shape, -1.0)
    for index in range(len(groups) - 1):
        count[index, index + 1] = count[index + 1, index] = 60 - index * 6
        zscore[index, index + 1] = zscore[index + 1, index] = 8.0 - index
    return coordinates, spot_ids, labels, image, list(groups), count, zscore, count * 54.0


def _spatial_components(
    path: str | Path | None,
    *,
    cluster_key: str,
) -> dict[str, object]:
    """Load the spatial arrays shared by the linked spatial demos."""

    repository = Path(__file__).resolve().parents[2]
    selected = _first_existing(
        path,
        "GUANACO_SPATIAL_H5AD",
        [
            repository / "data" / "anndata" / "visium_hne_spatial.h5ad",
            Path("/Users/xiyuanzhang/Documents/GUANACO_v2/data/visium_hne_spatial.h5ad"),
        ],
    )
    if selected is None:
        (
            coordinates, spot_ids, labels, tissue_image, categories,
            counts, enrichment, distance_sums,
        ) = _fallback_spatial()
        distances = np.linspace(1.0, 49.0, 49)
        occurrence = np.ones((len(categories), len(categories), len(distances)))
        for i in range(len(categories)):
            for j in range(len(categories)):
                if counts[i, j] > 0:
                    occurrence[i, j] += 0.45 * np.exp(-distances / (8 + i + j))
        palette = dict(zip(
            categories,
            ("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD"),
            strict=False,
        ))
        spatial_adata = ad.AnnData(
            X=np.zeros((len(spot_ids), 1), dtype=np.float32),
            obs=pd.DataFrame(
                {cluster_key: pd.Categorical(labels, categories=categories)},
                index=pd.Index(np.asarray(spot_ids, dtype=str)),
            ),
            var=pd.DataFrame(index=["spatial_placeholder"]),
        )
        spatial_adata.obsm["spatial"] = np.asarray(coordinates, dtype=np.float32)
        spatial_adata.uns["spatial"] = {
            "demo": {
                "images": {"lowres": np.asarray(tissue_image)},
                "scalefactors": {"tissue_lowres_scalef": 1.0},
            }
        }
        return {
            "adata": spatial_adata,
            "coordinates": coordinates,
            "spot_ids": spot_ids,
            "labels": labels,
            "tissue_image": tissue_image,
            "categories": categories,
            "counts": np.asarray(counts),
            "enrichment": np.asarray(enrichment),
            "distance_sums": np.asarray(distance_sums),
            "occurrence": occurrence,
            "distance": distances,
            "palette": palette,
            "selected": None,
        }

    spatial_adata = ad.read_h5ad(selected)
    neighborhood = spatial_adata.uns[f"{cluster_key}_nhood_enrichment"]
    occurrence_result = spatial_adata.uns.get(f"{cluster_key}_co_occurrence")
    library = next(iter(spatial_adata.uns["spatial"].values()))
    images, scalefactors = library["images"], library["scalefactors"]
    image_key = "lowres" if "lowres" in images else "hires"
    tissue_image = np.asarray(images[image_key], dtype=float)
    scale = float(scalefactors.get(f"tissue_{image_key}_scalef", 1.0))
    coordinates = np.asarray(spatial_adata.obsm["spatial"], dtype=float)[:, :2] * scale
    spot_ids = np.asarray(spatial_adata.obs_names.astype(str), dtype=object)
    groups = spatial_adata.obs[cluster_key]
    labels = groups.astype(str).to_numpy(dtype=object)
    categories = (
        list(map(str, groups.cat.categories))
        if isinstance(groups.dtype, pd.CategoricalDtype)
        else list(pd.unique(labels))
    )
    counts = np.asarray(neighborhood["count"])
    enrichment = np.asarray(neighborhood["zscore"])
    graph = spatial_adata.obsp["spatial_distances"].tocoo()
    codes = pd.Categorical(labels, categories=categories).codes.astype(np.int64)
    pair_codes = codes[graph.row] * len(categories) + codes[graph.col]
    distance_sums = np.bincount(
        pair_codes, weights=graph.data, minlength=len(categories) ** 2
    ).reshape(len(categories), len(categories))
    if (
        isinstance(occurrence_result, dict)
        and occurrence_result.get("occ") is not None
        and occurrence_result.get("interval") is not None
    ):
        occurrence = np.asarray(occurrence_result.get("occ"), dtype=float)
        distance = np.asarray(occurrence_result.get("interval"), dtype=float)[1:]
    else:
        occurrence = None
        distance = None
    colors = list(map(str, spatial_adata.uns.get(f"{cluster_key}_colors", ())))
    palette = dict(zip(categories, colors, strict=False))
    return {
        "adata": spatial_adata,
        "coordinates": coordinates,
        "spot_ids": spot_ids,
        "labels": labels,
        "tissue_image": tissue_image,
        "categories": categories,
        "counts": counts,
        "enrichment": enrichment,
        "distance_sums": distance_sums,
        "occurrence": occurrence,
        "distance": distance,
        "palette": palette,
        "selected": selected,
    }


def _spatial_relationship_tables(
    components: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create the pair overview and distance-level detail tables."""

    categories = list(components["categories"])
    counts = np.asarray(components["counts"])
    enrichment = np.asarray(components["enrichment"], dtype=float)
    distance_sums = np.asarray(components["distance_sums"], dtype=float)
    occurrence = components["occurrence"]
    distance = components["distance"]
    if occurrence is None or distance is None:
        raise ValueError("Spatial co-occurrence results are required for this demo.")
    occurrence = np.asarray(occurrence, dtype=float)
    distance = np.asarray(distance, dtype=float)
    if occurrence.shape[:2] != (len(categories), len(categories)):
        raise ValueError("Spatial co-occurrence groups do not match cluster groups.")
    if occurrence.shape[2] != len(distance):
        raise ValueError("Spatial co-occurrence distances do not match intervals.")

    pairs = []
    for i, source in enumerate(categories):
        for j, target in enumerate(categories):
            pair_id = f"{source}::{target}"
            pairs.append(
                {
                    "pair_id": pair_id,
                    "source_group": source,
                    "target_group": target,
                    "enrichment": float(enrichment[i, j]),
                    "neighbor_count": int(counts[i, j]),
                    "total_distance": float(distance_sums[i, j]),
                }
            )
    pair_table = pd.DataFrame(pairs).set_index("pair_id", drop=False)

    curves = []
    for pair in pair_table.itertuples(index=False):
        i = categories.index(pair.source_group)
        j = categories.index(pair.target_group)
        label = f"{pair.source_group} → {pair.target_group}"
        for index, value in enumerate(occurrence[i, j, :]):
            curves.append(
                {
                    "curve_id": f"{pair.pair_id}::{index}",
                    "pair_id": pair.pair_id,
                    "pair_label": label,
                    "conditional_group": pair.source_group,
                    "observed_group": pair.target_group,
                    "distance": float(distance[index]),
                    "co_occurrence": float(value),
                }
            )
    curve_table = pd.DataFrame(curves).set_index("curve_id")
    return pair_table, curve_table


def _spatial_pair_cells(
    components: dict[str, object], pair_table: pd.DataFrame
) -> pd.DataFrame:
    """Make heatmap marks carry the cell IDs consumed by ``gc.pl.spatial``."""

    spot_ids = np.asarray(components["spot_ids"], dtype=object).astype(str)
    labels = np.asarray(components["labels"], dtype=object).astype(str)
    rows: list[dict[str, object]] = []
    for pair in pair_table.itertuples(index=False):
        in_pair = (labels == pair.source_group) | (labels == pair.target_group)
        for position in np.flatnonzero(in_pair):
            rows.append(
                {
                    "membership_id": f"{pair.pair_id}::{spot_ids[position]}",
                    "pair_id": pair.pair_id,
                    "cell_id": spot_ids[position],
                    "source_group": pair.source_group,
                    "target_group": pair.target_group,
                    "enrichment": pair.enrichment,
                    "neighbor_count": pair.neighbor_count,
                    "total_distance": pair.total_distance,
                }
            )
    result = pd.DataFrame(rows).set_index("membership_id")
    return result


def load_spatial_relationship_demo(
    path: str | Path | None = None,
    *,
    cluster_key: str = "cluster",
) -> tuple[pd.DataFrame, ad.AnnData, pd.DataFrame, Path | None]:
    """Load a pair table, spatial AnnData, and co-occurrence detail table."""

    components = _spatial_components(path, cluster_key=cluster_key)
    pair_summary, curve_table = _spatial_relationship_tables(components)
    pair_cells = _spatial_pair_cells(components, pair_summary)
    return pair_cells, components["adata"], curve_table, components["selected"]


def load_pbmc_cd4_relationship(
    path: str | Path | None = None,
    *,
    mdata=None,
    max_cells: int = 8_000,
    seed: int = 930,
):
    """Build a shared-cell RNA CD4/protein CD4 relationship table."""

    if mdata is None:
        mdata, rna_feature, protein_feature, selected = load_pbmc_mudata_demo(
            path, max_cells=max_cells, seed=seed
        )
    else:
        selected = Path(path).expanduser().resolve() if path is not None else None
        protein_names = list(map(str, mdata.mod["protein"].var_names))
        rna_feature = next(
            (str(name) for name in mdata.mod["rna"].var_names if str(name) == "CD4"),
            str(mdata.mod["rna"].var_names[0]),
        )
        protein_feature = next(
            (
                feature
                for feature in protein_names
                if feature.lower().startswith("adt_cd4-")
                or feature.lower().startswith("cd4")
            ),
            protein_names[0],
        )
    rna = mdata.mod["rna"]
    protein = mdata.mod["protein"]
    protein_candidates = [
        str(feature)
        for feature in protein.var_names
        if str(feature).lower().startswith("adt_cd4-")
    ]
    if not protein_candidates:
        protein_candidates = [protein_feature]

    rna_values = _dense(rna[:, [rna_feature]].X).reshape(-1).astype(float)
    protein_values = np.mean(
        np.vstack(
            [
                _dense(protein[:, [feature]].X).reshape(-1).astype(float)
                for feature in protein_candidates
            ]
        ),
        axis=0,
    )
    group_key = "celltype_l1" if "celltype_l1" in rna.obs else "cell_type"
    level1 = rna.obs[group_key].astype(str).to_numpy()
    correlations: dict[str, float] = {}
    for group in pd.unique(level1):
        in_group = level1 == group
        rho = stats.spearmanr(
            rna_values[in_group], protein_values[in_group], nan_policy="omit"
        ).statistic
        correlations[str(group)] = float(0.0 if not np.isfinite(rho) else rho)
    cell_correlations = np.asarray([correlations[str(group)] for group in level1])
    labels = np.asarray(
        [f"{group} · ρ={correlations[str(group)]:.2f}" for group in level1],
        dtype=object,
    )
    cell_ids = pd.Index(rna.obs_names.astype(str), name="cell_id")
    relationship = pd.DataFrame(
        {
            "rna_cd4": rna_values,
            "protein_cd4": protein_values,
            "level1": level1,
            "level1_spearman": cell_correlations,
            "correlation_label": labels,
        },
        index=cell_ids,
    )
    embedding = rna.copy()
    embedding.obs["level1"] = pd.Categorical(level1)
    embedding.obs["level1_spearman"] = cell_correlations
    embedding.obs["correlation_label"] = pd.Categorical(labels)
    summary = (
        relationship.groupby("level1", as_index=False, sort=False)
        .agg(
            n_cells=("rna_cd4", "size"),
            spearman=("level1_spearman", "first"),
            mean_rna_cd4=("rna_cd4", "mean"),
            mean_protein_cd4=("protein_cd4", "mean"),
        )
        .sort_values("spearman", ascending=False, ignore_index=True)
    )
    protein_label = "mean(" + ", ".join(protein_candidates) + ")"
    return embedding, relationship, summary, selected, rna_feature, protein_label


__all__ = [
    "MARKERS",
    "load_liana_long",
    "load_pbmc_cd4_relationship",
    "load_pbmc_mudata_demo",
    "load_spatial_relationship_demo",
    "make_external_cell_table",
    "make_pathway_demo",
    "make_single_cell",
    "make_volcano_data",
    "significant_liana",
]
