"""Prepare the public CLL cohort as a patient-level GUANACO MuData example.

The source cohort contains four partially overlapping patient-level matrices:
mRNA expression, DNA methylation, ex-vivo drug viability, and genomic
alterations. Each AnnData modality keeps only patients with an actual
measurement; MuData-level ``obs`` stores the 200-patient union. The original
matrices are distributed as an RDS file, so this
script uses the system ``Rscript`` executable only to export those matrices to
temporary TSV files.  All annotation, statistics, embeddings, and H5MU writing
are performed in Python.

Source study:
    Dietrich et al., Drug-perturbation-based stratification of blood cancer,
    J Clin Invest (2018), PMID 29227286.

Processed data:
    Bioconductor MOFAdata / the public MOFA CLL tutorial resources.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import requests
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
from umap import UMAP


mu.set_options(pull_on_update=False)


BASE_URL = "https://usda-ree-ars.github.io/SEAStatsData/omics"
SOURCES = {
    "CLL_data.rds": f"{BASE_URL}/CLL_data.rds",
    "sample_metadata.txt": f"{BASE_URL}/sample_metadata.txt",
    "drugs.txt.gz": f"{BASE_URL}/drugs.txt.gz",
    "Hsapiens_genes_BioMart.87.txt.gz": (
        f"{BASE_URL}/Hsapiens_genes_BioMart.87.txt.gz"
    ),
}
SOURCE_STUDY = "https://doi.org/10.1172/JCI93801"
MOFA_TUTORIAL = (
    "https://usda-ree-ars.github.io/SEAStats/"
    "multiomics_demo/Demo2_MOFA_chronicleukemia.html"
)

def _download(url: str, destination: Path, *, force: bool) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with partial.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    partial.replace(destination)


def _export_rds(rds_path: Path, output_dir: Path) -> None:
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError(
            "Rscript is required to read the public CLL_data.rds source file."
        )
    helper = Path(__file__).with_name("prepare_cll_rds_export.R")
    if not helper.exists():
        raise RuntimeError(f"Missing RDS export helper: {helper}")
    subprocess.run([rscript, str(helper), str(rds_path), str(output_dir)], check=True)


def _read_view(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="gzip", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.T.apply(pd.to_numeric, errors="coerce")


def _status(values: pd.Series, present: str, absent: str) -> pd.Categorical:
    numeric = pd.to_numeric(values, errors="coerce")
    labels = np.full(len(numeric), "Unknown", dtype=object)
    labels[numeric.eq(1).fillna(False).to_numpy()] = present
    labels[numeric.eq(0).fillna(False).to_numpy()] = absent
    return pd.Categorical(labels, categories=[absent, present, "Unknown"])


def _yes_no_unknown(values: pd.Series) -> pd.Categorical:
    labels = np.full(len(values), "Unknown", dtype=object)
    normalized = values.astype("string").str.strip().str.lower()
    labels[normalized.isin(["true", "1", "yes"]).fillna(False).to_numpy()] = "Yes"
    labels[normalized.isin(["false", "0", "no"]).fillna(False).to_numpy()] = "No"
    return pd.Categorical(labels, categories=["No", "Yes", "Unknown"])


def _availability(frame: pd.DataFrame) -> pd.Categorical:
    available = frame.notna().any(axis=1)
    return pd.Categorical(
        np.where(available, "Available", "Missing"),
        categories=["Available", "Missing"],
    )


def _build_obs(
    metadata: pd.DataFrame,
    matrices: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    samples = matrices["rna"].index
    raw = metadata.set_index("sample").reindex(samples)
    obs = pd.DataFrame(index=pd.Index(samples, name="patient_id"))

    obs["sex"] = pd.Categorical(
        raw["Gender"].map({"m": "Male", "f": "Female"}).fillna("Unknown"),
        categories=["Female", "Male", "Unknown"],
    )
    obs["age_at_sampling"] = pd.to_numeric(raw["age"], errors="coerce")
    obs["age_group"] = pd.Categorical(
        pd.cut(
            obs["age_at_sampling"],
            bins=[-np.inf, 59, 69, 79, np.inf],
            labels=["<60", "60–69", "70–79", "≥80"],
        ).astype("string").fillna("Unknown"),
        categories=["<60", "60–69", "70–79", "≥80", "Unknown"],
        ordered=True,
    )
    obs["IGHV_status"] = _status(raw["IGHV"], "Mutated", "Unmutated")
    obs["trisomy12_status"] = _status(
        raw["trisomy12"], "Present", "Absent"
    )
    obs["treated_after_sampling"] = _yes_no_unknown(raw["treatedAfter"])
    died = _yes_no_unknown(raw["died"])
    obs["vital_status"] = pd.Categorical(
        pd.Series(died.astype(str), index=obs.index)
        .map({"Yes": "Died", "No": "Alive"})
        .fillna("Unknown"),
        categories=["Alive", "Died", "Unknown"],
    )
    obs["time_to_treatment_years"] = pd.to_numeric(raw["TTT"], errors="coerce")
    obs["time_to_death_or_followup_years"] = pd.to_numeric(
        raw["TTD"], errors="coerce"
    )

    mutation = matrices["mutation"]
    for feature, column in {
        "TP53": "TP53_status",
        "NOTCH1": "NOTCH1_status",
        "SF3B1": "SF3B1_status",
        "ATM": "ATM_status",
        "del17p13": "del17p13_status",
    }.items():
        obs[column] = _status(mutation[feature], "Altered", "Not altered")

    for key, frame in matrices.items():
        obs[f"{key}_data"] = _availability(frame)
    return obs


def _gene_var(
    raw_names: pd.Index,
    annotation: pd.DataFrame,
    mutation_features: set[str],
) -> pd.DataFrame:
    annotation = annotation.drop_duplicates("ens_id").set_index("ens_id")
    var = annotation.reindex(raw_names).copy()
    symbols = var["symbol"].astype("string").str.strip()
    symbols = symbols.mask(symbols.eq("") | symbols.isna(), raw_names.astype(str))
    symbols = symbols.map(
        lambda name: f"{name}_mRNA" if str(name) in mutation_features else str(name)
    )
    var["ensembl_id"] = raw_names.astype(str)
    var["feature_name"] = symbols.astype(str)
    var["feature_type"] = "gene expression"
    var.index = pd.Index(symbols.astype(str), name="feature")
    return var


def _drug_var(raw_names: pd.Index, annotation: pd.DataFrame) -> pd.DataFrame:
    annotation = annotation.drop_duplicates("drug_id").set_index("drug_id")
    records = []
    names = []
    pattern = re.compile(r"^(D_\d+)_(\d+)$")
    for raw_name in raw_names.astype(str):
        match = pattern.match(raw_name)
        drug_id = match.group(1) if match else raw_name
        dose_index = int(match.group(2)) if match else np.nan
        details = annotation.loc[drug_id] if drug_id in annotation.index else None
        drug = str(details["name"]) if details is not None else drug_id
        names.append(f"{drug}_{dose_index:g}" if np.isfinite(dose_index) else drug)
        records.append(
            {
                "source_feature": raw_name,
                "drug_id": drug_id,
                "drug": drug,
                "dose_index": dose_index,
                "main_targets": (
                    str(details["main_targets"]) if details is not None else ""
                ),
                "target_category": (
                    str(details["target_category"]) if details is not None else ""
                ),
                "pathway": str(details["pathway"]) if details is not None else "",
                "feature_type": "ex-vivo drug viability",
            }
        )
    return pd.DataFrame(records, index=pd.Index(names, name="feature"))


def _simple_var(raw_names: pd.Index, feature_type: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_feature": raw_names.astype(str),
            "feature_type": feature_type,
        },
        index=pd.Index(raw_names.astype(str), name="feature"),
    )


def _make_unique_names(var: pd.DataFrame) -> pd.DataFrame:
    holder = ad.AnnData(X=np.empty((0, len(var))), var=var.copy())
    holder.var_names_make_unique(join="-")
    return holder.var.copy()


def _imputed_scaled(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64).copy()
    medians = np.nanmedian(values, axis=0)
    medians[~np.isfinite(medians)] = 0.0
    missing_rows, missing_cols = np.where(~np.isfinite(values))
    values[missing_rows, missing_cols] = medians[missing_cols]
    varying = np.nanstd(values, axis=0) > np.finfo(float).eps
    if not np.any(varying):
        raise ValueError("A modality has no variable features after imputation.")
    return StandardScaler().fit_transform(values[:, varying])


def _embed(matrix: np.ndarray, *, seed: int) -> tuple[np.ndarray, np.ndarray]:
    scaled = _imputed_scaled(matrix)
    n_components = min(30, scaled.shape[0] - 1, scaled.shape[1])
    pca = PCA(n_components=n_components, random_state=seed).fit_transform(scaled)
    umap = UMAP(
        n_neighbors=15,
        min_dist=0.25,
        metric="euclidean",
        random_state=seed,
        transform_seed=seed,
    ).fit_transform(pca[:, : min(20, pca.shape[1])])
    return pca.astype(np.float32), np.asarray(umap, dtype=np.float32)


def _differential_entry(
    matrix: pd.DataFrame,
    group: pd.Series,
    group_a: str,
    group_b: str,
) -> dict[str, list]:
    labels = group.astype("string")
    values_a = matrix.loc[labels.eq(group_a).fillna(False)].to_numpy(dtype=float)
    values_b = matrix.loc[labels.eq(group_b).fillna(False)].to_numpy(dtype=float)
    result = ttest_ind(values_b, values_a, axis=0, equal_var=False, nan_policy="omit")
    score = np.asarray(result.statistic, dtype=float)
    pvalue = np.asarray(result.pvalue, dtype=float)
    valid_p = np.isfinite(pvalue)
    padj = np.ones_like(pvalue)
    if np.any(valid_p):
        padj[valid_p] = multipletests(pvalue[valid_p], method="fdr_bh")[1]
    mean_a = np.nanmean(values_a, axis=0)
    mean_b = np.nanmean(values_b, axis=0)
    return {
        "gene": matrix.columns.astype(str).tolist(),
        "logfoldchange": (mean_b - mean_a).tolist(),
        "score": score.tolist(),
        "padj": padj.tolist(),
        "mean_group_a": mean_a.tolist(),
        "mean_group_b": mean_b.tolist(),
        "group_a": [group_a] * matrix.shape[1],
        "group_b": [group_b] * matrix.shape[1],
    }


def _build_modalities(
    matrices: dict[str, pd.DataFrame],
    obs: pd.DataFrame,
    gene_annotation: pd.DataFrame,
    drug_annotation: pd.DataFrame,
    *,
    seed: int,
) -> dict[str, ad.AnnData]:
    gene_var = _make_unique_names(
        _gene_var(
            matrices["rna"].columns,
            gene_annotation,
            set(matrices["mutation"].columns.astype(str)),
        )
    )
    rna = matrices["rna"].copy()
    rna.columns = gene_var.index

    drug_var = _make_unique_names(_drug_var(matrices["drug"].columns, drug_annotation))
    drug = matrices["drug"].copy()
    drug.columns = drug_var.index

    prepared = {
        "rna": (rna, gene_var),
        "methylation": (
            matrices["methylation"],
            _simple_var(matrices["methylation"].columns, "CpG methylation M-value"),
        ),
        "drug": (drug, drug_var),
        "mutation": (
            matrices["mutation"],
            _simple_var(matrices["mutation"].columns, "binary genomic alteration"),
        ),
    }

    modalities: dict[str, ad.AnnData] = {}
    for offset, (name, (frame, var)) in enumerate(prepared.items()):
        measured_frame = frame.loc[frame.notna().any(axis=1)]
        adata = ad.AnnData(
            X=measured_frame.to_numpy(dtype=np.float32),
            obs=pd.DataFrame(index=measured_frame.index.copy()),
            var=var.copy(),
        )
        local_pca, local_umap = _embed(adata.X, seed=seed + offset)
        adata.obsm["X_pca"] = local_pca
        adata.obsm["X_umap"] = local_umap
        adata.uns["measurement"] = {
            "observation_type": "patient",
            "modality": name,
            "missing_values_preserved": True,
            "patients_in_cohort": int(len(obs)),
            "patients_measured": int(adata.n_obs),
        }
        if name == "drug":
            adata.uns["measurement"]["value"] = "ex-vivo viability"
            adata.uns["measurement"]["interpretation"] = (
                "Lower viability indicates stronger drug activity."
            )
        modalities[name] = adata

    rna_frame = prepared["rna"][0].dropna(axis=0, how="all")
    rna_obs = obs.reindex(rna_frame.index)
    modalities["rna"].uns["volcano"] = {
        "entries": {
            "IGHV: Mutated − Unmutated": _differential_entry(
                rna_frame,
                rna_obs["IGHV_status"],
                "Unmutated",
                "Mutated",
            ),
            "Trisomy 12: Present − Absent": _differential_entry(
                rna_frame,
                rna_obs["trisomy12_status"],
                "Absent",
                "Present",
            ),
        },
        "default_entry": "IGHV: Mutated − Unmutated",
        "effect_definition": "difference of mean normalized expression",
        "test": "Welch t-test with Benjamini-Hochberg correction",
    }
    return modalities


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest(mdata: mu.MuData, output: Path, cache_dir: Path) -> dict:
    return {
        "output": str(output.resolve()),
        "n_patients": int(mdata.n_obs),
        "modalities": {
            name: {
                "n_patients": int(adata.n_obs),
                "n_features": int(adata.n_vars),
                "missing_values": int(np.isnan(np.asarray(adata.X)).sum()),
                "embeddings": list(adata.obsm.keys()),
            }
            for name, adata in mdata.mod.items()
        },
        "clinical_columns": list(mdata.obs.columns),
        "sources": {
            filename: {
                "url": url,
                "path": str((cache_dir / filename).resolve()),
                "sha256": _sha256(cache_dir / filename),
            }
            for filename, url in SOURCES.items()
        },
        "source_study": SOURCE_STUDY,
        "tutorial": MOFA_TUTORIAL,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Persistent source-file cache (defaults beside the output file).",
    )
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()

    output = args.output.expanduser().resolve()
    cache_dir = (
        args.cache_dir.expanduser().resolve()
        if args.cache_dir
        else output.parent / "cll_mofa_source"
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    for filename, url in SOURCES.items():
        print(f"Fetching {filename}")
        _download(url, cache_dir / filename, force=args.force_download)

    with tempfile.TemporaryDirectory(prefix="guanaco-cll-") as temp_name:
        temp_dir = Path(temp_name)
        _export_rds(cache_dir / "CLL_data.rds", temp_dir)
        matrices = {
            "rna": _read_view(temp_dir / "mRNA.tsv.gz"),
            "methylation": _read_view(temp_dir / "Methylation.tsv.gz"),
            "drug": _read_view(temp_dir / "Drugs.tsv.gz"),
            "mutation": _read_view(temp_dir / "Mutations.tsv.gz"),
        }

    sample_order = matrices["rna"].index
    for name, frame in matrices.items():
        if set(frame.index) != set(sample_order):
            raise ValueError(f"Modality '{name}' does not contain the same patients.")
        matrices[name] = frame.reindex(sample_order)

    metadata = pd.read_csv(cache_dir / "sample_metadata.txt", sep="\t")
    gene_annotation = pd.read_csv(
        cache_dir / "Hsapiens_genes_BioMart.87.txt.gz",
        sep="\t",
        compression="gzip",
    )
    drug_annotation = pd.read_csv(cache_dir / "drugs.txt.gz", compression="gzip")
    obs = _build_obs(metadata, matrices)
    modalities = _build_modalities(
        matrices,
        obs,
        gene_annotation,
        drug_annotation,
        seed=args.seed,
    )

    mdata = mu.MuData(modalities)
    mdata.obs = obs.reindex(mdata.obs_names).copy()
    mdata.uns["guanaco_example"] = {
        "title": "CLL patient-level clinical multi-omics",
        "description": (
            "Partially overlapping patient-level mRNA, methylation, ex-vivo "
            "drug viability, and genomic alteration matrices with shared "
            "clinical metadata."
        ),
        "observation_type": "patient",
        "source_study": SOURCE_STUDY,
        "processed_data": "Bioconductor MOFAdata / public MOFA CLL tutorial",
        "tutorial": MOFA_TUTORIAL,
        "drug_value_interpretation": (
            "Lower ex-vivo viability indicates stronger drug activity; these are "
            "experimental measurements, not administered treatments."
        ),
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(output, compression="gzip")
    manifest = _manifest(mdata, output, cache_dir)
    manifest_path = output.with_name(f"{output.stem}_manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    print(f"Saved manifest to {manifest_path}")


if __name__ == "__main__":
    main()
