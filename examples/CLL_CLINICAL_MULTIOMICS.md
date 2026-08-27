# CLL patient-level clinical multi-omics example

This example demonstrates GUANACO on patient-level matrices rather than cells.
It combines four partially overlapping measurements from a cohort of 200
chronic lymphocytic leukaemia (CLL) patients. Each modality contains only
patients with an actual measurement in that matrix:

| Modality | Shape | Measurement |
| --- | ---: | --- |
| `rna` | 136 × 5,000 | Normalized mRNA expression |
| `methylation` | 196 × 4,248 | CpG methylation M-values |
| `drug` | 184 × 310 | Ex-vivo viability after drug-dose perturbations |
| `mutation` | 200 × 69 | Binary mutations and copy-number alterations |

Shared patient annotations include age, sex, IGHV status, trisomy 12, selected
driver alterations, treatment-after-sampling status, vital status, and
time-to-event fields. The MuData-level observation table stores the 200-patient
union; modality tabs do not create all-missing rows for unmeasured patients.

## Generate the data

The source matrices are distributed as an RDS file. `Rscript` is required only
for the temporary export step; the Python script performs the remaining
annotation, statistics, embeddings, and H5MU writing.

RNA symbols that also occur in the mutation panel receive an `_mRNA` suffix,
following the MOFA tutorial convention, so cross-modality feature names remain
unambiguous.

```bash
python examples/scripts/prepare_cll_clinical_multiomics.py \
  --output /Users/xiyuanzhang/Documents/GUANACO_v2/data/cll_clinical_multiomics.h5mu
```

The command also creates:

- `cll_clinical_multiomics_manifest.json`, containing shapes, missing-value
  counts, clinical columns, embeddings, and source links;
- `cll_mofa_source/`, a reusable cache of the four downloaded source files.

Launch the prepared example with:

```bash
guanaco -c examples/configs/cll_clinical_multiomics.json
```

## Suggested exploration

1. Use the RNA embedding and colour the 136 measured patients by `IGHV_status`.
2. Inspect the precomputed RNA contrast `IGHV: Mutated − Unmutated` in the
   volcano view, then validate features such as `LPL`, `ZAP70`, `CD38`, and
   `TCL1A` with heatmaps and comparative violins.
3. In the drug modality, compare `ibrutinib_1`, `ibrutinib_3`, and
   `ibrutinib_5` by IGHV status and trisomy 12 status.
4. In Omics comparison, compare RNA and drug-response profiles across shared
   clinical groups while retaining each modality's actual patient coverage.
5. Filter by driver alterations such as `TP53_status`, `NOTCH1_status`, or
   `del17p13_status` and repeat the molecular and drug-response comparisons.

Drug values are experimental ex-vivo viability measurements, not records of
administered treatment. Lower viability indicates stronger measured drug
activity. Dose suffixes preserve the source order and should not be interpreted
as increasing concentration without consulting the original study protocol.

## Provenance

- Dietrich et al. *Drug-perturbation-based stratification of blood cancer*.
  J Clin Invest (2018). <https://doi.org/10.1172/JCI93801>
- Bioconductor `MOFAdata` processed CLL matrices.
- Public MOFA/SEAStats tutorial resources:
  <https://usda-ree-ars.github.io/SEAStats/multiomics_demo/Demo2_MOFA_chronicleukemia.html>
