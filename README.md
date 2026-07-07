# Multi-Cohort Breast Cancer Survival Prediction

A machine-learning pipeline for breast cancer overall-survival prediction that integrates **11 public gene-expression cohorts (3,818 patients)** into a single harmonized survival-analysis workflow, benchmarks three gene-expression survival models against a clinical Cox baseline, evaluates cross-cohort generalization using leave-one-cohort-out validation, and validates the final signature on an independent external cohort.

**Headline result:** a sparse **284-gene Elastic-Net Cox signature** achieves **C-index 0.716** on a fully held-out external cohort, **SCAN-B / GSE96058**, generalizing across both cohort and expression-platform differences.

---

## Why This Project

Many breast-cancer prognostic signatures are developed and validated in limited cohort settings, often with incomplete evaluation of cross-cohort and cross-population transfer.

This project asks two main questions:

1. **Does a survival signature generalize across cohorts and platforms**, or does it silently overfit to one dataset?
2. **Does the prognostic signal transfer consistently across patient groups and study populations?**

Answering these questions requires more than fitting a model. It requires building a carefully harmonized multi-cohort dataset, controlling batch effects, using fair validation strategies, and testing the final frozen model on a fully external cohort.

The core contribution of this project is therefore the **end-to-end survival prediction pipeline**: data harmonization, batch correction, feature filtering, model benchmarking, leave-one-cohort-out validation, and external validation.

---

## Results

### Model benchmark

| Model | Same-split test C-index | LOCO event-weighted C-index | Notes |
|---|---:|---:|---|
| Clinical Cox baseline | 0.633 | — | Receptor-stratified clinical baseline |
| **Elastic-Net Cox** | **0.651** | **~0.638** | Final selected model |
| DeepSurv | 0.647 | 0.636 | Neural Cox benchmark |
| Random Survival Forest | 0.635 | 0.621 | Non-linear tree-based benchmark |

The flexible models, DeepSurv and Random Survival Forest, did not meaningfully outperform the regularized Elastic-Net Cox model. Elastic-Net Cox was retained as the final model because it was competitive, sparse, interpretable, robust under leave-one-cohort-out validation, and externally validated.

---

## External Validation on SCAN-B / GSE96058

The frozen Elastic-Net Cox model was validated on **SCAN-B / GSE96058**, an independent RNA-seq breast cancer cohort that was not used during model training or tuning.

### SCAN-B validation summary

| Metric | Result |
|---|---:|
| Final validation patients | 3,273 |
| Events | 336 |
| Censored patients | 2,937 |
| Frozen model genes | 284 |
| Frozen genes available in SCAN-B | 276 |
| External C-index | **0.716** |
| Log-rank p-value | ~1.1e-31 |

### Risk-group separation in SCAN-B

Patients were divided into low, intermediate, and high-risk groups using risk-score tertiles.

| Risk group | Patients | Events | Event rate |
|---|---:|---:|---:|
| Low | 1,091 | 34 | 3.1% |
| Intermediate | 1,091 | 98 | 9.0% |
| High | 1,091 | 204 | 18.7% |

The event rate increased clearly from low-risk to high-risk groups, showing that the frozen signature preserved prognostic separation in an independent external cohort.

---

## Clinical Comparison on SCAN-B

The frozen gene-expression signature was also compared against standard clinical variables in SCAN-B.

Clinical variables included:

- Age
- ER status
- PgR status
- HER2 status
- Nottingham grade

### Complete-case SCAN-B clinical comparison

| Model | C-index |
|---|---:|
| Signature only | 0.714 |
| Clinical only | 0.762 |
| Clinical + signature | **0.793** |

The combined **clinical + signature** model performed best. This suggests that the gene-expression signature does not replace clinical variables, but adds complementary prognostic information.

The signature remained independently prognostic after clinical adjustment:

| Measure | Result |
|---|---:|
| Adjusted HR | 3.38 |
| 95% CI | 2.41–4.74 |
| p-value | ~1.9e-12 |


---

## Pipeline

```text
Expression data from 11 cohorts
      │
      │  Platform-specific preprocessing
      │  Probe-to-gene mapping / collapsing where required
      │  Gene-symbol harmonization
      │  Common-gene merge
      ▼
Merged expression matrix
      │
      │  Sparse missing-value handling
      │  ComBat batch correction using cohort as batch
      ▼
Survival-matched expression matrix
      │
      │  3,818 patients with os_time / os_event
      │  MAD filter → top 5,000 genes
      ▼
Model benchmarking
      │
      ├─► Elastic-Net Cox
      ├─► Random Survival Forest
      ├─► DeepSurv
      └─► Clinical Cox baseline
              │
              ├─► Same-split test evaluation
              ├─► Leave-one-cohort-out validation
              └─► External validation on SCAN-B
```

---

## Reproducing the Analysis

Scripts are organized in numbered order. Set your local paths in `scripts/00_config_template.R` before running the pipeline.

| Step | Script | Purpose |
|---:|---|---|
| 0 | `00_config_template.R` | Path and parameter template |
| 1 | `01_merge_and_combat.R` | Merge cohorts and run ComBat batch correction |
| 2 | `02_clinical_harmonization_template.R` | Harmonize survival and clinical metadata |
| 3 | `03_mad_gene_filtering.R` | Select top MAD-ranked genes |
| 4 | `04_elastic_net_cox.R` | Train Elastic-Net Cox model |
| 5 | `05_univariate_cox_filter.R` | In-fold univariate Cox filtering |
| 6 | `06_random_survival_forest.R` | Random Survival Forest benchmark |
| 7 | `07_export_for_deepsurv.R` | Export expression and clinical data for DeepSurv |
| 8 | `08_deepsurv_train_test.py` | DeepSurv same-split benchmark |
| 9 | `09_elastic_net_loco_validation.R` | Elastic-Net leave-one-cohort-out validation |
| 10 | `10_scanb_external_validation.R` | External validation on SCAN-B |
| 11 | `11_clinical_model_comparison.R` | Signature vs clinical model comparison |

Environment setup files are available in:

```text
environment/r_packages.txt
environment/python_requirements.txt
```

Methodology and interpretation details are available in:

```text
docs/
```

---

## Cohorts

### Discovery cohorts

The integrated discovery dataset contains 11 breast cancer cohorts:

- METABRIC
- TCGA-BRCA
- CAL
- GSE1456
- GSE7390
- GSE20685
- GSE20711
- GSE42568
- GSE58812
- GSE88770
- GSE162228

### External validation cohort

- SCAN-B / GSE96058

---

## Data Availability

Raw expression matrices, clinical metadata, patient-level risk scores, and trained model objects are **not included** in this repository due to file size, dataset licensing, and patient-level privacy considerations.

All source cohorts are publicly accessible through their original repositories, including GEO, cBioPortal, UCSC Xena, and SCAN-B / GSE96058.

This repository contains:

- clean pipeline scripts
- documentation
- aggregate non-patient-level result summaries
- reproducibility templates

It does **not** contain:

- raw `.CEL` files
- expression matrices
- clinical metadata tables
- patient-level predictions
- patient-level risk scores
- trained `.rds` model objects

See `data/README.md` for data-access notes.

---

## Repository Structure

```text
scripts/                 Numbered R/Python pipeline scripts
docs/                    Methodology, model, and validation notes
results/summary_tables/  Aggregate non-patient-level results
data/README.md           Data access and privacy note
environment/             R and Python package lists
```

---

## Tech Stack

**R / Bioconductor:** `survival`, `glmnet`, `sva` / ComBat, `randomForestSRC`, `matrixStats`

**Python:** DeepSurv-style neural Cox modelling using `pycox`, `torchtuples`, and PyTorch

**Infrastructure:** Bash, SLURM / HPC

---

## Status and Roadmap

### Complete

- Multi-cohort expression preprocessing and integration
- Clinical survival metadata harmonization
- ComBat batch correction
- MAD-based feature filtering
- Elastic-Net Cox model development
- Random Survival Forest benchmark
- DeepSurv benchmark
- Leave-one-cohort-out validation
- External validation on SCAN-B
- Clinical vs signature comparison on SCAN-B

### In progress

- Higher-powered cross-cohort ancestry-transfer analysis
- Manuscript preparation

---

## Research Context

This project was conducted as part of a breast cancer survival modelling research project under **Prof. Nita Parekh, IIIT Hyderabad**.