# Multi-Cohort Breast Cancer Survival Prediction

This repository contains a clean, public-facing version of a machine-learning pipeline for breast cancer overall-survival prediction using multi-cohort gene-expression data.

## Overview

The project integrates public breast cancer expression cohorts into a harmonized survival-analysis workflow. The pipeline covers expression preprocessing, probe-to-gene mapping, common-gene merging, missing-value handling, ComBat batch correction, clinical metadata harmonization, feature filtering, survival-model training, leave-one-cohort-out validation, and external validation.

## Key Features

- Integrates 11 breast cancer cohorts into a common gene-expression space.
- Builds a survival-ready dataset with 10,435 genes and 3,818 patients.
- Applies sparse missing-value handling and ComBat batch correction.
- Uses MAD-based filtering to select high-variance genes.
- Benchmarks Elastic Net Cox, Random Survival Forest, and DeepSurv.
- Performs leave-one-cohort-out validation.
- Externally validates the frozen Elastic Net Cox model on SCAN-B / GSE96058.

## Models Compared

| Model | Type | Role |
|---|---|---|
| Elastic Net Cox | Regularized survival model | Main interpretable model |
| Random Survival Forest | Tree-based survival ensemble | Non-linear benchmark |
| DeepSurv | Neural Cox model | Deep-learning benchmark |

## Data Availability

Raw expression matrices, clinical metadata, patient-level risk scores, trained model objects, and private logs are not included in this public repository due to file size, licensing, and patient-level privacy considerations.

The repository contains clean scripts, documentation, and aggregate summary tables only.

## Repository Structure

```text
scripts/                 Clean R/Python scripts
utils/                   Reusable helper functions
docs/                    Methodology and interpretation notes
results/summary_tables/  Aggregate non-patient-level results
data/README.md           Data access and privacy note
environment/             Package lists
```

## Tech Stack

R, Bioconductor, `survival`, `glmnet`, `sva`, `randomForestSRC`, Python, PyTorch, `pycox`, Bash, SLURM/HPC.
