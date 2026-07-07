# Data

Raw and processed data are not included in this repository.

This project uses public breast cancer gene-expression and clinical datasets from GEO, ArrayExpress/BioStudies, cBioPortal, UCSC Xena, METABRIC, TCGA, and SCAN-B/GSE96058.

Large expression matrices, clinical metadata files, patient-level risk scores, and trained models are excluded because of file size, licensing, and patient-level privacy considerations.

To reproduce the pipeline, download the datasets from their original sources and update the paths in `scripts/00_config_template.R`.
