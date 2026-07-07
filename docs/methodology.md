# Methodology

## 1. Expression preprocessing

Raw microarray cohorts were normalized with RMA where raw CEL files were available. The resulting expression matrices were converted to gene-level matrices.

## 2. Probe-to-gene mapping

Microarray platforms measure probes, not genes. Probe identifiers were mapped to gene symbols using platform-specific Bioconductor annotation packages or platform annotation tables.

## 3. Probe collapsing

For genes represented by multiple probes, one gene-level expression value was retained. Affymetrix datasets used a MaxMean-style representative-probe strategy. In earlier METABRIC processing, Illumina probes were mapped using `illuminaHumanv3.db` and collapsed to gene-level expression.

## 4. Multi-cohort merging

Each cohort was converted to a gene x sample matrix. Gene symbols were cleaned and standardized, common genes across all cohorts were identified, and matrices were merged column-wise into a shared expression space.

## 5. Missing-value handling

Sparse missing values were handled with per-gene missingness thresholding followed by gene-wise median imputation. This avoided unnecessarily dropping genes affected by only a few missing entries.

## 6. Batch correction

ComBat was applied using cohort identity as the batch variable. PCA before and after ComBat was used to check whether strong cohort/platform effects were reduced.

## 7. Feature filtering

Gene-wise MAD was calculated across patients. The top 5,000 MAD-ranked genes were used as the primary modelling input, while top 2,500 and top 7,500 gene sets were kept for sensitivity analysis.

## 8. Survival modelling

Overall survival was modelled using Elastic Net Cox, Random Survival Forest, and DeepSurv. Model performance was evaluated using Harrell's C-index, leave-one-cohort-out validation, and external SCAN-B validation.
