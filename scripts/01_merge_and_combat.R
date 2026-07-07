#!/usr/bin/env Rscript
# 01_merge_and_combat.R
# Clean public version of the multi-cohort merge + missing-value + ComBat pipeline.
# This script assumes each cohort is already a gene x sample expression matrix.

suppressPackageStartupMessages({
  library(sva)
  library(matrixStats)
})

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit the paths first.")
}
source("config.R")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# User-edited input list
# -----------------------------
# Edit this named vector in your local copy. Names become cohort labels.
DATASET_FILES <- c(
  CAL       = file.path(DATA_DIR, "processed", "CAL_gene_mapped.csv"),
  GSE7390   = file.path(DATA_DIR, "processed", "GSE7390_gene_mapped.csv"),
  GSE20711  = file.path(DATA_DIR, "processed", "GSE20711_gene_mapped.csv"),
  GSE20685  = file.path(DATA_DIR, "processed", "GSE20685_gene_mapped.csv"),
  METABRIC  = file.path(DATA_DIR, "processed", "METABRIC_gene_expression.rds"),
  TCGA      = file.path(DATA_DIR, "processed", "TCGA_gene_expression.csv"),
  GSE1456   = file.path(DATA_DIR, "processed", "GSE1456_combined_gene_mapped.csv"),
  GSE42568  = file.path(DATA_DIR, "processed", "GSE42568_gene_mapped.csv"),
  GSE58812  = file.path(DATA_DIR, "processed", "GSE58812_gene_mapped.csv"),
  GSE88770  = file.path(DATA_DIR, "processed", "GSE88770_gene_mapped.csv"),
  GSE162228 = file.path(DATA_DIR, "processed", "GSE162228_gene_mapped.csv")
)

read_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    mat <- readRDS(path)
  } else {
    df <- read.csv(path, row.names = 1, check.names = FALSE)
    mat <- as.matrix(df)
  }
  storage.mode(mat) <- "numeric"
  mat
}

clean_gene_names <- function(mat) {
  genes <- trimws(toupper(rownames(mat)))
  keep <- !is.na(genes) & genes != ""
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  rownames(mat) <- genes
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  mat
}

message("Loading and cleaning cohort matrices...")
matrices <- lapply(names(DATASET_FILES), function(ds) {
  mat <- read_matrix(DATASET_FILES[[ds]])
  mat <- clean_gene_names(mat)
  colnames(mat) <- paste(ds, colnames(mat), sep = ".")
  message(ds, ": ", nrow(mat), " genes x ", ncol(mat), " samples")
  mat
})
names(matrices) <- names(DATASET_FILES)

message("Finding common genes...")
common_genes <- Reduce(intersect, lapply(matrices, rownames))
message("Common genes: ", length(common_genes))

matrices <- lapply(matrices, function(mat) mat[common_genes, , drop = FALSE])
stopifnot(all(vapply(matrices, function(x) identical(rownames(x), common_genes), logical(1))))

message("Merging matrices...")
merged <- do.call(cbind, matrices)
batch <- rep(names(matrices), times = vapply(matrices, ncol, integer(1)))
stopifnot(length(batch) == ncol(merged))

sample_manifest <- data.frame(
  merged_colname = colnames(merged),
  batch = batch,
  stringsAsFactors = FALSE
)
write.csv(sample_manifest, file.path(TABLE_DIR, "sample_manifest_TEMPLATE_DO_NOT_UPLOAD_IF_PATIENT_LEVEL.csv"), row.names = FALSE)

# -----------------------------
# Missing-value thresholding + median imputation
# -----------------------------
missing_fraction <- rowMeans(is.na(merged))
MISSINGNESS_THRESHOLD <- 0.20
keep_genes <- missing_fraction <= MISSINGNESS_THRESHOLD
filtered <- merged[keep_genes, , drop = FALSE]

message("Total NAs before imputation: ", sum(is.na(filtered)))
message("Genes removed by missingness threshold: ", sum(!keep_genes))

impute_row_median <- function(x) {
  if (anyNA(x)) x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}
imputed <- t(apply(filtered, 1, impute_row_median))
stopifnot(sum(is.na(imputed)) == 0)

# -----------------------------
# PCA helper
# -----------------------------
plot_pca <- function(mat, batch, out_png, title) {
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 2)
  png(out_png, width = 1200, height = 900, res = 150)
  plot(pca$x[, 1], pca$x[, 2],
       pch = 16,
       xlab = paste0("PC1 (", pct[1], "%)"),
       ylab = paste0("PC2 (", pct[2], "%)"),
       main = title,
       col = as.integer(as.factor(batch)))
  legend("topright", legend = levels(as.factor(batch)), col = seq_along(levels(as.factor(batch))), pch = 16, cex = 0.7)
  dev.off()
  invisible(pct)
}

pre_pct <- plot_pca(imputed, batch, file.path(FIG_DIR, "pca_before_combat.png"), "PCA before ComBat")

# -----------------------------
# ComBat
# -----------------------------
message("Running ComBat...")
combat <- ComBat(dat = imputed, batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
stopifnot(sum(is.na(combat)) == 0)
post_pct <- plot_pca(combat, batch, file.path(FIG_DIR, "pca_after_combat.png"), "PCA after ComBat")

summary_df <- data.frame(
  metric = c("common_genes", "samples", "missingness_threshold", "genes_removed_by_threshold", "nas_after_imputation", "nas_after_combat", "pca_before_pc1", "pca_before_pc2", "pca_after_pc1", "pca_after_pc2"),
  value = c(length(common_genes), ncol(merged), MISSINGNESS_THRESHOLD, sum(!keep_genes), sum(is.na(imputed)), sum(is.na(combat)), pre_pct[1], pre_pct[2], post_pct[1], post_pct[2])
)
write.csv(summary_df, file.path(TABLE_DIR, "combat_summary.csv"), row.names = FALSE)

# Keep this output local/private; do not commit .rds files.
saveRDS(combat, file.path(DATA_DIR, "processed", "merged_expression_ComBat_corrected_LOCAL_ONLY.rds"))
message("Done. ComBat matrix saved locally. Do not upload .rds files to GitHub.")
