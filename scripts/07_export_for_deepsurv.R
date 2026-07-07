#!/usr/bin/env Rscript
# 07_export_for_deepsurv.R
# Export R objects to CSV for Python DeepSurv training.
# CSVs generated here are patient-level; keep them local/private and do not upload.

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit paths first.")
}
source("config.R")

EXPORT_DIR <- file.path(DATA_DIR, "deepsurv_export_LOCAL_ONLY")
dir.create(EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)

expr <- readRDS(EXPR_MAD_TOP5000_RDS)  # genes x patients
clin <- readRDS(MASTER_CLINICAL_RDS)
stopifnot(identical(colnames(expr), clin$merged_colname))

# Patients x genes, with sample_id as first column.
expr_py <- data.frame(sample_id = colnames(expr), t(expr), check.names = FALSE)
clinical_py <- data.frame(
  sample_id = clin$merged_colname,
  cohort = clin$dataset,
  os_time = clin$os_time,
  os_event = clin$os_event,
  stringsAsFactors = FALSE
)

# Use existing split if provided; otherwise create a simple split in R.
if (file.exists(TRAIN_TEST_SPLIT_CSV)) {
  split_df <- read.csv(TRAIN_TEST_SPLIT_CSV, stringsAsFactors = FALSE)
  stopifnot(all(clin$merged_colname %in% split_df$sample_id))
  clinical_py$set <- split_df$set[match(clin$merged_colname, split_df$sample_id)]
} else {
  set.seed(2026)
  clinical_py$set <- "test"
  strata <- interaction(clin$dataset, clin$os_event, drop = TRUE)
  train_idx <- unlist(lapply(split(seq_len(nrow(clin)), strata), function(idx) {
    sample(idx, max(1, floor(length(idx) * TRAIN_FRACTION)))
  }), use.names = FALSE)
  clinical_py$set[train_idx] <- "train"
}

write.csv(expr_py, file.path(EXPORT_DIR, "deepsurv_expr_LOCAL_ONLY.csv"), row.names = FALSE)
write.csv(clinical_py, file.path(EXPORT_DIR, "deepsurv_clinical_LOCAL_ONLY.csv"), row.names = FALSE)
cat("Exported DeepSurv CSVs to ", EXPORT_DIR, "\n", sep = "")
cat("Do not upload these patient-level CSVs to GitHub.\n")
