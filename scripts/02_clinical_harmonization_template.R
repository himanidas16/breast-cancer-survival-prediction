#!/usr/bin/env Rscript
# 02_clinical_harmonization_template.R
# Template for converting dataset-specific clinical files into a standard survival table.
# This is intentionally a template because each cohort uses different metadata column names.

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit paths first.")
}
source("config.R")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

standardize_clinical <- function(df,
                                 dataset_name,
                                 sample_id_col,
                                 os_time_col,
                                 os_event_col,
                                 merged_colname_fun,
                                 time_unit = c("months", "days"),
                                 extra_cols = character()) {
  time_unit <- match.arg(time_unit)

  os_time <- as.numeric(as.character(df[[os_time_col]]))
  if (time_unit == "days") os_time <- os_time / 30.44

  out <- data.frame(
    sample_id = as.character(df[[sample_id_col]]),
    dataset = dataset_name,
    os_time = os_time,
    os_event = as.numeric(as.character(df[[os_event_col]])),
    stringsAsFactors = FALSE
  )

  out$merged_colname <- merged_colname_fun(out$sample_id, df)

  for (col in extra_cols) {
    if (col %in% colnames(df)) out[[col]] <- df[[col]]
  }

  out <- out[!is.na(out$os_time) & !is.na(out$os_event), , drop = FALSE]
  out <- out[out$os_time > 0, , drop = FALSE]
  out <- out[!duplicated(out$merged_colname), , drop = FALSE]
  out
}

# -----------------------------
# Example: GSE7390-style metadata
# -----------------------------
# pheno <- read.csv("data/raw/GSE7390_pheno.csv", check.names = FALSE)
# gse7390 <- standardize_clinical(
#   df = pheno,
#   dataset_name = "GSE7390",
#   sample_id_col = "geo_accession",
#   os_time_col = "t.os:ch1",
#   os_event_col = "e.os:ch1",
#   time_unit = "days",
#   merged_colname_fun = function(sample_id, df) paste0("GSE7390.", sample_id, ".cel.gz")
# )

# -----------------------------
# Combine already-cleaned cohort tables
# -----------------------------
# clean_files <- list.files("data/processed/clinical_clean", pattern = "\\.csv$", full.names = TRUE)
# master <- do.call(rbind, lapply(clean_files, read.csv, check.names = FALSE))
# stopifnot(!any(is.na(master$os_time)))
# stopifnot(!any(is.na(master$os_event)))
# stopifnot(!any(duplicated(master$merged_colname)))
# saveRDS(master, file.path(DATA_DIR, "processed", "master_clinical_LOCAL_ONLY.rds"))
# write.csv(aggregate(os_event ~ dataset, master, function(x) c(n = length(x), events = sum(x))),
#           file.path(TABLE_DIR, "clinical_event_summary.csv"), row.names = FALSE)

message("This is a template. Edit the dataset-specific column names before running.")
