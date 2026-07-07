#!/usr/bin/env Rscript
# 10_scanb_external_validation.R
# External validation template for applying a frozen Elastic Net Cox signature to SCAN-B.
# Keep SCAN-B expression/clinical files and patient-level risk scores private.

suppressPackageStartupMessages({
  library(survival)
})

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit paths first.")
}
source("config.R")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

cindex_from_risk <- function(time, event, risk) {
  as.numeric(survConcordance(Surv(time, event) ~ risk)$concordance)
}

# Expected coefficient file format:
# gene,coef
coef_df <- read.csv(FROZEN_COEF_CSV, stringsAsFactors = FALSE)
coef_df <- coef_df[coef_df$coef != 0, ]

# Optional discovery gene SD file format:
# gene,sd
sd_df <- read.csv(DISCOVERY_GENE_SD_CSV, stringsAsFactors = FALSE)

scanb_expr <- readRDS(SCANB_EXPR_RDS)  # genes x patients
scanb_clin <- read.csv(SCANB_CLINICAL_CSV, stringsAsFactors = FALSE)

# Required clinical columns: sample_id, os_time, os_event
stopifnot(all(c("sample_id", "os_time", "os_event") %in% colnames(scanb_clin)))
stopifnot(all(scanb_clin$sample_id %in% colnames(scanb_expr)))
scanb_expr <- scanb_expr[, scanb_clin$sample_id, drop = FALSE]

available <- intersect(coef_df$gene, rownames(scanb_expr))
missing <- setdiff(coef_df$gene, rownames(scanb_expr))
message("Frozen model genes: ", nrow(coef_df))
message("Available in SCAN-B: ", length(available))
message("Missing in SCAN-B: ", length(missing))

coef_use <- coef_df[match(available, coef_df$gene), ]
sd_use <- sd_df[match(available, sd_df$gene), ]
if (any(is.na(sd_use$sd))) stop("Missing discovery SD for some available genes.")

# Gene-wise z-score within SCAN-B.
expr_use <- scanb_expr[available, , drop = FALSE]
scanb_mean <- rowMeans(expr_use, na.rm = TRUE)
scanb_sd <- apply(expr_use, 1, sd, na.rm = TRUE)
scanb_sd[scanb_sd == 0 | is.na(scanb_sd)] <- 1
expr_z <- sweep(sweep(expr_use, 1, scanb_mean, "-"), 1, scanb_sd, "/")

# Convert coefficient to per-discovery-SD scale.
beta_per_sd <- coef_use$coef * sd_use$sd
risk <- as.numeric(crossprod(beta_per_sd, expr_z))

cind <- cindex_from_risk(scanb_clin$os_time, scanb_clin$os_event, risk)

# Equal tertile risk groups.
q <- quantile(risk, probs = c(1/3, 2/3), na.rm = TRUE)
risk_group <- cut(risk, breaks = c(-Inf, q[1], q[2], Inf), labels = c("Low", "Intermediate", "High"), include.lowest = TRUE)

group_summary <- aggregate(scanb_clin$os_event, by = list(risk_group = risk_group), FUN = function(x) c(patients = length(x), events = sum(x), event_rate = mean(x)))
group_summary <- do.call(data.frame, group_summary)
colnames(group_summary) <- c("risk_group", "patients", "events", "event_rate")

summary <- data.frame(
  item = c("external_dataset", "patients", "events", "frozen_model_genes", "available_genes", "missing_genes", "external_cindex"),
  value = c("SCAN-B / GSE96058", nrow(scanb_clin), sum(scanb_clin$os_event), nrow(coef_df), length(available), length(missing), cind)
)
write.csv(summary, file.path(TABLE_DIR, "scanb_validation_summary_generated.csv"), row.names = FALSE)
write.csv(group_summary, file.path(TABLE_DIR, "scanb_risk_group_summary_generated.csv"), row.names = FALSE)

if (isTRUE(SAVE_PATIENT_LEVEL_OUTPUTS)) {
  risk_df <- data.frame(sample_id = scanb_clin$sample_id, os_time = scanb_clin$os_time, os_event = scanb_clin$os_event, risk_score = risk, risk_group = risk_group)
  write.csv(risk_df, file.path(TABLE_DIR, "scanb_patient_risk_scores_LOCAL_ONLY.csv"), row.names = FALSE)
}

cat("SCAN-B external validation complete. Aggregate summaries saved.\n")
