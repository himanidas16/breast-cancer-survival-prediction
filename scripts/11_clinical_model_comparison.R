#!/usr/bin/env Rscript
# 11_clinical_model_comparison.R
# Compare clinical-only, signature-only, and clinical + signature Cox models.
# Intended for external validation cohorts where clinical variables are available.

suppressPackageStartupMessages(library(survival))

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit paths first.")
}
source("config.R")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Local/private input expected with one row per patient.
# Required columns: os_time, os_event, risk_score
# Optional clinical columns: age, er_status, pr_status, her2_status, grade, stage
INPUT_CSV <- file.path(DATA_DIR, "external", "external_clinical_with_signature_LOCAL_ONLY.csv")
if (!file.exists(INPUT_CSV)) stop("Missing local/private input: ", INPUT_CSV)

df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
required <- c("os_time", "os_event", "risk_score")
stopifnot(all(required %in% colnames(df)))

candidate_clinical <- intersect(c("age", "er_status", "pr_status", "her2_status", "grade", "stage"), colnames(df))
complete_cols <- c(required, candidate_clinical)
df <- df[complete.cases(df[, complete_cols, drop = FALSE]), ]

surv_obj <- Surv(df$os_time, df$os_event)
clinical_formula <- if (length(candidate_clinical) > 0) {
  as.formula(paste("surv_obj ~", paste(candidate_clinical, collapse = " + ")))
} else {
  NULL
}

fit_sig <- coxph(surv_obj ~ risk_score, data = df)

rows <- list(
  data.frame(model = "signature_only", n = nrow(df), cindex = summary(fit_sig)$concordance[1])
)

if (!is.null(clinical_formula)) {
  fit_clin <- coxph(clinical_formula, data = df)
  fit_combined <- coxph(as.formula(paste("surv_obj ~ risk_score +", paste(candidate_clinical, collapse = " + "))), data = df)
  lrt <- anova(fit_clin, fit_combined, test = "LRT")
  rows[[2]] <- data.frame(model = "clinical_only", n = nrow(df), cindex = summary(fit_clin)$concordance[1])
  rows[[3]] <- data.frame(model = "clinical_plus_signature", n = nrow(df), cindex = summary(fit_combined)$concordance[1])
  write.csv(lrt, file.path(TABLE_DIR, "clinical_signature_lrt_generated.csv"))
}

out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, "clinical_model_comparison_generated.csv"), row.names = FALSE)
cat("Clinical model comparison complete.\n")
