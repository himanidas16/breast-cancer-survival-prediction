#!/usr/bin/env Rscript
# 04_elastic_net_cox.R
# Clean Elastic Net Cox training/evaluation script.
# Saves aggregate summaries by default. Patient-level risk scores are optional and off by default.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit paths first.")
}
source("config.R")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

cindex_from_risk <- function(time, event, risk) {
  as.numeric(survConcordance(Surv(time, event) ~ risk)$concordance)
}

make_stratified_split <- function(cohort, event, train_fraction = 0.70, seed = 2026) {
  set.seed(seed)
  strata <- interaction(cohort, event, drop = TRUE)
  train_idx <- unlist(lapply(split(seq_along(event), strata), function(idx) {
    n_train <- max(1, floor(length(idx) * train_fraction))
    sample(idx, n_train)
  }), use.names = FALSE)
  train_idx <- sort(unique(train_idx))
  list(train = train_idx, test = setdiff(seq_along(event), train_idx))
}

expr <- readRDS(EXPR_MAD_TOP5000_RDS)  # genes x patients
clin <- readRDS(MASTER_CLINICAL_RDS)

stopifnot(identical(colnames(expr), clin$merged_colname))
X <- t(expr)  # patients x genes
time <- clin$os_time
event <- clin$os_event
cohort <- clin$dataset

split <- make_stratified_split(cohort, event, TRAIN_FRACTION, seed = 2026)
train_idx <- split$train
test_idx <- split$test

# Standardize using training data only.
train_mean <- colMeans(X[train_idx, , drop = FALSE])
train_sd <- apply(X[train_idx, , drop = FALSE], 2, sd)
train_sd[train_sd == 0 | is.na(train_sd)] <- 1
X_scaled <- sweep(sweep(X, 2, train_mean, "-"), 2, train_sd, "/")

Y_train <- Surv(time[train_idx], event[train_idx])

alpha_results <- lapply(ALPHA_GRID, function(a) {
  set.seed(2026)
  fit <- cv.glmnet(
    x = X_scaled[train_idx, , drop = FALSE],
    y = Y_train,
    family = "cox",
    alpha = a,
    nfolds = N_FOLDS,
    type.measure = "deviance"
  )
  coefs <- as.matrix(coef(fit, s = "lambda.min"))
  n_genes <- sum(coefs != 0)
  data.frame(alpha = a, lambda_min = fit$lambda.min, cvm_min = min(fit$cvm), genes_at_lambda_min = n_genes)
})
alpha_summary <- do.call(rbind, alpha_results)
write.csv(alpha_summary, file.path(TABLE_DIR, "elastic_net_alpha_tuning.csv"), row.names = FALSE)

best_alpha <- alpha_summary$alpha[which.min(alpha_summary$cvm_min)]
message("Best alpha: ", best_alpha)

set.seed(2026)
best_fit <- cv.glmnet(
  x = X_scaled[train_idx, , drop = FALSE],
  y = Y_train,
  family = "cox",
  alpha = best_alpha,
  nfolds = N_FOLDS,
  type.measure = "deviance"
)

risk_train <- as.numeric(predict(best_fit, newx = X_scaled[train_idx, , drop = FALSE], s = "lambda.min", type = "link"))
risk_test  <- as.numeric(predict(best_fit, newx = X_scaled[test_idx, , drop = FALSE], s = "lambda.min", type = "link"))

train_c <- cindex_from_risk(time[train_idx], event[train_idx], risk_train)
test_c  <- cindex_from_risk(time[test_idx], event[test_idx], risk_test)

coef_train <- as.matrix(coef(best_fit, s = "lambda.min"))
selected_genes <- rownames(coef_train)[coef_train[, 1] != 0]

summary <- data.frame(
  item = c("patients", "events", "train_patients", "test_patients", "train_events", "test_events", "best_alpha", "lambda_min", "selected_genes", "train_cindex", "test_cindex"),
  value = c(nrow(X), sum(event), length(train_idx), length(test_idx), sum(event[train_idx]), sum(event[test_idx]), best_alpha, best_fit$lambda.min, length(selected_genes), train_c, test_c)
)
write.csv(summary, file.path(TABLE_DIR, "elastic_net_summary_generated.csv"), row.names = FALSE)

if (isTRUE(SAVE_PATIENT_LEVEL_OUTPUTS)) {
  risk_df <- data.frame(
    sample_id = clin$merged_colname,
    set = ifelse(seq_len(nrow(clin)) %in% train_idx, "train", "test"),
    cohort = cohort,
    os_time = time,
    os_event = event,
    risk_score = NA_real_,
    stringsAsFactors = FALSE
  )
  risk_df$risk_score[train_idx] <- risk_train
  risk_df$risk_score[test_idx] <- risk_test
  write.csv(risk_df, file.path(TABLE_DIR, "elastic_net_patient_risk_scores_LOCAL_ONLY.csv"), row.names = FALSE)
}

# Save model locally only; .rds ignored by .gitignore.
saveRDS(best_fit, file.path(DATA_DIR, "processed", "elasticnet_cvfit_LOCAL_ONLY.rds"))
cat("Elastic Net Cox complete. Aggregate summary saved.\n")
