#!/usr/bin/env Rscript
# 09_elastic_net_loco_validation.R
# Leave-one-cohort-out validation for Elastic Net Cox.

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

expr <- readRDS(EXPR_MAD_TOP5000_RDS)
clin <- readRDS(MASTER_CLINICAL_RDS)
stopifnot(identical(colnames(expr), clin$merged_colname))

X <- t(expr)
time <- clin$os_time
event <- clin$os_event
cohort <- clin$dataset

results <- lapply(sort(unique(cohort)), function(held_out) {
  test_idx <- which(cohort == held_out)
  train_idx <- which(cohort != held_out)

  # Scale with training cohorts only.
  train_mean <- colMeans(X[train_idx, , drop = FALSE])
  train_sd <- apply(X[train_idx, , drop = FALSE], 2, sd)
  train_sd[train_sd == 0 | is.na(train_sd)] <- 1
  X_scaled <- sweep(sweep(X, 2, train_mean, "-"), 2, train_sd, "/")

  set.seed(2026)
  fit <- cv.glmnet(
    x = X_scaled[train_idx, , drop = FALSE],
    y = Surv(time[train_idx], event[train_idx]),
    family = "cox",
    alpha = 0.1,
    nfolds = N_FOLDS,
    type.measure = "deviance"
  )

  risk <- as.numeric(predict(fit, newx = X_scaled[test_idx, , drop = FALSE], s = "lambda.min", type = "link"))
  cind <- cindex_from_risk(time[test_idx], event[test_idx], risk)
  coefs <- as.matrix(coef(fit, s = "lambda.min"))

  data.frame(
    held_out_cohort = held_out,
    test_patients = length(test_idx),
    events = sum(event[test_idx]),
    lambda_min = fit$lambda.min,
    genes_at_lambda_min = sum(coefs != 0),
    cindex = cind,
    stringsAsFactors = FALSE
  )
})

loco <- do.call(rbind, results)
weighted <- weighted.mean(loco$cindex, w = loco$events)
unweighted <- mean(loco$cindex)

write.csv(loco, file.path(TABLE_DIR, "elastic_net_loco_validation_generated.csv"), row.names = FALSE)
write.csv(data.frame(metric = c("mean_unweighted", "mean_event_weighted"), value = c(unweighted, weighted)),
          file.path(TABLE_DIR, "elastic_net_loco_summary_generated.csv"), row.names = FALSE)
cat("Elastic Net LOCO validation complete.\n")
