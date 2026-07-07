#!/usr/bin/env Rscript
# 06_random_survival_forest.R
# RSF benchmark with in-fold univariate Cox filtering.

suppressPackageStartupMessages({
  library(survival)
  library(randomForestSRC)
})

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit paths first.")
}
source("config.R")
source("scripts/05_univariate_cox_filter.R")
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

selected <- univariate_cox_filter(
  X_train = X[train_idx, , drop = FALSE],
  time_train = time[train_idx],
  event_train = event[train_idx],
  p_threshold = UNIVAR_P,
  max_genes = UNIVAR_MAXN
)
message("RSF selected genes: ", length(selected))

train_df <- data.frame(os_time = time[train_idx], os_event = event[train_idx], X[train_idx, selected, drop = FALSE], check.names = FALSE)
test_df <- data.frame(X[test_idx, selected, drop = FALSE], check.names = FALSE)

set.seed(2026)
fit <- rfsrc(
  Surv(os_time, os_event) ~ .,
  data = train_df,
  ntree = RSF_NTREE,
  nodesize = RSF_NODESIZE,
  nsplit = RSF_NSPLIT,
  importance = FALSE,
  seed = -2026
)

pred <- predict(fit, newdata = test_df)
# For randomForestSRC survival, higher predicted mortality usually means higher risk.
risk_test <- pred$predicted
if (is.matrix(risk_test)) risk_test <- rowSums(risk_test)

test_c <- cindex_from_risk(time[test_idx], event[test_idx], risk_test)
summary <- data.frame(
  model = "Random Survival Forest",
  test_patients = length(test_idx),
  test_events = sum(event[test_idx]),
  selected_genes = length(selected),
  ntree = RSF_NTREE,
  nodesize = RSF_NODESIZE,
  nsplit = RSF_NSPLIT,
  test_cindex = test_c
)
write.csv(summary, file.path(TABLE_DIR, "rsf_summary_generated.csv"), row.names = FALSE)

saveRDS(fit, file.path(DATA_DIR, "processed", "rsf_model_LOCAL_ONLY.rds"))
cat("RSF complete. Aggregate summary saved.\n")
