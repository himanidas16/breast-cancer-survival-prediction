# 05_univariate_cox_filter.R
# Helper function used by RSF and optionally DeepSurv.
# Selects genes inside the training fold only to avoid information leakage.

suppressPackageStartupMessages(library(survival))

univariate_cox_filter <- function(X_train, time_train, event_train, p_threshold = 0.05, max_genes = 500, fallback_genes = 50) {
  # X_train must be patients x genes.
  stopifnot(nrow(X_train) == length(time_train), length(time_train) == length(event_train))

  pvals <- apply(X_train, 2, function(g) {
    if (sd(g, na.rm = TRUE) == 0 || all(is.na(g))) return(NA_real_)
    fit <- tryCatch(
      suppressWarnings(coxph(Surv(time_train, event_train) ~ g)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    tryCatch(summary(fit)$coefficients[1, "Pr(>|z|)"], error = function(e) NA_real_)
  })

  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0) stop("No valid univariate Cox p-values were computed.")

  keep <- names(which(pvals < p_threshold))
  if (length(keep) == 0) {
    keep <- names(sort(pvals))[seq_len(min(fallback_genes, length(pvals)))]
  }
  if (length(keep) > max_genes) {
    keep <- names(sort(pvals[keep]))[seq_len(max_genes)]
  }
  keep
}
