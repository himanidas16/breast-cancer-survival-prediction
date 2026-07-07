#!/usr/bin/env Rscript
# 03_mad_gene_filtering.R
# MAD-based gene filtering for a survival-matched ComBat-corrected expression matrix.

suppressPackageStartupMessages(library(matrixStats))

if (!file.exists("config.R")) {
  stop("Copy scripts/00_config_template.R to config.R and edit paths first.")
}
source("config.R")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

expr <- readRDS(EXPR_SURVIVAL_RDS)  # genes x patients
stopifnot(is.matrix(expr) || is.data.frame(expr))
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

cat("Input matrix checks:\n")
cat("Dimensions:", nrow(expr), "genes x", ncol(expr), "patients\n")
cat("NA count:", sum(is.na(expr)), "\n")
cat("Duplicate gene names:", sum(duplicated(rownames(expr))), "\n")
cat("Duplicate sample names:", sum(duplicated(colnames(expr))), "\n")

mad_vals <- rowMads(expr, na.rm = TRUE)
mad_table <- data.frame(gene = rownames(expr), MAD = mad_vals, stringsAsFactors = FALSE)
mad_table <- mad_table[order(mad_table$MAD, decreasing = TRUE), ]
mad_table$MAD_rank <- seq_len(nrow(mad_table))
mad_table$in_top2500 <- mad_table$MAD_rank <= 2500
mad_table$in_top5000 <- mad_table$MAD_rank <= 5000
mad_table$in_top7500 <- mad_table$MAD_rank <= 7500

thresholds <- data.frame(
  filter_set = c("Top 2500 MAD", "Top 5000 MAD", "Top 7500 MAD"),
  genes_retained = c(2500, 5000, 7500),
  mad_threshold = c(mad_table$MAD[2500], mad_table$MAD[5000], mad_table$MAD[7500]),
  use = c("Sensitivity analysis", "Primary modelling set", "Sensitivity analysis")
)
write.csv(thresholds, file.path(TABLE_DIR, "mad_thresholds.csv"), row.names = FALSE)

# For public GitHub, save only gene lists if advisor approves.
# Exact gene lists may be unpublished research outputs, so keep them local by default.
if (exists("SAVE_GENE_LISTS") && isTRUE(SAVE_GENE_LISTS)) {
  write.csv(mad_table, file.path(TABLE_DIR, "gene_MAD_ranked_table_LOCAL_REVIEW_BEFORE_UPLOAD.csv"), row.names = FALSE)
}

png(file.path(FIG_DIR, "mad_distribution.png"), width = 1200, height = 900, res = 150)
hist(mad_vals, breaks = 80, xlab = "Gene-wise MAD", main = "Distribution of Gene-wise MAD")
dev.off()

png(file.path(FIG_DIR, "ranked_mad_curve.png"), width = 1200, height = 900, res = 150)
plot(sort(mad_vals, decreasing = TRUE), type = "l", lwd = 2,
     xlab = "Gene rank by MAD", ylab = "MAD", main = "Ranked MAD curve")
dev.off()

# Save filtered matrices locally only; .rds is ignored by .gitignore.
for (k in c(2500, 5000, 7500)) {
  genes <- mad_table$gene[seq_len(k)]
  mat <- expr[genes, , drop = FALSE]
  saveRDS(mat, file.path(DATA_DIR, "processed", paste0("expr_MAD_top", k, "_LOCAL_ONLY.rds")))
}

cat("MAD filtering complete. Public-safe summary saved to results/summary_tables/mad_thresholds.csv\n")
