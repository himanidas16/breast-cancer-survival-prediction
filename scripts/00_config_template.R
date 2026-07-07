# 00_config_template.R
# Copy this file to config.R and edit paths for your machine.
# Do not commit config.R if it contains private paths.

set.seed(2026)

# -----------------------------
# Project directories
# -----------------------------
PROJECT_DIR <- "."
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
TABLE_DIR   <- file.path(RESULTS_DIR, "summary_tables")

# -----------------------------
# Input files
# -----------------------------
# These are placeholders. Keep actual files outside GitHub.
EXPR_SURVIVAL_RDS <- file.path(DATA_DIR, "processed", "expr_survival_matrix.rds")
MASTER_CLINICAL_RDS <- file.path(DATA_DIR, "processed", "master_clinical.rds")

# Primary modelling matrix after MAD filtering.
EXPR_MAD_TOP5000_RDS <- file.path(DATA_DIR, "processed", "expr_MAD_top5000.rds")

# Optional train/test split file from Elastic Net Cox.
TRAIN_TEST_SPLIT_CSV <- file.path(DATA_DIR, "processed", "elasticnet_train_test_split.csv")

# External validation placeholders.
SCANB_EXPR_RDS <- file.path(DATA_DIR, "external", "scanb_expression.rds")
SCANB_CLINICAL_CSV <- file.path(DATA_DIR, "external", "scanb_clinical.csv")
FROZEN_COEF_CSV <- file.path(DATA_DIR, "processed", "elasticnet_frozen_coefficients.csv")
DISCOVERY_GENE_SD_CSV <- file.path(DATA_DIR, "processed", "discovery_gene_sd.csv")

# -----------------------------
# Modelling settings
# -----------------------------
ALPHA_GRID <- c(0.1, 0.3, 0.5, 0.7, 0.9)
N_FOLDS <- 10
TRAIN_FRACTION <- 0.70

# RSF settings
RSF_NTREE <- 1000
RSF_NODESIZE <- 15
RSF_NSPLIT <- 10
UNIVAR_P <- 0.05
UNIVAR_MAXN <- 500

# Data safety settings
SAVE_PATIENT_LEVEL_OUTPUTS <- FALSE
