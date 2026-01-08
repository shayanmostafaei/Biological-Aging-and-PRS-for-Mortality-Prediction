# ============================================================
# scripts/00_setup.R
# Install/load packages + shared helpers
# ============================================================

quiet_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg, "\nInstall it and re-run.")
  }
}

pkgs <- c(
  "data.table","dplyr","tidyr","readr","stringr","purrr","tibble",
  "ggplot2","ggnewscale","scales","patchwork","BioAge",
  "survival","survminer","rms","timeROC","pROC",
  "SuperLearner","caret","xgboost","randomForest",
  "Matrix","glmnet"
)

invisible(lapply(pkgs, quiet_require))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

# -----------------------------
# Shared colors (match manuscript style)
# -----------------------------
# Flowchart/boxes (requested)
COL_TG      <- "#E6550D"  # TwinGene orange
COL_UKB     <- "#1F78B4"  # UK Biobank blue

# Heatmaps/density (used across your supplementary figures)
COL_TG_HEAT <- "#D55E00"  # TwinGene orange used in heatmaps/densities
COL_UKB_HEAT<- "#0072B2"  # UKB blue used in heatmaps

# -----------------------------
# Shared helpers
# -----------------------------
zscore <- function(x) as.numeric(scale(x))

age_residual <- function(x, age) {
  stats::residuals(stats::lm(x ~ age))
}

# Simple save helper
save_rds <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(x, path)
  invisible(path)
}
