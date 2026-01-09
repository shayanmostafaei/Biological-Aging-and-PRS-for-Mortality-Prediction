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
  "ggplot2","ggnewscale","scales","patchwork",
  "survival","survminer","rms","timeROC","pROC",
  "SuperLearner","caret","xgboost","randomForest",
  "Matrix","glmnet","BioAge"
)

invisible(lapply(pkgs, quiet_require))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# -----------------------------
# Colors used across figures (match manuscript style)
# -----------------------------
COL_TG_BOX   <- "#E6550D"  # TwinGene orange (flowchart boxes)
COL_UKB_BOX  <- "#1F78B4"  # UKB blue (flowchart boxes)

COL_TG_HEAT  <- "#D55E00"  # TwinGene orange (heatmaps/densities)
COL_UKB_HEAT <- "#0072B2"  # UKB blue (heatmaps)

# -----------------------------
# Utilities
# -----------------------------
zscore <- function(x) as.numeric(scale(x))

age_residual <- function(x, age) {
  stats::residuals(stats::lm(x ~ age))
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

stop_if_missing <- function(df, cols, df_name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(df_name, " is missing required columns: ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}
