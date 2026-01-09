# ============================================================
# TwinGene_analysis.R
# Main entrypoint for multicohort analyses:
# TwinGene discovery + UK Biobank external validation
# ============================================================

rm(list = ls())

# -----------------------------
# User configuration
# -----------------------------
output_dir <- "results_outputs"   # change as needed
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Provide local paths to your analytic datasets (NOT included in repo)
twingene_path <- "data/TwinGene_analytic.csv"
ukb_path      <- "data/UKB_analytic.csv"

# Seeds for reproducibility
set.seed(12345)

# -----------------------------
# Source modular scripts
# -----------------------------
source("scripts/00_setup.R")
source("scripts/01_prepare_inputs.R")
source("scripts/02_bioage_correlation.R")
source("scripts/03_models_cox_roc.R")
source("scripts/04_superlearner_cv.R")
source("scripts/05_subgroup_interactions.R")
source("scripts/06_figures_supp.R")

# -----------------------------
# Load inputs
# -----------------------------
tg_raw  <- read_input_twingene(twingene_path)
ukb_raw <- read_input_ukb(ukb_path)

# -----------------------------
# Harmonize / derive core fields
# -----------------------------
tg  <- prep_twingene(tg_raw)
ukb <- prep_ukb(ukb_raw)

# -----------------------------
# Run analyses
# -----------------------------
# 1) Correlation plot (Figure 2)
cor_out <- make_ba_correlation_plots(tg, ukb, outdir = output_dir)

# 2) Univariate ROC (secondary) + Cox models (primary)
roc_out <- run_univariate_roc(tg, ukb, outdir = output_dir)
cox_out <- run_cox_models(tg, ukb, outdir = output_dir)

# 3) Multivariable SuperLearner (10-fold CV; repeated CV for robustness)
sl_out <- run_superlearner_models(
  tg, ukb,
  outdir  = output_dir,
  K       = 10,
  repeats = 10,              # set repeats=1 for single CV run
  seeds   = 1001:1010
)

# 4) Subgroup AUCs at fixed horizons + interaction p-value heatmap 
s2_out  <- run_subgroup_followup_age(tg, ukb, outdir = output_dir,
                                     horizons = c(5,10,15),
                                     age_breaks = c(-Inf, 50, 60, Inf))
s3_out  <- run_interaction_heatmaps(tg, ukb, outdir = output_dir)

# 5) Supplementary biomarker density plots 
fig_s1 <- make_supp_figure1_density(tg, outdir = output_dir)

# -----------------------------
# Session info for reproducibility
# -----------------------------
sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("DONE. Outputs saved to: ", normalizePath(output_dir))
