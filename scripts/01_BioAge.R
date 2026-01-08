# ============================================================
# Correlation plots (BioAge::plot_baa) for TwinGene + UKB
# - Uses the SAME variable set and labeling across cohorts
# ============================================================

library(BioAge)

# ---------------------------
# 1) Define the target variables per cohort
# ---------------------------

# TwinGene columns 
vars_tg <- c(
  "PRS_Z",
  "age",
  "Phenoage_res",
  "KDM_res",
  "HD_res",
  "FI",
  "Telomere_Length"
)

# UKB columns 
vars_ukb <- c(
  "PRS_Z",
  "age",
  "Phenoage_res",
  "KDM_res",
  "HD_res",
  "FI",
  "Telomere_Length"
)

# ---------------------------
# 2) Axis types (BioAge uses "float" here)
# ---------------------------
axis_type_tg  <- setNames(rep("float", length(vars_tg)),  vars_tg)
axis_type_ukb <- setNames(rep("float", length(vars_ukb)), vars_ukb)

# ---------------------------
# 3) Labels (clean manuscript-consistent)
# ---------------------------
label_tg <- c(
  "PRS_Z"                 = "PRS",
  "age"                   = "CA",
  "Phenoage_res"          = "PhenoAge\n(residual)",
  "KDM_res"               = "KDM\n(residual)",
  "HD_res"                = "HD\n(residual)",
  "FI"                    = "FI\n(log-transformed)",
  "Telomere_Length"       = "TL"
)

label_ukb <- c(
  "PRS_Z"                   = "PRS",
  "age"                     = "CA",
  "Phenoage_res"       = "PhenoAge\n(residual)",
  "KDM_res"            = "KDM\n(residual)",
  "HD_res"             = "HD\n(residual)",
  "FI"                 = "FI\n(residual)",
  "Telomere_Length"    = "TL\n(residual)"
)

# ---------------------------
# 4) Small helper to check columns and plot
# ---------------------------
plot_baa_safe <- function(df, vars, labels, axis_type, cohort_name = "Cohort") {
  missing <- setdiff(vars, names(df))
  if (length(missing) > 0) {
    stop(
      cohort_name, ": missing required columns for plot_baa(): ",
      paste(missing, collapse = ", ")
    )
  }
  BioAge::plot_baa(df, vars, labels, axis_type)
}

# ---------------------------
# 5) Run plots
# ---------------------------
# TwinGene
p_tg <- plot_baa_safe(
  df = data_twingene,
  vars = vars_tg,
  labels = label_tg,
  axis_type = axis_type_tg,
  cohort_name = "TwinGene"
)

# UK Biobank
p_ukb <- plot_baa_safe(
  df = data_ukb,
  vars = vars_ukb,
  labels = label_ukb,
  axis_type = axis_type_ukb,
  cohort_name = "UK Biobank"
)

# Print to device (interactive)
print(p_tg)
print(p_ukb)

# ---------------------------
# 6) Optional: save (recommended for repo reproducibility)
# ---------------------------
# ggsave("figures/Figure2A_TwinGene_correlation.png", p_tg, width=8, height=7, dpi=300)
# ggsave("figures/Figure2B_UKB_correlation.png",     p_ukb, width=8, height=7, dpi=300)
