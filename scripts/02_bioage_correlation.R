# ============================================================
# scripts/02_bioage_correlation.R
# Correlation plot (BioAge::plot_baa) for TwinGene + UKB
# ============================================================

make_one_ba_corr <- function(df, cohort_name, outdir) {

  quiet_require("BioAge")

  agevar <- c("PRS_Z", "age", "phenoage_res", "kdm_res", "hd_res",
              "Frailty_Index_0.10_log", "Telomere")

  miss <- setdiff(agevar, names(df))
  if (length(miss) > 0) {
    stop("Missing variables for correlation plot in ", cohort_name, ": ",
         paste(miss, collapse = ", "))
  }

  axis_type <- setNames(rep("float", length(agevar)), agevar)

  label <- c(
    "PRS_Z"                 = "PRS",
    "age"                   = "Chronological age",
    "phenoage_res"          = "PhenoAge residual",
    "kdm_res"               = "KDM residual",
    "hd_res"                = "HD residual",
    "Frailty_Index_0.10_log"= "Frailty index (log)",
    "Telomere"              = "Telomere length"
  )

  png(file.path(outdir, paste0("Figure2_correlation_", gsub(" ", "_", cohort_name), ".png")),
      width = 2200, height = 1800, res = 300)
  BioAge::plot_baa(df, agevar, label, axis_type)
  dev.off()

  invisible(TRUE)
}

make_ba_correlation_plots <- function(tg, ukb, outdir) {
  make_one_ba_corr(tg,  "TwinGene",   outdir)
  make_one_ba_corr(ukb, "UK Biobank", outdir)
  invisible(TRUE)
}
