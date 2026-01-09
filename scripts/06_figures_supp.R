# ============================================================
# scripts/06_figures_supp.R
# Supplementary Figure X — TwinGene biomarker density plots
# ============================================================

make_supp_figure1_density <- function(tg, outdir) {

  tg_fill <- COL_TG_HEAT  

  # Detect CRP if present 
  crp_candidates <- names(tg)[stringr::str_detect(
    names(tg),
    regex("^crp$|crp|c_reactive|c-reactive", ignore_case = TRUE)
  )]
  crp_col <- if (length(crp_candidates) > 0) crp_candidates[1] else NA_character_

  # Reverse log if needed 
  maybe_exp <- function(x) {
    q95 <- suppressWarnings(stats::quantile(x, 0.95, na.rm = TRUE))
    if (is.finite(q95) && q95 < 20) exp(x) else x
  }

  # Convert mmol/L to mg/dL if needed (only if values look like mmol/L)
  mmol_to_mgdl_chol <- function(x) x * 38.67
  mmol_to_mgdl_trig <- function(x) x * 88.57

  maybe_convert_chol <- function(x) {
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (is.finite(med) && med < 30) mmol_to_mgdl_chol(x) else x
  }
  maybe_convert_trig <- function(x) {
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (is.finite(med) && med < 30) mmol_to_mgdl_trig(x) else x
  }

  # Build biomarker df 
  stop_if_missing(tg, c("glucose_mmol","hdl","ldl","totchol","trig","lnhba1c","lncreat_umol","cyst"),
                  "TwinGene (for Supplementary Figure 1)")

  tg_df <- tibble::tibble(
    glucose_mmol  = maybe_exp(tg$glucose_mmol),
    hdl_mgdl      = maybe_convert_chol(tg$hdl),
    ldl_mgdl      = maybe_convert_chol(tg$ldl),
    totchol_mgdl  = maybe_convert_chol(tg$totchol),
    trig_mgdl     = maybe_convert_trig(maybe_exp(tg$trig)),
    hba1c_pct     = maybe_exp(tg$lnhba1c),
    creat_umol    = maybe_exp(tg$lncreat_umol),
    cyst_mgL      = tg$cyst
  )

  if (!is.na(crp_col)) {
    tg_df <- tg_df %>% mutate(crp_mgdl = maybe_exp(tg[[crp_col]]))
  }

  long <- tg_df %>%
    tidyr::pivot_longer(everything(), names_to = "Biomarker", values_to = "Value") %>%
    dplyr::filter(is.finite(Value)) %>%
    dplyr::mutate(
      BiomarkerLabel = dplyr::recode(
        Biomarker,
        glucose_mmol = "Serum glucose (mmol/L)",
        hdl_mgdl     = "HDL cholesterol (mg/dL)",
        ldl_mgdl     = "LDL cholesterol (mg/dL)",
        totchol_mgdl = "Total cholesterol (mg/dL)",
        trig_mgdl    = "Triglyceride (mg/dL)",
        hba1c_pct    = "Glycohemoglobin (HbA1c, %)",
        creat_umol   = "Creatinine (µmol/L)",
        cyst_mgL     = "Cystatin C (mg/L)",
        crp_mgdl     = "C-reactive protein (mg/dL)"
      )
    )

  p <- ggplot(long, aes(x = Value)) +
    geom_density(fill = tg_fill, color = "black", linewidth = 0.45, alpha = 0.25) +
    facet_wrap(~ BiomarkerLabel, scales = "free", ncol = 4) +
    labs(
      title = "Supplementary Figure 1. Density plots of biomarkers in TwinGene",
      x = NULL,
      y = "Density"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "plain"),
      plot.title = element_text(face = "bold"),
      panel.spacing = grid::unit(1.0, "lines")
    )

  ggplot2::ggsave(
    filename = file.path(outdir, "Supplementary_Figure1_TwinGene_density.png"),
    plot = p, width = 10, height = 10, units = "in", dpi = 300
  )

  p
}
