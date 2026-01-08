# ============================================================
# scripts/05_figures_supp.R
# Supplementary figure 3
# ============================================================

make_supp_figure1_density <- function(tg, outdir) {

  # Expect TwinGene columns (as you shared):
  # glucose_mmol, hdl, ldl, totchol, trig, lnhba1c, lncreat_umol, cyst
  # CRP is optional and may have different column names.
  #
  # Output: TwinGene-only overlay-style density panels,
  # with TwinGene fill = "#D55E00".

  tg_fill <- "#D55E00"

  # --- Detect CRP if present
  crp_candidates <- names(tg)[stringr::str_detect(
    names(tg),
    regex("^crp$|crp|c_reactive|c-reactive", ignore_case = TRUE)
  )]
  crp_col <- if (length(crp_candidates) > 0) crp_candidates[1] else NA_character_

  # --- Reverse log if variable looks logged 
  maybe_exp <- function(x) {
    q95 <- suppressWarnings(stats::quantile(x, 0.95, na.rm = TRUE))
    if (is.finite(q95) && q95 < 20) exp(x) else x
  }

  # --- Unit heuristics 
  mmol_to_mgdl_chol <- function(x) x * 38.67
  mmol_to_mgdl_trig <- function(x) x * 88.57

  maybe_convert_chol_to_mgdl <- function(x) {
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (is.finite(med) && med < 30) mmol_to_mgdl_chol(x) else x
  }
  maybe_convert_trig_to_mgdl <- function(x) {
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (is.finite(med) && med < 30) mmol_to_mgdl_trig(x) else x
  }

  # ---------------------------
  # Build plotting data in publication-style units
  # ---------------------------
  tg_df <- tibble::tibble(
    glucose_mmol  = maybe_exp(tg$glucose_mmol),           # mmol/L
    hba1c_pct     = maybe_exp(tg$lnhba1c),                # %
    creat_umol    = maybe_exp(tg$lncreat_umol),           # µmol/L
    cyst_mgL      = tg$cyst,                              # mg/L
    trig_mgdl     = maybe_convert_trig_to_mgdl(maybe_exp(tg$trig)),   # mg/dL
    totchol_mgdl  = maybe_convert_chol_to_mgdl(tg$totchol),          # mg/dL
    ldl_mgdl      = maybe_convert_chol_to_mgdl(tg$ldl),              # mg/dL
    hdl_mgdl      = maybe_convert_chol_to_mgdl(tg$hdl)               # mg/dL
  )

  # Add CRP panel if present 
  if (!is.na(crp_col)) {
    tg_df <- tg_df %>%
      dplyr::mutate(crp_mgdl = maybe_exp(tg[[crp_col]]))  # mg/dL
  }

  # ---------------------------
  # Long format + labels
  # ---------------------------
  long <- tg_df %>%
    tidyr::pivot_longer(everything(), names_to = "Biomarker", values_to = "Value") %>%
    dplyr::filter(is.finite(Value)) %>%
    dplyr::mutate(
      Biomarker = factor(
        Biomarker,
        levels = c(
          "glucose_mmol","hba1c_pct","creat_umol","cyst_mgL",
          "crp_mgdl","trig_mgdl","totchol_mgdl","ldl_mgdl","hdl_mgdl"
        )
      ),
      BiomarkerLabel = dplyr::recode(
        Biomarker,
        glucose_mmol = "Serum glucose (mmol/L)",
        hba1c_pct    = "Glycohemoglobin (HbA1c, %)",
        creat_umol   = "Creatinine (µmol/L)",
        cyst_mgL     = "Cystatin C (mg/L)",
        crp_mgdl     = "C-reactive protein (mg/dL)",
        trig_mgdl    = "Triglyceride (mg/dL)",
        totchol_mgdl = "Total cholesterol (mg/dL)",
        ldl_mgdl     = "LDL cholesterol (mg/dL)",
        hdl_mgdl     = "HDL cholesterol (mg/dL)"
      )
    ) %>%
    # If CRP is missing, drop that facet cleanly:
    dplyr::filter(!(Biomarker == "crp_mgdl" & is.na(crp_col)))

  # ---------------------------
  # Plot
  # ---------------------------
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

  # Save
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(
    filename = file.path(outdir, "Supplementary_Figure1_TwinGene_density.png"),
    plot = p, width = 13.5, height = 10.5, units = "in", dpi = 300
  )

  return(p)
}
