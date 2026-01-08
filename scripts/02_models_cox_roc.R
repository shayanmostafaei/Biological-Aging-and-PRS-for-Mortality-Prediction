# ============================================================
# scripts/02_models_cox_roc.R
# Univariate ROC (secondary) + Cox PH models (primary)
# TwinGene discovery + UK Biobank external validation
# ============================================================

run_univariate_roc <- function(tg, ukb, outdir) {

  quiet_require("pROC")
  quiet_require("dplyr")
  quiet_require("readr")

  # ---------------------------
  # Helper: resolve a variable name from candidates
  # ---------------------------
  resolve_var <- function(df, candidates) {
    candidates <- candidates[!is.na(candidates)]
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }

  # ---------------------------
  # Predictor sets (labels -> candidate column names)
  # ---------------------------
  preds_tg <- list(
    "CA"               = c("age"),
    "PhenoAge_residual"= c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual"),
    "KDM_residual"     = c("kdm_res", "KDM.Residual"),
    "HD_residual"      = c("hd_res", "HD.Residual"),
    "FI"               = c("Frailty_Index_0.10_log", "FI"),
    "TL"               = c("Telomere", "Telomeres"),
    "PRS"              = c("PRS_Z")
  )

  preds_ukb <- list(
    "CA"               = c("age"),
    "PhenoAge_residual"= c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual"),
    "KDM_residual"     = c("kdm_res", "KDM.Residual"),
    "HD_residual"      = c("hd_res", "HD.Residual"),
    "FI"               = c("FI", "fi_10percent"),
    "TL"               = c("Telomeres", "Telomere", "Telomere_Length"),
    "PRS"              = c("PRS_Z")
  )

  # ---------------------------
  # Core ROC worker
  # ---------------------------
  roc_one <- function(df, status_col, pred_map, cohort_name) {

    out <- lapply(names(pred_map), function(lbl) {

      v <- resolve_var(df, pred_map[[lbl]])
      if (is.na(v)) return(NULL)

      x <- df[[v]]
      y <- df[[status_col]]

      ok <- is.finite(x) & !is.na(y)
      if (sum(ok) < 50) return(NULL)

      # direction="auto" (default) avoids AUC<0.5 for protective predictors
      r  <- pROC::roc(y[ok], x[ok], quiet = TRUE)
      ci <- pROC::ci.auc(r)

      data.frame(
        Cohort        = cohort_name,
        Predictor     = lbl,
        Predictor_col = v,
        AUC           = as.numeric(pROC::auc(r)),
        L             = as.numeric(ci[1]),
        U             = as.numeric(ci[3]),
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(out)
  }

  tg_roc  <- roc_one(tg,  status_col = "status",       pred_map = preds_tg,  cohort_name = "TwinGene")
  ukb_roc <- roc_one(ukb, status_col = "death_status", pred_map = preds_ukb, cohort_name = "UK Biobank")

  roc_tab <- dplyr::bind_rows(tg_roc, ukb_roc)
  readr::write_csv(roc_tab, file.path(outdir, "univariate_ROC_AUCs.csv"))

  roc_tab
}


run_cox_models <- function(tg, ukb, outdir) {

  quiet_require("survival")
  quiet_require("dplyr")
  quiet_require("readr")
  quiet_require("stringr")

  # ---------------------------
  # Helper: resolve variable name from candidates
  # ---------------------------
  resolve_var <- function(df, candidates) {
    candidates <- candidates[!is.na(candidates)]
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }

  # ---------------------------
  # Predictors (same “labels” as ROC)
  # ---------------------------
  preds_tg <- list(
    "CA"               = c("age"),
    "PhenoAge_residual"= c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual"),
    "KDM_residual"     = c("kdm_res", "KDM.Residual"),
    "HD_residual"      = c("hd_res", "HD.Residual"),
    "FI"               = c("Frailty_Index_0.10_log", "FI"),
    "TL"               = c("Telomere", "Telomeres"),
    "PRS"              = c("PRS_Z")
  )

  preds_ukb <- list(
    "CA"               = c("age"),
    "PhenoAge_residual"= c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual"),
    "KDM_residual"     = c("kdm_res", "KDM.Residual"),
    "HD_residual"      = c("hd_res", "HD.Residual"),
    "FI"               = c("FI", "fi_10percent"),
    "TL"               = c("Telomeres", "Telomere", "Telomere_Length"),
    "PRS"              = c("PRS_Z")
  )

  # ---------------------------
  # Identify PCs if present (PC1..PC10, or any PC* then keep first 10)
  # ---------------------------
  find_pcs <- function(df) {
    pc <- names(df)[stringr::str_detect(names(df), "^PC\\d+$")]
    if (length(pc) == 0) pc <- names(df)[stringr::str_detect(names(df), "^PC")]
    if (length(pc) == 0) return(character(0))
    pc <- pc[order(pc)]
    head(pc, 10)
  }

  pcs_tg  <- find_pcs(tg)
  pcs_ukb <- find_pcs(ukb)

  # ---------------------------
  # Univariate Cox (time-since-baseline)
  # ---------------------------
  cox_uni <- function(df, time_col, status_col, pred_label, pred_candidates, cohort_name) {

    v <- resolve_var(df, pred_candidates)
    if (is.na(v)) return(NULL)

    f <- stats::as.formula(paste0("survival::Surv(", time_col, ",", status_col, ") ~ ", v))
    fit <- survival::coxph(f, data = df)
    s <- summary(fit)

    data.frame(
      Cohort         = cohort_name,
      Model          = "Univariate",
      Predictor      = pred_label,
      Predictor_col  = v,
      HR             = unname(s$coef[,"exp(coef)"]),
      L              = unname(s$conf.int[,"lower .95"]),
      U              = unname(s$conf.int[,"upper .95"]),
      P              = unname(s$coef[,"Pr(>|z|)"]),
      stringsAsFactors = FALSE
    )
  }

  tg_uni <- dplyr::bind_rows(lapply(names(preds_tg), function(lbl) {
    cox_uni(tg, "time_to_death", "death_status", lbl, preds_tg[[lbl]], "TwinGene")
  }))

  ukb_uni <- dplyr::bind_rows(lapply(names(preds_ukb), function(lbl) {
    cox_uni(ukb, "time_to_death", "death_status", lbl, preds_ukb[[lbl]], "UK Biobank")
  }))

  # ---------------------------
  # Multivariable Cox (cohort-available covariates)
  # TwinGene: sex, BMI, smoking, education (+ PCs if present)
  # UKB: sex, BMI, smoking, alcohol, education (+ PCs if present)
  # Plus CA + PRS + BA measures used in the manuscript
  # ---------------------------
  tg_cov_candidates <- list(
    "gender"          = c("gender", "sex"),
    "bmi"             = c("bmi", "BMI"),
    "smoker"          = c("smoker", "smoking"),
    "education_years" = c("education_years", "education")
  )

  ukb_cov_candidates <- list(
    "sex"       = c("sex", "gender"),
    "bmi"       = c("bmi", "BMI"),
    "smoking"   = c("smoking", "smoker"),
    "alcohol"   = c("alcohol", "alcohol_consumption"),
    "education" = c("education", "education_years")
  )

  resolve_terms <- function(df, term_candidates_list) {
    out <- character(0)
    for (nm in names(term_candidates_list)) {
      v <- resolve_var(df, term_candidates_list[[nm]])
      if (!is.na(v)) out <- c(out, v)
    }
    out
  }

  # Predictors to include in multivariable model 
  tg_core_preds <- c(
    resolve_var(tg, c("age")),
    resolve_var(tg, c("PRS_Z")),
    resolve_var(tg, c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual")),
    resolve_var(tg, c("Frailty_Index_0.10_log", "FI")),
    resolve_var(tg, c("Telomere", "Telomeres"))
  )

  ukb_core_preds <- c(
    resolve_var(ukb, c("age")),
    resolve_var(ukb, c("PRS_Z")),
    resolve_var(ukb, c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual")),
    resolve_var(ukb, c("FI", "fi_10percent")),
    resolve_var(ukb, c("Telomeres", "Telomere", "Telomere_Length"))
  )

  tg_terms  <- unique(na.omit(c(resolve_terms(tg,  tg_cov_candidates),  pcs_tg,  tg_core_preds)))
  ukb_terms <- unique(na.omit(c(resolve_terms(ukb, ukb_cov_candidates), pcs_ukb, ukb_core_preds)))

  fit_multi <- function(df, time_col, status_col, terms, cohort_name) {

    # keep only terms that exist
    terms <- terms[terms %in% names(df)]

    f <- stats::as.formula(
      paste0("survival::Surv(", time_col, ",", status_col, ") ~ ", paste(terms, collapse = " + "))
    )

    fit <- survival::coxph(f, data = df)
    s <- summary(fit)

    data.frame(
      Cohort        = cohort_name,
      Model         = "Multivariable",
      Predictor     = rownames(s$coef),
      HR            = unname(s$coef[,"exp(coef)"]),
      L             = unname(s$conf.int[,"lower .95"]),
      U             = unname(s$conf.int[,"upper .95"]),
      P             = unname(s$coef[,"Pr(>|z|)"]),
      stringsAsFactors = FALSE
    )
  }

  tg_mv  <- fit_multi(tg,  "time_to_death", "death_status",   tg_terms,  "TwinGene")
  ukb_mv <- fit_multi(ukb, "time_to_death", "death_status", ukb_terms, "UK Biobank")

  out <- dplyr::bind_rows(tg_uni, ukb_uni, tg_mv, ukb_mv)
  readr::write_csv(out, file.path(outdir, "cox_models_time_since_baseline.csv"))

  out
}
