# ============================================================
# scripts/03_models_cox_roc.R
# Univariate ROC (secondary) + Cox PH models (primary)
# Uses standardized columns:
#   time, event, cohort, predictors
# ============================================================

run_univariate_roc <- function(tg, ukb, outdir) {

  quiet_require("pROC")

  pred_list <- c("age","phenoage_res","kdm_res","hd_res",
                 "Frailty_Index_0.10_log","Telomere","PRS_Z")

  roc_one <- function(df, cohort_name) {
    out <- lapply(pred_list, function(v) {
      if (!v %in% names(df)) return(NULL)
      x <- df[[v]]
      y <- df$event
      ok <- is.finite(x) & !is.na(y)
      if (sum(ok) < 50) return(NULL)
      r  <- pROC::roc(y[ok], x[ok], quiet = TRUE, direction = "<")
      ci <- pROC::ci.auc(r)
      data.frame(
        Cohort = cohort_name,
        Predictor = v,
        AUC = as.numeric(pROC::auc(r)),
        L = as.numeric(ci[1]),
        U = as.numeric(ci[3])
      )
    })
    dplyr::bind_rows(out)
  }

  tg_roc  <- roc_one(tg,  "TwinGene")
  ukb_roc <- roc_one(ukb, "UK Biobank")

  roc_tab <- dplyr::bind_rows(tg_roc, ukb_roc)
  readr::write_csv(roc_tab, file.path(outdir, "univariate_ROC_AUCs.csv"))

  roc_tab
}

run_cox_models <- function(tg, ukb, outdir) {

  quiet_require("survival")

  preds <- c("age","phenoage_res","kdm_res","hd_res",
             "Frailty_Index_0.10_log","Telomere","PRS_Z")

  cox_uni <- function(df, pred, cohort_name) {
    if (!pred %in% names(df)) return(NULL)
    f <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", pred))
    fit <- survival::coxph(f, data = df)
    s <- summary(fit)
    data.frame(
      Cohort = cohort_name,
      Model = "Univariate (time-since-baseline)",
      Predictor = pred,
      HR = unname(s$coef[,"exp(coef)"]),
      L  = unname(s$conf.int[,"lower .95"]),
      U  = unname(s$conf.int[,"upper .95"]),
      P  = unname(s$coef[,"Pr(>|z|)"])
    )
  }

  tg_uni  <- bind_rows(lapply(preds, cox_uni, df = tg,  cohort_name = "TwinGene"))
  ukb_uni <- bind_rows(lapply(preds, cox_uni, df = ukb, cohort_name = "UK Biobank"))

  # Multivariable covariates (cohort-available)
  tg_covars  <- c("age","sex","bmi","smoking","education")
  ukb_covars <- c("age","sex","bmi","smoking","alcohol","education")

  tg_terms  <- c(tg_covars,  "PRS_Z","phenoage_res","Frailty_Index_0.10_log","Telomere")
  ukb_terms <- c(ukb_covars, "PRS_Z","phenoage_res","Frailty_Index_0.10_log","Telomere")

  fit_multi <- function(df, terms, cohort_name) {
    terms <- terms[terms %in% names(df)]
    f <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", paste(terms, collapse = " + ")))
    fit <- survival::coxph(f, data = df)
    s <- summary(fit)
    data.frame(
      Cohort = cohort_name,
      Model = "Multivariable (time-since-baseline)",
      Predictor = rownames(s$coef),
      HR = unname(s$coef[,"exp(coef)"]),
      L  = unname(s$conf.int[,"lower .95"]),
      U  = unname(s$conf.int[,"upper .95"]),
      P  = unname(s$coef[,"Pr(>|z|)"])
    )
  }

  tg_mv  <- fit_multi(tg,  tg_terms,  "TwinGene")
  ukb_mv <- fit_multi(ukb, ukb_terms, "UK Biobank")

  out <- bind_rows(tg_uni, ukb_uni, tg_mv, ukb_mv)
  readr::write_csv(out, file.path(outdir, "cox_models_time_since_baseline.csv"))

  out
}
