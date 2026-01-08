# ============================================================
# scripts/03_superlearner_cv.R
# 10-fold CV SuperLearner + optional repeats with different seeds
# ============================================================

run_superlearner_models <- function(tg, ukb, outdir, K = 10, repeats = 10, seeds = 1001:1010) {

  quiet_require("SuperLearner")
  quiet_require("pROC")
  quiet_require("dplyr")
  quiet_require("readr")
  quiet_require("stringr")

  # ---------------------------
  # Helpers
  # ---------------------------
  resolve_var <- function(df, candidates) {
    candidates <- candidates[!is.na(candidates)]
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }

  find_pcs <- function(df) {
    pc <- names(df)[stringr::str_detect(names(df), "^PC\\d+$")]
    if (length(pc) == 0) pc <- names(df)[stringr::str_detect(names(df), "^PC")]
    if (length(pc) == 0) return(character(0))
    pc <- pc[order(pc)]
    head(pc, 10)
  }

  get_auc <- function(pred, y) {
    # direction="auto" by default, avoids AUC flipping for protective markers
    r <- pROC::roc(y, pred, quiet = TRUE)
    as.numeric(pROC::auc(r))
  }

  # ---------------------------
  # Outcome columns
  # ---------------------------
  y_tg  <- "death_status"
  y_ukb <- "death_status"

  if (!y_tg %in% names(tg))  stop("TwinGene outcome column not found: ", y_tg)
  if (!y_ukb %in% names(ukb)) stop("UKB outcome column not found: ", y_ukb)

  # Ensure 0/1 integer outcomes
  tg_y  <- as.integer(tg[[y_tg]])
  ukb_y <- as.integer(ukb[[y_ukb]])

  # ---------------------------
  # Resolve core predictors and covariates per cohort
  # ---------------------------
  # CA
  tg_age  <- resolve_var(tg,  c("age"))
  ukb_age <- resolve_var(ukb, c("age"))

  # PRS
  tg_prs  <- resolve_var(tg,  c("PRS_Z"))
  ukb_prs <- resolve_var(ukb, c("PRS_Z"))

  # BA components 
  tg_pheno <- resolve_var(tg,  c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual"))
  ukb_pheno <- resolve_var(ukb, c("phenoage_res", "Phenoage.Residual", "PhenoAge.Residual"))

  tg_fi <- resolve_var(tg,  c("Frailty_Index_0.10_log", "FI"))              
  ukb_fi <- resolve_var(ukb, c("FI", "fi_10percent"))       

  tg_tl <- resolve_var(tg,  c("Telomere", "Telomeres", "Telomere_Length.Residual"))
  ukb_tl <- resolve_var(ukb, c("Telomeres", "Telomere", "Telomere_Length"))

  # Covariates (cohort-available)
  tg_cov <- c(
    resolve_var(tg, c("gender", "sex")),
    resolve_var(tg, c("bmi", "BMI")),
    resolve_var(tg, c("smoker", "smoking")),
    resolve_var(tg, c("education_years", "education"))
  )
  tg_cov <- tg_cov[!is.na(tg_cov)]

  ukb_cov <- c(
    resolve_var(ukb, c("sex", "gender")),
    resolve_var(ukb, c("bmi", "BMI")),
    resolve_var(ukb, c("smoking", "smoker")),
    resolve_var(ukb, c("alcohol", "alcohol_consumption")),
    resolve_var(ukb, c("education", "education_years"))
  )
  ukb_cov <- ukb_cov[!is.na(ukb_cov)]

  # PCs if present
  tg_pcs  <- find_pcs(tg)
  ukb_pcs <- find_pcs(ukb)

  # ---------------------------
  # Model sets 
  # ---------------------------
  pred_sets_tg <- list(
    CA         = c(tg_age),
    CA_cov     = c(tg_age, tg_cov, tg_pcs),
    CA_cov_PRS = c(tg_age, tg_cov, tg_pcs, tg_prs),
    Full       = c(tg_age, tg_cov, tg_pcs, tg_prs, tg_pheno, tg_fi, tg_tl)
  )

  pred_sets_ukb <- list(
    CA         = c(ukb_age),
    CA_cov     = c(ukb_age, ukb_cov, ukb_pcs),
    CA_cov_PRS = c(ukb_age, ukb_cov, ukb_pcs, ukb_prs),
    Full       = c(ukb_age, ukb_cov, ukb_pcs, ukb_prs, ukb_pheno, ukb_fi, ukb_tl)
  )

  # Drop NAs and keep only columns that exist
  clean_set <- function(df, v) unique(v[!is.na(v) & v %in% names(df)])
  pred_sets_tg  <- lapply(pred_sets_tg,  clean_set, df = tg)
  pred_sets_ukb <- lapply(pred_sets_ukb, clean_set, df = ukb)

  # ---------------------------
  # Learners 
  # ---------------------------
  learners <- c(
    "SL.glm",
    "SL.glmnet",
    "SL.randomForest",
    "SL.xgboost"
  )

  # ---------------------------
  # Fit one cohort
  # ---------------------------
  fit_one <- function(df, y_vec, pred_sets, cohort_name, seeds_local) {

    run_set_one_seed <- function(pred_names, seed) {

      X <- df[, pred_names, drop = FALSE]
      ok <- stats::complete.cases(X) & !is.na(y_vec)
      X2 <- X[ok, , drop = FALSE]
      y2 <- y_vec[ok]

      # guardrails
      if (nrow(X2) < 200) return(list(auc = NA_real_))

      set.seed(seed)
      sl <- SuperLearner::CV.SuperLearner(
        Y = y2,
        X = X2,
        family = stats::binomial(),
        SL.library = learners,
        V = K
      )

      pred_oof <- as.numeric(sl$SL.predict)
      auc <- get_auc(pred_oof, y2)

      list(auc = auc, n = length(y2))
    }

    res <- lapply(seeds_local, function(sd) {
      one_seed <- lapply(names(pred_sets), function(m) {
        r <- run_set_one_seed(pred_sets[[m]], sd)
        data.frame(
          Cohort = cohort_name,
          Model  = m,
          Seed   = sd,
          AUC    = r$auc,
          N      = r$n,
          stringsAsFactors = FALSE
        )
      })
      dplyr::bind_rows(one_seed)
    }) %>% dplyr::bind_rows()

    res
  }

  # ---------------------------
  # Seeds handling
  # ---------------------------
  if (repeats <= 1) {
    seeds_use <- seeds[1]
  } else {
    seeds_use <- seeds[seq_len(min(repeats, length(seeds)))]
  }

  # ---------------------------
  # Run cohorts
  # ---------------------------
  tg_res  <- fit_one(tg,  tg_y,  pred_sets_tg,  "TwinGene",  seeds_use)
  ukb_res <- fit_one(ukb, ukb_y, pred_sets_ukb, "UK Biobank", seeds_use)

  all_res <- dplyr::bind_rows(tg_res, ukb_res)

  # ---------------------------
  # Summaries for Table 3 + robustness reporting
  # ---------------------------
  summ <- all_res %>%
    dplyr::group_by(Cohort, Model) %>%
    dplyr::summarise(
      AUC_mean = mean(AUC, na.rm = TRUE),
      AUC_sd   = sd(AUC, na.rm = TRUE),
      AUC_min  = min(AUC, na.rm = TRUE),
      AUC_max  = max(AUC, na.rm = TRUE),
      Repeats  = sum(is.finite(AUC)),
      N_used   = round(mean(N, na.rm = TRUE)),
      .groups  = "drop"
    )

  readr::write_csv(all_res, file.path(outdir, "superlearner_cv_all_seeds.csv"))
  readr::write_csv(summ,    file.path(outdir, "superlearner_cv_summary.csv"))

  # Also write the exact variable sets used (reproducibility)
  used_sets <- dplyr::bind_rows(
    dplyr::bind_rows(lapply(names(pred_sets_tg), function(m) {
      data.frame(Cohort="TwinGene", Model=m, Variables=paste(pred_sets_tg[[m]], collapse="; "))
    })),
    dplyr::bind_rows(lapply(names(pred_sets_ukb), function(m) {
      data.frame(Cohort="UK Biobank", Model=m, Variables=paste(pred_sets_ukb[[m]], collapse="; "))
    }))
  )
  readr::write_csv(used_sets, file.path(outdir, "superlearner_variable_sets.csv"))

  list(all = all_res, summary = summ, variable_sets = used_sets)
}
