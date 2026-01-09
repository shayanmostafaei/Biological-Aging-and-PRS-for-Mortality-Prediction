# ============================================================
# scripts/04_superlearner_cv.R
# 10-fold CV SuperLearner + optional repeats with different seeds
# ============================================================

run_superlearner_models <- function(tg, ukb, outdir, K = 10, repeats = 10, seeds = 1001:1010) {

  quiet_require("SuperLearner")
  quiet_require("pROC")

  # Predictor sets 
  cov_common <- c("sex","bmi","smoking","education")
  cov_ukb_extra <- c("alcohol")

  pred_sets_tg <- list(
    CA         = c("age"),
    CA_cov     = c("age", cov_common),
    CA_cov_PRS = c("age", cov_common, "PRS_Z"),
    Full       = c("age", cov_common, "PRS_Z", "phenoage_res", "Frailty_Index_0.10_log", "Telomere")
  )

  pred_sets_ukb <- list(
    CA         = c("age"),
    CA_cov     = c("age", cov_common, cov_ukb_extra),
    CA_cov_PRS = c("age", cov_common, cov_ukb_extra, "PRS_Z"),
    Full       = c("age", cov_common, cov_ukb_extra, "PRS_Z", "phenoage_res", "Frailty_Index_0.10_log", "Telomere")
  )

  learners <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.xgboost")

  get_auc <- function(pred, y) {
    r <- pROC::roc(y, pred, quiet = TRUE, direction = "<")
    as.numeric(pROC::auc(r))
  }

  run_one_seed <- function(df, pred_names, seed) {
    pred_names <- pred_names[pred_names %in% names(df)]
    X <- df[, pred_names, drop = FALSE]
    y <- as.integer(df$event)

    ok <- complete.cases(X) & !is.na(y)
    X <- X[ok, , drop = FALSE]
    y <- y[ok]

    set.seed(seed)
    sl <- SuperLearner::CV.SuperLearner(
      Y = y,
      X = X,
      family = stats::binomial(),
      SL.library = learners,
      V = K
    )
    pred_oof <- sl$SL.predict
    auc <- get_auc(pred_oof, y)
    list(auc = auc, n = length(y))
  }

  if (repeats <= 1) {
    seeds_use <- seeds[1]
  } else {
    seeds_use <- seeds[seq_len(min(repeats, length(seeds)))]
  }

  fit_cohort <- function(df, pred_sets, cohort_name) {
    res <- purrr::map_dfr(seeds_use, function(sd) {
      purrr::map_dfr(names(pred_sets), function(m) {
        r <- run_one_seed(df, pred_sets[[m]], sd)
        data.frame(
          Cohort = cohort_name,
          Model  = m,
          Seed   = sd,
          AUC    = r$auc,
          N_used = r$n
        )
      })
    })
    res
  }

  tg_res  <- fit_cohort(tg,  pred_sets_tg,  "TwinGene")
  ukb_res <- fit_cohort(ukb, pred_sets_ukb, "UK Biobank")

  all_res <- bind_rows(tg_res, ukb_res)

  summ <- all_res %>%
    group_by(Cohort, Model) %>%
    summarise(
      AUC_mean = mean(AUC),
      AUC_sd   = sd(AUC),
      Seeds    = n_distinct(Seed),
      .groups  = "drop"
    )

  readr::write_csv(all_res, file.path(outdir, "superlearner_cv_all_seeds.csv"))
  readr::write_csv(summ,    file.path(outdir, "superlearner_cv_summary.csv"))

  list(all = all_res, summary = summ)
}
