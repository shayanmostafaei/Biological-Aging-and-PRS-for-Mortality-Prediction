# ============================================================
# scripts/05_subgroup_interactions.R
# AUC by Age group × Horizon (5/10/15) faceted by Predictor and Cohort
# Interaction p-value significance map (formatted with <0.001)
# ============================================================

# ---- AUC at fixed horizon as binary outcome:
auc_horizon_binary <- function(df, horizon_years, pred) {
  # time in years
  t <- df$time
  e <- df$event

  y <- ifelse(e == 1 & t <= horizon_years, 1,
              ifelse(t > horizon_years, 0, NA))
  x <- df[[pred]]

  ok <- is.finite(x) & !is.na(y)
  if (sum(ok) < 30) return(NA_real_)

  r <- pROC::roc(y[ok], x[ok], quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(r))
}

run_subgroup_followup_age <- function(tg, ukb, outdir,
                                      horizons = c(5,10,15),
                                      age_breaks = c(-Inf, 50, 60, Inf)) {

  quiet_require("pROC")
  quiet_require("ggnewscale")

  predictors <- c("PRS_Z","phenoage_res","kdm_res","hd_res",
                  "Frailty_Index_0.10_log","Telomere")

  age_labels <- c("<50","50–60",">60")

  make_age_group <- function(age) {
    cut(age, breaks = age_breaks, labels = age_labels, right = FALSE)
  }

  make_grid <- function(df, cohort_name) {
    df <- df %>% mutate(AgeGroup = make_age_group(age))
    expand.grid(
      Cohort = cohort_name,
      Predictor = predictors,
      AgeGroup = age_labels,
      FollowUpYears = horizons,
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      rowwise() %>%
      mutate(
        AUC = {
          dsub <- df %>% filter(AgeGroup == .data$AgeGroup)
          if (nrow(dsub) < 30) NA_real_ else auc_horizon_binary(dsub, .data$FollowUpYears, .data$Predictor)
        }
      ) %>%
      ungroup()
  }

  s2_tg  <- make_grid(tg,  "TwinGene")
  s2_ukb <- make_grid(ukb, "UK Biobank")
  s2_all <- bind_rows(s2_tg, s2_ukb)

  readr::write_csv(s2_all, file.path(outdir, "TableS2_AUC_by_age_followup.csv"))

  tg_fill <- COL_TG_HEAT
  uk_fill <- COL_UKB_HEAT

  df_tg  <- s2_all %>% filter(Cohort == "TwinGene")
  df_ukb <- s2_all %>% filter(Cohort == "UK Biobank")

  df_tg <- df_tg %>%
    mutate(
      FollowUpYears = factor(FollowUpYears, levels = horizons),
      AgeGroup      = factor(AgeGroup, levels = age_labels),
      Predictor     = factor(Predictor, levels = c("FI","HD","KDM","PhenoAge","PRS","TL")) %>% as.character()
    )

  # Ensure Predictor order/labels (short names)
  pred_levels <- c("FI","HD","KDM","PhenoAge","PRS","TL")
  pred_lab <- c(
    "Frailty_Index_0.10_log" = "FI",
    "hd_res" = "HD",
    "kdm_res" = "KDM",
    "phenoage_res" = "PhenoAge",
    "PRS_Z" = "PRS",
    "Telomere" = "TL"
  )

  prep_plot <- function(df) {
    df %>%
      mutate(
        Predictor = recode(Predictor, !!!pred_lab),
        Predictor = factor(Predictor, levels = pred_levels),
        FollowUpYears = factor(FollowUpYears, levels = horizons),
        AgeGroup = factor(AgeGroup, levels = age_labels),
        label_auc = ifelse(is.na(AUC), "", sprintf("%.3f", AUC))
      )
  }

  df_tg  <- prep_plot(df_tg)
  df_ukb <- prep_plot(df_ukb)

  p <- ggplot() +
    geom_tile(data = df_tg, aes(FollowUpYears, AgeGroup, fill = AUC),
              color = "white", linewidth = 0.7) +
    geom_text(data = df_tg, aes(FollowUpYears, AgeGroup, label = label_auc), size = 3.6) +
    scale_fill_gradient(
      name = "AUC (TwinGene)",
      low  = "#FCE0D6", high = tg_fill,
      limits = c(0.45, 1.00), oob = squish
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = df_ukb, aes(FollowUpYears, AgeGroup, fill = AUC),
              color = "white", linewidth = 0.7) +
    geom_text(data = df_ukb, aes(FollowUpYears, AgeGroup, label = label_auc), size = 3.6) +
    scale_fill_gradient(
      name = "AUC (UK Biobank)",
      low  = "#D9ECFF", high = uk_fill,
      limits = c(0.45, 1.00), oob = squish
    ) +
    facet_grid(Cohort ~ Predictor, switch = "y") +
    labs(
      title = "Figure S2. Univariate AUCs by baseline age group and follow-up horizon",
      x = "Follow-up duration (years)",
      y = "Age group (years)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )

  ggplot2::ggsave(file.path(outdir, "Figure_S2_AUC_heatmap.png"), p,
                  width = 13.5, height = 7.5, units = "in", dpi = 300)

  list(table = s2_all, plot = p)
}

# ---- Interaction p-values (Cox LRT comparing models with vs without interaction)
interaction_pvalue_lrt <- function(df, modifier, ba, covars) {
  # base model
  terms0 <- c(covars, ba, modifier)
  terms0 <- terms0[terms0 %in% names(df)]
  f0 <- as.formula(paste0("survival::Surv(time, event) ~ ", paste(terms0, collapse = " + ")))

  # interaction model
  f1 <- as.formula(paste0("survival::Surv(time, event) ~ ",
                          paste(covars, ba, modifier, paste0(ba, ":", modifier)), collapse = " + ")))

  # Ensure required columns exist
  needed <- unique(c("time","event", covars, ba, modifier))
  needed <- needed[needed %in% names(df)]
  dsub <- df %>% dplyr::select(all_of(needed)) %>% na.omit()

  if (nrow(dsub) < 200) return(NA_real_)

  fit0 <- survival::coxph(f0, data = dsub)
  fit1 <- survival::coxph(f1, data = dsub)

  an <- anova(fit0, fit1, test = "LRT")
  as.numeric(an[2, "P(>|Chi|)"])
}

run_interaction_heatmaps <- function(tg, ukb, outdir) {

  quiet_require("survival")
  quiet_require("ggnewscale")

  modifiers <- c("bmi","sex","smoking","education")
  mod_labels <- c("bmi"="BMI", "sex"="Sex", "smoking"="Smoking", "education"="Education")

  predictors <- c("phenoage_res","kdm_res","hd_res","Frailty_Index_0.10_log","Telomere","PRS_Z")
  pred_labels <- c(
    "phenoage_res"="PhenoAge",
    "kdm_res"="KDM",
    "hd_res"="HD",
    "Frailty_Index_0.10_log"="FI",
    "Telomere"="TL",
    "PRS_Z"="PRS"
  )

  # Cohort covariates (exclude modifier itself from covars when testing that modifier)
  cov_base <- c("age","sex","bmi","smoking","education","alcohol")

  one_cohort <- function(df, cohort_name) {
    out <- expand.grid(
      Cohort = cohort_name,
      InteractionFactor = modifiers,
      Predictor = predictors,
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      rowwise() %>%
      mutate(
        P_value = {
          mod <- InteractionFactor
          covars <- setdiff(cov_base, mod)
          covars <- covars[covars %in% names(df)]
          interaction_pvalue_lrt(df, mod, Predictor, covars)
        }
      ) %>%
      ungroup()

    out
  }

  s3 <- bind_rows(
    one_cohort(tg,  "TwinGene"),
    one_cohort(ukb, "UK Biobank")
  ) %>%
    mutate(
      InteractionFactor = recode(InteractionFactor, !!!mod_labels),
      Predictor = recode(Predictor, !!!pred_labels),
      minuslog10p = -log10(P_value),
      p_label = format_p(P_value),
      InteractionFactor = factor(InteractionFactor, levels = c("BMI","Sex","Smoking","Education")),
      Predictor = factor(Predictor, levels = c("PhenoAge","KDM","HD","FI","TL","PRS"))
    )

  readr::write_csv(s3, file.path(outdir, "TableS3_interaction_pvalues.csv"))

  tg_fill <- COL_TG_HEAT
  uk_fill <- COL_UKB_HEAT

  df_tg  <- s3 %>% filter(Cohort == "TwinGene")
  df_ukb <- s3 %>% filter(Cohort == "UK Biobank")

  p <- ggplot() +
    geom_tile(data = df_tg, aes(InteractionFactor, Predictor, fill = minuslog10p),
              color = "white", linewidth = 0.6) +
    geom_text(data = df_tg, aes(InteractionFactor, Predictor, label = p_label),
              size = 3.6) +
    scale_fill_gradient(
      name = expression(-log[10](P)~"(TwinGene)"),
      low  = "#FCE0D6", high = tg_fill
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = df_ukb, aes(InteractionFactor, Predictor, fill = minuslog10p),
              color = "white", linewidth = 0.6) +
    geom_text(data = df_ukb, aes(InteractionFactor, Predictor, label = p_label),
              size = 3.6) +
    scale_fill_gradient(
      name = expression(-log[10](P)~"(UK Biobank)"),
      low  = "#D9ECFF", high = uk_fill
    ) +
    facet_grid(Cohort ~ ., switch = "y") +
    labs(
      title = "Figure S3. Interaction significance map for lifestyle and demographic modifiers",
      x = "Effect modifier",
      y = "Predictor"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )

  ggplot2::ggsave(file.path(outdir, "Figure_S3_interactions_heatmap.png"), p,
                  width = 10.5, height = 7.5, units = "in", dpi = 300)

  list(table = s3, plot = p)
}
