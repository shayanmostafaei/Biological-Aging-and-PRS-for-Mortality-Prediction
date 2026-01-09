# ============================================================
# scripts/01_prepare_inputs.R
# Read + harmonize TwinGene and UKB to a shared schema
# ============================================================

read_input_twingene <- function(path) {
  if (!file.exists(path)) stop("TwinGene file not found: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

read_input_ukb <- function(path) {
  if (!file.exists(path)) stop("UKB file not found: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

# ---- Standardized names used downstream:
# time, event, sex, smoking, education, alcohol, bmi
# predictors: age, phenoage_res, kdm_res, hd_res, Frailty_Index_0.10_log, Telomere, PRS_Z
prep_twingene <- function(df) {

  # Required for TwinGene
  req <- c("age","time_to_event","status",
           "phenoage_res","kdm_res","hd_res",
           "Frailty_Index_0.10_log","Telomere","PRS_Z",
           "bmi","smoker","gender","education_years")
  stop_if_missing(df, req, "TwinGene")

  out <- df %>%
    mutate(
      cohort = "TwinGene",
      time   = as.numeric(time_to_event),
      event  = as.integer(status),
      sex    = as.integer(gender),    
      smoking = as.integer(smoker),
      education = as.numeric(education_years),
      bmi = as.numeric(bmi)
    ) %>%
    mutate(
      age = as.numeric(age),
      phenoage_res = as.numeric(phenoage_res),
      kdm_res      = as.numeric(kdm_res),
      hd_res       = as.numeric(hd_res),
      Frailty_Index_0.10_log = as.numeric(Frailty_Index_0.10_log),
      Telomere = as.numeric(Telomere),
      PRS_Z    = as.numeric(PRS_Z)
    )

  # sanity on event
  out$event[out$event != 0 & out$event != 1] <- NA_integer_

  out
}

prep_ukb <- function(df) {

  req_core <- c("age","time_to_death","death_status","bmi","sex","smoking","education")
  stop_if_missing(df, req_core, "UK Biobank")

  # Resolve residual names (support both styles)
  # phenoage_res priority:
  if (!"phenoage_res" %in% names(df)) {
    if ("Phenoage.Residual" %in% names(df)) df$phenoage_res <- df$Phenoage.Residual
    if ("Phenoage.Residual" %in% names(df) == FALSE && "Phenoage.Residual" %in% names(df) == FALSE) {
      # also support "Phenoage.Residual" only; if missing, try "Phenoage.Residual" etc.
    }
  }
  if (!"kdm_res" %in% names(df)) {
    if ("KDM.Residual" %in% names(df)) df$kdm_res <- df$KDM.Residual
  }
  if (!"hd_res" %in% names(df)) {
    if ("HD.Residual" %in% names(df)) df$hd_res <- df$HD.Residual
  }
  if (!"PRS_Z" %in% names(df)) {
    if ("PRS_Z" %in% names(df)) { } else if ("PRS_Z" %in% names(df) == FALSE && "PRS_Z" %in% names(df) == FALSE) {
      if ("PRS_Z" %in% names(df)) { }  
    }
  }
  if (!"Telomere" %in% names(df)) {
    if ("Telomere_Length.Residual" %in% names(df)) df$Telomere <- df$Telomere_Length.Residual
    if ("Telomeres" %in% names(df) && !"Telomere" %in% names(df)) df$Telomere <- df$Telomeres
  }
  if (!"Frailty_Index_0.10_log" %in% names(df)) {
    # If you have FI residual or FI, choose FI (or FI.Residual) depending on your porpuse. 
    if ("FI" %in% names(df)) df$Frailty_Index_0.10_log <- df$FI
    if ("FI.Residual" %in% names(df)) df$Frailty_Index_0.10_log <- df$FI.Residual
  }

  if (!"alcohol" %in% names(df)) {
    if ("alcohol" %in% names(df)) df$alcohol <- df$alcohol
  }

  req_preds <- c("phenoage_res","kdm_res","hd_res","Frailty_Index_0.10_log","Telomere","PRS_Z")
  stop_if_missing(df, req_preds, "UK Biobank (after residual mapping)")

  out <- df %>%
    mutate(
      cohort = "UK Biobank",
      time   = as.numeric(time_to_death),
      event  = as.integer(death_status),
      age    = as.numeric(age),
      sex    = as.integer(sex),
      smoking = as.integer(smoking),
      education = as.numeric(education),
      bmi = as.numeric(bmi),
      alcohol = if ("alcohol" %in% names(df)) as.numeric(alcohol) else NA_real_
    ) %>%
    mutate(
      phenoage_res = as.numeric(phenoage_res),
      kdm_res      = as.numeric(kdm_res),
      hd_res       = as.numeric(hd_res),
      Frailty_Index_0.10_log = as.numeric(Frailty_Index_0.10_log),
      Telomere = as.numeric(Telomere),
      PRS_Z    = as.numeric(PRS_Z)
    )

  out$event[out$event != 0 & out$event != 1] <- NA_integer_

  out
}
