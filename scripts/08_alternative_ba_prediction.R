# ============================================================
# 08_alternative_ba_prediction.R
#
# Conventional
# Conventional + PhenoAge residual + FI + TL
# Conventional + KDM residual + FI + TL
# Conventional + log-HD + FI + TL
# ============================================================

source(
    "scripts/00_setup.R"
)


cohorts <-
    load_analysis_cohorts()


tg <-
    cohorts$TwinGene


ukb <-
    cohorts$UKB


# Reuse EXACT primary folds
tg_fold <-
    readRDS(
        file.path(
            PRIVATE_OUTPUT_DIR,
            "TwinGene_10fold_partition.rds"
        )
    )


ukb_fold <-
    readRDS(
        file.path(
            PRIVATE_OUTPUT_DIR,
            "UKB_314064_10fold_partition.rds"
        )
    )


tg$fold <-
    tg_fold$fold


stopifnot(
    identical(
        as.character(
            ukb$id
        ),
        as.character(
            ukb_fold$id
        )
    )
)


ukb$fold <-
    ukb_fold$fold


formula_conventional <-
    Surv(
        time,
        event
    ) ~
        splines::ns(
            age,
            df = 4
        ) +
        sex_model +
        splines::ns(
            bmi,
            df = 3
        ) +
        smoke_model +
        edu_model


formula_pheno <-
    update(
        formula_conventional,
        . ~
            . +
            pheno_z_cv +
            fi_z_cv +
            tl_z_cv
    )


formula_kdm <-
    update(
        formula_conventional,
        . ~
            . +
            kdm_z_cv +
            fi_z_cv +
            tl_z_cv
    )


formula_hd <-
    update(
        formula_conventional,
        . ~
            . +
            hd_z_cv +
            fi_z_cv +
            tl_z_cv
    )


formulas <- list(

    Conventional =
        formula_conventional,

    PhenoAge =
        formula_pheno,

    KDM =
        formula_kdm,

    HD =
        formula_hd
)


tg_oof <-
    cv_cox_models(
        tg,
        formulas
    )


ukb_oof <-
    cv_cox_models(
        ukb,
        formulas
    )


performance <-
    bind_rows(

        evaluate_oof(
            tg,
            tg_oof,
            "TwinGene"
        ),

        evaluate_oof(
            ukb,
            ukb_oof,
            "UKB"
        )
    )


ref <- performance %>%

    filter(
        model ==
            "Conventional"
    ) %>%

    select(
        cohort,
        ref_C =
            C_index,
        ref_AUC5 =
            AUC_5,
        ref_AUC10 =
            AUC_10,
        ref_Brier5 =
            Brier_5,
        ref_Brier10 =
            Brier_10
    )


performance <-
    performance %>%

    left_join(
        ref,
        by = "cohort"
    ) %>%

    mutate(

        Delta_C =
            C_index -
            ref_C,

        Delta_AUC5 =
            AUC_5 -
            ref_AUC5,

        Delta_AUC10 =
            AUC_10 -
            ref_AUC10,

        Delta_Brier5 =
            Brier_5 -
            ref_Brier5,

        Delta_Brier10 =
            Brier_10 -
            ref_Brier10
    )


print(
    performance,
    n = Inf,
    width = Inf
)


write_csv(
    performance,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Alternative_BA_prediction_models.csv"
    )
)


# OOF predictions remain private
saveRDS(
    tg_oof,
    file.path(
        PRIVATE_OUTPUT_DIR,
        "TwinGene_alternative_BA_OOF.rds"
    )
)


saveRDS(
    ukb_oof,
    file.path(
        PRIVATE_OUTPUT_DIR,
        "UKB_alternative_BA_OOF.rds"
    )
)
