# ============================================================
# 04_prediction_cv.R
# Fixed 10-fold out-of-fold mortality prediction
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


# ============================================================
# Generate / save fixed folds
# ============================================================

tg_fold_file <-
    file.path(
        PRIVATE_OUTPUT_DIR,
        "TwinGene_10fold_partition.rds"
    )


ukb_fold_file <-
    file.path(
        PRIVATE_OUTPUT_DIR,
        "UKB_314064_10fold_partition.rds"
    )


if (
    file.exists(
        tg_fold_file
    )
) {

    fold <-
        readRDS(
            tg_fold_file
        )

    tg$fold <-
        fold$fold

} else {

    tg <-
        make_twingene_folds(
            tg
        )

    saveRDS(
        list(
            id = tg$id,
            pair_id = tg$pair_id,
            fold = tg$fold,
            seed = RANDOM_SEED
        ),
        tg_fold_file
    )
}


if (
    file.exists(
        ukb_fold_file
    )
) {

    fold <-
        readRDS(
            ukb_fold_file
        )

    stopifnot(
        identical(
            as.character(
                fold$id
            ),
            as.character(
                ukb$id
            )
        )
    )

    ukb$fold <-
        fold$fold

} else {

    ukb <-
        make_ukb_folds(
            ukb
        )

    saveRDS(
        list(
            id = ukb$id,
            fold = ukb$fold,
            seed = RANDOM_SEED
        ),
        ukb_fold_file
    )
}


# ============================================================
# Fold QC
# ============================================================

fold_qc <-
    bind_rows(

        tg %>%
            group_by(
                fold
            ) %>%
            summarise(
                Cohort = "TwinGene",
                N = n(),
                Deaths = sum(event),
                .groups = "drop"
            ),

        ukb %>%
            group_by(
                fold
            ) %>%
            summarise(
                Cohort = "UKB",
                N = n(),
                Deaths = sum(event),
                .groups = "drop"
            )
    )


write_csv(
    fold_qc,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Prediction_fold_QC.csv"
    )
)


# ============================================================
# Primary prediction models
# ============================================================

formula_age <-
    Surv(
        time,
        event
    ) ~
        splines::ns(
            age,
            df = 4
        )


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


formula_enhanced <-
    update(
        formula_conventional,
        . ~
            . +
            pheno_z_cv +
            fi_z_cv +
            tl_z_cv
    )


formulas <- list(

    Age =
        formula_age,

    Conventional =
        formula_conventional,

    Enhanced =
        formula_enhanced
)


# ============================================================
# Run CV
# ============================================================

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


# Participant-level OOF predictions are PRIVATE
saveRDS(
    tg_oof,
    file.path(
        PRIVATE_OUTPUT_DIR,
        "TwinGene_primary_OOF.rds"
    )
)


saveRDS(
    ukb_oof,
    file.path(
        PRIVATE_OUTPUT_DIR,
        "UKB_primary_OOF.rds"
    )
)


# ============================================================
# Performance
# ============================================================

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


print(
    performance,
    n = Inf,
    width = Inf
)


write_csv(
    performance,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Table3_prediction_performance.csv"
    )
)
