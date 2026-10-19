# ============================================================
# 00_setup.R
# Shared setup and analysis functions
# ============================================================


# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

required_packages <- c(
    "dplyr",
    "tidyr",
    "tibble",
    "readr",
    "ggplot2",
    "survival",
    "splines",
    "riskRegression",
    "emmeans",
    "BioAge",
    "dcurves"
)


missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (length(missing_packages) > 0) {

    stop(
        paste(
            "Install the following packages before running:",
            paste(
                missing_packages,
                collapse = ", "
            )
        )
    )
}


suppressPackageStartupMessages({

    library(dplyr)
    library(tidyr)
    library(tibble)
    library(readr)
    library(ggplot2)
    library(survival)
    library(splines)
})


# ------------------------------------------------------------
# Local configuration
# ------------------------------------------------------------

if (!file.exists("config.R")) {

    stop(
        paste0(
            "config.R was not found.\n",
            "Copy config.example.R to config.R and edit the local paths."
        )
    )
}


source("config.R")


PUBLIC_OUTPUT_DIR <-
    file.path(
        OUTPUT_DIR,
        "public"
    )


PRIVATE_OUTPUT_DIR <-
    file.path(
        OUTPUT_DIR,
        "private"
    )


dir.create(
    PUBLIC_OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


dir.create(
    PRIVATE_OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# Generic utilities
# ============================================================

assert_columns <- function(
        data,
        required,
        name = "data"
) {

    missing <- setdiff(
        required,
        names(data)
    )

    if (length(missing) > 0) {

        stop(
            paste0(
                name,
                " is missing: ",
                paste(
                    missing,
                    collapse = ", "
                )
            )
        )
    }

    invisible(TRUE)
}


zscore <- function(x) {

    as.numeric(
        scale(x)
    )
}


# ============================================================
# Harmonize TwinGene
# ============================================================

harmonize_twingene <- function(x) {

    required <- c(
        "id",
        "pair_id",
        "age",
        "sex",
        "bmi",
        "ever_smoker",
        "education_hi_int",
        "phenoage_raw",
        "kdm_raw",
        "hd_raw",
        "fi",
        "telomere",
        "time",
        "event"
    )

    assert_columns(
        x,
        required,
        "TwinGene"
    )


    tibble(

        id =
            as.character(
                x$id
            ),

        pair_id =
            as.character(
                x$pair_id
            ),

        age =
            as.numeric(
                x$age
            ),

        sex_model =
            factor(
                x$sex
            ),

        bmi =
            as.numeric(
                x$bmi
            ),

        smoke_model =
            factor(
                x$ever_smoker
            ),

        edu_model =
            factor(
                x$education_hi_int
            ),

        phenoage_raw =
            as.numeric(
                x$phenoage_raw
            ),

        kdm_raw =
            as.numeric(
                x$kdm_raw
            ),

        hd_raw =
            as.numeric(
                x$hd_raw
            ),

        fi =
            as.numeric(
                x$fi
            ),

        telomere =
            as.numeric(
                x$telomere
            ),

        time =
            as.numeric(
                x$time
            ),

        event =
            as.integer(
                x$event
            )
    )
}


# ============================================================
# Harmonize UKB
# ============================================================

harmonize_ukb <- function(x) {

    required <- c(
        "id",
        "age",
        "sex",
        "bmi",
        "smoking",
        "education",
        "phenoage",
        "kdm",
        "hd",
        "fi",
        "telomere",
        "time",
        "event"
    )

    assert_columns(
        x,
        required,
        "UKB"
    )


    x %>%

        transmute(

            id =
                as.character(id),

            age =
                as.numeric(age),

            sex_model =
                factor(
                    ifelse(
                        as.numeric(sex) == 1,
                        "Men",
                        "Women"
                    ),
                    levels =
                        c(
                            "Women",
                            "Men"
                        )
                ),

            bmi =
                as.numeric(bmi),

            smoke_model =
                factor(
                    ifelse(
                        as.numeric(smoking) == 0,
                        "Never",
                        "Ever"
                    ),
                    levels =
                        c(
                            "Never",
                            "Ever"
                        )
                ),

            edu_model =
                factor(
                    ifelse(
                        as.numeric(education) == 3,
                        "Low",
                        "High/Intermediate"
                    ),
                    levels =
                        c(
                            "Low",
                            "High/Intermediate"
                        )
                ),

            phenoage_raw =
                as.numeric(phenoage),

            kdm_raw =
                as.numeric(kdm),

            hd_raw =
                as.numeric(hd),

            fi =
                as.numeric(fi),

            telomere =
                as.numeric(telomere),

            time =
                as.numeric(time),

            event =
                as.integer(event)
        )
}


# ============================================================
# Load both final analysis cohorts
# ============================================================

load_analysis_cohorts <- function() {

    tg_source <-
        readRDS(
            TWIN_GENE_FILE
        )


    ukb_source <-
        readRDS(
            UKB_ANALYSIS_FILE
        )


    tg <-
        harmonize_twingene(
            tg_source
        )


    ukb <-
        harmonize_ukb(
            ukb_source
        )


    stopifnot(
        nrow(tg) == 9617
    )


    stopifnot(
        nrow(ukb) == 314064
    )


    list(
        TwinGene = tg,
        UKB = ukb
    )
}


# ============================================================
# Full-cohort preprocessing for association models
# ============================================================

prepare_full_cohort <- function(data) {

    stopifnot(
        all(
            data$hd_raw > 0
        )
    )


    pheno_fit <-
        lm(
            phenoage_raw ~
                splines::ns(
                    age,
                    df = 3
                ),
            data = data
        )


    kdm_fit <-
        lm(
            kdm_raw ~
                splines::ns(
                    age,
                    df = 3
                ),
            data = data
        )


    data %>%

        mutate(

            pheno_res =
                residuals(
                    pheno_fit
                ),

            kdm_res =
                residuals(
                    kdm_fit
                ),

            hd_log =
                log(
                    hd_raw
                ),

            pheno_z =
                zscore(
                    pheno_res
                ),

            kdm_z =
                zscore(
                    kdm_res
                ),

            hd_z =
                zscore(
                    hd_log
                ),

            fi_z =
                zscore(
                    fi
                ),

            tl_z =
                zscore(
                    telomere
                )
        )
}


# ============================================================
# Fold-safe preprocessing
# ============================================================

standardize_train_test <- function(
        train_x,
        test_x
) {

    mu <-
        mean(
            train_x,
            na.rm = TRUE
        )


    sigma <-
        sd(
            train_x,
            na.rm = TRUE
        )


    if (
        !is.finite(sigma) ||
        sigma <= 0
    ) {

        stop(
            "Invalid training-fold SD."
        )
    }


    list(

        train =
            (train_x - mu) /
            sigma,

        test =
            (test_x - mu) /
            sigma
    )
}


prepare_ba_fold <- function(
        train,
        test
) {

    pheno_fit <-
        lm(
            phenoage_raw ~
                splines::ns(
                    age,
                    df = 3
                ),
            data = train
        )


    train$pheno_res_cv <-
        residuals(
            pheno_fit
        )


    test$pheno_res_cv <-
        test$phenoage_raw -
        predict(
            pheno_fit,
            newdata = test
        )


    kdm_fit <-
        lm(
            kdm_raw ~
                splines::ns(
                    age,
                    df = 3
                ),
            data = train
        )


    train$kdm_res_cv <-
        residuals(
            kdm_fit
        )


    test$kdm_res_cv <-
        test$kdm_raw -
        predict(
            kdm_fit,
            newdata = test
        )


    train$hd_log_cv <-
        log(
            train$hd_raw
        )


    test$hd_log_cv <-
        log(
            test$hd_raw
        )


    variables <- list(

        pheno_z_cv =
            c(
                "pheno_res_cv"
            ),

        kdm_z_cv =
            c(
                "kdm_res_cv"
            ),

        hd_z_cv =
            c(
                "hd_log_cv"
            ),

        fi_z_cv =
            c(
                "fi"
            ),

        tl_z_cv =
            c(
                "telomere"
            )
    )


    for (
        new_name in
        names(
            variables
        )
    ) {

        old_name <-
            variables[
                [new_name]
            ]


        z <-
            standardize_train_test(

                train[
                    [old_name]
                ],

                test[
                    [old_name]
                ]
            )


        train[
            [new_name]
        ] <-
            z$train


        test[
            [new_name]
        ] <-
            z$test
    }


    list(
        train = train,
        test = test
    )
}


# ============================================================
# Cross-validation folds
# ============================================================

make_twingene_folds <- function(
        data,
        K = 10,
        seed = RANDOM_SEED
) {

    set.seed(seed)


    pair_table <-
        data %>%

        group_by(
            pair_id
        ) %>%

        summarise(

            pair_n =
                n(),

            pair_events =
                sum(event),

            .groups =
                "drop"
        )


    strata <-
        interaction(
            pair_table$pair_n,
            pair_table$pair_events,
            drop = TRUE
        )


    pair_table$fold <-
        NA_integer_


    for (
        s in
        levels(
            strata
        )
    ) {

        index <-
            which(
                strata == s
            )


        index <-
            sample(
                index
            )


        pair_table$fold[
            index
        ] <-
            rep(
                seq_len(K),
                length.out =
                    length(index)
            )
    }


    data %>%

        left_join(
            pair_table %>%
                select(
                    pair_id,
                    fold
                ),
            by = "pair_id"
        )
}


make_ukb_folds <- function(
        data,
        K = 10,
        seed = RANDOM_SEED
) {

    set.seed(seed)


    data$fold <-
        NA_integer_


    for (
        event_value in
        c(0, 1)
    ) {

        index <-
            which(
                data$event ==
                    event_value
            )


        index <-
            sample(
                index
            )


        data$fold[
            index
        ] <-
            rep(
                seq_len(K),
                length.out =
                    length(index)
            )
    }


    data
}


# ============================================================
# Absolute-risk helper
# ============================================================

get_h0_at_time <- function(
        fit,
        time_value
) {

    bh <-
        survival::basehaz(
            fit,
            centered = TRUE
        )


    index <-
        which(
            bh$time <=
                time_value
        )


    if (
        length(index) == 0
    ) {

        return(0)
    }


    bh$hazard[
        max(index)
    ]
}


# ============================================================
# Generic CV Cox engine
# ============================================================

cv_cox_models <- function(
        data,
        formulas,
        times = c(5, 10)
) {

    N <-
        nrow(data)


    model_names <-
        names(
            formulas
        )


    LP <-
        setNames(

            lapply(
                model_names,
                function(x)
                    rep(
                        NA_real_,
                        N
                    )
            ),

            model_names
        )


    RISK <-
        setNames(

            lapply(

                model_names,

                function(x)

                    matrix(
                        NA_real_,
                        nrow = N,
                        ncol =
                            length(times),
                        dimnames =
                            list(
                                NULL,
                                as.character(times)
                            )
                    )
            ),

            model_names
        )


    for (
        current_fold in
        sort(
            unique(
                data$fold
            )
        )
    ) {

        message(
            "Running fold ",
            current_fold
        )


        train_index <-
            which(
                data$fold !=
                    current_fold
            )


        test_index <-
            which(
                data$fold ==
                    current_fold
            )


        train <-
            data[
                train_index,
                ,
                drop = FALSE
            ]


        test <-
            data[
                test_index,
                ,
                drop = FALSE
            ]


        prep <-
            prepare_ba_fold(
                train,
                test
            )


        train <-
            prep$train


        test <-
            prep$test


        for (
            model_name in
            model_names
        ) {

            fit <-
                survival::coxph(

                    formula =
                        formulas[
                            [model_name]
                        ],

                    data =
                        train,

                    ties =
                        "efron",

                    x = TRUE,

                    model = TRUE
                )


            lp_test <-
                predict(

                    fit,

                    newdata =
                        test,

                    type =
                        "lp",

                    reference =
                        "sample"
                )


            LP[
                [model_name]
            ][
                test_index
            ] <-
                as.numeric(
                    lp_test
                )


            for (
                j in
                seq_along(times)
            ) {

                h0 <-
                    get_h0_at_time(
                        fit,
                        times[j]
                    )


                RISK[
                    [model_name]
                ][
                    test_index,
                    j
                ] <-

                    1 -
                    exp(
                        -h0 *
                        exp(
                            lp_test
                        )
                    )
            }
        }
    }


    list(
        lp = LP,
        risk = RISK,
        times = times,
        folds = data$fold
    )
}


# ============================================================
# Evaluate OOF predictions
# ============================================================

evaluate_oof <- function(
        data,
        oof,
        cohort
) {

    model_names <-
        names(
            oof$lp
        )


    cindex <-
        sapply(

            model_names,

            function(model_name) {

                fit <-
                    survival::concordance(

                        Surv(
                            data$time,
                            data$event
                        ) ~
                            oof$lp[
                                [model_name]
                            ],

                        reverse = TRUE,

                        timewt = "n"
                    )


                as.numeric(
                    fit$concordance
                )
            }
        )


    score <-
        riskRegression::Score(

            object =
                oof$risk,

            formula =
                Surv(
                    time,
                    event
                ) ~ 1,

            data =
                data,

            metrics =
                c(
                    "AUC",
                    "Brier"
                ),

            times =
                oof$times,

            null.model =
                FALSE,

            se.fit =
                FALSE,

            conf.int =
                FALSE,

            contrasts =
                FALSE,

            cens.method =
                "ipcw",

            verbose = 0
        )


    auc <-
        as.data.frame(
            score$AUC$score
        ) %>%

        transmute(

            model =
                as.character(model),

            time =
                times,

            value =
                AUC
        ) %>%

        pivot_wider(

            names_from =
                time,

            values_from =
                value,

            names_prefix =
                "AUC_"
        )


    brier <-
        as.data.frame(
            score$Brier$score
        ) %>%

        transmute(

            model =
                as.character(model),

            time =
                times,

            value =
                Brier
        ) %>%

        pivot_wider(

            names_from =
                time,

            values_from =
                value,

            names_prefix =
                "Brier_"
        )


    tibble(

        cohort =
            cohort,

        model =
            model_names,

        C_index =
            as.numeric(
                cindex
            )
    ) %>%

        left_join(
            auc,
            by = "model"
        ) %>%

        left_join(
            brier,
            by = "model"
        )
}
