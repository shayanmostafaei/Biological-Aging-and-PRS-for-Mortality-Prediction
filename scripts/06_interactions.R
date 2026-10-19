# ============================================================
# 06_interactions.R
# Prespecified BA x modifier interactions
# ============================================================

source(
    "scripts/00_setup.R"
)


cohorts <-
    load_analysis_cohorts()


tg <-
    prepare_full_cohort(
        cohorts$TwinGene
    )


ukb <-
    prepare_full_cohort(
        cohorts$UKB
    )


measures <- c(

    PhenoAge =
        "pheno_z",

    KDM =
        "kdm_z",

    HD =
        "hd_z",

    FI =
        "fi_z",

    TL =
        "tl_z"
)


modifiers <- c(

    Sex =
        "sex_model",

    BMI =
        "bmi",

    Smoking =
        "smoke_model",

    Education =
        "edu_model"
)


interaction_test <- function(
        data,
        ba_var,
        modifier,
        cohort,
        cluster_pair = FALSE
) {

    d <-
        data


    d$x <-
        d[
            [ba_var]
        ]


    interaction_term <-
        paste0(
            "x:",
            modifier
        )


    formula <-
        as.formula(

            paste0(

                "Surv(time,event) ~ ",
                "splines::ns(age,4) + ",
                "sex_model + ",
                "splines::ns(bmi,3) + ",
                "smoke_model + ",
                "edu_model + ",
                "x + ",
                interaction_term
            )
        )


    if (
        cluster_pair
    ) {

        fit <-
            coxph(
                formula,
                data = d,
                cluster =
                    d$pair_id,
                robust = TRUE,
                ties = "efron"
            )

    } else {

        fit <-
            coxph(
                formula,
                data = d,
                ties = "efron"
            )
    }


    beta_names <-
        names(
            coef(
                fit
            )
        )


    index <-
        grep(
            "^x:|:x",
            beta_names
        )


    beta <-
        coef(
            fit
        )[
            index
        ]


    V <-
        vcov(
            fit
        )[
            index,
            index,
            drop = FALSE
        ]


    statistic <-
        as.numeric(
            t(beta) %*%
                solve(V) %*%
                beta
        )


    p <-
        pchisq(
            statistic,
            df =
                length(index),
            lower.tail =
                FALSE
        )


    tibble(

        Cohort =
            cohort,

        Measure =
            names(
                measures
            )[
                measures ==
                    ba_var
            ][1],

        Modifier =
            names(
                modifiers
            )[
                modifiers ==
                    modifier
            ][1],

        df =
            length(index),

        Wald =
            statistic,

        P =
            p
    )
}


run_grid <- function(
        data,
        cohort,
        cluster_pair
) {

    bind_rows(

        lapply(

            names(measures),

            function(m) {

                bind_rows(

                    lapply(

                        names(modifiers),

                        function(z)

                            interaction_test(

                                data,

                                measures[
                                    [m]
                                ],

                                modifiers[
                                    [z]
                                ],

                                cohort,

                                cluster_pair
                            )
                    )
                )
            }
        )
    )
}


interaction_results <-
    bind_rows(

        run_grid(
            tg,
            "TwinGene",
            TRUE
        ),

        run_grid(
            ukb,
            "UKB",
            FALSE
        )
    ) %>%

    group_by(
        Cohort
    ) %>%

    mutate(
        P_FDR =
            p.adjust(
                P,
                method =
                    "BH"
            )
    ) %>%

    ungroup()


print(
    interaction_results,
    n = Inf
)


write_csv(
    interaction_results,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Interaction_tests_20_per_cohort.csv"
    )
)
