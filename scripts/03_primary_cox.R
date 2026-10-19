# ============================================================
# 03_primary_cox.R
# Primary five biological-aging mortality associations
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


run_one <- function(
        data,
        variable,
        cohort,
        cluster_pair = FALSE
) {

    d <-
        data


    d$x <-
        d[
            [variable]
        ]


    formula <-
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
            edu_model +
            x


    if (
        cluster_pair
    ) {

        fit <-
            survival::coxph(
                formula,
                data = d,
                cluster =
                    d$pair_id,
                robust = TRUE,
                ties = "efron"
            )

    } else {

        fit <-
            survival::coxph(
                formula,
                data = d,
                ties = "efron"
            )
    }


    s <-
        summary(
            fit
        )


    row <-
        which(
            rownames(
                s$coefficients
            ) ==
                "x"
        )


    tibble(

        Cohort =
            cohort,

        Measure =
            names(
                measures
            )[
                measures ==
                    variable
            ][1],

        HR =
            s$coefficients[
                row,
                "exp(coef)"
            ],

        Lower95 =
            s$conf.int[
                row,
                "lower .95"
            ],

        Upper95 =
            s$conf.int[
                row,
                "upper .95"
            ],

        P =
            s$coefficients[
                row,
                "Pr(>|z|)"
            ]
    )
}


results <-
    bind_rows(

        lapply(

            measures,

            function(x)
                run_one(
                    tg,
                    x,
                    "TwinGene",
                    TRUE
                )
        ),

        lapply(

            measures,

            function(x)
                run_one(
                    ukb,
                    x,
                    "UKB",
                    FALSE
                )
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
    results,
    n = Inf
)


write_csv(
    results,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Table2_primary_Cox_models.csv"
    )
)
