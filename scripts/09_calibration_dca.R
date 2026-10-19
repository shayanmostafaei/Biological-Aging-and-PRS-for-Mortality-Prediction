# ============================================================
# 09_calibration_dca.R
# Calibration + exploratory decision curve analysis
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


tg_oof <-
    readRDS(
        file.path(
            PRIVATE_OUTPUT_DIR,
            "TwinGene_primary_OOF.rds"
        )
    )


ukb_oof <-
    readRDS(
        file.path(
            PRIVATE_OUTPUT_DIR,
            "UKB_primary_OOF.rds"
        )
    )


km_risk <- function(
        data,
        horizon
) {

    fit <-
        survfit(
            Surv(
                time,
                event
            ) ~ 1,
            data = data
        )


    s <-
        summary(
            fit,
            times =
                horizon,
            extend =
                TRUE
        )


    1 -
    as.numeric(
        s$surv
    )
}


calibration_table <- function(
        data,
        risks,
        cohort,
        model,
        horizon,
        groups = 10
) {

    d <-
        data %>%

        mutate(

            predicted =
                as.numeric(
                    risks[
                        ,
                        as.character(
                            horizon
                        )
                    ]
                ),

            group =
                dplyr::ntile(
                    predicted,
                    groups
                )
        )


    bind_rows(

        lapply(

            sort(
                unique(
                    d$group
                )
            ),

            function(g) {

                x <-
                    d[
                        d$group ==
                            g,
                    ]


                tibble(

                    Cohort =
                        cohort,

                    Model =
                        model,

                    Horizon =
                        horizon,

                    Group =
                        g,

                    N =
                        nrow(x),

                    Predicted =
                        mean(
                            x$predicted
                        ),

                    Observed =
                        km_risk(
                            x,
                            horizon
                        )
                )
            }
        )
    )
}


make_all_calibration <- function(
        data,
        oof,
        cohort
) {

    bind_rows(

        lapply(

            c(
                5,
                10
            ),

            function(horizon)

                bind_rows(

                    calibration_table(
                        data,
                        oof$risk$Conventional,
                        cohort,
                        "Conventional",
                        horizon
                    ),

                    calibration_table(
                        data,
                        oof$risk$Enhanced,
                        cohort,
                        "Enhanced",
                        horizon
                    )
                )
        )
    )
}


calibration <-
    bind_rows(

        make_all_calibration(
            tg,
            tg_oof,
            "TwinGene"
        ),

        make_all_calibration(
            ukb,
            ukb_oof,
            "UKB"
        )
    )


write_csv(
    calibration,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Calibration_5y_10y.csv"
    )
)


# ============================================================
# Calibration plot
# ============================================================

p <-
    ggplot(
        calibration,
        aes(
            x =
                Predicted,
            y =
                Observed,
            linetype =
                factor(
                    Horizon
                ),
            group =
                interaction(
                    Model,
                    Horizon
                )
        )
    ) +

    geom_abline(
        intercept = 0,
        slope = 1,
        linewidth = 0.5
    ) +

    geom_point() +

    geom_line() +

    facet_grid(
        Cohort ~ Model
    ) +

    labs(
        x =
            "Predicted mortality risk",

        y =
            "Observed mortality risk",

        linetype =
            "Horizon (years)"
    ) +

    theme_classic()


ggsave(
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Calibration_5y_10y.png"
    ),
    p,
    width = 9,
    height = 6,
    dpi = 300
)


# ============================================================
# Exploratory DCA
# ============================================================

run_dca <- function(
        data,
        oof,
        cohort,
        horizon
) {

    d <-
        data %>%

        mutate(

            conventional =
                oof$risk$Conventional[
                    ,
                    as.character(
                        horizon
                    )
                ],

            enhanced =
                oof$risk$Enhanced[
                    ,
                    as.character(
                        horizon
                    )
                ]
        )


    thresholds <-
        if (
            horizon == 5
        ) {

            seq(
                0.005,
                0.20,
                by =
                    0.005
            )

        } else {

            seq(
                0.01,
                0.30,
                by =
                    0.01
            )
        }


    fit <-
        dcurves::dca(

            Surv(
                time,
                event
            ) ~
                conventional +
                enhanced,

            data =
                d,

            time =
                horizon,

            thresholds =
                thresholds
        )


    plot_object <-
        plot(
            fit
        )


    ggsave(
        file.path(
            PUBLIC_OUTPUT_DIR,
            paste0(
                "DCA_",
                cohort,
                "_",
                horizon,
                "y.png"
            )
        ),
        plot_object,
        width = 7,
        height = 5,
        dpi = 300
    )


    invisible(
        fit
    )
}


run_dca(
    tg,
    tg_oof,
    "TwinGene",
    5
)


run_dca(
    tg,
    tg_oof,
    "TwinGene",
    10
)


run_dca(
    ukb,
    ukb_oof,
    "UKB",
    5
)


run_dca(
    ukb,
    ukb_oof,
    "UKB",
    10
)
