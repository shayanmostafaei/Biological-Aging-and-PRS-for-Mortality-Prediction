# ============================================================
# 02_descriptives_correlations.R
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
# Descriptive summary
# ============================================================

summary_one <- function(
        data,
        cohort
) {

    fu <-
        quantile(
            data$time,
            c(
                0.25,
                0.50,
                0.75
            )
        )


    tibble(

        Cohort =
            cohort,

        N =
            nrow(data),

        Deaths =
            sum(
                data$event
            ),

        Death_percent =
            100 *
            mean(
                data$event
            ),

        Age_mean =
            mean(
                data$age
            ),

        Age_SD =
            sd(
                data$age
            ),

        BMI_mean =
            mean(
                data$bmi
            ),

        BMI_SD =
            sd(
                data$bmi
            ),

        PhenoAge_mean =
            mean(
                data$phenoage_raw
            ),

        PhenoAge_SD =
            sd(
                data$phenoage_raw
            ),

        KDM_mean =
            mean(
                data$kdm_raw
            ),

        KDM_SD =
            sd(
                data$kdm_raw
            ),

        HD_mean =
            mean(
                data$hd_raw
            ),

        HD_SD =
            sd(
                data$hd_raw
            ),

        FI_mean =
            mean(
                data$fi
            ),

        FI_SD =
            sd(
                data$fi
            ),

        TL_mean =
            mean(
                data$telomere
            ),

        TL_SD =
            sd(
                data$telomere
            ),

        Followup_Q1 =
            fu[1],

        Followup_median =
            fu[2],

        Followup_Q3 =
            fu[3]
    )
}


table1 <-
    bind_rows(

        summary_one(
            tg,
            "TwinGene"
        ),

        summary_one(
            ukb,
            "UKB"
        )
    )


write_csv(
    table1,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Table1_descriptive_summary.csv"
    )
)


print(
    table1,
    width = Inf
)


# ============================================================
# Raw correlations
# ============================================================

cor_vars <- c(
    "age",
    "phenoage_raw",
    "kdm_raw",
    "hd_raw",
    "fi",
    "telomere"
)


cor_tg <-
    cor(
        tg[
            cor_vars
        ],
        use =
            "pairwise.complete.obs",
        method =
            "pearson"
    )


cor_ukb <-
    cor(
        ukb[
            cor_vars
        ],
        use =
            "pairwise.complete.obs",
        method =
            "pearson"
    )


write_csv(
    as.data.frame(
        cor_tg
    ) %>%
        rownames_to_column(
            "Variable"
        ),
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Figure2_TwinGene_correlations.csv"
    )
)


write_csv(
    as.data.frame(
        cor_ukb
    ) %>%
        rownames_to_column(
            "Variable"
        ),
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Figure2_UKB_correlations.csv"
    )
)


# ============================================================
# BioAge-style correlation figures
# ============================================================

agevar <-
    cor_vars


axis_type <-
    setNames(
        rep(
            "float",
            length(agevar)
        ),
        agevar
    )


labels <- c(

    age =
        "Chronological\nage",

    phenoage_raw =
        "PhenoAge",

    kdm_raw =
        "KDM\nbiological age",

    hd_raw =
        "Homeostatic\ndysregulation",

    fi =
        "Frailty\nindex",

    telomere =
        "Telomere\nlength"
)


png(
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Figure2A_TwinGene_correlations.png"
    ),
    width = 2400,
    height = 1800,
    res = 300
)

BioAge::plot_baa(
    tg,
    agevar,
    labels,
    axis_type
)

dev.off()


png(
    file.path(
        PUBLIC_OUTPUT_DIR,
        "Figure2B_UKB_correlations.png"
    ),
    width = 2400,
    height = 1800,
    res = 300
)

BioAge::plot_baa(
    ukb,
    agevar,
    labels,
    axis_type
)

dev.off()
