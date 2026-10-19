# ============================================================
# 01_build_ukb_complete_case.R
#
# Reconstruct phenotype-only UKB complete-case sample
# ============================================================

source(
    "scripts/00_setup.R"
)


bio <-
    readRDS(
        UKB_BIOAGE_FILE
    )


frailty <-
    readRDS(
        UKB_FRAILTY_FILE
    )


death <-
    readRDS(
        UKB_DEATH_FILE
    )


telomere <-
    readRDS(
        UKB_TELOMERE_FILE
    )


detect_id <- function(
        x,
        object_name
) {

    candidates <-
        intersect(
            c(
                "eid",
                "f.eid",
                "id"
            ),
            names(x)
        )


    if (
        length(candidates) != 1
    ) {

        stop(
            paste(
                "Unable to identify one ID field in",
                object_name
            )
        )
    }


    candidates[1]
}


bio_id <-
    detect_id(
        bio,
        "UKB_BioAge"
    )


frailty_id <-
    detect_id(
        frailty,
        "frailty"
    )


death_id <-
    detect_id(
        death,
        "mortality"
    )


tl_id <-
    detect_id(
        telomere,
        "telomere"
    )


tl_field <-
    "f.22191.0.0"


stopifnot(
    tl_field %in%
        names(
            telomere
        )
)


bioage_inputs <- c(

    "waist",
    "fev_1000",
    "albumin_gL",
    "alp",
    "bun",
    "creat_umol",
    "glucose_mmol",
    "uap",
    "lymph",
    "mcv",
    "rbc",
    "rdw",
    "crp",
    "dbp",
    "sbp",
    "hba1c",
    "trig",
    "totchol"
)


assert_columns(
    bio,
    c(
        bioage_inputs,
        "phenoage",
        "kdm",
        "hd"
    ),
    "UKB_BioAge"
)


bio <-
    bio %>%

    mutate(

        id =
            sprintf(
                "%.0f",
                as.numeric(
                    .data[
                        [bio_id]
                    ]
                )
            )
    )


bio_matrix <-
    as.matrix(
        bio[
            ,
            bioage_inputs,
            drop = FALSE
        ]
    )


bio$inputs_complete <-

    complete.cases(
        bio_matrix
    ) &

    rowSums(
        !is.finite(
            bio_matrix
        )
    ) == 0


bio$ba_complete <-

    is.finite(
        bio$phenoage
    ) &

    is.finite(
        bio$kdm
    ) &

    is.finite(
        bio$hd
    )


eligible_ba <-

    bio$inputs_complete &

    bio$ba_complete


bio$pheno_qc_z <-
    NA_real_


bio$kdm_qc_z <-
    NA_real_


bio$pheno_qc_z[
    eligible_ba
] <-

    as.numeric(
        scale(
            bio$phenoage[
                eligible_ba
            ]
        )
    )


bio$kdm_qc_z[
    eligible_ba
] <-

    as.numeric(
        scale(
            bio$kdm[
                eligible_ba
            ]
        )
    )


bio$extreme_ba <-

    eligible_ba &

    (
        abs(
            bio$pheno_qc_z
        ) > 5 |

        abs(
            bio$kdm_qc_z
        ) > 5
    )


bio$bioage_ok <-

    eligible_ba &

    !bio$extreme_ba


frailty2 <-
    frailty %>%

    transmute(

        id =
            sprintf(
                "%.0f",
                as.numeric(
                    .data[
                        [frailty_id]
                    ]
                )
            ),

        age =
            as.numeric(age),

        sex =
            as.numeric(sex),

        bmi =
            as.numeric(bmi),

        smoking =
            as.numeric(smoking),

        education =
            as.numeric(education),

        fi =
            as.numeric(fi)
    )


death2 <-
    death %>%

    transmute(

        id =
            sprintf(
                "%.0f",
                as.numeric(
                    .data[
                        [death_id]
                    ]
                )
            ),

        date_baseline =
            as.Date(
                date_baseline
            ),

        date_death =
            as.Date(
                date_death
            )
    )


tl2 <-
    telomere %>%

    transmute(

        id =
            sprintf(
                "%.0f",
                as.numeric(
                    .data[
                        [tl_id]
                    ]
                )
            ),

        telomere =
            as.numeric(
                .data[
                    [tl_field]
                ]
            )
    )


ukb <-
    bio %>%

    left_join(
        frailty2,
        by = "id"
    ) %>%

    left_join(
        tl2,
        by = "id"
    ) %>%

    left_join(
        death2,
        by = "id"
    )


admin_censor <-
    as.Date(
        "2022-12-30"
    )


ukb <-
    ukb %>%

    mutate(

        mortality_ok =

            !is.na(
                date_baseline
            ) &

            date_baseline <=
                admin_censor &

            (
                is.na(
                    date_death
                ) |

                date_death >=
                    date_baseline
            ),

        fi_ok =

            is.finite(fi) &

            fi >= 0 &

            fi <= 1,

        tl_ok =

            is.finite(
                telomere
            ) &

            telomere > 0,

        covariates_ok =

            is.finite(age) &

            sex %in%
                c(0, 1) &

            is.finite(bmi) &

            smoking %in%
                c(0, 1, 2) &

            education %in%
                c(1, 2, 3),

        complete_case =

            bioage_ok &

            mortality_ok &

            fi_ok &

            tl_ok &

            covariates_ok
    )


N1 <-
    nrow(ukb)


N2 <-
    sum(
        ukb$inputs_complete &
        ukb$ba_complete
    )


N3 <-
    sum(
        ukb$bioage_ok
    )


N4 <-
    sum(
        ukb$bioage_ok &
        ukb$mortality_ok
    )


N5 <-
    sum(
        ukb$bioage_ok &
        ukb$mortality_ok &
        ukb$fi_ok
    )


N6 <-
    sum(
        ukb$bioage_ok &
        ukb$mortality_ok &
        ukb$fi_ok &
        ukb$tl_ok
    )


N7 <-
    sum(
        ukb$complete_case
    )


sample_flow <-
    tibble(

        stage = c(

            "UKB BioAge source",

            "Complete biomarker inputs and BA estimates",

            "After >5-SD age-scale BA QC",

            "Plus valid mortality follow-up",

            "Plus complete FI",

            "Plus complete TL",

            "Plus complete age, sex, BMI, smoking and education"
        ),

        N =
            c(
                N1,
                N2,
                N3,
                N4,
                N5,
                N6,
                N7
            )
    ) %>%

    mutate(

        excluded_from_previous =
            c(
                NA_integer_,
                N[
                    -length(N)
                ] -
                N[
                    -1
                ]
            )
    )


print(
    sample_flow,
    n = Inf
)


stopifnot(
    N7 == 314064
)


ukb_final <-
    ukb %>%

    filter(
        complete_case
    ) %>%

    mutate(

        event =
            as.integer(

                !is.na(
                    date_death
                ) &

                date_death <=
                    admin_censor
            ),

        end_date =
            if_else(
                event == 1L,
                date_death,
                admin_censor
            ),

        time =
            as.numeric(
                end_date -
                date_baseline
            ) /
            365.25
    )


stopifnot(
    sum(
        ukb_final$event
    ) ==
    25084
)


observed_followup <-
    quantile(
        ukb_final$time,
        c(
            0.25,
            0.50,
            0.75
        )
    )


cat(
    sprintf(
        "\nFinal UKB N = %s\nDeaths = %s\nObserved follow-up = %.2f [%.2f-%.2f]\n",
        format(
            nrow(
                ukb_final
            ),
            big.mark = ","
        ),
        format(
            sum(
                ukb_final$event
            ),
            big.mark = ","
        ),
        observed_followup[2],
        observed_followup[1],
        observed_followup[3]
    )
)


write_csv(
    sample_flow,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "UKB_sample_flow.csv"
    )
)


# Participant-level file: PRIVATE only
saveRDS(
    ukb_final,
    UKB_ANALYSIS_FILE
)
