# ============================================================
# 05_prediction_uncertainty.R
# Paired uncertainty for incremental prediction performance
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


# ============================================================
# UKB: paired C-index covariance
# ============================================================

ukb_c <- data.frame(

    time =
        ukb$time,

    event =
        ukb$event,

    Conventional =
        ukb_oof$lp$Conventional,

    Enhanced =
        ukb_oof$lp$Enhanced
)


cc <-
    survival::concordance(

        Surv(
            time,
            event
        ) ~
            Conventional +
            Enhanced,

        data =
            ukb_c,

        reverse =
            TRUE,

        timewt =
            "n"
    )


c_est <-
    as.numeric(
        coef(
            cc
        )
    )


names(
    c_est
) <-
    c(
        "Conventional",
        "Enhanced"
    )


V <-
    vcov(
        cc
    )


delta_c <-
    c_est[
        "Enhanced"
    ] -
    c_est[
        "Conventional"
    ]


se_delta_c <-
    sqrt(
        V[2, 2] +
        V[1, 1] -
        2 *
        V[1, 2]
    )


ukb_c_result <-
    tibble(

        Metric =
            "C-index",

        Delta =
            delta_c,

        SE =
            se_delta_c,

        Lower95 =
            delta_c -
            1.96 *
            se_delta_c,

        Upper95 =
            delta_c +
            1.96 *
            se_delta_c
    )


# ============================================================
# UKB: paired AUC/Brier inference
# ============================================================

ukb_score <-
    riskRegression::Score(

        object = list(

            Conventional =
                ukb_oof$risk$Conventional,

            Enhanced =
                ukb_oof$risk$Enhanced
        ),

        formula =
            Surv(
                time,
                event
            ) ~ 1,

        data =
            ukb,

        metrics =
            c(
                "AUC",
                "Brier"
            ),

        times =
            c(
                5,
                10
            ),

        null.model =
            FALSE,

        se.fit =
            TRUE,

        conf.int =
            0.95,

        contrasts =
            TRUE,

        cens.method =
            "ipcw",

        verbose =
            0
    )


write_csv(
    ukb_c_result,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "UKB_paired_Cindex_difference.csv"
    )
)


write_csv(
    as.data.frame(
        ukb_score$AUC$contrasts
    ),
    file.path(
        PUBLIC_OUTPUT_DIR,
        "UKB_paired_AUC_differences.csv"
    )
)


write_csv(
    as.data.frame(
        ukb_score$Brier$contrasts
    ),
    file.path(
        PUBLIC_OUTPUT_DIR,
        "UKB_paired_Brier_differences.csv"
    )
)


# ============================================================
# TwinGene pair-cluster bootstrap
# ============================================================

metric_difference <- function(
        data,
        oof,
        index
) {

    d <-
        data[
            index,
            ,
            drop = FALSE
        ]


    conv_lp <-
        oof$lp$Conventional[
            index
        ]


    enh_lp <-
        oof$lp$Enhanced[
            index
        ]


    c_conv <-
        survival::concordance(

            Surv(
                d$time,
                d$event
            ) ~
                conv_lp,

            reverse =
                TRUE,

            timewt =
                "n"
        )$concordance


    c_enh <-
        survival::concordance(

            Surv(
                d$time,
                d$event
            ) ~
                enh_lp,

            reverse =
                TRUE,

            timewt =
                "n"
        )$concordance


    risk_list <- list(

        Conventional =
            oof$risk$Conventional[
                index,
                ,
                drop = FALSE
            ],

        Enhanced =
            oof$risk$Enhanced[
                index,
                ,
                drop = FALSE
            ]
    )


    score <-
        riskRegression::Score(

            object =
                risk_list,

            formula =
                Surv(
                    time,
                    event
                ) ~ 1,

            data =
                d,

            metrics =
                c(
                    "AUC",
                    "Brier"
                ),

            times =
                c(
                    5,
                    10
                ),

            null.model =
                FALSE,

            se.fit =
                FALSE,

            contrasts =
                FALSE,

            cens.method =
                "ipcw",

            verbose =
                0
        )


    auc <-
        as.data.frame(
            score$AUC$score
        )


    brier <-
        as.data.frame(
            score$Brier$score
        )


    get_delta <- function(
            tab,
            value_name,
            time_value
    ) {

        e <-
            tab[
                tab$model ==
                    "Enhanced" &
                tab$times ==
                    time_value,
                value_name
            ]


        c <-
            tab[
                tab$model ==
                    "Conventional" &
                tab$times ==
                    time_value,
                value_name
            ]


        e - c
    }


    c(

        Delta_C =
            c_enh -
            c_conv,

        Delta_AUC5 =
            get_delta(
                auc,
                "AUC",
                5
            ),

        Delta_AUC10 =
            get_delta(
                auc,
                "AUC",
                10
            ),

        Delta_Brier5 =
            get_delta(
                brier,
                "Brier",
                5
            ),

        Delta_Brier10 =
            get_delta(
                brier,
                "Brier",
                10
            )
    )
}


set.seed(
    RANDOM_SEED
)


pairs <-
    unique(
        tg$pair_id
    )


B <-
    1000


boot_results <-
    matrix(
        NA_real_,
        nrow = B,
        ncol = 5
    )


colnames(
    boot_results
) <-
    c(
        "Delta_C",
        "Delta_AUC5",
        "Delta_AUC10",
        "Delta_Brier5",
        "Delta_Brier10"
    )


for (
    b in
    seq_len(B)
) {

    sampled_pairs <-
        sample(
            pairs,
            length(pairs),
            replace = TRUE
        )


    index <-
        unlist(

            lapply(

                sampled_pairs,

                function(p)
                    which(
                        tg$pair_id ==
                            p
                    )
            ),

            use.names =
                FALSE
        )


    boot_results[
        b,
    ] <-
        metric_difference(
            tg,
            tg_oof,
            index
        )


    if (
        b %% 50 == 0
    ) {

        message(
            "Bootstrap ",
            b,
            "/",
            B
        )
    }
}


tg_uncertainty <-
    tibble(

        Metric =
            colnames(
                boot_results
            ),

        Estimate =
            apply(
                boot_results,
                2,
                mean,
                na.rm = TRUE
            ),

        Bootstrap_SE =
            apply(
                boot_results,
                2,
                sd,
                na.rm = TRUE
            ),

        Lower95 =
            apply(
                boot_results,
                2,
                quantile,
                probs =
                    0.025,
                na.rm =
                    TRUE
            ),

        Upper95 =
            apply(
                boot_results,
                2,
                quantile,
                probs =
                    0.975,
                na.rm =
                    TRUE
            )
    )


write_csv(
    tg_uncertainty,
    file.path(
        PUBLIC_OUTPUT_DIR,
        "TwinGene_pair_bootstrap_prediction_differences.csv"
    )
)
