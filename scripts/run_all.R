# ============================================================
# run_all.R
#
# Reproduce revised mortality analyses from authorized
# analysis data.
# ============================================================


source(
    "scripts/01_build_ukb_complete_case.R"
)


source(
    "scripts/02_descriptives_correlations.R"
)


source(
    "scripts/03_primary_cox.R"
)


source(
    "scripts/04_prediction_cv.R"
)


source(
    "scripts/05_prediction_uncertainty.R"
)


source(
    "scripts/06_interactions.R"
)


source(
    "scripts/07_survival_sensitivity.R"
)


source(
    "scripts/08_alternative_ba_prediction.R"
)


source(
    "scripts/09_calibration_dca.R"
)


capture.output(

    sessionInfo(),

    file =
        file.path(
            PUBLIC_OUTPUT_DIR,
            "sessionInfo.txt"
        )
)


cat(
    "\nAll analyses completed.\n"
)
