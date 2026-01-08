# Predicting All-Cause Mortality Using Genetic Data and Biological Aging Markers (TwinGene discovery + UK Biobank validation)

## Overview

This project evaluates the predictive performance of **chronological age (CA)**, **blood-based biological aging (BA) measures**, **frailty**, **telomere length**, and an aging-related **polygenic risk score (PRS; mvAge)** for **all-cause mortality** in a multicohort design:

- **Discovery cohort:** Swedish **TwinGene** (analytic n = 9,617; median follow-up 16.70 years)
- **External validation cohort:** **UK Biobank (UKB)** (analytic n = 179,504; median follow-up 11.83 years)

The analysis includes:
- Calculation and/or use of BA measures: **PhenoAge**, **KDM**, **homeostatic dysregulation (HD)**, **frailty index (FI)**, and **leukocyte telomere length (TL)**
- Construction/use of an aging-related PRS (**mvAge**) using GWAS summary statistics  
  https://www.nature.com/articles/s43587-023-00455-5
- **Univariate ROC/AUC** comparisons (secondary; ignores censoring)
- **Univariate and multivariable Cox proportional hazards models** for time-to-death (primary)
- **Ensemble machine learning** using **SuperLearner** with cross-validated out-of-fold predictions
- **Subgroup analyses** by baseline age group and fixed follow-up horizons (5/10/15 years)
- **Interaction analyses** testing BA-by-(BMI, sex, smoking, education) effect modification

### Key analytic choice: age-adjusted (“age-gap”) BA residuals
Because PhenoAge, KDM, and HD are strongly correlated with CA, analyses that include CA use **age-adjusted residuals** within each cohort:
- `BA_res = residuals(lm(BA ~ age))`
This improves interpretability by representing deviation from expected BA at a given CA.

### Note on contemporary comparator BA measures (e.g., GOLD BioAge)
We considered adding contemporary clinical BA comparators (e.g., **GOLD BioAge**), but required inputs (e.g., albumin, ALP, GGT, WBC, RDW, MCV, lymphocytes) are not available/harmonized in the current analytic sets for TwinGene and UKB; therefore, GOLD BioAge is not computed in the main pipeline.

---

## Directory Structure

├── README.md
├── LICENSE
├── Requirements.txt
├── TwinGene_analysis.R           # Main runner script (multicohort pipeline)
└── scripts/
├── 00_setup.R                    # Package loading + shared helpers
├── 01_prepare_inputs.R           # Input schemas + preprocessing helpers
├── 02_models_cox_roc.R           # Univariate ROC + Cox models
├── 03_superlearner_cv.R          # 10-fold CV SuperLearner (+ optional repeats)
├── 04_subgroup_interactions.R    # Subgroup scaffolding + interaction hooks
└── 05_figures_supp.R             # Supplementary figures (e.g., biomarker densities)

> Cohort data files are **not** included in the repository due to access restrictions.

---

## Requirements

See `Requirements.txt` for the list of R packages used in the analysis.

---

## License

This project is licensed under the MIT License — see the `LICENSE` file for details.

---

## Contact

For questions or contributions, please contact:
- Dr. Shayan Mostafaei (shayan.mostafaei@ki.se)
- Dr. Sara Hägg (sara.hagg@ki.se)
