# Predicting All-Cause Mortality Using Genetic Data and Biological Aging Markers (TwinGene discovery + UK Biobank validation)

## Overview

This project evaluates the predictive performance of **chronological age (CA)**, **blood-based biological aging (BA) measures**, **frailty**, **telomere length**, and an aging-related **polygenic risk score (PRS; mvAge)** for **all-cause mortality** in a multicohort design:

- **Discovery cohort:** Swedish **TwinGene** (analytic n = 9,617; median follow-up 16.70 years)
- **External validation cohort:** **UK Biobank (UKB)** (analytic n = 179,504; median follow-up 11.83 years)

The analysis includes:
- Calculation and/or use of BA measures: **PhenoAge**, **KDM**, **homeostatic dysregulation (HD)**, **frailty index (FI)**, and **leukocyte telomere length (TL)**
- Construction/use of an aging-related PRS (**mvAge**) using GWAS summary statistics:  
  https://www.nature.com/articles/s43587-023-00455-5
- **Univariate ROC/AUC** comparisons *(secondary; ignores censoring)*
- **Univariate and multivariable Cox proportional hazards models** *(primary; time-to-death)*
- **Ensemble machine learning** using **SuperLearner** with cross-validated out-of-fold predictions  
  *(10-fold CV; optional repeated CV across multiple random seeds for robustness)*
- **Subgroup analyses** by baseline age group and fixed follow-up horizons *(5/10/15 years)*
- **Interaction analyses** testing BA-by-(BMI, sex, smoking, education) effect modification

---

## Key analytic choice: age-adjusted (“age-gap”) BA residuals

Because PhenoAge, KDM, and HD are strongly correlated with CA, analyses that include CA use **age-adjusted residuals** within each cohort:

- `BA_res = residuals(lm(BA ~ age))`

This improves interpretability by representing deviation from expected biological aging at a given chronological age. Residual versions are used in ROC, Cox, and prediction models wherever CA is included.

---

## Note on contemporary comparator BA measures (e.g., GOLD BioAge)

We considered adding contemporary clinical BA comparators (e.g., **GOLD BioAge**) to reflect newer clinically oriented BA algorithms. However, required inputs (e.g., albumin, ALP, GGT, WBC, RDW, MCV, lymphocytes) are **not available and/or not harmonized** in the current analytic sets for **both** TwinGene and UKB. Therefore, GOLD BioAge is **not computed** in the main pipeline and this limitation is documented in the manuscript.

---

## Directory structure

├── README.md
├── LICENSE.txt
├── Requirements.txt
├── TwinGene_analysis.R             # Main runner (TwinGene discovery + UKB validation)
└── scripts/
├── 00_setup.R                      # Package loading + shared helpers + shared colors
├── 01_prepare_inputs.R             # Read/clean/harmonize analytic datasets
├── 01_BioAge.R                     # Correlation plotting utilities (TwinGene + UKB)
├── 02_models_cox_roc.R             # Univariate ROC + Cox models (time-to-event primary)
├── 03_superlearner_cv.R            # 10-fold CV SuperLearner (+ optional repeats)
├── 04_subgroup_interactions.R      # Subgroup analyses (age x follow-up) + interaction maps
└── 05_figures_supp.R               # Supplementary figures (S3 densities; S1/S2 heatmaps)

---

> Cohort data files are **not** included in the repository due to access restrictions.

---

## Inputs (not included)

You must provide local analytic datasets (paths can be edited in `TwinGene_analysis.R`), e.g.:

- `data/TwinGene_analytic.csv`
- `data/UKB_analytic.csv`

Scripts assume the datasets include outcome/time variables and harmonized predictor/covariate columns (see comments inside `scripts/01_prepare_inputs.R`).

---

## How to run

From R:

```r
source("TwinGene_analysis.R")

Key outputs are written to the configured output directory (default in runner), including:
	•	Univariate ROC AUC tables
	•	Cox model results tables
	•	SuperLearner CV results (all seeds + summary mean±SD when repeats enabled)
	•	Supplementary figures

A sessionInfo.txt file is saved for reproducibility.

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
