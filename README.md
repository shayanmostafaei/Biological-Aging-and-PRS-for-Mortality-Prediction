# Biological aging measures for all-cause mortality prediction in two population-based cohorts

Analysis code accompanying the manuscript:

**Mostafaei S, Lindh J, Qin C, Mak JKL, Hägg S. Biological aging measures for all-cause mortality prediction in two population-based cohorts.**

**Manuscript status:** Under review at *npj Aging*.

## Overview

This study evaluates five complementary biological aging measures in relation to all-cause mortality:

- PhenoAge
- Klemera–Doubal method (KDM)
- Homeostatic dysregulation (HD)
- Frailty index (FI)
- Leukocyte telomere length (TL)

The same analytical framework was applied independently in TwinGene and UK Biobank (UKB).

The repository retains its historical name, which includes “PRS.” Polygenic risk scores are **not included in the revised manuscript analyses**. Earlier PRS-based analyses are preserved only under `legacy/` for provenance and should not be used to reproduce the revised manuscript.

## Study cohorts

| Characteristic | TwinGene | UKB |
|---|---:|---:|
| Analytic sample | 9,617 | 314,064 |
| Deaths | 3,090 | 25,084 |
| Mortality proportion | 32.1% | 8.0% |
| Median observed follow-up, years | 16.76 | 13.79 |

The primary outcome was all-cause mortality.

Time since baseline was used as the primary survival timescale.

TwinGene mortality follow-up extended through 1 February 2024.

For UKB, follow-up was calculated from baseline assessment to death or administrative censoring on 30 December 2022, whichever occurred first.

## UKB phenotype-only analytic sample

The revised UKB sample was reconstructed independently of genotype, ancestry, PRS, and unrelatedness criteria.

Starting from the biomarker-based aging dataset, eligibility required:

1. valid biomarker inputs and PhenoAge, KDM, and HD estimates;
2. removal of extreme age-scale biological-age values (>5 SD for PhenoAge or KDM);
3. valid mortality follow-up;
4. complete frailty index;
5. complete leukocyte telomere length;
6. complete age, sex, BMI, smoking, and education.

The resulting phenotype-only complete-case sample contained 314,064 participants.

The exact reconstruction is implemented in:

`scripts/01_build_ukb_complete_case.R`

## Biological aging measures

The clinical biomarker-based aging measures used cohort-specific implementations.

TwinGene and UKB differed in the number and availability of clinical biomarkers, with eight clinical biomarkers used in the TwinGene implementation and 18 in the UKB implementation.

For descriptive analyses and correlations, the raw cohort-specific measures were used.

For mortality association analyses:

- PhenoAge was residualized on chronological age using a natural spline with 3 degrees of freedom.
- KDM was residualized on chronological age using the same approach.
- HD was natural-log transformed.
- PhenoAge residual, KDM residual, log-HD, FI, and TL were standardized within cohort.

For prediction analyses, residualization and standardization were estimated within each training fold and then applied to the held-out fold to avoid information leakage.

## Primary mortality association analyses

Each biological aging measure was evaluated separately in an adjusted Cox proportional hazards model.

The biological aging measures were **not mutually adjusted** in the primary association models.

All models included:

- chronological age, modeled using a natural spline with 4 degrees of freedom;
- sex;
- BMI, modeled using a natural spline with 3 degrees of freedom;
- smoking;
- education.

TwinGene models additionally used twin-pair-clustered robust variance estimation.

Hazard ratios were expressed per 1-SD higher biological aging measure.

Benjamini-Hochberg false-discovery-rate correction was applied separately within each cohort across the five primary biological aging tests.

## Correlation analyses

Pairwise Pearson correlations were calculated among:

- chronological age;
- raw PhenoAge;
- raw KDM;
- raw HD;
- FI;
- TL.

The correlation analyses are descriptive and are not used for model selection.

## Cross-validated mortality prediction

Three nested Cox prediction models were evaluated independently within each cohort:

1. **Age model:** chronological age only.
2. **Conventional model:** chronological age, sex, BMI, smoking, and education.
3. **Biological-aging-enhanced model:** conventional predictors plus PhenoAge residual, FI, and TL.

KDM residual and HD were not entered simultaneously with PhenoAge in the primary enhanced model because the biomarker-derived aging measures contain overlapping information.

Prediction performance was estimated from a fixed 10-fold cross-validation partition.

The same folds were used for all compared prediction models.

In TwinGene, co-twins were retained within the same cross-validation fold.

For UKB, folds were balanced with respect to mortality status.

All transformations required for prediction—including PhenoAge/KDM residualization and standardization—were estimated in the training data and applied to held-out participants.

Performance metrics included:

- Harrell's C-index;
- censoring-aware 5-year AUC;
- censoring-aware 10-year AUC;
- 5-year Brier score;
- 10-year Brier score.

Participant-level out-of-fold predictions are not distributed because they derive from restricted cohort data.

## Prediction uncertainty

Incremental performance of the biological-aging-enhanced model relative to the conventional model was evaluated using paired out-of-fold predictions.

For TwinGene, uncertainty was estimated using 1,000 twin-pair cluster-bootstrap replicates, retaining co-twins together during resampling.

For UKB:

- paired C-index differences were evaluated using the covariance of concordance estimates;
- AUC and Brier-score differences were evaluated using paired influence-function inference.

Uncertainty therefore reflects performance conditional on the fixed cross-validation partitions.

## Alternative biomarker-based aging models

As a sensitivity analysis, the primary enhanced model was compared with otherwise identical models replacing PhenoAge residual with:

- KDM residual; or
- log-transformed HD.

FI and TL were retained in every enhanced model.

The same cross-validation folds and fold-specific preprocessing procedures were used for all models.

These analyses were intended to assess robustness to the choice of biomarker-derived aging measure rather than to perform formal data-driven model selection.

In UKB, the PhenoAge-based enhanced model showed the strongest overall discrimination among the three alternatives. In TwinGene, performance differences between the PhenoAge- and HD-based models were small, supporting the interpretation that the biomarker-derived aging measures contain overlapping prognostic information.

## Effect modification

Interactions between each of the five biological aging measures and the following modifiers were evaluated:

- sex;
- BMI;
- smoking;
- education.

This produced 20 interaction tests per cohort.

One-degree-of-freedom Wald tests were used for the prespecified interaction terms and Benjamini-Hochberg correction was applied separately within each cohort across the 20 tests.

TwinGene interaction models retained twin-pair-clustered robust variance estimation.

## Diagnostics and sensitivity analyses

Sensitivity analyses included:

- Schoenfeld residual tests for proportional hazards;
- interval-specific associations during 0-5, 5-10, and >10 years;
- a 2-year landmark analysis;
- attained age as the survival timescale;
- alternative PhenoAge/KDM/HD prediction models;
- calibration of out-of-fold absolute mortality risk;
- exploratory decision-curve analysis.

Decision-curve analysis was considered exploratory because no clinical intervention or action threshold was prespecified.

## Repository structure

`config.example.R`  
Template containing paths to restricted cohort data.

`scripts/00_setup.R`  
Shared functions, data harmonization, preprocessing, cross-validation, and performance utilities.

`scripts/01_build_ukb_complete_case.R`  
Reconstructs the phenotype-only UKB complete-case cohort and corrected survival follow-up.

`scripts/02_descriptives_correlations.R`  
Produces descriptive statistics and raw biological-aging correlation analyses.

`scripts/03_primary_cox.R`  
Runs the five primary adjusted biological-aging mortality models.

`scripts/04_prediction_cv.R`  
Runs the fixed 10-fold Cox prediction analysis.

`scripts/05_prediction_uncertainty.R`  
Calculates paired uncertainty for incremental prediction performance.

`scripts/06_interactions.R`  
Runs prespecified biological-aging effect-modification analyses.

`scripts/07_survival_sensitivity.R`  
Runs PH diagnostics, interval-specific models, landmark analyses, and attained-age analyses.

`scripts/08_alternative_ba_prediction.R`  
Compares PhenoAge-, KDM-, and HD-based enhanced prediction models.

`scripts/09_calibration_dca.R`  
Produces calibration summaries and exploratory decision-curve analyses from out-of-fold absolute-risk predictions.

`legacy/`  
Historical analyses from an earlier PRS-based version of the project. These scripts do not reproduce the revised manuscript.

## Software

Analyses were conducted in R 4.5.1.

Principal packages include:

- `survival`
- `riskRegression`
- `emmeans`
- `BioAge`
- `dcurves`
- `ggplot2`
- `dplyr`
- `tidyr`
- `readr`

See `Requirements.txt` for the analysis package inventory.

## Reproducibility

Participant-level cohort data are not included in this public repository.

To run the analysis with authorized cohort data:

1. copy `config.example.R` to `config.R`;
2. edit the local paths in `config.R`;
3. run the scripts in numerical order, or use `scripts/run_all.R`.

`config.R`, participant-level datasets, fold assignments, and participant-level out-of-fold predictions are excluded from version control.

Aggregate tables generated by the scripts may be exported for manuscript verification.

## Data availability

UK Biobank participant data are available to approved researchers through the UK Biobank access process.

TwinGene data are part of the Swedish Twin Registry and cannot be distributed publicly because of ethical and legal restrictions. Access requires appropriate ethical and data-access approvals.

For information about Swedish Twin Registry data access, contact:

str-research@meb.ki.se

## Citation

Mostafaei S, Lindh J, Qin C, Mak JKL, Hägg S. Biological aging measures for all-cause mortality prediction in two population-based cohorts. Manuscript under review at *npj Aging*.

When using this repository, please also cite the repository release or commit corresponding to the manuscript version.

## License

Code is distributed under the MIT License; see `LICENSE.txt`.

The software license does not provide access to, or permission to redistribute, the underlying participant-level cohort data.

## Contact

Shayan Mostafaei  
Karolinska Institutet  
shayan.mostafaei@ki.se

Sara Hägg  
Karolinska Institutet  
sara.hagg@ki.se
