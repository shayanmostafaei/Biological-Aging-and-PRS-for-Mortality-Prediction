# Biological aging measures for all-cause mortality prediction in two population-based cohorts

Analysis materials accompanying the manuscript by Shayan Mostafaei, Jakob Lindh, Chenxi Qin, Jonathan K. L. Mak, and Sara Hägg.

**Manuscript status:** Under review at *npj Aging*.

## Overview

This study evaluates associations of five complementary biological aging measures with all-cause mortality and their incremental prognostic contribution beyond chronological age and conventional risk factors:

- PhenoAge
- Klemera–Doubal method (KDM)
- Homeostatic dysregulation (HD)
- Frailty index (FI)
- Leukocyte telomere length (TL)

A common analytical framework was applied independently in TwinGene and UK Biobank (UKB). Prediction models were fitted and evaluated separately within each cohort.

The repository retains its historical name, which includes “PRS.” Polygenic risk scores are not included in the current analyses.

## Code status

The manuscript methods are documented below. Earlier analysis scripts are retained in `legacy/` for historical reference and do not reproduce the revised manuscript. The final scripts used to generate the revised manuscript results have not yet been incorporated into the current analysis workflow.

## Study cohorts

| Characteristic | TwinGene | UK Biobank |
|---|---:|---:|
| Analytic sample | 9,617 | 314,064 |
| Deaths | 3,090 | 25,084 |
| Mortality proportion | 32.1% | 8.0% |
| Median follow-up, years | 16.76 | 13.79 |

The outcome was all-cause mortality. Time since baseline was the primary survival timescale.

TwinGene mortality follow-up extended through 1 February 2024. UKB administrative censoring was 30 December 2022.

## Biological aging measures and preprocessing

The clinical biomarker-based measures used cohort-specific implementations, incorporating eight biomarkers in TwinGene and 18 in UKB.

For descriptive summaries and correlations, analyses used raw PhenoAge, raw KDM, raw HD, FI on its deficit-proportion scale, and cohort-specific adjusted TL.

For mortality analyses:

- PhenoAge and KDM were separately regressed on chronological age using natural splines with 3 degrees of freedom.
- The resulting residuals represented deviation from the biological-age value expected for chronological age.
- HD was natural-log transformed and was not residualized on chronological age.
- PhenoAge residual, KDM residual, log HD, FI, and TL were standardized to mean 0 and SD 1.

For association analyses, preprocessing was estimated within each full analytic cohort. For prediction analyses, residualization models and standardization parameters were estimated within training folds and applied to held-out participants.

## Mortality association analyses

Each biological aging measure was evaluated individually in a covariate-adjusted Cox proportional hazards model. The biological aging measures were not mutually adjusted.

Adjustment included:

- Chronological age: natural spline with 4 degrees of freedom
- Sex
- BMI: natural spline with 3 degrees of freedom
- Smoking: ever versus never
- Education: high/intermediate versus low

TwinGene models used twin-pair-clustered robust variance estimation. Hazard ratios were reported per 1-SD higher biological aging measure.

Benjamini–Hochberg correction was applied separately within each cohort across the five primary biological aging associations.

## Cross-validated mortality prediction

Three nested Cox prediction models were evaluated:

1. **Age model:** chronological age.
2. **Conventional model:** chronological age, sex, BMI, smoking, and education.
3. **Biological aging-enhanced model:** conventional predictors plus PhenoAge residual, FI, and TL.

Age and BMI used the spline specifications described above.

Performance was evaluated using a fixed 10-fold partition within each cohort. The same folds were reused across prediction models and alternative biomarker-based sensitivity models.

In TwinGene, co-twins were assigned to the same fold. UKB used a mortality-balanced partition generated with a fixed random seed.

Residualization and standardization were estimated within each training fold before application to the corresponding held-out fold.

## Prediction performance

Evaluation used out-of-fold predictions and included:

- Harrell’s C-index
- Censoring-aware time-dependent AUC at 5 and 10 years
- Censoring-aware Brier scores at 5 and 10 years
- Calibration across deciles of predicted mortality risk
- Paired uncertainty analyses for incremental performance
- Exploratory survival decision-curve analysis

| Cohort | Conventional model C-index | Biological aging-enhanced model C-index |
|---|---:|---:|
| TwinGene | 0.780 | 0.787 |
| UKB | 0.720 | 0.746 |

The corresponding increases in C-index were 0.0073 (95% CI 0.0052–0.0094) in TwinGene and 0.0263 (95% CI 0.0240–0.0286) in UKB.

Uncertainty analyses used paired out-of-fold predictions. TwinGene used 1,000 twin-pair cluster-bootstrap replicates. UKB used paired concordance covariance for C-index differences and paired influence-function inference for AUC and Brier-score differences.

These uncertainty estimates condition on the fixed folds and out-of-fold predictions. 

Decision-curve analyses were exploratory and do not establish clinical utility.

## Effect modification

Interactions with sex, BMI, smoking, and education were evaluated for each of the five biological aging measures.

This yielded 20 interaction tests per cohort. One-degree-of-freedom Wald tests were used, with Benjamini–Hochberg correction applied separately within each cohort across the 20 tests.

TwinGene interaction models retained twin-pair-clustered robust variance estimation.

Subgroup-specific associations were estimated for interactions surviving correction. BMI-related estimates were evaluated at the cohort-specific mean and at −1 and +1 SD.

## Diagnostics and sensitivity analyses

The revised study included:

- Schoenfeld residual diagnostics.
- Interval-specific associations during 0–5, 5–10, and >10 years of follow-up.
- A 2-year landmark analysis.
- Attained age as the underlying survival timescale.
- Alternative prediction models replacing PhenoAge residual with KDM residual or log HD while retaining FI and TL.

Alternative prediction models used the same cross-validation folds, conventional covariates, training-fold preprocessing approach, and performance metrics as the primary enhanced model.

## Software

The manuscript analyses were conducted in R 4.5.1.

Principal packages reported in the manuscript were:

- `survival`
- `riskRegression`
- `emmeans`
- `BioAge`
- `dcurves`
- `ggplot2`

Exact package versions, installation sources, and execution instructions should accompany the final analysis-code release.

## Data availability

Participant-level data are not distributed through this repository.

UK Biobank data are available to registered researchers through the UK Biobank access process:
https://www.ukbiobank.ac.uk/

TwinGene data are part of the Swedish Twin Registry and are not publicly available because of ethical and legal restrictions. Access requires approval from the Swedish Ethical Review Authority and the Steering Committee of the Swedish Twin Registry.

For information about TwinGene data access, contact:
str-research@meb.ki.se

## Citation

Mostafaei S, Lindh J, Qin C, Mak JKL, Hägg S. Biological aging measures for all-cause mortality prediction in two population-based cohorts. Manuscript under review at npj Aging.

Please also identify the repository release or commit used when citing these analysis materials.

## License

Repository code is distributed under the MIT License; see `LICENSE.txt`.

The software license does not grant access to, or permission to redistribute, restricted cohort data.

## Contact

Shayan Mostafaei: shayan.mostafaei@ki.se

Sara Hägg: sara.hagg@ki.se
