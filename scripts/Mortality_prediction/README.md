# Five-year mortality prediction

This module evaluates whether current affective-functional status and the history-deviation representation improve five-year mortality discrimination beyond demographic covariates.

## Scientific design

- Population: participants aged at least 50 years at the landmark.
- Landmark: each participant's earliest eligible direct joint visit at or after age 50 years.
- Outcome: all-cause mortality within five years.
- Development: participant-level 70% HRS sample.
- Internal validation: remaining 30% of HRS participants.
- External evaluation: MHAS and SHARE.
- Sensitivity evidence: ELSA.
- Support-only evaluation: CHARLS.
- Model 0: sex, standardized education and calendar year; attained age is represented by the Cox baseline hazard.
- Model 1: Model 0 plus current affective and functional scores.
- Model 2: Model 0 plus affective history, affective deviation, functional history and functional deviation.
- Primary comparison: Model 2 minus Model 1.

HRS coefficients and cumulative baseline hazards are frozen for every internal and external evaluation. Five-year AUC, Brier score, calibration slope, calibration intercept and grouped calibration are estimated with inverse-probability-of-censoring weighting. Differences in AUC use 1,000 paired participant-level bootstrap repetitions.

## Input

The module reuses the locked dynamic candidate produced by `scripts/Mortality/01_build_dynamic_candidate.R`:

```text
results/mortality_bridge/candidate/dynamic_candidate_exact.rds
results/mortality_bridge/candidate/dynamic_candidate_lock.txt
results/mortality_bridge/candidate/dynamic_candidate_gates.csv
results/mortality_bridge/candidate/fixed_scaling.csv
```

No new participant-level input is required.

## Run

From the repository root:

```r
source("config/config.R", encoding = "UTF-8")
source("run_mortality_prediction.R", encoding = "UTF-8")
```

The global dispatcher can run the mortality bridge and prediction module in sequence:

```r
Sys.setenv(
  DCV_RUN_PIPELINE = paste(
    "mortality_bridge",
    "mortality_prediction",
    sep = ","
  )
)
source("run_all.R", encoding = "UTF-8")
```

For a short execution test, set `MORTALITY_PREDICTION_BOOTSTRAP_B=100` in the process environment. Manuscript results use 1,000 repetitions.

## Stages

### 01_build_age50_landmark_dataset.R

Builds and locks one age-eligible landmark per participant. It verifies event terminality, a single landmark per participant, temporal ordering, age eligibility, history depth and the current-score identity.

### 02_fit_and_validate_five_year_models.R

Fits the three HRS development models, freezes coefficients and cumulative baseline hazards, evaluates all cohorts, performs nested Wald and likelihood-ratio tests, and runs the paired participant bootstrap. Cox cumulative baseline hazards are evaluated as right-continuous step functions.

### 03_build_submission_outputs.R

Reads locked aggregate results and produces:

- Table 2 in CSV, XLSX and DOCX formats;
- Supplementary Fig. 34 in SVG, PDF, 600 dpi LZW TIFF and PNG preview formats;
- clean figure source data;
- Supplementary Tables S54-S57 in XLSX and CSV formats;
- an integrated `Supplementary_Tables.xlsx` when an existing S1-S53 workbook is configured;
- figure legend, reporting summary, quality-assurance record and MD5 output manifest.

No model is refitted in the submission-output stage.

## Output structure

```text
results/mortality_prediction/
├── landmark/
├── models/
└── submission/
```

`results/` is excluded from Git. Manuscript-facing tables, figures and source data are generated locally under `results/mortality_prediction/submission/` and are not required in the code repository.

## Supplementary-table mapping

- S54: age-50 landmark sample, follow-up, baseline-hazard support and history depth.
- S55: HRS development coefficients and nested model comparisons.
- S56: five-year AUC, Brier score and calibration across all three models.
- S57: paired bootstrap differences in AUC and analytic audit.
- S58-S59: ADNI validation outputs, if included in the reporting package.

## Reproducibility controls

- fixed age threshold, horizon, HRS split seed and default bootstrap count;
- input and output MD5 manifests;
- stage locks and session information;
- no participant identifiers in submission outputs;
- replicate-level bootstrap results contain iteration numbers and aggregate metric differences only;
- analysis gates must pass before locks or submission outputs are written.
