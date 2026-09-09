[Uploading README.md…]()
# Affective-functional dynamics across ageing cohorts and genetic architectures

This repository contains the reproducible analysis workflow for a five-cohort study of affective burden, functional vulnerability, mortality, dementia-related outcomes and genetic architecture. The workflow separates population-average, between-person and within-person information, distinguishes historical burden from current deviation, evaluates five-year mortality risk stratification and performs complementary genomic and molecular analyses.

The cohorts are CHARLS, ELSA, HRS, MHAS and SHARE. Participant-level data, mortality records, dementia-related outcomes, genotypes, controlled summary statistics, molecular resources and reference panels remain with their original providers and are accessed locally under their applicable data-use conditions.

## Scientific workflow

The repository contains the following analysis modules:

1. cohort harmonisation for affective, functional and cognition domains;
2. longitudinal mixed models, lagged associations and within-between decompositions;
3. exact-date attained-age history-deviation mortality models;
4. five-year mortality prediction with internal and external evaluation;
5. clinical translation using continuous history-deviation models in adults aged 50–64 years with at least one prior functional measurement;
6. external calibration and cohort-specific intercept recalibration for transportability assessment;
7. dementia-related cause-specific analyses with death retained as a competing event;
8. LDSC, GenomicSEM, factor GWAS, QSNP, FUMA and MAGMA analyses;
9. local shared architecture, gene-set, external omics and molecular follow-up analyses;
10. external endpoint mapping and leave-one-chromosome or leave-one-locus analyses;
11. optional individual-level ADNI factor-PGS validation.

## Repository structure

```text
config/                            local path template and input manifests
data/README.md                     governed input-data contract
scripts/Track_A/                   harmonisation and longitudinal models
scripts/Longitudinal_robustness/   timing, lagged and within-between analyses
scripts/Mortality/                 dynamic mortality analysis
scripts/Mortality_prediction/      five-year prediction and clinical translation
scripts/Dementia/                  dementia-related clinical-event analysis
scripts/Track_B/                   LDSC, GenomicSEM, factor GWAS and FUMA/MAGMA
scripts/Track_C/                   local architecture, omics and molecular analyses
scripts/Track_D/                   external endpoint mapping and influence analyses
scripts/ADNI/                      optional factor-PGS validation
scripts/SCRIPT_MANIFEST.csv        ordered script inventory
software/README.md                 software-environment recording contract
run_*.R                            stage-level entry points
DESCRIPTION                        R dependency inventory
DATA_ACCESS.md                     provider access routes and redistribution rules
CITATION.cff                       repository citation metadata
```

Generated participant-level and intermediate results are written under `results/`, which is excluded from Git.

## Requirements

- R 4.4 or later;
- packages listed in `DESCRIPTION`;
- GenomicSEM, LAVA, cfdr.pleio and ldscr for the genomic workflow;
- LDSC, FUMA, MAGMA, PLINK and SMR for their corresponding external stages.

Each analysis script checks direct R dependencies before reading controlled data. Locked stages write session information and input or output lineage files beside their results.

## Configuration

Run commands from the repository root. Create the local configuration file:

```r
file.copy("config/config.example.R", "config/config.R")
source("config/config.R", encoding = "UTF-8")
```

Edit local paths and executable locations in `config/config.R`. This file is excluded from Git.

The configuration includes the clinical translation paths:

```r
CLINICAL_LANDMARK_INPUT
CLINICAL_PRIOR1_OUTPUT_ROOT
CLINICAL_RECALIBRATION_OUTPUT_ROOT
```

The default example uses:

```text
results/mortality_prediction/landmark/age50_landmark_dataset.rds
results/clinical_risk_prior1_continuous/
results/clinical_risk_recalibration/
```

The local architecture workflow also requires copies of the two manifest templates:

```r
file.copy("config/trackC1_trait_manifest.example.csv", "config/trackC1_trait_manifest.csv")
file.copy("config/trackC1_pair_manifest.example.csv", "config/trackC1_pair_manifest.csv")
```

Required input schemas are documented in `data/README.md`. Provider access routes are listed in `DATA_ACCESS.md`.

## Cohort and clinical-event workflow

```r
source("config/config.R", encoding = "UTF-8")
source("run_track_A.R", encoding = "UTF-8")
source("run_longitudinal_robustness.R", encoding = "UTF-8")
source("run_mortality_bridge.R", encoding = "UTF-8")
source("run_mortality_prediction.R", encoding = "UTF-8")
source("run_clinical_translation.R", encoding = "UTF-8")
source("run_dementia_bridge.R", encoding = "UTF-8")
```

The global dispatcher accepts the same modules:

```r
Sys.setenv(
  DCV_RUN_PIPELINE = paste(
    "track_a",
    "longitudinal_robustness",
    "mortality_bridge",
    "mortality_prediction",
    "clinical_translation",
    "dementia_bridge",
    sep = ","
  )
)
source("run_all.R", encoding = "UTF-8")
```

## Five-year mortality prediction

The prediction module selects the earliest eligible landmark at or after age 50 years for each participant and compares three nested models:

- Model 0: sex, standardized education and calendar year, with attained age represented by the Cox baseline hazard;
- Model 1: Model 0 plus current affective and functional scores;
- Model 2: Model 0 plus affective history, affective deviation, functional history and functional deviation.

HRS provides a 70% development sample and a 30% internal-validation sample. Model coefficients and cumulative baseline hazards are frozen for evaluation in MHAS and SHARE. ELSA provides sensitivity evidence and CHARLS contributes support-only results. The module reports five-year AUC, Brier score, calibration, grouped calibration, nested model comparisons and paired participant-level bootstrap differences.

The main prediction stages are:

```text
scripts/Mortality_prediction/01_build_age50_landmark_dataset.R
scripts/Mortality_prediction/02_fit_and_validate_five_year_models.R
scripts/Mortality_prediction/03_build_submission_outputs.R
```

## Clinical translation and external calibration

The clinical translation module evaluates adults aged 50–64 years with `prior_n >= 1`. Continuous history and deviation terms provide the primary model. High historical functional burden based on the cohort-specific H_F P75 is retained for descriptive risk-gradient displays.

Run both clinical stages together:

```r
source("run_clinical_translation.R", encoding = "UTF-8")
```

The clinical stages are:

```text
scripts/Mortality_prediction/05_clinical_risk_prior1_continuous.R
scripts/Mortality_prediction/06_external_recalibration.R
```

The first stage audits event support and evaluates the continuous history-deviation model in HRS development, HRS internal validation and external cohorts. The second stage freezes the HRS development coefficients and updates only a cohort-specific five-year calibration intercept:

```text
logit(P_recalibrated) = logit(P_original) + cohort-specific intercept
```

The recalibration stage reports calibration-in-the-large, calibration slope, observed risk, predicted risk before and after recalibration, Brier score before and after recalibration, decile-level calibration data and publication-ready calibration plots. The recalibrated output is a transportability analysis with fixed predictor coefficients.

Clinical translation outputs are written locally to:

```text
results/clinical_risk_prior1_continuous/
results/clinical_risk_recalibration/
```

Important output files include:

```text
prior1_event_audit.csv
prior1_event_gates.csv
prior1_continuous_model_metrics.csv
prior1_continuous_auc_bootstrap.csv
prior1_preserved_current_HF_gradient.csv
prior1_continuous_decision_curve.csv
external_recalibration_summary.csv
external_recalibration_decile_data.csv
Prior1_Continuous_HF_Gradient.tiff
External_Recalibration_Calibration.tiff
```

## Dementia-related outcomes

The dementia workflow reuses governed dynamic intervals and links ELSA and SHARE outcome records. First survey-reported dementia or Alzheimer disease events are analysed with cause-specific models that retain death as a competing event when it occurs within a valid risk interval.

## Genomic and molecular workflow

Track B contains a FUMA checkpoint. Run the pre-FUMA stage, complete the documented FUMA submission, place downloaded outputs at the configured location and continue with the post-FUMA stages.

```r
source("config/config.R", encoding = "UTF-8")

Sys.setenv(DCV_RUN_PIPELINE = "track_b_pre_fuma")
source("run_all.R", encoding = "UTF-8")

Sys.setenv(
  DCV_RUN_PIPELINE = paste(
    "track_b_post_fuma",
    "track_d",
    "track_c_architecture",
    "track_c_gene_sets",
    "track_c_omics",
    "track_c_molecular",
    "track_b_audit",
    sep = ","
  )
)
source("run_all.R", encoding = "UTF-8")
```

## ADNI factor-PGS validation

After configuring `ADNI_CLINICAL_FILE`, `ADNI_FACTOR_PGS_FILE` and `ADNI_OUTPUT_ROOT`:

```r
source("config/config.R", encoding = "UTF-8")
source("run_adni_validation.R", encoding = "UTF-8")
```

## Output and data policy

```text
results/
├── mortality_bridge/
├── mortality_prediction/
├── clinical_risk_prior1_continuous/
├── clinical_risk_recalibration/
├── dementia_bridge/
├── adni_validation/
└── ...
```

The repository contains code, configuration templates, input schemas and documentation. It does not contain participant-level records, identifiers, exact dates, genotypes, controlled summary statistics, licensed reference panels or manuscript result binaries. Local results are generated by the analysis stages and remain excluded through `.gitignore`.

## Reproducibility controls

- repository-relative paths and one local configuration file;
- ordered script manifest and stage-level runners;
- explicit input schemas and provider access routes;
- fixed seeds for bootstrap procedures;
- frozen HRS development coefficients for external evaluation;
- calibration intercept updating isolated to the recalibration stage;
- stage locks and session information;
- input and output lineage records where applicable;
- syntax-checkable R scripts without controlled data.

Full numerical reproduction requires authorised access to the underlying cohort and reference resources.

## Data and code availability

Data access conditions are listed in `DATA_ACCESS.md`. The repository does not grant access to third-party data or change provider terms. Local paths, controlled data and generated results are excluded through `.gitignore`.

## Citation and licence

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. Data remain subject to the terms of their original providers.
