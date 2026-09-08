# Affective-functional dynamics across ageing cohorts and genetic architectures

This repository contains the complete analysis workflow for a five-cohort study of affective burden and functional vulnerability. The code separates population-average, between-person and within-person associations; distinguishes historical burden from current deviation; links these estimands to mortality and dementia-related events; evaluates five-year mortality discrimination; and resolves correlated affective and functional genetic architectures.

The study uses CHARLS, ELSA, HRS, MHAS and SHARE. Participant-level cohort data, mortality records, dementia-related outcomes, genotypes, controlled summary statistics, molecular resources and reference panels remain with their original providers and are not redistributed.

## Scientific workflow

The repository implements the following linked analyses:

1. cohort-specific harmonisation of affective, functional and cognition measures;
2. longitudinal mixed models, lagged associations, participant-plus-wave fixed effects and within-between decomposition;
3. exact-date attained-age history-deviation models for all-cause mortality;
4. HRS-trained five-year mortality prediction with internal validation and external evaluation;
5. cohort-specific cause-specific models for dementia-related events in ELSA and SHARE, with death handled as a competing event;
6. LDSC, genomic structural equation modelling, factor GWAS, QSNP, FUMA and MAGMA analyses;
7. local shared architecture, gene-set definition, external omics and molecular follow-up;
8. TPMI and FinnGen endpoint mapping with chromosome and locus influence analysis;
9. optional individual-level ADNI factor-PGS validation.

## Repository structure

```text
config/                            local path template and input manifests
data/README.md                     governed input-data contract
scripts/Track_A/                   harmonisation and longitudinal models
scripts/Longitudinal_robustness/   timing, lagged and within-between analyses
scripts/Mortality/                 dynamic mortality analysis
scripts/Mortality_prediction/      five-year mortality prediction and submission outputs
scripts/Dementia/                  dementia-related clinical-event analysis
scripts/Track_B/                   LDSC, GenomicSEM, factor GWAS and FUMA/MAGMA
scripts/Track_C/                   local architecture, omics and molecular analyses
scripts/Track_D/                   external endpoint mapping and influence analysis
scripts/ADNI/                      optional individual-level factor-PGS validation
scripts/SCRIPT_MANIFEST.csv        ordered script inventory and runner mapping
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

Each script checks its direct R dependencies before reading controlled data. Every locked analysis stage writes a session record and MD5 lineage files beside its results.

## Configuration

Run all commands from the repository root. Create the local configuration file:

```r
file.copy("config/config.example.R", "config/config.R")
source("config/config.R", encoding = "UTF-8")
```

Edit only paths and executable locations in `config/config.R`. The file is excluded from Git.

The local architecture workflow also requires local copies of the two manifest templates:

```r
file.copy(
  "config/trackC1_trait_manifest.example.csv",
  "config/trackC1_trait_manifest.csv"
)
file.copy(
  "config/trackC1_pair_manifest.example.csv",
  "config/trackC1_pair_manifest.csv"
)
```

Required input schemas are defined in `data/README.md`. Provider access routes are listed in `DATA_ACCESS.md`.

## Cohort and clinical-event workflow

```r
source("config/config.R", encoding = "UTF-8")
source("run_track_A.R", encoding = "UTF-8")
source("run_longitudinal_robustness.R", encoding = "UTF-8")
source("run_mortality_bridge.R", encoding = "UTF-8")
source("run_mortality_prediction.R", encoding = "UTF-8")
source("run_dementia_bridge.R", encoding = "UTF-8")
```

The same stages can be selected through the global dispatcher:

```r
Sys.setenv(
  DCV_RUN_PIPELINE = paste(
    "track_a",
    "longitudinal_robustness",
    "mortality_bridge",
    "mortality_prediction",
    "dementia_bridge",
    sep = ","
  )
)
source("run_all.R", encoding = "UTF-8")
```

### Longitudinal and dynamic mortality analyses

Track A prepares the cohort-specific measures and longitudinal models. The robustness module evaluates timing, lagged, fixed-effect and within-between specifications. The mortality bridge constructs exact-date attained-age intervals and estimates age-varying associations for affective history, affective deviation, functional history and functional deviation.

### Five-year mortality prediction

The prediction module reuses the locked mortality candidate and adds no new participant-level source file. It selects each participant's earliest eligible landmark at or after age 50 years and compares three nested models:

- Model 0: sex, standardized education and calendar year, with attained age represented by the Cox baseline hazard;
- Model 1: Model 0 plus current affective and functional scores;
- Model 2: Model 0 plus affective history, affective deviation, functional history and functional deviation.

HRS provides a 70% development sample and 30% internal-validation sample. Coefficients and cumulative baseline hazards are frozen for evaluation in MHAS and SHARE, with ELSA providing sensitivity evidence and CHARLS contributing support-only results. The module reports five-year AUC, Brier score, calibration, grouped calibration, nested model tests and 1,000 paired participant-level bootstrap repetitions.

It produces:

- main Table 2 in CSV, XLSX and DOCX formats;
- Supplementary Fig. 34 in SVG, PDF, 600 dpi TIFF and PNG formats;
- Supplementary Tables S54-S57;
- source data, figure legend, quality-assurance record and MD5 manifest.

Full design and output details are in `scripts/Mortality_prediction/README.md`.

### Dementia-related outcomes

The dementia workflow reuses the dynamic intervals and links governed ELSA and SHARE outcome records. ELSA contributes first survey-reported dementia. SHARE contributes first survey-reported Alzheimer disease or dementia. Event dates represent first survey detection. Cause-specific models retain death as a competing event when it occurs within a valid risk interval.

## Genomic workflow

Track B contains a FUMA checkpoint. Run the pre-FUMA stage, complete the documented FUMA submission, place the downloaded outputs at the configured path and continue with the remaining stages.

```r
source("config/config.R", encoding = "UTF-8")

Sys.setenv(DCV_RUN_PIPELINE = "track_b_pre_fuma")
source("run_all.R", encoding = "UTF-8")

# Complete the FUMA submission using the settings written by the pre-FUMA stage.

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

The ADNI module generates aggregate Figure 6 assets and Supplementary Tables S58-S59. Its input and output contract is documented in `scripts/ADNI/README.md`.

## Output structure

```text
results/
├── mortality_bridge/
├── mortality_prediction/
│   ├── landmark/
│   ├── models/
│   └── submission/
├── dementia_bridge/
├── adni_validation/
└── ...
```

`results/` remains local. Do not copy landmark RDS files, participant identifiers, dates, genotypes, controlled summary statistics or licensed resources into the repository. Tables, figures and other manuscript outputs are generated locally by the submission-output stage and are not required for running the analysis.

## Supplementary numbering

- S1-S53: longitudinal, mortality, dementia, genomic and molecular analyses;
- S54-S57: five-year mortality prediction;
- S58-S59: optional ADNI validation.

## Reproducibility controls

- repository-relative paths and one local configuration file;
- ordered script manifest and stage-level runners;
- explicit input schemas and provider access routes;
- fixed seeds and locked analysis protocols;
- input and output MD5 manifests;
- stage gates that stop before a lock is written;
- session information for every major stage;
- vector, raster and tabular source-data exports for manuscript figures are written locally under `results/`;
- no participant-level results or manuscript binaries are required in the code repository.

All R scripts can be syntax-checked without controlled data. Full numerical reproduction requires authorised access to the underlying resources.

## Data and code availability

Data access conditions are listed in `DATA_ACCESS.md`. The repository does not grant access to third-party data or alter provider terms. Local paths, controlled data and results are excluded through `.gitignore`.

## Citation and licence

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. Data remain subject to the terms of their original providers.
