[README.md](https://github.com/user-attachments/files/31404634/README.md)
# Affective-functional dynamics across ageing cohorts and genetic architectures

This repository contains the analysis code for a five-cohort study of affective burden and functional vulnerability. The workflow separates population-average, between-person and within-person associations; relates history-anchored burden and current deviation to mortality and dementia-related clinical events; and evaluates parallel affective and functional genetic architectures.

Participant-level cohort data, mortality and dementia records, licensed genome-wide association study summary statistics, molecular resources and reference panels are controlled by their original providers and are not redistributed.

## Analyses

The repository implements five linked components:

1. harmonisation and longitudinal modelling in CHARLS, ELSA, HRS, MHAS and SHARE;
2. timing, fixed-effect and within-between sensitivity analyses;
3. attained-age history-deviation models for all-cause mortality;
4. cohort-specific cause-specific models for subsequent dementia-related outcomes in ELSA and SHARE, with death handled as a competing event;
5. linkage disequilibrium score regression, genomic structural equation modelling, factor genome-wide association analyses, local architecture, molecular follow-up and external endpoint mapping.

## Repository structure

```text
config/                            local path template and input schemas
data/README.md                     governed input-data contract
scripts/Track_A/                   cohort harmonisation and longitudinal models
scripts/Longitudinal_robustness/   estimand and timing sensitivity analyses
scripts/Mortality/                 dynamic mortality analysis
scripts/Dementia/                  dementia-related clinical-event analysis
scripts/Track_B/                   LDSC, GenomicSEM and factor GWAS
scripts/Track_C/                   local architecture and molecular analyses
scripts/Track_D/                   external endpoint mapping
scripts/ADNI/                      individual-level ADNI factor-PGS clinical validation
software/                          recorded analysis environments
run_*.R                            stage-level entry points
```

The analysis scripts are listed in `scripts/SCRIPT_MANIFEST.csv` in execution order and mapped to their reporting role.

## Requirements

Use R 4.4 or later. Required R packages are listed in `DESCRIPTION`. Recorded environments for the longitudinal, mortality and figure workflows are provided in `software/`.

The genomic workflow also requires GenomicSEM, LDSC, LAVA, FUMA, MAGMA, PLINK and SMR. Install external software from its official distribution and comply with its licence. Local executables and reference resources are declared in the configuration file.

## Configuration

Copy the configuration template and replace the example paths with approved local resources:

```r
file.copy("config/config.example.R", "config/config.R")
source("config/config.R", encoding = "UTF-8")
```

`config/config.R` is excluded from Git.

The local architecture analysis also requires copies of the two manifest templates:

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

The mortality and dementia-related outcome input schemas are documented in `data/README.md`. Provider access routes are listed in `DATA_ACCESS.md`.

## Cohort and clinical-event workflow

Run R from the repository root:

```r
source("config/config.R", encoding = "UTF-8")
source("run_track_A.R", encoding = "UTF-8")
source("run_longitudinal_robustness.R", encoding = "UTF-8")
source("run_mortality_bridge.R", encoding = "UTF-8")
source("run_dementia_bridge.R", encoding = "UTF-8")
```

Track A prepares the cohort-specific domains and longitudinal models. The robustness workflow performs the fixed-effect, timing and within-between analyses. The mortality workflow constructs exact-date attained-age intervals and estimates age-varying history and deviation associations. The dementia workflow reuses those dynamic intervals, links governed ELSA and SHARE outcome records, fits cohort-specific cause-specific models and creates aggregate source data.

The two dementia-related endpoints retain their cohort wording. ELSA contributes first survey-reported dementia; SHARE contributes first survey-reported Alzheimer disease or dementia. The recorded date represents first survey detection rather than exact clinical onset.

Individual stages can also be selected with `run_all.R`:

```r
Sys.setenv(
  DCV_RUN_PIPELINE = paste(
    "track_a",
    "longitudinal_robustness",
    "mortality_bridge",
    "dementia_bridge",
    sep = ","
  )
)
source("run_all.R", encoding = "UTF-8")
```

## Genomic workflow

Track B contains a FUMA checkpoint. Run the pre-FUMA stage, submit the generated files using the recorded settings, place the downloaded outputs at the configured path and continue with the remaining stages.

```r
source("config/config.R", encoding = "UTF-8")

Sys.setenv(DCV_RUN_PIPELINE = "track_b_pre_fuma")
source("run_all.R", encoding = "UTF-8")

# Complete the documented FUMA submission before continuing.

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

## ADNI factor-PGS clinical validation

The ADNI module is an independent individual-level validation of the genetic domains. It uses frozen affective-factor and functional-factor polygenic scores with longitudinal Clinical Dementia Rating-Sum of Boxes (CDR-SB) measurements. It is distinct from the TPMI and FinnGen gene-set endpoint-mapping workflow.

After setting the three `ADNI_*` paths in `config/config.R`, run:

```r
source("config/config.R", encoding = "UTF-8")
source("run_adni_validation.R", encoding = "UTF-8")
```

Alternatively, select the module through the global dispatcher:

```r
Sys.setenv(DCV_RUN_PIPELINE = "adni_validation")
source("run_all.R", encoding = "UTF-8")
```

The module writes a locked analysis dataset, model summaries, Figure 6 assets, aggregate source data, residual diagnostics and supplementary tables under `results/adni_validation/`. The estimates are longitudinal associations; they do not establish causality, clinical utility or a prediction tool. Full input, execution and interpretation details are in `scripts/ADNI/README.md`.

## Outputs

Generated files are written under `results/`, which is excluded from Git. The clinical workflows create model summaries, vector and raster figures, aggregate source data and software-session records. Only aggregate outputs permitted by the applicable data-use agreements may be shared.

Do not upload participant identifiers, dates, row-level observations, genotype files, controlled summary statistics, credentials or licensed reference resources.

Because the principal inputs are access controlled, the complete workflow requires separate approval from the relevant providers. The repository supplies the custom analysis logic, input schemas, execution order and public access routes needed to reproduce the analyses in an authorised environment.

## Data and code availability

Data access conditions and provider links are given in `DATA_ACCESS.md`. The repository does not grant access to third-party data or alter their terms of use.

## Citation and licence

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. Data remain subject to the terms of their original providers.
