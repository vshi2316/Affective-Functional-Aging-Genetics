# Affective-functional dynamics across ageing cohorts and genetic architectures

Analysis code for a five-cohort longitudinal study of affective burden and functional vulnerability, a history-anchored mortality analysis, and parallel multivariate genomic analyses.

The repository contains the custom code used for the reported analyses. Participant-level cohort data, mortality records, licensed genome-wide association study summary statistics, molecular resources and reference panels are governed by their original providers and are not redistributed.

## Repository structure

```text
config/                            local path template and Track C input schemas
data/README.md                     governed-input contract
scripts/Track_A/                   cohort harmonisation and longitudinal models
scripts/Longitudinal_robustness/   estimand and timing sensitivity analyses
scripts/Mortality/                 dynamic mortality analysis
scripts/Track_B/                   LDSC, GenomicSEM and factor GWAS
scripts/Track_C/                   local architecture and molecular analyses
scripts/Track_D/                   external endpoint mapping
software/                          session information from the analyses
run_*.R                            stage-level entry points
```

The 35 analysis scripts are listed in `scripts/SCRIPT_MANIFEST.csv` in execution order and mapped to their reporting role.

## Requirements

Use R 4.4 or later. Required R packages are listed in `DESCRIPTION`. Recorded package versions are provided in:

- `software/longitudinal_sessionInfo.txt`
- `software/mortality_sessionInfo.txt`
- `software/figure_sessionInfo.txt`

The genomic workflow also uses GenomicSEM, LDSC, LAVA, FUMA, MAGMA, PLINK, SMR, coloc and susieR. Install third-party tools from their official distributions and comply with their licences. Paths to local executables and reference resources are set in the configuration file.

## Configuration

Copy the configuration template, edit local paths, and load it before an analysis stage:

```r
file.copy("config/config.example.R", "config/config.R")
source("config/config.R", encoding = "UTF-8")
```

`config/config.R` is excluded from Git. Track C local architecture additionally requires copies of the two manifest templates:

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

Replace the example paths with the approved local inputs. The required mortality input columns and interview-date rules are documented in `data/README.md`. Data providers, catalogue identifiers and access routes are listed in `DATA_ACCESS.md`.

## Execution order

### Cohort and mortality analyses

```r
source("config/config.R", encoding = "UTF-8")
source("run_track_A.R", encoding = "UTF-8")
source("run_longitudinal_robustness.R", encoding = "UTF-8")
source("run_mortality_bridge.R", encoding = "UTF-8")
```

Track A harmonises CHARLS, ELSA, HRS, MHAS and SHARE before cohort-specific models. The robustness stage performs the within-between decomposition and timing analyses. The mortality stage constructs attained-age start-stop data from quality-controlled visit and endpoint files, fits the prespecified models, and generates aggregate figure data.

### Genomic analyses

Track B contains a FUMA checkpoint. Run the pre-FUMA stage, submit the generated files with the settings recorded by the scripts, download the required FUMA outputs, and then continue with the remaining stages.

```r
source("config/config.R", encoding = "UTF-8")

Sys.setenv(DCV_RUN_PIPELINE = "track_b_pre_fuma")
source("run_all.R", encoding = "UTF-8")

# Pause here for the FUMA submission and download.

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

Individual stage entry points are retained for users who have access to only part of the governed input set.

## Outputs

Generated files are written under `results/`, which is excluded from Git. Only aggregate tables and figure source data permitted by the original data-use agreements may be shared. Do not upload participant identifiers, dates, row-level records, genotype files, controlled summary statistics, credentials or licensed reference resources.

Because the primary inputs are access controlled, the public workflow cannot execute without separately approved data access. The repository provides the complete custom analysis logic, required input schemas, execution order, software environments and public data-access routes.

## Citation and licence

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. Data remain subject to the terms of their original providers.
