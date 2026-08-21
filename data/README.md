[README.md](https://github.com/user-attachments/files/31285564/README.md)
# Governed input-data contract

No participant-level data or licensed summary statistics are distributed with this repository. All controlled files must remain in an access-controlled local environment.

Provider names, catalogue identifiers and public access routes are listed in `DATA_ACCESS.md`.

## Dynamic mortality inputs

The mortality workflow requires five quality-controlled inputs. Paths are set in `config/config.R`.

### Harmonised longitudinal input

Format: RDS data frame.

Required columns:

- `cohort`: CHARLS, ELSA, HRS, MHAS or SHARE;
- `id`: cohort-native linkage identifier;
- `subject_id`: analysis identifier used by the harmonised cohort workflow;
- `wave`: numeric wave identifier;
- `affective_harmonized_v1`: harmonised affective score, higher values indicating greater burden;
- `functional_harmonized_v1`: harmonised functional score, higher values indicating greater limitation.

### Exact visit-date input

Format: CSV.

Required columns:

- `cohort`, `id_norm`, `wave`;
- `interview_year_raw`, `interview_month_raw`.

Month is set to July only when a valid interview year is available but month is unavailable. Conflicting subject-wave dates are prohibited.

### Wave interview-status input

Format: CSV.

Required columns:

- `cohort`, `id_norm`, `wave`;
- `proxy_status`: `direct`, `proxy_or_exit` or `unknown`.

Only direct interviews contribute time-updated affective and functional exposure measurements. CHARLS wave 1 is classified as direct under the accepted questionnaire-mode audit.

### Main mortality-endpoint input

Format: RDS data frame containing HRS, MHAS, SHARE and ELSA.

Required columns:

- `cohort`, `subject_id`;
- `baseline_date_direct`, `baseline_age_direct`;
- `primary_exit_date`, `primary_event`, `primary_eligible`;
- `sex_baseline`, `education_baseline`.

### Support-only mortality-endpoint input

Format: RDS data frame containing CHARLS with the same columns as the main endpoint input.

Mortality and censoring dates must reflect the adjudicated endpoint hierarchy documented for each cohort. A current visit must precede the endpoint date. Exposure information from proxy or end-of-life interviews is not used.

## Dementia-related outcome inputs

The dementia-related outcome workflow reuses `dynamic_candidate_exact.rds` from the mortality workflow and requires five controlled inputs.

### Identifier crosswalk

Format: RDS data frame.

Required columns:

- `cohort`: ELSA or SHARE;
- `candidate_id`: analysis identifier in `dynamic_candidate_exact.rds`;
- `provider_id`: identifier in the corresponding harmonized provider file.

The pair `cohort`, `candidate_id` and the pair `cohort`, `provider_id` must each be unique. The crosswalk is governed participant-level information and must not be committed.

### ELSA harmonized outcome panel

Format: Stata file.

Required identifier: `idauniqc`.

For waves 1-10, the workflow reads the harmonized variables `r[wave]demene`, `r[wave]iwstat`, `r[wave]iwy`, `r[wave]iwm` and `r[wave]proxy`. Dementia-related status must use 0 for no and 1 for yes. Proxy status uses 1 for a proxy interview.

### SHARE harmonized outcome panel

Format: Stata file.

Required identifier and country variables: `mergeid`, `country`.

For waves 2 and 4-9, the workflow reads `r[wave]alzdeme`, `r[wave]iwstat`, `r[wave]iwy`, `r[wave]iwm` and `r[wave]proxy`. Alzheimer disease or dementia status must use 0 for no and 1 for yes. Proxy values 1 and 2 are both classified as proxy interviews.

### Cohort cognition files

Format: one RDS long-format data frame per cohort.

Required columns:

- `id`: provider identifier corresponding to `provider_id` in the crosswalk;
- `wave`: cohort wave;
- `cognition_score`: cohort-native observed cognition score.

The cognition sensitivity uses the dynamic-entry wave or immediately preceding cohort wave and standardizes observed scores within the eligible cohort sample.

### Event-time rules

The primary analysis requires an explicit negative dementia-related status at dynamic entry. First subsequent positive status defines the survey-detected event date. ELSA and SHARE endpoints are analysed separately. Known death is retained as a competing event when supported by a dynamic risk interval. A sensitivity analysis censors follow-up at the last observed dementia-related assessment.

## Other governed resources

The full repository additionally requires:

- harmonised source files for CHARLS, ELSA, HRS, MHAS and SHARE;
- five formatted GWAS summary-statistic files for major depression, loneliness, frailty, self-rated health and walking pace;
- LDSC HapMap3 and European LD-score references;
- PLINK-format LD reference data for MAGMA and regional analyses;
- the documented FUMA job downloads;
- TPMI and FinnGen endpoint summary resources;
- prepared external omics and regional QTL resources used by the frozen analyses.

Each resource remains subject to its provider's access, acknowledgement and redistribution conditions. Scripts stop when a required file, column, checksum or stage marker is missing.

## Prohibited repository content

Do not commit participant identifiers, dates, row-level observations, genotype files, restricted GWAS/QTL files, credentials, access tokens or data-provider download links that expose controlled content.
