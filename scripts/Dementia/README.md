[Uploading README.md…]()
# Dynamic dementia-related outcome workflow

This directory contains the analysis of history-anchored affective and functional measures in relation to subsequent dementia-related outcomes in ELSA and SHARE. Death is handled as a competing event in cohort-specific cause-specific models.

## Script order

1. `01_build_event_dataset.R`
   - reads the dynamic history-deviation intervals created by the mortality workflow;
   - links provider-governed ELSA and SHARE outcome records through a controlled identifier crosswalk;
   - requires an explicit negative dementia-related status at dynamic entry;
   - identifies the first subsequent dementia-related outcome and known death;
   - retains known death when supported by a dynamic interval;
   - constructs strict outcome-observation and baseline-status sensitivity datasets;
   - records event counts, proxy status, later negative reports and observation support.

2. `02_fit_competing_event_models.R`
   - fits cohort-specific attained-age cause-specific Cox models;
   - compares current affective and functional scores with separate historical burden and current deviation components;
   - uses participant-clustered robust standard errors;
   - stratifies SHARE baseline hazards by country;
   - tests the two history-deviation equality constraints jointly;
   - evaluates cognition adjustment, early events, proxy reports, later negative reports, baseline status and observation censoring;
   - fits a complementary log-log interval model as a timing sensitivity analysis.

3. `03_build_figure_and_source_data.R`
   - reads model outputs without refitting;
   - builds the four-panel dementia-related clinical-event figure;
   - exports SVG, PDF, 600-dpi TIFF and PNG preview files;
   - writes aggregate source data and submission tables with descriptive worksheet names for insertion at their final citation positions.

## Outcome definitions

ELSA contributes the first survey-reported dementia outcome. SHARE contributes the first survey-reported Alzheimer disease or dementia outcome. These endpoints are analysed separately because their wording differs. The recorded wave date represents first survey detection rather than exact clinical onset.

## Reproduction

Complete the controlled input paths in `config/config.R`. Run the mortality workflow before this module because its dynamic interval dataset supplies the history and deviation measures.

```r
source("run_mortality_bridge.R", encoding = "UTF-8")
source("run_dementia_bridge.R", encoding = "UTF-8")
```

Use an empty `results/dementia_bridge/` directory for a full reproduction.

## Interpretation

The models estimate cause-specific prognostic associations. They do not estimate exact dementia onset, absolute risk, diagnostic accuracy, treatment effects or causal effects.
