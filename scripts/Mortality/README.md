# Dynamic mortality workflow

This directory contains the complete public implementation of the history-anchored estimand–mortality analysis.

## Script order

1. `01_build_dynamic_candidate.R`
   - reads only quality-controlled governed inputs;
   - keeps direct visits with jointly observed affective and functional scores;
   - enforces visit-before-endpoint ordering;
   - constructs attained-age start–stop intervals;
   - calculates prior history means and current deviations without using future observations;
   - applies common cohort-specific first-visit scaling;
   - audits events and age-band support;
   - writes a candidate lock without estimating an association.

2. `02_fit_dynamic_models.R`
   - reads the locked candidate and feasibility gates;
   - fits four exposure components simultaneously;
   - permits exposure coefficients to vary across four attained-age bands;
   - uses subject-cluster robust Cox variance;
   - applies the frozen primary, replication, sensitivity and support-only cohort roles;
   - calculates the four primary global Wald tests and Holm-adjusted P values;
   - reports random-effects meta-analysis, heterogeneity and prediction intervals;
   - runs only the prespecified delayed-entry, lagged-visit, history-stability and interaction sensitivities;
   - writes a final model lock.

3. `03_build_figure_and_source_data.R`
   - does not refit a model;
   - creates cohort-specific, meta-analytic and history–deviation contrast panels;
   - exports PDF, SVG and 600-dpi LZW-compressed TIFF;
   - writes aggregate figure Source Data and checksums.

## Interpretation boundary

`H_A` and `H_F` are history-anchored prior burden, not pure time-invariant traits. `D_A` and `D_F` are current deviations from each participant's observed prior history. Associations are prognostic and do not establish causation or validated individual prediction.

## Reproduction

Set the governed input paths in `config/config.R`, start R in the repository root and run:

```r
source("run_mortality_bridge.R", encoding = "UTF-8")
```

Existing locks intentionally stop an attempted rerun. Remove a lock only when reproducing the workflow in a new, empty output directory and never after inspecting results.
