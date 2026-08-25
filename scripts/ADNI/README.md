[Uploading README.md…]()
# ADNI factor-PGS CDR-SB validation

This module provides a prespecified individual-level clinical validation of the affective and functional genetic domains in the Alzheimer's Disease Neuroimaging Initiative (ADNI). It uses frozen affective-factor and functional-factor polygenic scores (factor-PGSs) and longitudinal Clinical Dementia Rating-Sum of Boxes (CDR-SB) measurements.

This module is distinct from the TPMI and FinnGen gene-set endpoint-mapping analysis. The ADNI exposure is an individual-level factor-PGS, and the outcome is longitudinal CDR-SB. It must not be described as an endpoint-mapping analysis, a causal analysis, or a clinical prediction tool.

## Script order

1. `01_build_locked_cdrsb_dataset.R`
   - links governed ADNIMERGE clinical data to the frozen factor-PGS file;
   - retains the primary ancestry-proxy sample;
   - establishes participant baseline, time since baseline and frozen covariates;
   - writes a locked analysis dataset and sample-flow audit.

2. `02_fit_factor_pgs_cdrsb_models.R`
   - fits the single-score primary models and applies BH correction to the two prespecified PGS-by-time tests;
   - fits the joint-score model to assess whether the two associations persist after simultaneous adjustment for both scores.

3. `03_build_adni_validation_outputs.R`
   - fits the confirmatory model allowing baseline diagnosis-specific CDR-SB time slopes;
   - creates Figure 6 assets, aggregate source data, a residual diagnostic figure and submission-ready supplementary tables;
   - exports no participant-level data outside the governed local results directory.

## Inputs

Set these local paths in `config/config.R`; do not commit that file.

```r
Sys.setenv(
  ADNI_CLINICAL_FILE = "/approved/path/ADNIMERGE.csv",
  ADNI_FACTOR_PGS_FILE = "/approved/path/ADNI_factor_PGS_scores.tsv",
  ADNI_OUTPUT_ROOT = file.path(project_dir, "results", "adni_validation")
)
```

The clinical CSV requires `PTID`, `EXAMDATE`, `CDRSB`, `AGE`, `PTGENDER`, `PTEDUCAT`, `DX_bl` and `APOE4`. The frozen PGS TSV requires `IID` or `PTID`, `PGS_affective_z`, `PGS_functional_z`, `primary_ancestry_proxy` and `PC1` to `PC10`.

## Reproduction

Run R from the repository root after creating `config/config.R`.

```r
source("config/config.R", encoding = "UTF-8")
source("run_adni_validation.R", encoding = "UTF-8")
```

Or use the global pipeline dispatcher:

```r
Sys.setenv(DCV_RUN_PIPELINE = "adni_validation")
source("run_all.R", encoding = "UTF-8")
```

## Primary model and interpretation

The confirmatory model estimates both factor-PGS-by-time terms jointly, adjusts for baseline age, sex, education, APOE4, PC1-PC10, baseline diagnosis and diagnosis-specific CDR-SB time slopes, and includes a participant-specific random intercept. BH correction is applied across the two prespecified PGS-by-time terms.

The resulting estimates are longitudinal associations in an independent clinical cohort. They do not establish causal effects, neuroprotection, dementia risk prediction, clinical utility or mechanisms. Do not re-tune PGS weights, select SNP thresholds, switch outcomes or choose follow-up windows using ADNI results.

## Outputs

All generated outputs are written under `results/adni_validation/`, which is excluded from Git. The module creates a locked analysis dataset, model summaries, Figure 6 in SVG/PDF/TIFF/PNG formats, figure source data, residual diagnostics and `ADNI_CDRSB_Submission_Tables.xlsx`. Share only aggregate outputs permitted by the ADNI data-use agreement.
