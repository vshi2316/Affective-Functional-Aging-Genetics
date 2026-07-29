project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
scripts <- file.path(project_dir, "scripts", "Track_A", sprintf("%02d_%s.R", 1:6, c(
  "prepare_cohort_inputs", "derive_and_harmonize_domains", "measurement_invariance",
  "longitudinal_models", "meta_analysis_and_robustness", "build_trackA_evidence_package"
)))
for (script in scripts) source(script, encoding = "UTF-8", chdir = FALSE)
