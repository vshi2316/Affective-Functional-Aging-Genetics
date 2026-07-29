## Frozen longitudinal robustness workflow.
## Run from the repository root after configuring DCV_BASE_DIR.

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")

if (!nzchar(base_dir) || !dir.exists(base_dir)) {
  stop("Set DCV_BASE_DIR to the controlled five-cohort data root before running the longitudinal robustness workflow.")
}

scripts <- file.path(
  project_dir,
  "scripts",
  "Longitudinal_robustness",
  c(
    "01_freeze_and_audit.R",
    "02_longitudinal_robustness.R",
    "03_genomic_boundary_audit.R",
    "04_lagged_and_sample_reconciliation.R",
    "05_within_between_synthesis.R"
  )
)

missing <- scripts[!file.exists(scripts)]
if (length(missing)) stop("Missing longitudinal robustness scripts: ", paste(missing, collapse = "; "))

for (script in scripts) {
  message("Running ", basename(script))
  source(script, encoding = "UTF-8", echo = FALSE)
}

message("Frozen longitudinal robustness workflow completed. No additional model family is authorised.")
