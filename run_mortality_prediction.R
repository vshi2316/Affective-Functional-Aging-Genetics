#!/usr/bin/env Rscript
# Dispatcher for the age-50 five-year mortality prediction workflow.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(DCV_PROJECT_DIR = project_dir)

config_file <- file.path(project_dir, "config", "config.R")
if (file.exists(config_file)) source(config_file, encoding = "UTF-8")

scripts <- file.path(
  project_dir,
  "scripts",
  "Mortality_prediction",
  c(
    "01_build_age50_landmark_dataset.R",
    "02_fit_and_validate_five_year_models.R",
    "03_build_submission_outputs.R"
  )
)
missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts)) {
  stop("Missing mortality-prediction script(s): ", paste(basename(missing_scripts), collapse = ", "))
}

Sys.setenv(
  MORTALITY_PREDICTION_LANDMARK_SCRIPT = scripts[[1]],
  MORTALITY_PREDICTION_MODEL_SCRIPT = scripts[[2]],
  MORTALITY_PREDICTION_PRESENTATION_SCRIPT = scripts[[3]]
)
for (script in scripts) {
  message("Running ", basename(script))
  source(script, encoding = "UTF-8", chdir = FALSE)
}

message("Five-year mortality prediction workflow completed.")
