#!/usr/bin/env Rscript
# Run the prespecified clinical translation and external recalibration modules.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- normalizePath(
  Sys.getenv("DCV_PROJECT_DIR", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
Sys.setenv(DCV_PROJECT_DIR = project_dir)

config_file <- file.path(project_dir, "config", "config.R")
if (file.exists(config_file)) source(config_file, encoding = "UTF-8")

scripts <- file.path(
  project_dir,
  "scripts",
  "Mortality_prediction",
  c("05_clinical_risk_prior1_continuous.R", "06_external_recalibration.R")
)
missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts)) {
  stop("Missing clinical-translation script(s): ", paste(basename(missing_scripts), collapse = ", "))
}

for (script in scripts) {
  message("Running ", basename(script))
  source(script, encoding = "UTF-8", chdir = FALSE)
}

message("Clinical translation and external recalibration completed.")
