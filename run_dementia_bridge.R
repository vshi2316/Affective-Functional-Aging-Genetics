#!/usr/bin/env Rscript
# Entry point for the dynamic dementia-related outcome analysis.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(DCV_PROJECT_DIR = project_dir)

config_file <- file.path(project_dir, "config", "config.R")
if (file.exists(config_file)) {
  source(config_file, encoding = "UTF-8")
}

scripts <- file.path(
  project_dir,
  "scripts",
  "Dementia",
  c(
    "01_build_event_dataset.R",
    "02_fit_competing_event_models.R",
    "03_build_figure_and_source_data.R"
  )
)

missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts)) {
  stop("Missing dementia-analysis script(s): ", paste(basename(missing_scripts), collapse = ", "))
}

for (script in scripts) {
  message("Running ", basename(script))
  source(script, encoding = "UTF-8", chdir = FALSE)
}

output_root <- Sys.getenv(
  "DEMENTIA_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "dementia_bridge")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "dementia_sessionInfo.txt"),
  useBytes = TRUE
)

message("Dynamic dementia-related outcome analysis completed.")
