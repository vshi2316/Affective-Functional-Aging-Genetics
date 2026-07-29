#!/usr/bin/env Rscript
# One-way dispatcher for the frozen dynamic mortality bridge.
# Set controlled-access input paths in config/config.R before running.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(DCV_PROJECT_DIR = project_dir)

config_file <- file.path(project_dir, "config", "config.R")
if (file.exists(config_file)) {
  source(config_file, encoding = "UTF-8")
}

scripts <- c(
  file.path(project_dir, "scripts", "Mortality", "01_build_dynamic_candidate.R"),
  file.path(project_dir, "scripts", "Mortality", "02_fit_dynamic_models.R"),
  file.path(project_dir, "scripts", "Mortality", "03_build_figure_and_source_data.R")
)
missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts)) {
  stop("Missing mortality script(s): ", paste(basename(missing_scripts), collapse = ", "))
}

for (script in scripts) {
  message("Running ", basename(script))
  source(script, encoding = "UTF-8", chdir = FALSE)
}

message("Dynamic mortality bridge completed. Frozen locks prevent result-driven reruns.")
