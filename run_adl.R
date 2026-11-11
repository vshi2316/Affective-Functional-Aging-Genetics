#!/usr/bin/env Rscript
# Entry point for first-recorded basic ADL difficulty at the next visit.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- normalizePath(
  Sys.getenv("DCV_PROJECT_DIR", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
Sys.setenv(DCV_PROJECT_DIR = project_dir)

config_file <- file.path(project_dir, "config", "config.R")
if (file.exists(config_file)) source(config_file, encoding = "UTF-8")

script <- file.path(project_dir, "scripts", "ADL", "01_build_first_recorded_adl_evidence.R")
if (!file.exists(script)) stop("Missing ADL script: ", script)

message("Running ", basename(script))
source(script, encoding = "UTF-8", chdir = FALSE)
message("First-recorded ADL analysis completed.")
