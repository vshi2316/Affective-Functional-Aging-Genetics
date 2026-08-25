#!/usr/bin/env Rscript
# Entry point for the prespecified ADNI factor-PGS CDR-SB validation workflow.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(DCV_PROJECT_DIR = project_dir)

config_file <- file.path(project_dir, "config", "config.R")
if (file.exists(config_file)) source(config_file, encoding = "UTF-8")

scripts <- file.path(
  project_dir, "scripts", "ADNI",
  c(
    "01_build_locked_cdrsb_dataset.R",
    "02_fit_factor_pgs_cdrsb_models.R",
    "03_build_adni_validation_outputs.R"
  )
)
missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts)) {
  stop("Missing ADNI-analysis script(s): ", paste(basename(missing_scripts), collapse = ", "))
}

for (script in scripts) {
  message("Running ", basename(script))
  source(script, encoding = "UTF-8", chdir = FALSE)
}

output_root <- Sys.getenv(
  "ADNI_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "adni_validation")
)
writeLines(capture.output(sessionInfo()), file.path(output_root, "adni_sessionInfo.txt"), useBytes = TRUE)
message("ADNI factor-PGS CDR-SB validation completed.")
