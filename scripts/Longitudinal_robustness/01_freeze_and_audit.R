## Frozen input and interval audit.
## Run first. Creates a lineage manifest and audits current five-cohort inputs.

options(stringsAsFactors = FALSE)
base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir)) stop("Set DCV_BASE_DIR before running this script.")
if (!dir.exists(base_dir)) stop("DCV_BASE_DIR does not exist: ", base_dir)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
out_dir <- file.path(base_dir, "robustness_analysis", "00_freeze_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("dplyr", "readr", "tibble", "purrr")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install before running: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble); library(purrr)})

input_files <- c(
  longitudinal_rds = file.path(base_dir, "analysis_ready_core", "longitudinal_models_v1", "affective_functional_longitudinal_input_v1.rds"),
  longitudinal_effects = file.path(base_dir, "analysis_ready_core", "longitudinal_models_v1", "affective_functional_longitudinal_effects_v1.csv"),
  genomic_covstruc = file.path(base_dir, "analysis_ready_core", "trackB_v1", "factor_gwas_qsnp_v1", "trackB_afffunc_ldsc_covstruc_v1.rds"),
  genomic_sumstats = file.path(base_dir, "analysis_ready_core", "trackB_v1", "factor_gwas_qsnp_v1", "trackB_afffunc_sumstats_v1.rds"),
  genomic_model_fit = file.path(base_dir, "analysis_ready_core", "trackB_v1", "genomicsem_competing_models_v4", "trackB_competing_model_fit_v4.csv")
)
manifest <- tibble(input_name = names(input_files), path = unname(input_files), exists = file.exists(input_files),
  size_bytes = ifelse(file.exists(input_files), file.info(input_files)$size, NA_real_),
  modified = ifelse(file.exists(input_files), as.character(file.info(input_files)$mtime), NA_character_),
  md5 = vapply(input_files, function(x) if (!file.exists(x)) NA_character_ else unname(tools::md5sum(x)), character(1)))
write_csv(manifest, file.path(out_dir, "input_manifest.csv"), na = "")

dat <- readRDS(input_files[["longitudinal_rds"]])
required <- c("cohort", "subject_id", "year", "wave", "affective_harmonized_v1", "functional_harmonized_v1",
              "age_baseline", "sex_baseline", "education_baseline")
stopifnot(all(required %in% names(dat)))
dat <- dat %>% mutate(joint = !is.na(affective_harmonized_v1) & !is.na(functional_harmonized_v1))
repeat_joint <- dat %>% filter(joint) %>% count(cohort, subject_id, name = "n_joint") %>%
  group_by(cohort) %>% summarise(persons_ge2_joint = sum(n_joint >= 2), .groups = "drop")
sample_flow <- dat %>% group_by(cohort) %>% summarise(rows = n(), persons = n_distinct(subject_id),
  affective_rows = sum(!is.na(affective_harmonized_v1)), functional_rows = sum(!is.na(functional_harmonized_v1)),
  joint_rows = sum(joint), joint_persons = n_distinct(subject_id[joint]), .groups = "drop") %>%
  left_join(repeat_joint, by = "cohort")
write_csv(sample_flow, file.path(out_dir, "sample_flow_by_cohort.csv"), na = "")

intervals <- dat %>% arrange(subject_id, year, wave) %>% group_by(subject_id) %>%
  mutate(delta_year = year - lag(year)) %>% ungroup() %>% filter(is.finite(delta_year), delta_year > 0)
interval_audit <- intervals %>% group_by(cohort) %>% summarise(n_transitions = n(), n_persons = n_distinct(subject_id),
  median_interval = median(delta_year, na.rm = TRUE), q025 = quantile(delta_year, .025, na.rm = TRUE),
  q975 = quantile(delta_year, .975, na.rm = TRUE), pct_lt_1 = mean(delta_year < 1, na.rm = TRUE),
  pct_1to3 = mean(delta_year >= 1 & delta_year <= 3, na.rm = TRUE), pct_gt_3 = mean(delta_year > 3, na.rm = TRUE), .groups = "drop")
write_csv(interval_audit, file.path(out_dir, "interval_audit.csv"), na = "")

writeLines(c(
  "Retain this manifest unchanged before interpreting any result.",
  "Primary estimand: association in observed repeated measures, conditional on baseline age, sex and education.",
  "Death is not treated as ordinary missingness; no mortality/IPCW claim is allowed without a verified death variable.",
  "The external GWAS input family is fixed by the recorded manifest.",
  "Qualification-grid P values are not new primary hypotheses."
), file.path(out_dir, "FREEZE_README.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Freeze/audit complete: ", out_dir)
