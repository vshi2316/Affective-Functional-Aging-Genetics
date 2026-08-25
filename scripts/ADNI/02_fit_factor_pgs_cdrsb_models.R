#!/usr/bin/env Rscript
# Fit prespecified primary and joint factor-PGS models for longitudinal CDR-SB.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

required <- c("dplyr", "readr", "tibble", "lme4", "lmerTest", "broom.mixed")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))

project_dir <- normalizePath(Sys.getenv("DCV_PROJECT_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
output_root <- Sys.getenv("ADNI_OUTPUT_ROOT", unset = file.path(project_dir, "results", "adni_validation"))
input_file <- file.path(output_root, "dataset", "03_analysis_dataset_locked.csv")
out_dir <- file.path(output_root, "models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Run 01_build_locked_cdrsb_dataset.R first.")

d <- readr::read_csv(input_file, show_col_types = FALSE) |>
  dplyr::filter(!is.na(.data$baseline_APOE4)) |>
  dplyr::mutate(PTID = as.character(.data$PTID), sex = factor(.data$sex), baseline_diagnosis = factor(.data$baseline_diagnosis)) |>
  dplyr::filter(dplyr::if_all(c("PGS_affective_z", "PGS_functional_z", "time_years",
                                "baseline_age_c", "education_c", dplyr::starts_with("PC")), is.finite))
if (dplyr::n_distinct(d$PTID) < 200L || nrow(d) < 500L) stop("Insufficient complete APOE4-adjusted sample.")

fit_lmm <- function(formula, model_name) {
  fit <- lmerTest::lmer(
    formula, data = d, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
  )
  fixed <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) |>
    dplyr::mutate(
      model = model_name, n_rows = nrow(d), n_people = dplyr::n_distinct(d$PTID),
      n_people_2plus_visits = sum(table(d$PTID) >= 2L), singular = lme4::isSingular(fit, tol = 1e-4)
    )
  list(fit = fit, fixed = fixed)
}

covariates <- "+ baseline_age_c + sex + education_c + baseline_diagnosis + baseline_APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (1 | PTID)"
affective <- fit_lmm(stats::as.formula(paste("CDRSB ~ PGS_affective_z * time_years", covariates)), "single_affective")
functional <- fit_lmm(stats::as.formula(paste("CDRSB ~ PGS_functional_z * time_years", covariates)), "single_functional")
joint <- fit_lmm(stats::as.formula(paste("CDRSB ~ PGS_affective_z * time_years + PGS_functional_z * time_years", covariates)), "joint_scores")

single <- dplyr::bind_rows(affective$fixed, functional$fixed)
primary <- single |>
  dplyr::filter(.data$term %in% c("PGS_affective_z:time_years", "PGS_functional_z:time_years")) |>
  dplyr::mutate(p_fdr_two_primary_terms = stats::p.adjust(.data$p.value, method = "BH"))
joint_primary <- joint$fixed |>
  dplyr::filter(.data$term %in% c("PGS_affective_z:time_years", "PGS_functional_z:time_years"))
readr::write_csv(single, file.path(out_dir, "01_single_score_model_coefficients.csv"), na = "")
readr::write_csv(primary, file.path(out_dir, "02_primary_interaction_summary.csv"), na = "")
readr::write_csv(joint$fixed, file.path(out_dir, "03_joint_score_model_coefficients.csv"), na = "")
readr::write_csv(joint_primary, file.path(out_dir, "04_joint_interaction_summary.csv"), na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
message("Primary and joint ADNI CDR-SB models written to: ", out_dir)
