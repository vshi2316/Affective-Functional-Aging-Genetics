# ==============================================================================
# Track A compact analysis script 4
# Final longitudinal analysis components.
# ==============================================================================

# ------------------------------------------------------------------------------
# Component: run_longitudinal_models_affective_functional.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it before running this script.")
  }
}

check_pkg("nlme")

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "harmonized_affective_functional")
out_dir <- file.path(base_dir, "analysis_ready_core", "longitudinal_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "pooled_affective_functional_harmonized.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

na_codes_to_na <- function(x, extra = c(-1, -3, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

prepare_data <- function(df) {
  df %>%
    mutate(
      cohort = as.character(cohort),
      subject_id = paste(cohort, id, sep = "::"),
      across(
        c(
          year, age, sex_female, education_level_clean,
          affective_harmonized, functional_harmonized
        ),
        na_codes_to_na
      )
    ) %>%
    arrange(cohort, subject_id, year, wave) %>%
    group_by(cohort, subject_id) %>%
    mutate(
      baseline_year = min(year, na.rm = TRUE),
      baseline_year = ifelse(is.infinite(baseline_year), NA_real_, baseline_year),
      time_since_baseline_years = year - baseline_year,
      age_baseline = first(age[!is.na(age)]),
      education_baseline = first(education_level_clean[!is.na(education_level_clean)]),
      sex_baseline = first(sex_female[!is.na(sex_female)])
    ) %>%
    ungroup()
}

extract_fixed_effects <- function(fit, cohort_name, model_id, outcome_name, n_obs, n_subjects) {
  coef_tbl <- as.data.frame(summary(fit)$tTable)
  coef_tbl$term <- rownames(coef_tbl)
  rownames(coef_tbl) <- NULL
  names(coef_tbl) <- sub(" ", "_", names(coef_tbl), fixed = TRUE)

  tibble(
    cohort = cohort_name,
    model_id = model_id,
    outcome = outcome_name,
    term = coef_tbl$term,
    estimate = coef_tbl$Value,
    std.error = coef_tbl$Std.Error,
    statistic = coef_tbl$t.value,
    p.value = coef_tbl$p.value,
    n_obs = n_obs,
    n_subjects = n_subjects,
    status = "ok",
    error_message = NA_character_
  )
}

fit_lmm <- function(data, formula_text) {
  tryCatch(
    nlme::lme(
      fixed = as.formula(formula_text),
      random = ~ 1 | subject_id,
      data = data,
      method = "ML",
      na.action = na.omit,
      control = nlme::lmeControl(
        returnObject = TRUE,
        msMaxIter = 200,
        opt = "optim"
      )
    ),
    error = function(e) e
  )
}

run_model_for_cohort <- function(df, cohort_name, model_id, outcome_name, formula_text, outcome_vars) {
  dat <- df %>%
    filter(cohort == cohort_name) %>%
    filter(if_all(all_of(c("time_since_baseline_years", "age_baseline", "sex_baseline", "education_baseline")), ~ !is.na(.x))) %>%
    filter(if_all(all_of(outcome_vars), ~ !is.na(.x)))

  subject_counts <- dat %>%
    count(subject_id, name = "n_rows")

  dat <- dat %>%
    inner_join(subject_counts %>% filter(n_rows >= 2), by = "subject_id")

  n_obs <- nrow(dat)
  n_subjects <- dplyr::n_distinct(dat$subject_id)

  if (n_obs < 200 || n_subjects < 100) {
    return(tibble(
      cohort = cohort_name,
      model_id = model_id,
      outcome = outcome_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_obs = n_obs,
      n_subjects = n_subjects,
      status = "skipped_low_sample",
      error_message = "Insufficient repeated-measures sample after filtering."
    ))
  }

  fit <- fit_lmm(dat, formula_text)

  if (inherits(fit, "error")) {
    return(tibble(
      cohort = cohort_name,
      model_id = model_id,
      outcome = outcome_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_obs = n_obs,
      n_subjects = n_subjects,
      status = "failed",
      error_message = conditionMessage(fit)
    ))
  }

  extract_fixed_effects(fit, cohort_name, model_id, outcome_name, n_obs, n_subjects)
}

model_specs <- tribble(
  ~model_id, ~outcome, ~formula_text, ~required_vars,
  "affective_trajectory", "affective_harmonized",
  "affective_harmonized ~ time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "affective_harmonized",
  "functional_trajectory", "functional_harmonized",
  "functional_harmonized ~ time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "functional_harmonized",
  "functional_on_affective", "functional_harmonized",
  "functional_harmonized ~ affective_harmonized + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "functional_harmonized|affective_harmonized",
  "affective_on_functional", "affective_harmonized",
  "affective_harmonized ~ functional_harmonized + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "affective_harmonized|functional_harmonized"
)

dat <- prepare_data(readRDS(input_file))
cohorts <- sort(unique(dat$cohort))

effects_tbl <- purrr::pmap_dfr(model_specs, function(model_id, outcome, formula_text, required_vars) {
  req <- strsplit(required_vars, "\\|")[[1]]
  bind_rows(lapply(cohorts, function(cc) {
    run_model_for_cohort(dat, cc, model_id, outcome, formula_text, req)
  }))
})

cohort_model_summary <- effects_tbl %>%
  group_by(cohort, model_id, outcome, status) %>%
  summarise(
    n_terms = sum(!is.na(term)),
    n_obs = dplyr::first(n_obs),
    n_subjects = dplyr::first(n_subjects),
    error_message = dplyr::first(error_message),
    .groups = "drop"
  )

write_csv(effects_tbl, file.path(out_dir, "affective_functional_longitudinal_effects.csv"), na = "")
write_csv(cohort_model_summary, file.path(out_dir, "affective_functional_longitudinal_model_summary.csv"), na = "")
saveRDS(dat, file.path(out_dir, "affective_functional_longitudinal_input.rds"))

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: run_longitudinal_models_basic_cognition.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it before running this script.")
  }
}

check_pkg("nlme")

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "harmonization_input")
out_dir <- file.path(base_dir, "analysis_ready_core", "longitudinal_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "pooled_basic_cognition_harmonization_input.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

na_codes_to_na <- function(x, extra = c(-1, -3, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

prepare_data <- function(df) {
  df %>%
    mutate(
      cohort = as.character(cohort),
      subject_id = paste(cohort, id, sep = "::"),
      across(
        c(
          year, age, sex_female, education_level_clean,
          cognition_harmonized
        ),
        na_codes_to_na
      )
    ) %>%
    arrange(cohort, subject_id, year, wave) %>%
    group_by(cohort, subject_id) %>%
    mutate(
      baseline_year = min(year, na.rm = TRUE),
      baseline_year = ifelse(is.infinite(baseline_year), NA_real_, baseline_year),
      time_since_baseline_years = year - baseline_year,
      age_baseline = first(age[!is.na(age)]),
      education_baseline = first(education_level_clean[!is.na(education_level_clean)]),
      sex_baseline = first(sex_female[!is.na(sex_female)])
    ) %>%
    ungroup()
}

extract_fixed_effects <- function(fit, cohort_name, model_id, outcome_name, n_obs, n_subjects) {
  coef_tbl <- as.data.frame(summary(fit)$tTable)
  coef_tbl$term <- rownames(coef_tbl)
  rownames(coef_tbl) <- NULL
  names(coef_tbl) <- sub(" ", "_", names(coef_tbl), fixed = TRUE)

  tibble(
    cohort = cohort_name,
    model_id = model_id,
    outcome = outcome_name,
    term = coef_tbl$term,
    estimate = coef_tbl$Value,
    std.error = coef_tbl$Std.Error,
    statistic = coef_tbl$t.value,
    p.value = coef_tbl$p.value,
    n_obs = n_obs,
    n_subjects = n_subjects,
    status = "ok",
    error_message = NA_character_
  )
}

fit_lmm <- function(data, formula_text) {
  tryCatch(
    nlme::lme(
      fixed = as.formula(formula_text),
      random = ~ 1 | subject_id,
      data = data,
      method = "ML",
      na.action = na.omit,
      control = nlme::lmeControl(
        returnObject = TRUE,
        msMaxIter = 200,
        opt = "optim"
      )
    ),
    error = function(e) e
  )
}

run_model_for_cohort <- function(df, cohort_name, model_id, outcome_name, formula_text, outcome_vars) {
  dat <- df %>%
    filter(cohort == cohort_name) %>%
    filter(if_all(all_of(c("time_since_baseline_years", "age_baseline", "sex_baseline", "education_baseline")), ~ !is.na(.x))) %>%
    filter(if_all(all_of(outcome_vars), ~ !is.na(.x)))

  subject_counts <- dat %>%
    count(subject_id, name = "n_rows")

  dat <- dat %>%
    inner_join(subject_counts %>% filter(n_rows >= 2), by = "subject_id")

  n_obs <- nrow(dat)
  n_subjects <- dplyr::n_distinct(dat$subject_id)

  if (n_obs < 200 || n_subjects < 100) {
    return(tibble(
      cohort = cohort_name,
      model_id = model_id,
      outcome = outcome_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_obs = n_obs,
      n_subjects = n_subjects,
      status = "skipped_low_sample",
      error_message = "Insufficient repeated-measures sample after filtering."
    ))
  }

  fit <- fit_lmm(dat, formula_text)

  if (inherits(fit, "error")) {
    return(tibble(
      cohort = cohort_name,
      model_id = model_id,
      outcome = outcome_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_obs = n_obs,
      n_subjects = n_subjects,
      status = "failed",
      error_message = conditionMessage(fit)
    ))
  }

  extract_fixed_effects(fit, cohort_name, model_id, outcome_name, n_obs, n_subjects)
}

model_specs <- tribble(
  ~model_id, ~outcome, ~formula_text, ~required_vars,
  "cognition_trajectory", "cognition_harmonized",
  "cognition_harmonized ~ time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "cognition_harmonized"
)

dat <- prepare_data(readRDS(input_file))
cohorts <- sort(unique(dat$cohort))

effects_tbl <- purrr::pmap_dfr(model_specs, function(model_id, outcome, formula_text, required_vars) {
  req <- strsplit(required_vars, "\\|")[[1]]
  bind_rows(lapply(cohorts, function(cc) {
    run_model_for_cohort(dat, cc, model_id, outcome, formula_text, req)
  }))
})

cohort_model_summary <- effects_tbl %>%
  group_by(cohort, model_id, outcome, status) %>%
  summarise(
    n_terms = sum(!is.na(term)),
    n_obs = dplyr::first(n_obs),
    n_subjects = dplyr::first(n_subjects),
    error_message = dplyr::first(error_message),
    .groups = "drop"
  )

write_csv(effects_tbl, file.path(out_dir, "basic_cognition_longitudinal_effects.csv"), na = "")
write_csv(cohort_model_summary, file.path(out_dir, "basic_cognition_longitudinal_model_summary.csv"), na = "")
saveRDS(dat, file.path(out_dir, "basic_cognition_longitudinal_input.rds"))

message("Done.")
message("Outputs written to: ", out_dir)

})


