# ==============================================================================
# Track A compact analysis script 3
# Final longitudinal analysis components.
# ==============================================================================

# ------------------------------------------------------------------------------
# Component: run_measurement_invariance_affective_functional.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(lavaan)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

fit_stats <- function(fit, model_name) {
  idx <- c("cfi", "tli", "rmsea", "srmr", "aic", "bic", "chisq.scaled", "df.scaled", "pvalue.scaled")
  vals <- lavaan::fitMeasures(fit, idx)
  tibble(
    model = model_name,
    cfi = unname(vals["cfi"]),
    tli = unname(vals["tli"]),
    rmsea = unname(vals["rmsea"]),
    srmr = unname(vals["srmr"]),
    aic = unname(vals["aic"]),
    bic = unname(vals["bic"]),
    chisq = unname(vals["chisq.scaled"]),
    df = unname(vals["df.scaled"]),
    pvalue = unname(vals["pvalue.scaled"])
  )
}

compare_fits <- function(fits_tbl) {
  fits_tbl %>%
    mutate(
      delta_cfi = c(NA_real_, diff(cfi)),
      delta_rmsea = c(NA_real_, diff(rmsea)),
      invariance_pass = case_when(
        is.na(delta_cfi) ~ "baseline",
        abs(delta_cfi) <= 0.01 & abs(delta_rmsea) <= 0.015 ~ "pass",
        TRUE ~ "fail"
      )
    )
}

na_codes_to_na <- function(x, extra = c(-1, -3, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

safe_cfa <- function(model, data, group, group.equal = NULL) {
  tryCatch(
    lavaan::cfa(
      model,
      data = data,
      group = group,
      group.equal = group.equal,
      estimator = "MLR",
      std.lv = TRUE,
      missing = "listwise"
    ),
    error = function(e) e
  )
}

is_successful_fit <- function(x) inherits(x, "lavaan")

summarise_fit_result <- function(fit_obj, model_name) {
  if (is_successful_fit(fit_obj)) {
    fit_stats(fit_obj, model_name) %>%
      mutate(status = "ok", error_message = NA_character_)
  } else {
    tibble(
      model = model_name,
      cfi = NA_real_,
      tli = NA_real_,
      rmsea = NA_real_,
      srmr = NA_real_,
      aic = NA_real_,
      bic = NA_real_,
      chisq = NA_real_,
      df = NA_real_,
      pvalue = NA_real_,
      status = "failed",
      error_message = conditionMessage(fit_obj)
    )
  }
}

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "harmonized_affective_functional")
out_dir <- file.path(base_dir, "analysis_ready_core", "measurement_invariance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "pooled_affective_functional_harmonized.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

aff_func <- readRDS(input_file)

analysis_df <- aff_func %>%
  mutate(
    cohort = as.character(cohort),
    subject_id = paste(cohort, id, sep = "::"),
    across(
      c(depression_score, loneliness_score, adl_count, iadl_count, self_rated_health),
      na_codes_to_na
    ),
    indicator_non_missing_n =
      rowSums(!is.na(cbind(depression_score, loneliness_score, adl_count, iadl_count, self_rated_health))),
    affective_indicator_n =
      rowSums(!is.na(cbind(depression_score, loneliness_score))),
    functional_indicator_n =
      rowSums(!is.na(cbind(adl_count, iadl_count, self_rated_health)))
  ) %>%
  filter(harmonization_eligibility == 1) %>%
  arrange(cohort, subject_id, desc(wave), desc(indicator_non_missing_n)) %>%
  group_by(cohort, subject_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    complete_indicator_case = ifelse(
      affective_indicator_n == 2 & functional_indicator_n == 3,
      1,
      0
    )
  )

cohort_snapshot_summary <- analysis_df %>%
  group_by(cohort) %>%
  summarise(
    n_subjects = n(),
    n_complete_indicator_cases = sum(complete_indicator_case == 1, na.rm = TRUE),
    depression_non_missing_n = sum(!is.na(depression_score)),
    loneliness_non_missing_n = sum(!is.na(loneliness_score)),
    adl_non_missing_n = sum(!is.na(adl_count)),
    iadl_non_missing_n = sum(!is.na(iadl_count)),
    srh_non_missing_n = sum(!is.na(self_rated_health)),
    depression_var = stats::var(depression_score, na.rm = TRUE),
    loneliness_var = stats::var(loneliness_score, na.rm = TRUE),
    adl_var = stats::var(adl_count, na.rm = TRUE),
    iadl_var = stats::var(iadl_count, na.rm = TRUE),
    srh_var = stats::var(self_rated_health, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cohort_pass_screen =
      n_complete_indicator_cases >= 100 &
      depression_non_missing_n >= 100 &
      loneliness_non_missing_n >= 100 &
      adl_non_missing_n >= 100 &
      iadl_non_missing_n >= 100 &
      srh_non_missing_n >= 100 &
      !is.na(depression_var) & depression_var > 0 &
      !is.na(loneliness_var) & loneliness_var > 0 &
      !is.na(adl_var) & adl_var > 0 &
      !is.na(iadl_var) & iadl_var > 0 &
      !is.na(srh_var) & srh_var > 0
  )

included_cohorts <- cohort_snapshot_summary %>%
  filter(cohort_pass_screen) %>%
  pull(cohort)

analysis_df_mi <- analysis_df %>%
  filter(
    cohort %in% included_cohorts,
    complete_indicator_case == 1
  )

if (length(included_cohorts) < 2) {
  write_csv(cohort_snapshot_summary, file.path(out_dir, "affective_functional_mi_snapshot_summary.csv"), na = "")
  stop(
    "Fewer than two cohorts passed the measurement invariance screening. ",
    "See affective_functional_mi_snapshot_summary.csv."
  )
}

model_syntax <- "
  affective =~ depression_score + loneliness_score
  functional =~ adl_count + iadl_count + self_rated_health
  affective ~~ functional
"

fit_configural <- safe_cfa(model_syntax, analysis_df_mi, "cohort")
fit_metric <- safe_cfa(model_syntax, analysis_df_mi, "cohort", group.equal = c("loadings"))
fit_scalar <- safe_cfa(model_syntax, analysis_df_mi, "cohort", group.equal = c("loadings", "intercepts"))
fit_strict <- safe_cfa(model_syntax, analysis_df_mi, "cohort", group.equal = c("loadings", "intercepts", "residuals"))

fits_tbl <- bind_rows(
  summarise_fit_result(fit_configural, "configural"),
  summarise_fit_result(fit_metric, "metric"),
  summarise_fit_result(fit_scalar, "scalar"),
  summarise_fit_result(fit_strict, "strict")
)

fits_tbl <- if (all(fits_tbl$status == "ok")) compare_fits(fits_tbl) else fits_tbl %>%
  mutate(delta_cfi = NA_real_, delta_rmsea = NA_real_, invariance_pass = ifelse(status == "ok", "baseline_or_partial", "failed"))

parameter_summary <- bind_rows(
  if (is_successful_fit(fit_configural)) parameterEstimates(fit_configural, standardized = TRUE) %>% mutate(model = "configural") else tibble(),
  if (is_successful_fit(fit_metric)) parameterEstimates(fit_metric, standardized = TRUE) %>% mutate(model = "metric") else tibble(),
  if (is_successful_fit(fit_scalar)) parameterEstimates(fit_scalar, standardized = TRUE) %>% mutate(model = "scalar") else tibble(),
  if (is_successful_fit(fit_strict)) parameterEstimates(fit_strict, standardized = TRUE) %>% mutate(model = "strict") else tibble()
)

saveRDS(analysis_df, file.path(out_dir, "affective_functional_mi_snapshot.rds"))
saveRDS(analysis_df_mi, file.path(out_dir, "affective_functional_mi_analysis_subset.rds"))
saveRDS(
  list(
    configural = fit_configural,
    metric = fit_metric,
    scalar = fit_scalar,
    strict = fit_strict,
    included_cohorts = included_cohorts
  ),
  file.path(out_dir, "affective_functional_mi_fits.rds")
)

writeLines(model_syntax, file.path(out_dir, "affective_functional_mi_model.txt"))
write_csv(cohort_snapshot_summary, file.path(out_dir, "affective_functional_mi_snapshot_summary.csv"), na = "")
write_csv(fits_tbl, file.path(out_dir, "affective_functional_mi_fit_summary.csv"), na = "")
write_csv(parameter_summary, file.path(out_dir, "affective_functional_mi_parameter_summary.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: run_measurement_invariance_basic_cognition.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "harmonization_input")
out_dir <- file.path(base_dir, "analysis_ready_core", "measurement_invariance")
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

basic_cog <- readRDS(input_file) %>%
  mutate(
    cohort = as.character(cohort),
    subject_id = paste(cohort, id, sep = "::"),
    across(
      c(age, sex_female, education_years, education_level_clean,
        marital_status_clean, cognition_score, cognition_harmonized),
      na_codes_to_na
    ),
    completeness_n =
      rowSums(!is.na(cbind(
        cognition_score, age, sex_female,
        education_level_clean, marital_status_clean
      )))
  ) %>%
  arrange(cohort, subject_id, desc(wave), desc(completeness_n)) %>%
  group_by(cohort, subject_id) %>%
  slice(1) %>%
  ungroup()

snapshot_summary <- basic_cog %>%
  group_by(cohort) %>%
  summarise(
    n_subjects = n(),
    cognition_non_missing_n = sum(!is.na(cognition_score)),
    harmonized_non_missing_n = sum(!is.na(cognition_harmonized)),
    eligibility_n = sum(cognition_harmonization_eligibility == 1, na.rm = TRUE),
    complete_case_n = sum(cognition_complete_case == 1, na.rm = TRUE),
    age_non_missing_n = sum(!is.na(age)),
    sex_non_missing_n = sum(!is.na(sex_female)),
    education_non_missing_n = sum(!is.na(education_level_clean)),
    marital_non_missing_n = sum(!is.na(marital_status_clean)),
    .groups = "drop"
  )

observed_score_summary <- basic_cog %>%
  group_by(cohort) %>%
  summarise(
    n_subjects = n(),
    cognition_mean = mean(cognition_score, na.rm = TRUE),
    cognition_sd = stats::sd(cognition_score, na.rm = TRUE),
    cognition_min = min(cognition_score, na.rm = TRUE),
    cognition_p25 = stats::quantile(cognition_score, 0.25, na.rm = TRUE),
    cognition_median = stats::median(cognition_score, na.rm = TRUE),
    cognition_p75 = stats::quantile(cognition_score, 0.75, na.rm = TRUE),
    cognition_max = max(cognition_score, na.rm = TRUE),
    harmonized_mean = mean(cognition_harmonized, na.rm = TRUE),
    harmonized_sd = stats::sd(cognition_harmonized, na.rm = TRUE),
    n_unique_raw = dplyr::n_distinct(cognition_score[!is.na(cognition_score)]),
    source_profile_example = dplyr::first(na.omit(cognition_source_profile)),
    .groups = "drop"
  )

source_profile_summary <- basic_cog %>%
  count(cohort, cognition_source_profile, name = "n_rows") %>%
  arrange(cohort, desc(n_rows))

mi_status <- tibble(
  method = "single_indicator_observed_score_fallback",
  mg_cfa_run = FALSE,
  reason = "Current v1 cognition input contains one harmonized cognition score per row rather than a shared multi-indicator item set across all cohorts.",
  recommendation = "Use this output as comparability diagnostics only. Run formal multi-group CFA after constructing a shared item-level cognition battery or a cohort-aligned component set."
)

writeLines(
  c(
    "basic_cognition measurement invariance v1",
    "status: fallback diagnostics only",
    "formal multi-indicator MG-CFA not run",
    "reason: the current harmonization input contains one observed cognition score rather than a shared multi-indicator battery across cohorts",
    "next step: use cohort-aligned raw cognition components if a formal invariance model is required"
  ),
  file.path(out_dir, "basic_cognition_mi_method_note.txt")
)

saveRDS(basic_cog, file.path(out_dir, "basic_cognition_mi_snapshot.rds"))
write_csv(snapshot_summary, file.path(out_dir, "basic_cognition_mi_snapshot_summary.csv"), na = "")
write_csv(observed_score_summary, file.path(out_dir, "basic_cognition_observed_score_summary.csv"), na = "")
write_csv(source_profile_summary, file.path(out_dir, "basic_cognition_source_profile_summary.csv"), na = "")
write_csv(mi_status, file.path(out_dir, "basic_cognition_mi_status.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})


