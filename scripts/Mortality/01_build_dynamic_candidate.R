#!/usr/bin/env Rscript
# Build the exact-date, history-anchored dynamic mortality candidate.
# This script performs deterministic data assembly and feasibility gates only.
# It does not estimate exposure–mortality associations.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

long_path <- Sys.getenv(
  "MORTALITY_LONGITUDINAL_INPUT",
  unset = file.path(project_dir, "data", "controlled", "longitudinal_input.rds")
)
date_path <- Sys.getenv(
  "MORTALITY_EXACT_VISIT_DATE_INPUT",
  unset = file.path(project_dir, "data", "controlled", "exact_visit_dates.csv")
)
proxy_path <- Sys.getenv(
  "MORTALITY_WAVE_STATUS_INPUT",
  unset = file.path(project_dir, "data", "controlled", "wave_interview_status.csv")
)
mort_main_path <- Sys.getenv(
  "MORTALITY_ENDPOINT_MAIN_INPUT",
  unset = file.path(project_dir, "data", "controlled", "mortality_endpoints_main.rds")
)
mort_support_path <- Sys.getenv(
  "MORTALITY_ENDPOINT_SUPPORT_INPUT",
  unset = file.path(project_dir, "data", "controlled", "mortality_endpoints_support.rds")
)

output_root <- Sys.getenv(
  "MORTALITY_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_bridge")
)
out_dir <- Sys.getenv(
  "MORTALITY_CANDIDATE_OUT",
  unset = file.path(output_root, "candidate")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lock <- file.path(out_dir, "dynamic_candidate_lock.txt")
if (file.exists(lock)) {
  stop("The dynamic candidate is already locked; do not overwrite it.")
}

required_inputs <- c(
  longitudinal = long_path,
  exact_visit_dates = date_path,
  wave_interview_status = proxy_path,
  mortality_endpoints_main = mort_main_path,
  mortality_endpoints_support = mort_support_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required controlled-access input(s): ",
       paste(names(missing_inputs), collapse = ", "))
}

pkgs <- c("dplyr", "tidyr", "readr", "tibble", "survival")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install before running: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(survival)
})

as_num <- function(x) suppressWarnings(as.numeric(x))
as_chr <- function(x) {
  y <- trimws(as.character(x)); y[y %in% c("", "NA", "NaN", "-1", "-8", "-9")] <- NA_character_; y
}
normalize_id <- function(x) {
  y <- toupper(gsub("\\s+", "", as_chr(x))); y <- sub("\\.0+$", "", y)
  z <- !is.na(y) & grepl("^[0-9]+$", y); y[z] <- sub("^0+", "", y[z]); y[y == ""] <- NA_character_; y
}
month_date <- function(year, month) {
  yy <- as_num(year); mm <- as_num(month)
  yy[yy < 1900 | yy > 2100] <- NA_real_; mm[mm < 1 | mm > 12] <- NA_real_
  mm[is.na(mm) & !is.na(yy)] <- 7
  out <- rep(as.Date(NA), length(yy)); ok <- is.finite(yy) & is.finite(mm)
  out[ok] <- as.Date(sprintf("%04d-%02d-15", as.integer(yy[ok]), as.integer(mm[ok]))); out
}

roles <- tribble(
  ~cohort, ~cohort_role,
  "HRS", "primary", "MHAS", "replication", "SHARE", "replication",
  "CHARLS", "supportive", "ELSA", "sensitivity_only"
)
main_cohorts <- c("HRS", "MHAS", "SHARE")

dates <- read_csv(date_path, show_col_types = FALSE, progress = FALSE) %>% transmute(
  cohort = as.character(cohort), id_norm = normalize_id(id_norm), wave = as_num(wave),
  interview_year = as_num(interview_year_raw), interview_month = as_num(interview_month_raw)
) %>% group_by(cohort, id_norm, wave) %>% summarise(
  n_year = n_distinct(interview_year, na.rm = TRUE), n_month = n_distinct(interview_month, na.rm = TRUE),
  interview_year = first(interview_year[is.finite(interview_year)], default = NA_real_),
  interview_month = first(interview_month[is.finite(interview_month)], default = NA_real_),
  .groups = "drop"
) %>% mutate(interview_date = month_date(interview_year, interview_month))
if (any(dates$n_year > 1 | dates$n_month > 1)) stop("Conflicting exact visit dates.")

proxy <- read_csv(proxy_path, show_col_types = FALSE, progress = FALSE) %>% transmute(
  cohort = as.character(cohort), id_norm = normalize_id(id_norm), wave = as_num(wave),
  proxy_status = as.character(proxy_status)
) %>% group_by(cohort, id_norm, wave) %>% summarise(
  proxy_status = case_when(
    any(proxy_status == "proxy_or_exit") ~ "proxy_or_exit",
    any(proxy_status == "direct") ~ "direct", TRUE ~ "unknown"
  ), .groups = "drop"
)

mort <- bind_rows(
  readRDS(mort_main_path) %>% filter(as.character(cohort) != "CHARLS"),
  readRDS(mort_support_path)
) %>% transmute(
  cohort = as.character(cohort), subject_id = as.character(subject_id),
  baseline_date = as.Date(baseline_date_direct), baseline_age = as_num(baseline_age_direct),
  exit_date = as.Date(primary_exit_date), death_event = as.logical(primary_event),
  endpoint_eligible = as.logical(primary_eligible),
  sex_baseline = as_num(sex_baseline), education_baseline = as_num(education_baseline)
)
if (anyDuplicated(mort %>% select(cohort, subject_id))) stop("Duplicate endpoint subjects.")

long0 <- readRDS(long_path)
need <- c("cohort", "id", "subject_id", "wave", "affective_harmonized_v1", "functional_harmonized_v1")
if (!all(need %in% names(long0))) stop("Longitudinal input lacks required fields.")

visits <- long0 %>% transmute(
  cohort = as.character(cohort), id_norm = normalize_id(id), subject_id = as.character(subject_id),
  wave = as_num(wave), affective_raw = as_num(affective_harmonized_v1),
  functional_raw = as_num(functional_harmonized_v1)
) %>%
  inner_join(mort, by = c("cohort", "subject_id")) %>%
  left_join(dates %>% select(cohort, id_norm, wave, interview_date),
            by = c("cohort", "id_norm", "wave")) %>%
  left_join(proxy, by = c("cohort", "id_norm", "wave")) %>%
  mutate(
    proxy_status = if_else(cohort == "CHARLS" & wave == 1, "direct", proxy_status),
    explicit_direct = proxy_status == "direct",
    joint_complete = is.finite(affective_raw) & is.finite(functional_raw)
  ) %>%
  filter(
    endpoint_eligible, explicit_direct, joint_complete,
    !is.na(interview_date), !is.na(exit_date), interview_date < exit_date,
    is.finite(baseline_age), !is.na(baseline_date),
    is.finite(sex_baseline), is.finite(education_baseline)
  ) %>%
  left_join(roles, by = "cohort") %>%
  arrange(cohort, subject_id, interview_date, wave)

dup <- visits %>% count(cohort, subject_id, interview_date, name = "n") %>% filter(n > 1)
write_csv(dup, file.path(out_dir, "duplicate_exact_visit_dates.csv"), na = "")
if (nrow(dup)) stop("Duplicate exact subject-visit dates; no candidate written.")

visits <- visits %>% group_by(cohort, subject_id) %>%
  mutate(joint_visit_number = row_number()) %>% ungroup()

scaling <- visits %>% filter(joint_visit_number == 1) %>% group_by(cohort) %>% summarise(
  first_visit_subjects = n(),
  affective_first_mean = mean(affective_raw), affective_first_sd = sd(affective_raw),
  functional_first_mean = mean(functional_raw), functional_first_sd = sd(functional_raw),
  education_mean = mean(education_baseline), education_sd = sd(education_baseline),
  dynamic_entry_year_mean = mean(as.numeric(format(interview_date, "%Y"))),
  .groups = "drop"
)
if (any(!is.finite(unlist(scaling %>% select(ends_with("_sd"))))) ||
    any(unlist(scaling %>% select(ends_with("_sd"))) <= 0)) stop("Invalid scaling SD.")
write_csv(scaling, file.path(out_dir, "fixed_scaling.csv"), na = "")

v <- visits %>% left_join(scaling, by = "cohort") %>%
  group_by(cohort, subject_id) %>% arrange(interview_date, wave, .by_group = TRUE) %>%
  mutate(
    prior_n = row_number() - 1L,
    H_A_raw = lag(cummean(affective_raw)), H_F_raw = lag(cummean(functional_raw)),
    D_A_raw = affective_raw - H_A_raw, D_F_raw = functional_raw - H_F_raw,
    H_A = (H_A_raw - affective_first_mean) / affective_first_sd,
    D_A = D_A_raw / affective_first_sd,
    H_F = (H_F_raw - functional_first_mean) / functional_first_sd,
    D_F = D_F_raw / functional_first_sd,
    next_visit_date = lead(interview_date),
    lag1_H_A = lag(H_A), lag1_D_A = lag(D_A), lag1_H_F = lag(H_F), lag1_D_F = lag(D_F)
  ) %>% ungroup() %>% filter(prior_n >= 1L) %>%
  mutate(
    interval_start_date = interview_date,
    interval_end_date = pmin(next_visit_date, exit_date, na.rm = TRUE),
    interval_event = death_event & (is.na(next_visit_date) | exit_date < next_visit_date),
    interval_days = as.numeric(interval_end_date - interval_start_date),
    interval_start_age = baseline_age + as.numeric(interval_start_date - baseline_date) / 365.25,
    interval_end_age = baseline_age + as.numeric(interval_end_date - baseline_date) / 365.25,
    calendar_year_c = as.numeric(format(interval_start_date, "%Y")) - dynamic_entry_year_mean,
    education_z = (education_baseline - education_mean) / education_sd,
    age_band_start = factor(case_when(
      interval_start_age < 65 ~ "lt65", interval_start_age < 75 ~ "65_74",
      interval_start_age < 85 ~ "75_84", TRUE ~ "ge85"
    ), levels = c("lt65", "65_74", "75_84", "ge85")),
    death_within_12m = interval_event & interval_days <= 365.25,
    lag1_available = if_all(c(lag1_H_A, lag1_D_A, lag1_H_F, lag1_D_F), is.finite)
  )

candidate <- v %>% filter(
  interval_days > 0, interval_end_date <= exit_date,
  if_all(c(H_A, D_A, H_F, D_F, interval_start_age, interval_end_age,
           education_z, calendar_year_c), is.finite)
) %>% select(
  cohort, cohort_role, subject_id, wave, prior_n,
  interval_start_date, interval_end_date, interval_start_age, interval_end_age,
  interval_days, interval_event, age_band_start, death_within_12m,
  H_A, D_A, H_F, D_F, lag1_H_A, lag1_D_A, lag1_H_F, lag1_D_F, lag1_available,
  sex_baseline, education_z, calendar_year_c
)

if (anyDuplicated(candidate %>% select(cohort, subject_id, interval_start_date))) stop("Duplicate intervals.")
if (any(candidate$interval_end_date <= candidate$interval_start_date)) stop("Non-positive interval.")
if (any(candidate$interval_event & candidate$interval_end_date < candidate$interval_start_date)) stop("Invalid event interval.")

audit <- candidate %>% group_by(cohort, cohort_role) %>% summarise(
  subjects = n_distinct(subject_id), intervals = n(), deaths = sum(interval_event),
  deaths_within_12m = sum(death_within_12m),
  lag1_subjects = n_distinct(subject_id[lag1_available]),
  lag1_intervals = sum(lag1_available), lag1_deaths = sum(interval_event & lag1_available),
  person_years = sum(interval_days) / 365.25,
  D_A_sd = sd(D_A), D_F_sd = sd(D_F),
  .groups = "drop"
)
write_csv(audit, file.path(out_dir, "exact_candidate_audit.csv"), na = "")

age_support <- survSplit(
  Surv(interval_start_age, interval_end_age, interval_event) ~ .,
  data = candidate %>% select(cohort, subject_id, interval_start_age, interval_end_age, interval_event),
  cut = c(65, 75, 85), episode = "age_episode"
) %>% mutate(
  age_band = factor(case_when(
    interval_start_age < 65 ~ "lt65", interval_start_age < 75 ~ "65_74",
    interval_start_age < 85 ~ "75_84", TRUE ~ "ge85"
  ), levels = c("lt65", "65_74", "75_84", "ge85")),
  py = interval_end_age - interval_start_age
) %>% group_by(cohort, age_band) %>% summarise(
  subjects = n_distinct(subject_id), intervals = n(), deaths = sum(interval_event),
  person_years = sum(py), .groups = "drop"
)
write_csv(age_support, file.path(out_dir, "exact_age_band_support.csv"), na = "")

main_audit <- audit %>% filter(cohort %in% main_cohorts)
main_age <- age_support %>% filter(cohort %in% main_cohorts)
gates <- tribble(
  ~gate, ~value, ~pass,
  "zero_duplicate_exact_visits", nrow(dup), nrow(dup) == 0,
  "three_main_cohorts_present", sum(main_cohorts %in% main_audit$cohort), setequal(main_audit$cohort, main_cohorts),
  "main_cohort_deaths_ge500", min(main_audit$deaths), all(main_audit$deaths >= 500),
  "main_age_bands_deaths_ge50", min(main_age$deaths), all(main_age$deaths >= 50),
  "all_components_finite", sum(!complete.cases(candidate[, c("H_A", "D_A", "H_F", "D_F")])),
    all(complete.cases(candidate[, c("H_A", "D_A", "H_F", "D_F")])),
  "no_association_model_fitted", 1, TRUE
)
write_csv(gates, file.path(out_dir, "dynamic_candidate_gates.csv"), na = "")
if (!all(gates$pass %in% TRUE)) stop("Exact dynamic candidate failed a gate; no lock written.")

saveRDS(candidate, file.path(out_dir, "dynamic_candidate_exact.rds"))
writeLines(c(
  "No association model or effect estimate was produced.",
  paste0("Main-cohort deaths: ", paste(paste(main_audit$cohort, main_audit$deaths, sep = "="), collapse = "; "), ".")
), file.path(out_dir, "CANDIDATE_DECISION.txt"), useBytes = TRUE)
writeLines(c(capture.output(sessionInfo())), file.path(out_dir, "candidate_sessionInfo.txt"), useBytes = TRUE)
writeLines(c(
  paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  paste0("Candidate MD5: ", unname(tools::md5sum(file.path(out_dir, "dynamic_candidate_exact.rds")))),
  "No association model was fitted. Do not rerun."
), lock, useBytes = TRUE)

