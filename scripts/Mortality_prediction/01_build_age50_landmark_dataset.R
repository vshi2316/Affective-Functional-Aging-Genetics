#!/usr/bin/env Rscript
# Build and lock one age-eligible five-year mortality landmark per participant.
# This stage prepares data only and fits no association or prediction model.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
mortality_root <- Sys.getenv(
  "MORTALITY_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_bridge")
)
candidate_path <- Sys.getenv(
  "MORTALITY_PREDICTION_CANDIDATE_INPUT",
  unset = file.path(mortality_root, "candidate", "dynamic_candidate_exact.rds")
)
output_root <- Sys.getenv(
  "MORTALITY_PREDICTION_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_prediction")
)
out_dir <- file.path(output_root, "landmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lock_path <- file.path(out_dir, "landmark_lock.txt")
if (file.exists(lock_path)) {
  stop("Age-50 landmark stage is already locked: ", lock_path)
}
if (!file.exists(candidate_path)) stop("Missing dynamic mortality candidate: ", candidate_path)

pkgs <- c("dplyr", "readr", "tibble")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install before running: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

min_landmark_age <- 50
horizon_years <- 5
split_seed <- 20260828

protocol <- tibble(
  field = c(
    "target", "population_eligibility", "landmark", "horizon_years",
    "development", "external_evaluation", "seed", "future_information"
  ),
  value = c(
    "all-cause mortality after the first age-eligible landmark",
    "attained age at landmark >=50 years",
    "one earliest eligible landmark at or after age 50 per participant",
    as.character(horizon_years),
    "HRS participant-level 70% development split",
    "HRS 30% internal validation; MHAS and SHARE external; ELSA sensitivity; CHARLS support-only",
    as.character(split_seed),
    "predictors restricted to information observed at or before the landmark"
  )
)
write_csv(protocol, file.path(out_dir, "landmark_protocol.csv"), na = "")

candidate_dir <- dirname(candidate_path)
manifest_paths <- c(
  dynamic_candidate = candidate_path,
  candidate_lock = file.path(candidate_dir, "dynamic_candidate_lock.txt"),
  candidate_gates = file.path(candidate_dir, "dynamic_candidate_gates.csv"),
  fixed_scaling = file.path(candidate_dir, "fixed_scaling.csv")
)
input_manifest <- tibble(
  input = names(manifest_paths),
  path = unname(manifest_paths),
  exists = file.exists(path),
  bytes = ifelse(exists, file.info(path)$size, NA_real_),
  md5 = ifelse(exists, unname(tools::md5sum(path)), NA_character_)
)
write_csv(input_manifest, file.path(out_dir, "input_manifest.csv"), na = "")

dat0 <- readRDS(candidate_path)
needed <- c(
  "cohort", "subject_id", "prior_n", "interval_start_age", "interval_end_age",
  "interval_event", "H_A", "D_A", "H_F", "D_F", "sex_baseline",
  "education_z", "calendar_year_c"
)
if (!all(needed %in% names(dat0))) {
  stop("Candidate missing required columns: ", paste(setdiff(needed, names(dat0)), collapse = ", "))
}

dat <- dat0 %>%
  transmute(
    cohort = as.character(cohort),
    subject_id = as.character(subject_id),
    prior_n = as.integer(prior_n),
    start_age = as.numeric(interval_start_age),
    end_age = as.numeric(interval_end_age),
    event = as.integer(as.logical(interval_event)),
    H_A = as.numeric(H_A), D_A = as.numeric(D_A),
    H_F = as.numeric(H_F), D_F = as.numeric(D_F),
    sex = factor(sex_baseline),
    education = as.numeric(education_z),
    calendar_year = as.numeric(calendar_year_c)
  ) %>%
  filter(cohort %in% c("HRS", "MHAS", "SHARE", "ELSA", "CHARLS")) %>%
  filter(
    if_all(c(start_age, end_age, H_A, D_A, H_F, D_F, education, calendar_year), is.finite),
    !is.na(sex), end_age > start_age, prior_n >= 1L
  ) %>%
  arrange(cohort, subject_id, start_age, end_age)

endpoint_person <- dat %>%
  group_by(cohort, subject_id) %>%
  summarise(
    event_rows = sum(event == 1L),
    endpoint_age = max(end_age),
    endpoint_event = as.integer(any(event == 1L)),
    event_end_age = if (any(event == 1L)) max(end_age[event == 1L]) else NA_real_,
    .groups = "drop"
  )

endpoint_audit <- endpoint_person %>%
  group_by(cohort) %>%
  summarise(
    subjects = n(),
    event_subjects = sum(endpoint_event == 1L),
    subjects_with_multiple_event_rows = sum(event_rows > 1L),
    nonterminal_event_subjects = sum(
      !is.na(event_end_age) & abs(event_end_age - endpoint_age) > 1e-8
    ),
    .groups = "drop"
  )
write_csv(endpoint_audit, file.path(out_dir, "endpoint_integrity_audit.csv"), na = "")

eligible_dat <- dat %>% filter(start_age >= min_landmark_age)
if (!nrow(eligible_dat)) stop("No candidate rows satisfy landmark age >=50 years.")

landmark <- eligible_dat %>%
  group_by(cohort, subject_id) %>%
  slice_min(start_age, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-end_age, -event) %>%
  left_join(
    endpoint_person %>% select(cohort, subject_id, endpoint_age, endpoint_event),
    by = c("cohort", "subject_id")
  ) %>%
  transmute(
    cohort, subject_id, prior_n,
    start_age, end_age = endpoint_age, event = endpoint_event,
    H_A, D_A, H_F, D_F,
    current_A = H_A + D_A,
    current_F = H_F + D_F,
    sex, education, calendar_year,
    followup = end_age - start_age,
    event5 = as.integer(event == 1L & followup <= horizon_years),
    observed5 = as.integer((event == 1L & followup <= horizon_years) | followup >= horizon_years),
    time5 = pmin(followup, horizon_years)
  )

if (anyDuplicated(landmark %>% select(cohort, subject_id))) {
  stop("Landmark construction did not yield one row per participant.")
}
if (any(landmark$followup <= 0, na.rm = TRUE)) {
  stop("At least one landmark has non-positive follow-up.")
}

hrs_ids <- landmark %>% filter(cohort == "HRS") %>% distinct(subject_id)
if (nrow(hrs_ids) < 1000) stop("Insufficient HRS participants for development and validation.")
set.seed(split_seed)
development_ids <- sample(hrs_ids$subject_id, floor(0.70 * nrow(hrs_ids)))
landmark <- landmark %>%
  mutate(
    split = case_when(
      cohort != "HRS" ~ "external",
      subject_id %in% development_ids ~ "development",
      TRUE ~ "internal_validation"
    )
  )

participant_split <- landmark %>% count(cohort, split, name = "subjects")
write_csv(participant_split, file.path(out_dir, "participant_split_audit.csv"), na = "")

landmark_flow <- landmark %>%
  group_by(cohort, split) %>%
  summarise(
    subjects = n_distinct(subject_id),
    events_5y = sum(event5),
    observed_5y = sum(observed5),
    censored_before_5y = sum(observed5 == 0L),
    median_followup = median(followup),
    min_landmark_age = min(start_age),
    median_landmark_age = median(start_age),
    max_landmark_age = max(start_age),
    median_prior_measurements = median(prior_n),
    .groups = "drop"
  )
write_csv(landmark_flow, file.path(out_dir, "landmark_flow.csv"), na = "")

history_depth <- landmark %>%
  group_by(cohort, split) %>%
  summarise(
    subjects = n_distinct(subject_id),
    min_prior_measurements = min(prior_n),
    q1_prior_measurements = quantile(prior_n, 0.25),
    median_prior_measurements = median(prior_n),
    q3_prior_measurements = quantile(prior_n, 0.75),
    max_prior_measurements = max(prior_n),
    one_prior_measurement = sum(prior_n == 1L),
    one_prior_measurement_pct = 100 * mean(prior_n == 1L),
    two_or_more_prior_measurements = sum(prior_n >= 2L),
    two_or_more_prior_measurements_pct = 100 * mean(prior_n >= 2L),
    .groups = "drop"
  )
write_csv(history_depth, file.path(out_dir, "landmark_history_depth.csv"), na = "")

integrity <- tibble(
  check = c(
    "at_most_one_event_row_per_participant",
    "all_event_rows_terminal",
    "one_landmark_per_participant",
    "landmark_age_at_least_50",
    "history_uses_at_least_one_prior_measurement",
    "current_affective_identity",
    "current_functional_identity",
    "endpoint_after_landmark",
    "five_cohorts_present"
  ),
  value = c(
    max(endpoint_person$event_rows),
    sum(!is.na(endpoint_person$event_end_age) &
          abs(endpoint_person$event_end_age - endpoint_person$endpoint_age) > 1e-8),
    nrow(landmark) - n_distinct(landmark$cohort, landmark$subject_id),
    min(landmark$start_age),
    min(landmark$prior_n),
    max(abs(landmark$current_A - (landmark$H_A + landmark$D_A))),
    max(abs(landmark$current_F - (landmark$H_F + landmark$D_F))),
    min(landmark$end_age - landmark$start_age),
    n_distinct(landmark$cohort)
  ),
  pass = c(
    all(endpoint_person$event_rows <= 1L),
    all(is.na(endpoint_person$event_end_age) |
          abs(endpoint_person$event_end_age - endpoint_person$endpoint_age) <= 1e-8),
    nrow(landmark) == n_distinct(landmark$cohort, landmark$subject_id),
    all(landmark$start_age >= min_landmark_age),
    all(landmark$prior_n >= 1L),
    max(abs(landmark$current_A - (landmark$H_A + landmark$D_A))) < 1e-10,
    max(abs(landmark$current_F - (landmark$H_F + landmark$D_F))) < 1e-10,
    all(landmark$end_age > landmark$start_age),
    setequal(unique(landmark$cohort), c("HRS", "MHAS", "SHARE", "ELSA", "CHARLS"))
  )
)
write_csv(integrity, file.path(out_dir, "landmark_integrity_audit.csv"), na = "")
if (!all(integrity$pass %in% TRUE)) {
  stop("Age-50 landmark integrity gate failed; inspect landmark_integrity_audit.csv.")
}

landmark_path <- file.path(out_dir, "age50_landmark_dataset.rds")
saveRDS(landmark, landmark_path)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
script_path <- Sys.getenv(
  "MORTALITY_PREDICTION_LANDMARK_SCRIPT",
  unset = file.path(project_dir, "scripts", "Mortality_prediction", "01_build_age50_landmark_dataset.R")
)
writeLines(
  c(
    paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
    paste0("Candidate MD5: ", unname(tools::md5sum(candidate_path))),
    paste0("Landmark MD5: ", unname(tools::md5sum(landmark_path))),
    paste0("Landmark script MD5: ", if (file.exists(script_path)) unname(tools::md5sum(script_path)) else NA_character_),
    paste0("Participants: ", nrow(landmark)),
    "No prediction model was fitted in this stage."
  ),
  lock_path,
  useBytes = TRUE
)

message("Age-50 mortality landmark completed: ", out_dir)
