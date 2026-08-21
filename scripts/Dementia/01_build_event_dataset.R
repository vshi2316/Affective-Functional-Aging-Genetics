#!/usr/bin/env Rscript
# Build attained-age risk intervals for first detected dementia-related outcomes.
# This stage derives cohort-specific outcomes and does not fit association models.

options(stringsAsFactors = FALSE)
required <- c("haven", "dplyr", "readr", "tibble", "tidyr", "stringr")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
candidate_path <- Sys.getenv(
  "DEMENTIA_DYNAMIC_CANDIDATE_INPUT",
  unset = file.path(project_dir, "results", "mortality_bridge", "candidate", "dynamic_candidate_exact.rds")
)
crosswalk_path <- Sys.getenv(
  "DEMENTIA_ID_CROSSWALK_INPUT",
  unset = file.path(project_dir, "data", "controlled", "dementia_id_crosswalk.rds")
)
elsa_outcome_path <- Sys.getenv(
  "DEMENTIA_ELSA_OUTCOME_INPUT",
  unset = file.path(project_dir, "data", "controlled", "elsa_harmonized_panel.dta")
)
share_outcome_path <- Sys.getenv(
  "DEMENTIA_SHARE_OUTCOME_INPUT",
  unset = file.path(project_dir, "data", "controlled", "share_harmonized_panel.dta")
)
elsa_cognition_path <- Sys.getenv(
  "DEMENTIA_ELSA_COGNITION_INPUT",
  unset = file.path(project_dir, "data", "controlled", "elsa_analysis_ready_long.rds")
)
share_cognition_path <- Sys.getenv(
  "DEMENTIA_SHARE_COGNITION_INPUT",
  unset = file.path(project_dir, "data", "controlled", "share_analysis_ready_long.rds")
)
out_dir <- Sys.getenv(
  "DEMENTIA_EVENT_DATA_OUT",
  unset = file.path(project_dir, "results", "dementia_bridge", "event_dataset")
)

required_inputs <- c(
  dynamic_candidate = candidate_path,
  id_crosswalk = crosswalk_path,
  ELSA_outcomes = elsa_outcome_path,
  SHARE_outcomes = share_outcome_path,
  ELSA_cognition = elsa_cognition_path,
  SHARE_cognition = share_cognition_path
)
for (p in required_inputs) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lock_path <- file.path(out_dir, "event_dataset_lock.txt")
if (file.exists(lock_path)) stop("The dementia event dataset is already locked; use an empty output directory for reproduction.")

norm_chr <- function(x) {
  y <- trimws(as.character(haven::zap_labels(x)))
  y[y %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  y
}
plain_num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
make_date <- function(year, month) {
  yy <- plain_num(year)
  mm <- plain_num(month)
  mm[is.na(mm) | mm < 1 | mm > 12] <- 7
  out <- as.Date(rep(NA_character_, length(yy)))
  ok <- !is.na(yy) & yy >= 1900 & yy <= 2100
  out[ok] <- as.Date(sprintf("%04d-%02d-15", round(yy[ok]), round(mm[ok])))
  out
}
labels_to_table <- function(x, cohort, wave, variable, role) {
  labs <- attr(x, "labels")
  if (is.null(labs) || !length(labs)) return(tibble::tibble(
    cohort = cohort, wave = wave, variable = variable, role = role,
    value = numeric(), value_label = character()
  ))
  tibble::tibble(
    cohort = cohort, wave = wave, variable = variable, role = role,
    value = as.numeric(unname(labs)), value_label = names(labs)
  )
}

analysis_specification <- tibble::tribble(
  ~item, ~decision,
  "cohorts", "ELSA and SHARE analysed separately",
  "outcome_ELSA", "First survey detection of reported dementia",
  "outcome_SHARE", "First survey detection of reported Alzheimer's disease or dementia",
  "outcome_interpretation", "Wave-level first detection; not exact clinical onset",
  "main_baseline_negative", "Dementia-related outcome explicitly coded 0 at the first eligible dynamic anchor wave",
  "sensitivity_baseline_negative", "Outcome coded 0 at the anchor wave or immediately preceding available cohort wave, with no earlier positive record",
  "competing_event", "Exact all-cause death date from the accepted dynamic mortality endpoint",
  "same_date_rule", "Death takes precedence in the main analysis; exclude same-date records in sensitivity analysis",
  "outcome_observation_censoring", "Primary rule retains a known death as a competing event even when it follows the last dementia assessment; a strict last-assessment censoring rule is analysed as sensitivity",
  "main_exposures", "H_A, D_A, H_F and D_F entered simultaneously",
  "current_score_comparator", "current_affective=H_A+D_A and current_functional=H_F+D_F",
  "main_covariates", "Sex, education and calendar year; SHARE stratified by country",
  "variance", "Participant-clustered robust standard errors",
  "primary_model", "Cause-specific attained-age start-stop Cox model; death censors dementia-related outcome follow-up",
  "interval_censoring_sensitivity", "Wave-interval complementary log-log model",
  "cognition_sensitivity", "Cognition from the anchor wave or immediately preceding wave only; standardized within cohort among eligible participants",
  "prodromal_sensitivity", "Exclude dementia-related outcomes first detected at the first post-anchor observed wave",
  "proxy_reporting", "Retain proxy-reported outcome in the main analysis; SHARE proxy codes 1 and 2 are both classified as proxy interviews; direct-interview-only endpoint is a sensitivity analysis",
  "outcome_stability_sensitivity", "Exclude first-positive events followed by a later explicit negative dementia-related assessment",
  "model_comparison", "M1 current scores versus M2 history-deviation components; test beta_HA=beta_DA and beta_HF=beta_DF"
)
readr::write_csv(analysis_specification, file.path(out_dir, "analysis_specification.csv"), na = "")

candidate <- readRDS(candidate_path) |>
  as.data.frame() |>
  dplyr::mutate(
    cohort = toupper(as.character(.data$cohort)),
    candidate_id = norm_chr(.data$subject_id),
    wave = as.integer(.data$wave),
    interval_start_date = as.Date(.data$interval_start_date),
    interval_end_date = as.Date(.data$interval_end_date),
    interval_start_age = as.numeric(.data$interval_start_age),
    interval_end_age = as.numeric(.data$interval_end_age),
    interval_event = as.integer(.data$interval_event),
    current_affective = .data$H_A + .data$D_A,
    current_functional = .data$H_F + .data$D_F,
    dynamic_row_complete = stats::complete.cases(
      .data$H_A, .data$D_A, .data$H_F, .data$D_F,
      .data$sex_baseline, .data$education_z, .data$calendar_year_c,
      .data$interval_start_date, .data$interval_end_date,
      .data$interval_start_age, .data$interval_end_age
    ) & .data$prior_n >= 1
  ) |>
  dplyr::filter(.data$cohort %in% c("ELSA", "SHARE"))

crosswalk <- readRDS(crosswalk_path) |>
  as.data.frame() |>
  dplyr::transmute(
    cohort = toupper(as.character(.data$cohort)),
    candidate_id = norm_chr(.data$candidate_id),
    raw_id = norm_chr(.data$provider_id)
  ) |>
  dplyr::distinct()
candidate <- candidate |>
  dplyr::left_join(crosswalk, by = c("cohort", "candidate_id"))
if (any(is.na(candidate$raw_id))) stop("Crosswalk is incomplete.")

spec <- tibble::tibble(
  cohort = c("ELSA", "SHARE"),
  id_col = c("idauniqc", "mergeid"),
  dementia_stem = c("demene", "alzdeme"),
  waves = list(as.integer(1:10), as.integer(c(2, 4:9))),
  resolved_path = c(elsa_outcome_path, share_outcome_path)
)

raw_long_rows <- list()
label_rows <- list()
country_rows <- list()

for (i in seq_len(nrow(spec))) {
  coh <- spec$cohort[i]
  waves <- as.integer(unlist(spec$waves[[i]], recursive = TRUE, use.names = FALSE))
  id_col <- spec$id_col[i]
  dvars <- paste0("r", waves, spec$dementia_stem[i])
  svars <- paste0("r", waves, "iwstat")
  yvars <- paste0("r", waves, "iwy")
  mvars <- paste0("r", waves, "iwm")
  pvars <- paste0("r", waves, "proxy")
  country_var <- if (coh == "SHARE") "country" else character()
  requested <- c(id_col, dvars, svars, yvars, mvars, pvars, country_var)
  header <- haven::read_dta(spec$resolved_path[i], n_max = 0)
  missing_vars <- setdiff(requested, names(header))
  if (length(missing_vars)) stop(coh, " raw file lacks: ", paste(missing_vars, collapse = ", "))
  raw <- haven::read_dta(spec$resolved_path[i], col_select = tidyselect::all_of(requested)) |>
    as.data.frame()
  raw$raw_id <- norm_chr(raw[[id_col]])
  if (anyDuplicated(raw$raw_id[!is.na(raw$raw_id)])) stop(coh, " raw ID is not unique.")
  if (coh == "SHARE") {
    country_rows[[coh]] <- tibble::tibble(
      cohort = coh, raw_id = raw$raw_id,
      country = as.factor(haven::as_factor(raw$country, levels = "values"))
    )
  } else {
    country_rows[[coh]] <- tibble::tibble(
      cohort = coh, raw_id = raw$raw_id, country = factor("ELSA")
    )
  }
  for (j in seq_along(waves)) {
    w <- waves[j]
    label_rows[[length(label_rows) + 1L]] <- dplyr::bind_rows(
      labels_to_table(raw[[dvars[j]]], coh, w, dvars[j], "dementia_related_outcome"),
      labels_to_table(raw[[pvars[j]]], coh, w, pvars[j], "proxy_interview"),
      labels_to_table(raw[[svars[j]]], coh, w, svars[j], "interview_status")
    )
    raw_long_rows[[length(raw_long_rows) + 1L]] <- tibble::tibble(
      cohort = coh,
      raw_id = raw$raw_id,
      wave = w,
      outcome_value = plain_num(raw[[dvars[j]]]),
      proxy_value = plain_num(raw[[pvars[j]]]),
      interview_status = plain_num(raw[[svars[j]]]),
      assessment_date = make_date(raw[[yvars[j]]], raw[[mvars[j]]])
    )
  }
}

raw_long <- dplyr::bind_rows(raw_long_rows) |>
  dplyr::filter(!is.na(.data$raw_id)) |>
  dplyr::mutate(
    outcome_observed = .data$outcome_value %in% c(0, 1),
    outcome_negative = .data$outcome_value == 0,
    outcome_positive = .data$outcome_value == 1,
    proxy_interview = dplyr::case_when(
      .data$cohort == "SHARE" ~ .data$proxy_value %in% c(1, 2),
      .data$cohort == "ELSA" ~ .data$proxy_value == 1,
      TRUE ~ FALSE
    )
  )
country <- dplyr::bind_rows(country_rows) |>
  dplyr::filter(!is.na(.data$raw_id)) |>
  dplyr::distinct(.data$cohort, .data$raw_id, .keep_all = TRUE)
readr::write_csv(
  dplyr::bind_rows(label_rows) |> dplyr::distinct(),
  file.path(out_dir, "outcome_value_labels.csv"), na = ""
)

anchors <- candidate |>
  dplyr::filter(.data$dynamic_row_complete) |>
  dplyr::arrange(.data$cohort, .data$candidate_id, .data$interval_start_date, .data$wave) |>
  dplyr::group_by(.data$cohort, .data$candidate_id, .data$raw_id) |>
  dplyr::slice(1L) |>
  dplyr::ungroup() |>
  dplyr::select(
    .data$cohort, .data$candidate_id, .data$raw_id,
    anchor_wave = .data$wave, anchor_date = .data$interval_start_date
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    previous_scheduled_wave = {
      cohort_waves <- as.integer(unlist(
        spec$waves[[match(.data$cohort, spec$cohort)]],
        recursive = TRUE, use.names = FALSE
      ))
      earlier <- cohort_waves[cohort_waves < .data$anchor_wave]
      if (length(earlier)) max(earlier) else NA_integer_
    }
  ) |>
  dplyr::ungroup()

anchor_history <- anchors |>
  dplyr::left_join(raw_long, by = c("cohort", "raw_id")) |>
  dplyr::arrange(.data$cohort, .data$candidate_id, .data$wave) |>
  dplyr::group_by(.data$cohort, .data$candidate_id, .data$raw_id,
                  .data$anchor_wave, .data$anchor_date, .data$previous_scheduled_wave) |>
  dplyr::summarise(
    exact_anchor_value = {
      z <- .data$outcome_value[.data$wave == .data$anchor_wave]
      if (length(z)) z[1L] else NA_real_
    },
    previous_available_wave = {
      eligible <- !is.na(.data$previous_scheduled_wave) &
        !is.na(.data$wave) &
        .data$wave == .data$previous_scheduled_wave &
        .data$outcome_observed
      if (any(eligible, na.rm = TRUE)) {
        .data$wave[which(eligible)[1L]]
      } else {
        NA_integer_
      }
    },
    previous_available_value = {
      eligible <- !is.na(.data$previous_scheduled_wave) &
        !is.na(.data$wave) &
        .data$wave == .data$previous_scheduled_wave &
        .data$outcome_observed
      if (any(eligible, na.rm = TRUE)) {
        .data$outcome_value[which(eligible)[1L]]
      } else {
        NA_real_
      }
    },
    any_positive_before_anchor = any(.data$wave < .data$anchor_wave & .data$outcome_positive, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    main_baseline_eligible = .data$exact_anchor_value == 0 & !.data$any_positive_before_anchor,
    sensitivity_baseline_eligible = !.data$any_positive_before_anchor & (
      .data$exact_anchor_value == 0 |
        (is.na(.data$exact_anchor_value) & .data$previous_available_value == 0)
    ),
    previous_wave_gap = .data$anchor_wave - .data$previous_available_wave
  )

death_endpoint <- candidate |>
  dplyr::filter(.data$interval_event == 1L, !is.na(.data$interval_end_date)) |>
  dplyr::group_by(.data$cohort, .data$candidate_id, .data$raw_id) |>
  dplyr::summarise(exact_death_date = min(.data$interval_end_date), .groups = "drop")

event_ledger <- anchors |>
  dplyr::left_join(
    anchor_history,
    by = c("cohort", "candidate_id", "raw_id", "anchor_wave", "anchor_date", "previous_scheduled_wave")
  ) |>
  dplyr::left_join(death_endpoint, by = c("cohort", "candidate_id", "raw_id")) |>
  dplyr::left_join(country, by = c("cohort", "raw_id")) |>
  dplyr::left_join(raw_long, by = c("cohort", "raw_id")) |>
  dplyr::filter(.data$wave > .data$anchor_wave, .data$assessment_date > .data$anchor_date) |>
  dplyr::group_by(.data$cohort, .data$candidate_id, .data$raw_id,
                  .data$anchor_wave, .data$anchor_date,
                  .data$main_baseline_eligible, .data$sensitivity_baseline_eligible,
                  .data$exact_death_date, .data$country) |>
  dplyr::summarise(
    post_anchor_observed_waves = sum(.data$outcome_observed, na.rm = TRUE),
    last_observed_dementia_assessment_date = ifelse(
      any(.data$outcome_observed),
      as.character(max(.data$assessment_date[.data$outcome_observed], na.rm = TRUE)),
      NA_character_
    ),
    first_post_anchor_observed_wave = ifelse(
      any(.data$outcome_observed), min(.data$wave[.data$outcome_observed]), NA_real_
    ),
    first_positive_wave = ifelse(
      any(.data$outcome_positive), min(.data$wave[.data$outcome_positive]), NA_real_
    ),
    first_positive_date = ifelse(
      any(.data$outcome_positive),
      as.character(min(.data$assessment_date[.data$outcome_positive], na.rm = TRUE)),
      NA_character_
    ),
    first_positive_proxy = ifelse(
      any(.data$outcome_positive),
      any(.data$proxy_interview[.data$wave == min(.data$wave[.data$outcome_positive])], na.rm = TRUE),
      NA
    ),
    negative_after_first_positive = ifelse(
      any(.data$outcome_positive),
      any(.data$outcome_negative & .data$wave > min(.data$wave[.data$outcome_positive]), na.rm = TRUE),
      FALSE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    first_positive_date = as.Date(.data$first_positive_date),
    last_observed_dementia_assessment_date = as.Date(.data$last_observed_dementia_assessment_date),
    same_date_event = !is.na(.data$first_positive_date) & !is.na(.data$exact_death_date) &
      .data$first_positive_date == .data$exact_death_date,
    outcome_before_death = !is.na(.data$first_positive_date) &
      (is.na(.data$exact_death_date) | .data$first_positive_date < .data$exact_death_date),
    death_before_or_same_date_outcome = !is.na(.data$exact_death_date) &
      (is.na(.data$first_positive_date) | .data$exact_death_date <= .data$first_positive_date),
    first_interval_positive = !is.na(.data$first_positive_wave) &
      .data$first_positive_wave == .data$first_post_anchor_observed_wave
  )

## Add people with no post-anchor raw record so their loss is explicit.
event_ledger <- anchors |>
  dplyr::select(.data$cohort, .data$candidate_id, .data$raw_id,
                .data$anchor_wave, .data$anchor_date) |>
  dplyr::left_join(event_ledger,
    by = c("cohort", "candidate_id", "raw_id", "anchor_wave", "anchor_date")
  ) |>
  dplyr::left_join(anchor_history |>
    dplyr::select(.data$cohort, .data$candidate_id,
                  baseline_main = .data$main_baseline_eligible,
                  baseline_sensitivity = .data$sensitivity_baseline_eligible),
    by = c("cohort", "candidate_id")
  ) |>
  dplyr::mutate(
    main_baseline_eligible = dplyr::coalesce(.data$main_baseline_eligible, .data$baseline_main, FALSE),
    sensitivity_baseline_eligible = dplyr::coalesce(.data$sensitivity_baseline_eligible, .data$baseline_sensitivity, FALSE),
    post_anchor_observed_waves = dplyr::coalesce(.data$post_anchor_observed_waves, 0L)
  ) |>
  dplyr::select(-.data$baseline_main, -.data$baseline_sensitivity)

## Cognition at the anchor wave or immediately preceding cohort wave only.
cognition_spec <- tibble::tribble(
  ~cohort, ~path,
  "ELSA", elsa_cognition_path,
  "SHARE", share_cognition_path
)
if (any(!file.exists(cognition_spec$path))) stop("Cognition source is unavailable.")
cognition_rows <- list()
for (i in seq_len(nrow(cognition_spec))) {
  x <- readRDS(cognition_spec$path[i]) |> as.data.frame()
  cognition_rows[[i]] <- x |>
    dplyr::transmute(
      cohort = cognition_spec$cohort[i], raw_id = norm_chr(.data$id),
      wave = as.integer(plain_num(.data$wave)), cognition_raw = plain_num(.data$cognition_score)
    ) |>
    dplyr::mutate(cognition = ifelse(is.finite(.data$cognition_raw) & .data$cognition_raw >= 0,
                                     .data$cognition_raw, NA_real_)) |>
    dplyr::filter(!is.na(.data$raw_id), !is.na(.data$wave)) |>
    dplyr::select(-.data$cognition_raw)
}
cognition <- dplyr::bind_rows(cognition_rows)

cognition_at_anchor <- anchors |>
  dplyr::left_join(cognition, by = c("cohort", "raw_id")) |>
  dplyr::filter(.data$wave <= .data$anchor_wave,
                .data$wave >= .data$anchor_wave - 1L,
                is.finite(.data$cognition)) |>
  dplyr::arrange(.data$cohort, .data$candidate_id, dplyr::desc(.data$wave)) |>
  dplyr::group_by(.data$cohort, .data$candidate_id) |>
  dplyr::slice(1L) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    cohort, candidate_id, cognition_wave = wave,
    cognition_wave_gap = anchor_wave - wave,
    cognition
  )

event_ledger <- event_ledger |>
  dplyr::left_join(cognition_at_anchor, by = c("cohort", "candidate_id")) |>
  dplyr::group_by(.data$cohort) |>
  dplyr::mutate(
    cognition_mean_eligible = mean(.data$cognition[.data$main_baseline_eligible], na.rm = TRUE),
    cognition_sd_eligible = stats::sd(.data$cognition[.data$main_baseline_eligible], na.rm = TRUE),
    cognition_z = ifelse(
      is.finite(.data$cognition) & is.finite(.data$cognition_sd_eligible) & .data$cognition_sd_eligible > 0,
      (.data$cognition - .data$cognition_mean_eligible) / .data$cognition_sd_eligible,
      NA_real_
    )
  ) |>
  dplyr::ungroup()

flow <- event_ledger |>
  dplyr::group_by(.data$cohort) |>
  dplyr::summarise(
    dynamic_anchor_people = dplyr::n(),
    main_same_wave_negative = sum(.data$main_baseline_eligible, na.rm = TRUE),
    sensitivity_same_or_previous_negative = sum(.data$sensitivity_baseline_eligible, na.rm = TRUE),
    main_with_post_anchor_outcome_observation = sum(
      .data$main_baseline_eligible & .data$post_anchor_observed_waves > 0, na.rm = TRUE
    ),
    main_first_detected_outcome_before_death = sum(
      .data$main_baseline_eligible & .data$outcome_before_death, na.rm = TRUE
    ),
    main_competing_deaths = sum(
      .data$main_baseline_eligible & .data$death_before_or_same_date_outcome, na.rm = TRUE
    ),
    same_date_events = sum(.data$main_baseline_eligible & .data$same_date_event, na.rm = TRUE),
    first_interval_outcomes = sum(
      .data$main_baseline_eligible & .data$outcome_before_death & .data$first_interval_positive,
      na.rm = TRUE
    ),
    outcome_records_from_proxy = sum(
      .data$main_baseline_eligible & .data$outcome_before_death & .data$first_positive_proxy,
      na.rm = TRUE
    ),
    negative_after_positive_records = sum(
      .data$main_baseline_eligible & .data$negative_after_first_positive, na.rm = TRUE
    ),
    main_with_anchor_or_previous_wave_cognition = sum(
      .data$main_baseline_eligible & is.finite(.data$cognition_z), na.rm = TRUE
    ),
    outcome_events_with_cognition = sum(
      .data$main_baseline_eligible & .data$outcome_before_death & is.finite(.data$cognition_z),
      na.rm = TRUE
    ),
    .groups = "drop"
  )
readr::write_csv(flow, file.path(out_dir, "event_flow.csv"), na = "")

quality <- event_ledger |>
  dplyr::filter(.data$main_baseline_eligible) |>
  dplyr::group_by(.data$cohort) |>
  dplyr::summarise(
    same_date_events = sum(.data$same_date_event, na.rm = TRUE),
    proxy_outcome_events = sum(.data$outcome_before_death & .data$first_positive_proxy, na.rm = TRUE),
    direct_outcome_events = sum(.data$outcome_before_death & !.data$first_positive_proxy, na.rm = TRUE),
    outcome_with_later_negative = sum(.data$negative_after_first_positive, na.rm = TRUE),
    no_post_anchor_outcome_observation = sum(.data$post_anchor_observed_waves == 0, na.rm = TRUE),
    last_observed_assessment_missing = sum(is.na(.data$last_observed_dementia_assessment_date)),
    .groups = "drop"
  )
readr::write_csv(quality, file.path(out_dir, "endpoint_quality.csv"), na = "")

candidate_with_endpoint <- candidate |>
  dplyr::filter(.data$dynamic_row_complete) |>
  dplyr::left_join(
    event_ledger |>
      dplyr::select(
        .data$cohort, .data$candidate_id, .data$main_baseline_eligible,
        .data$sensitivity_baseline_eligible, .data$country,
        .data$first_positive_date, .data$first_positive_wave,
        .data$last_observed_dementia_assessment_date,
        .data$first_post_anchor_observed_wave, .data$first_positive_proxy,
        .data$negative_after_first_positive, .data$exact_death_date,
        .data$same_date_event, .data$outcome_before_death,
        .data$death_before_or_same_date_outcome,
        .data$cognition_wave_gap, .data$cognition_z
      ),
    by = c("cohort", "candidate_id")
  )

build_model_candidate <- function(x, eligibility_column, baseline_definition, observation_rule) {
  out <- x |>
    dplyr::mutate(
      baseline_definition = baseline_definition,
      observation_rule = observation_rule,
      analysis_person = dplyr::coalesce(.data[[eligibility_column]], FALSE),
      at_risk_start = dplyr::case_when(
        observation_rule == "strict_last_assessment_censoring" ~
          .data$analysis_person &
          !is.na(.data$last_observed_dementia_assessment_date) &
          .data$interval_start_date < .data$last_observed_dementia_assessment_date,
        TRUE ~
          .data$analysis_person &
          (
            !is.na(.data$first_positive_date) |
              !is.na(.data$exact_death_date) |
              !is.na(.data$last_observed_dementia_assessment_date)
          )
      ),
      observation_end_date = dplyr::if_else(
        observation_rule == "strict_last_assessment_censoring",
        pmin(
          .data$interval_end_date,
          dplyr::coalesce(.data$first_positive_date, as.Date("9999-12-31")),
          dplyr::coalesce(.data$exact_death_date, as.Date("9999-12-31")),
          dplyr::coalesce(.data$last_observed_dementia_assessment_date, as.Date("9999-12-31")),
          na.rm = TRUE
        ),
        pmin(
          .data$interval_end_date,
          dplyr::coalesce(.data$first_positive_date, as.Date("9999-12-31")),
          dplyr::coalesce(.data$exact_death_date, as.Date("9999-12-31")),
          dplyr::if_else(
            is.na(.data$exact_death_date),
            dplyr::coalesce(.data$last_observed_dementia_assessment_date, as.Date("9999-12-31")),
            as.Date("9999-12-31")
          ),
          na.rm = TRUE
        )
      ),
      analysis_end_date = .data$observation_end_date,
      dementia_event = as.integer(
        .data$at_risk_start & !is.na(.data$first_positive_date) &
          .data$first_positive_date > .data$interval_start_date &
          .data$first_positive_date == .data$analysis_end_date &
          (is.na(.data$exact_death_date) | .data$first_positive_date < .data$exact_death_date)
      ),
      competing_death = as.integer(
        .data$at_risk_start & !is.na(.data$exact_death_date) &
          .data$exact_death_date > .data$interval_start_date &
          .data$exact_death_date == .data$analysis_end_date &
          (is.na(.data$first_positive_date) | .data$exact_death_date <= .data$first_positive_date)
      ),
      analysis_end_age = .data$interval_start_age +
        as.numeric(.data$analysis_end_date - .data$interval_start_date) / 365.25,
      positive_at_first_followup_wave = !is.na(.data$first_positive_wave) &
        .data$first_positive_wave == .data$first_post_anchor_observed_wave
    ) |>
    dplyr::filter(
      .data$analysis_person, .data$at_risk_start,
      !is.na(.data$analysis_end_date),
      .data$analysis_end_date > .data$interval_start_date,
      .data$analysis_end_age > .data$interval_start_age
    ) |>
    dplyr::group_by(.data$cohort, .data$candidate_id) |>
    dplyr::mutate(
      terminal_end_age = max(.data$analysis_end_age, na.rm = TRUE),
      endpoint_not_terminal = (.data$dementia_event == 1L | .data$competing_death == 1L) &
        abs(.data$analysis_end_age - .data$terminal_end_age) > 1e-8
    ) |>
    dplyr::ungroup()
  if (any(out$endpoint_not_terminal, na.rm = TRUE)) {
    stop("At least one event is not on the terminal retained interval for ", baseline_definition, ".")
  }
  if (any(
    out$observation_rule == "strict_last_assessment_censoring" &
      out$analysis_end_date > out$last_observed_dementia_assessment_date,
    na.rm = TRUE
  )) stop("Strict observation-censoring rule was violated.")
  out |>
    dplyr::select(-.data$terminal_end_age, -.data$endpoint_not_terminal)
}

candidate_model <- build_model_candidate(
  candidate_with_endpoint, "main_baseline_eligible", "same_wave_negative", "retain_known_death"
)
candidate_baseline_sensitivity <- build_model_candidate(
  candidate_with_endpoint, "sensitivity_baseline_eligible", "same_or_previous_wave_negative", "retain_known_death"
)
candidate_strict_observation <- build_model_candidate(
  candidate_with_endpoint, "main_baseline_eligible", "same_wave_negative", "strict_last_assessment_censoring"
)
candidate_strict_observation_baseline_sensitivity <- build_model_candidate(
  candidate_with_endpoint, "sensitivity_baseline_eligible", "same_or_previous_wave_negative", "strict_last_assessment_censoring"
)

model_support <- dplyr::bind_rows(
  candidate_model,
  candidate_baseline_sensitivity,
  candidate_strict_observation,
  candidate_strict_observation_baseline_sensitivity
) |>
  dplyr::group_by(.data$observation_rule, .data$baseline_definition, .data$cohort) |>
  dplyr::summarise(
    people = dplyr::n_distinct(.data$candidate_id),
    intervals = dplyr::n(),
    dementia_related_events = sum(.data$dementia_event),
    competing_deaths = sum(.data$competing_death),
    proxy_outcome_events = sum(.data$dementia_event == 1L & .data$first_positive_proxy, na.rm = TRUE),
    first_followup_wave_events = sum(.data$dementia_event == 1L & .data$positive_at_first_followup_wave, na.rm = TRUE),
    people_with_cognition = dplyr::n_distinct(.data$candidate_id[is.finite(.data$cognition_z)]),
    events_with_cognition = sum(.data$dementia_event == 1L & is.finite(.data$cognition_z)),
    .groups = "drop"
  )
readr::write_csv(model_support, file.path(out_dir, "model_support.csv"), na = "")
saveRDS(candidate_model, file.path(out_dir, "dementia_event_main.rds"))
saveRDS(candidate_baseline_sensitivity, file.path(out_dir, "dementia_event_baseline_sensitivity.rds"))
saveRDS(candidate_strict_observation, file.path(out_dir, "dementia_event_strict_observation.rds"))
saveRDS(
  candidate_strict_observation_baseline_sensitivity,
  file.path(out_dir, "dementia_event_strict_observation_baseline_sensitivity.rds")
)

censoring_audit <- event_ledger |>
  dplyr::filter(.data$main_baseline_eligible) |>
  dplyr::group_by(.data$cohort) |>
  dplyr::summarise(
    baseline_eligible_people = dplyr::n(),
    no_post_anchor_outcome_assessment = sum(
      is.na(.data$last_observed_dementia_assessment_date)
    ),
    with_post_anchor_outcome_assessment = sum(
      !is.na(.data$last_observed_dementia_assessment_date)
    ),
    death_after_last_outcome_assessment = sum(
      !is.na(.data$exact_death_date) &
        !is.na(.data$last_observed_dementia_assessment_date) &
        .data$exact_death_date > .data$last_observed_dementia_assessment_date
    ),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    dplyr::bind_rows(candidate_model, candidate_strict_observation) |>
      dplyr::group_by(.data$observation_rule, .data$cohort) |>
      dplyr::summarise(
        retained_people = dplyr::n_distinct(.data$candidate_id),
        retained_intervals = dplyr::n(),
        intervals_censored_at_last_outcome_assessment = sum(
          .data$dementia_event == 0L & .data$competing_death == 0L &
            .data$analysis_end_date == .data$last_observed_dementia_assessment_date,
          na.rm = TRUE
        ),
        retained_competing_deaths = sum(.data$competing_death),
        retained_dementia_events = sum(.data$dementia_event),
        .groups = "drop"
      ),
    by = "cohort"
  )
readr::write_csv(censoring_audit, file.path(out_dir, "observation_censoring_summary.csv"), na = "")

death_after_last_assessment_audit <- event_ledger |>
  dplyr::filter(
    .data$main_baseline_eligible,
    !is.na(.data$exact_death_date),
    !is.na(.data$last_observed_dementia_assessment_date),
    .data$exact_death_date > .data$last_observed_dementia_assessment_date
  ) |>
  dplyr::select(
    .data$cohort, .data$candidate_id, .data$raw_id,
    .data$anchor_date, .data$last_observed_dementia_assessment_date,
    .data$exact_death_date, .data$first_positive_date
  ) |>
  dplyr::left_join(
    candidate_with_endpoint |>
      dplyr::group_by(.data$cohort, .data$candidate_id) |>
      dplyr::summarise(
        last_dynamic_interval_end = max(.data$interval_end_date, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("cohort", "candidate_id")
  ) |>
  dplyr::left_join(
    candidate_model |>
      dplyr::group_by(.data$cohort, .data$candidate_id) |>
      dplyr::summarise(
        retained_in_primary = TRUE,
        primary_competing_death = any(.data$competing_death == 1L),
        primary_terminal_date = max(.data$analysis_end_date),
        .groups = "drop"
      ),
    by = c("cohort", "candidate_id")
  ) |>
  dplyr::left_join(
    candidate_strict_observation |>
      dplyr::group_by(.data$cohort, .data$candidate_id) |>
      dplyr::summarise(
        retained_in_strict_sensitivity = TRUE,
        strict_competing_death = any(.data$competing_death == 1L),
        strict_terminal_date = max(.data$analysis_end_date),
        .groups = "drop"
      ),
    by = c("cohort", "candidate_id")
  ) |>
  dplyr::mutate(
    retained_in_primary = dplyr::coalesce(.data$retained_in_primary, FALSE),
    primary_competing_death = dplyr::coalesce(.data$primary_competing_death, FALSE),
    retained_in_strict_sensitivity = dplyr::coalesce(.data$retained_in_strict_sensitivity, FALSE),
    strict_competing_death = dplyr::coalesce(.data$strict_competing_death, FALSE),
    death_within_dynamic_support = .data$exact_death_date <= .data$last_dynamic_interval_end
  )
readr::write_csv(
  death_after_last_assessment_audit |>
    dplyr::group_by(.data$cohort) |>
    dplyr::summarise(
      deaths_after_last_assessment = dplyr::n(),
      deaths_within_dynamic_support = sum(.data$death_within_dynamic_support, na.rm = TRUE),
      retained_as_primary_competing_death = sum(.data$primary_competing_death, na.rm = TRUE),
      retained_as_strict_competing_death = sum(.data$strict_competing_death, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(out_dir, "death_after_last_assessment_summary.csv"),
  na = ""
)

manifest <- tibble::tibble(
  role = names(required_inputs),
  path = unname(required_inputs),
  md5 = unname(tools::md5sum(required_inputs))
)
readr::write_csv(manifest, file.path(out_dir, "input_manifest.csv"), na = "")

writeLines(c(
  "Dementia-related event datasets completed.",
  "No association model was fitted in this stage.",
  paste0("Completed at: ", format(Sys.time(), tz = "Asia/Shanghai"))
), lock_path)
message("Dementia event datasets completed: ", out_dir)
