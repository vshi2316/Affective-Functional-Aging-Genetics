#!/usr/bin/env Rscript
# Fit cause-specific models for first detected dementia-related outcomes.
# ELSA and SHARE are analysed separately because their outcome wording differs.

options(stringsAsFactors = FALSE)
required <- c("survival", "dplyr", "readr", "tibble", "sandwich")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
output_root <- Sys.getenv(
  "DEMENTIA_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "dementia_bridge")
)
event_dir <- Sys.getenv(
  "DEMENTIA_EVENT_DATA_OUT",
  unset = file.path(output_root, "event_dataset")
)
input_path <- file.path(event_dir, "dementia_event_main.rds")
baseline_sensitivity_path <- file.path(event_dir, "dementia_event_baseline_sensitivity.rds")
strict_observation_path <- file.path(event_dir, "dementia_event_strict_observation.rds")
strict_observation_baseline_sensitivity_path <- file.path(
  event_dir, "dementia_event_strict_observation_baseline_sensitivity.rds"
)
protocol_path <- file.path(event_dir, "analysis_specification.csv")
out_dir <- Sys.getenv(
  "DEMENTIA_MODEL_OUT",
  unset = file.path(output_root, "models")
)
if (!file.exists(input_path)) stop("Missing main event dataset: ", input_path)
if (!file.exists(baseline_sensitivity_path)) stop("Missing baseline-sensitivity input: ", baseline_sensitivity_path)
if (!file.exists(strict_observation_path)) stop("Missing strict-observation input: ", strict_observation_path)
if (!file.exists(strict_observation_baseline_sensitivity_path)) stop("Missing strict-observation baseline input: ", strict_observation_baseline_sensitivity_path)
if (!file.exists(protocol_path)) stop("Missing analysis specification: ", protocol_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lock_path <- file.path(out_dir, "model_results_lock.txt")
if (file.exists(lock_path)) stop("Model results already exist; use an empty output directory for reproduction.")

inference_specification <- tibble::tribble(
  ~item, ~decision,
  "primary_reporting_unit", "Cohort-specific estimates; ELSA and SHARE are not assumed to measure an identical endpoint",
  "main_exposure_family", "Four M2 exposure coefficients within each cohort; Holm correction within cohort",
  "estimand_comparison_family", "Two history-versus-deviation contrasts within each cohort; Holm correction within cohort",
  "primary_M1_M2_test", "Two-degree-of-freedom robust Wald test of beta_HA=beta_DA and beta_HF=beta_DF",
  "partial_likelihood_test", "Supportive nested-model comparison only",
  "sensitivity_inference", "Direction, magnitude and interval stability; no significance-based selection",
  "cross_cohort_rule", "Report both cohorts separately; no forced pooled estimate because endpoint wording differs"
)
readr::write_csv(
  inference_specification,
  file.path(out_dir, "00_formal_inference_specification.csv"),
  na = ""
)

prepare_dat <- function(path) {
  readRDS(path) |>
    as.data.frame() |>
  dplyr::mutate(
    cohort = toupper(as.character(.data$cohort)),
    candidate_id = as.character(.data$candidate_id),
    country = as.factor(.data$country),
    interval_years = .data$analysis_end_age - .data$interval_start_age,
    age_band = cut(
      .data$interval_start_age,
      breaks = c(-Inf, 65, 75, 85, Inf),
      labels = c("younger_than_65", "65_to_74", "75_to_84", "85_or_older"),
      right = FALSE
    )
    ) |>
    dplyr::filter(
    .data$cohort %in% c("ELSA", "SHARE"),
    is.finite(.data$interval_start_age),
    is.finite(.data$analysis_end_age),
    .data$analysis_end_age > .data$interval_start_age,
    .data$dementia_event %in% c(0L, 1L),
    .data$competing_death %in% c(0L, 1L)
    )
}

dat <- prepare_dat(input_path)
dat_baseline_sensitivity <- prepare_dat(baseline_sensitivity_path)
dat_strict_observation <- prepare_dat(strict_observation_path)
dat_strict_observation_baseline_sensitivity <- prepare_dat(strict_observation_baseline_sensitivity_path)

if (any(dat$dementia_event == 1L & dat$competing_death == 1L)) {
  stop("Dementia-related outcome and competing death coexist on a retained row.")
}
if (any(dat_strict_observation$analysis_end_date > dat_strict_observation$last_observed_dementia_assessment_date, na.rm = TRUE)) {
  stop("Strict-observation input extends beyond the last observed dementia-related assessment.")
}
if (any(dat_strict_observation_baseline_sensitivity$analysis_end_date > dat_strict_observation_baseline_sensitivity$last_observed_dementia_assessment_date, na.rm = TRUE)) {
  stop("Strict-observation baseline sensitivity input extends beyond the last observed dementia-related assessment.")
}
terminal_audit <- dat |>
  dplyr::group_by(.data$cohort, .data$candidate_id) |>
  dplyr::mutate(
    terminal_age = max(.data$analysis_end_age),
    nonterminal_event = (.data$dementia_event == 1L | .data$competing_death == 1L) &
      abs(.data$analysis_end_age - .data$terminal_age) > 1e-8
  ) |>
  dplyr::ungroup()
if (any(terminal_audit$nonterminal_event)) stop("At least one event is not terminal.")
dat <- terminal_audit |> dplyr::select(-.data$terminal_age, -.data$nonterminal_event)

flow <- dat |>
  dplyr::group_by(.data$cohort) |>
  dplyr::summarise(
    people = dplyr::n_distinct(.data$candidate_id),
    intervals = dplyr::n(),
    dementia_related_events = sum(.data$dementia_event),
    competing_deaths = sum(.data$competing_death),
    censored_people = dplyr::n_distinct(.data$candidate_id) -
      dplyr::n_distinct(.data$candidate_id[.data$dementia_event == 1L | .data$competing_death == 1L]),
    person_years = sum(.data$interval_years),
    proxy_outcome_events = sum(.data$dementia_event == 1L & .data$first_positive_proxy, na.rm = TRUE),
    first_followup_wave_events = sum(.data$dementia_event == 1L & .data$positive_at_first_followup_wave, na.rm = TRUE),
    outcome_after_later_negative_flag = sum(.data$dementia_event == 1L & .data$negative_after_first_positive, na.rm = TRUE),
    cognition_complete_people = dplyr::n_distinct(.data$candidate_id[is.finite(.data$cognition_z)]),
    cognition_complete_events = sum(.data$dementia_event == 1L & is.finite(.data$cognition_z)),
    .groups = "drop"
  )
readr::write_csv(flow, file.path(out_dir, "01_analysis_flow.csv"), na = "")

covariates <- c("sex_baseline", "education_z", "calendar_year_c")
exposure_m1 <- c("current_affective", "current_functional")
exposure_m2 <- c("H_A", "D_A", "H_F", "D_F")

make_rhs <- function(cohort, terms) {
  rhs <- paste(c(covariates, terms), collapse = " + ")
  if (cohort == "SHARE") rhs <- paste(rhs, "+ strata(country)")
  paste(rhs, "+ cluster(candidate_id)")
}

fit_cox <- function(d, cohort, terms) {
  f <- stats::as.formula(paste(
    "survival::Surv(interval_start_age, analysis_end_age, dementia_event) ~",
    make_rhs(cohort, terms)
  ))
  fit_warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(f, data = d, ties = "efron", x = TRUE, model = TRUE),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(fit, "captured_warnings") <- unique(fit_warnings)
  fit
}

cox_fit_status <- function(fit) {
  warnings <- attr(fit, "captured_warnings")
  warnings <- if (is.null(warnings)) character() else warnings
  numerical_warning <- any(grepl(
    "converg|infinite|singular|did not|failed",
    warnings, ignore.case = TRUE
  ))
  b <- stats::coef(fit)
  V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  valid <- length(b) > 0L && all(is.finite(b)) &&
    !is.null(V) && all(is.finite(V)) &&
    length(fit$loglik) >= 2L && all(is.finite(fit$loglik)) &&
    !numerical_warning
  list(
    converged = isTRUE(valid),
    warning_count = length(warnings),
    warning_message = paste(warnings, collapse = " | ")
  )
}

extract_cox <- function(fit, d, cohort, model_name, analysis_name) {
  s <- summary(fit)
  z <- as.data.frame(s$coefficients)
  z$term <- rownames(z)
  rownames(z) <- NULL
  robust_col <- if ("robust se" %in% names(z)) "robust se" else "se(coef)"
  p_col <- grep("Pr\\(", names(z), value = TRUE)[1L]
  out <- tibble::tibble(
    cohort = cohort,
    analysis = analysis_name,
    model = model_name,
    term = z$term,
    beta = z$coef,
    robust_se = z[[robust_col]],
    hazard_ratio = exp(z$coef),
    conf_low_95 = exp(z$coef - 1.96 * z[[robust_col]]),
    conf_high_95 = exp(z$coef + 1.96 * z[[robust_col]]),
    z_value = z$coef / z[[robust_col]],
    p_value = z[[p_col]],
    people = dplyr::n_distinct(d$candidate_id),
    intervals = nrow(d),
    events = s$nevent,
    concordance = unname(s$concordance[1L])
  )
  exposure_terms <- out$term %in% c(exposure_m1, exposure_m2)
  out$p_holm_exposure_family <- NA_real_
  out$p_holm_exposure_family[exposure_terms] <- stats::p.adjust(
    out$p_value[exposure_terms], method = "holm"
  )
  out
}

contrast_test <- function(fit, cohort, analysis_name) {
  b <- stats::coef(fit)
  V <- stats::vcov(fit)
  needed <- c("H_A", "D_A", "H_F", "D_F")
  if (!all(needed %in% names(b))) stop("M2 coefficients are incomplete in ", cohort, " ", analysis_name)
  one_contrast <- function(a, d, label) {
    cvec <- setNames(rep(0, length(b)), names(b))
    cvec[a] <- 1
    cvec[d] <- -1
    est <- sum(cvec * b)
    se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    tibble::tibble(
      cohort = cohort, analysis = analysis_name, contrast = label,
      beta_difference = est, robust_se = se,
      z_value = est / se,
      p_value = 2 * stats::pnorm(abs(est / se), lower.tail = FALSE),
      ratio_of_hazard_ratios = exp(est),
      conf_low_95 = exp(est - 1.96 * se),
      conf_high_95 = exp(est + 1.96 * se)
    )
  }
  individual <- dplyr::bind_rows(
    one_contrast("H_A", "D_A", "affective_history_equals_deviation"),
    one_contrast("H_F", "D_F", "functional_history_equals_deviation")
  ) |>
    dplyr::mutate(p_holm_contrast_family = stats::p.adjust(.data$p_value, method = "holm"))
  C <- matrix(0, nrow = 2L, ncol = length(b), dimnames = list(
    c("affective", "functional"), names(b)
  ))
  C[1L, "H_A"] <- 1; C[1L, "D_A"] <- -1
  C[2L, "H_F"] <- 1; C[2L, "D_F"] <- -1
  delta <- as.numeric(C %*% b)
  middle <- C %*% V %*% t(C)
  joint_stat <- as.numeric(t(delta) %*% solve(middle, delta))
  joint <- tibble::tibble(
    cohort = cohort, analysis = analysis_name,
    test = "joint_history_deviation_constraints",
    statistic = joint_stat, df = 2L,
    p_value = stats::pchisq(joint_stat, df = 2L, lower.tail = FALSE)
  )
  list(individual = individual, joint = joint)
}

partial_likelihood_comparison <- function(fit1, fit2, cohort, analysis_name) {
  statistic <- 2 * (fit2$loglik[2L] - fit1$loglik[2L])
  tibble::tibble(
    cohort = cohort, analysis = analysis_name,
    comparison = "M2_history_deviation_vs_M1_current_scores",
    statistic = statistic, df = 2L,
    p_value = stats::pchisq(statistic, df = 2L, lower.tail = FALSE),
    interpretation = "Supportive partial-likelihood comparison; robust two-constraint Wald test is primary"
  )
}

main_coefficients <- list()
main_contrasts <- list()
main_joint <- list()
model_comparisons <- list()
model_status <- list()

for (coh in c("ELSA", "SHARE")) {
  d <- dat |> dplyr::filter(.data$cohort == coh)
  fits <- list(
    M0_demographic = fit_cox(d, coh, character()),
    M1_current_scores = fit_cox(d, coh, exposure_m1),
    M2_history_deviation = fit_cox(d, coh, exposure_m2)
  )
  for (nm in names(fits)) {
    fs <- cox_fit_status(fits[[nm]])
    main_coefficients[[paste(coh, nm)]] <- extract_cox(
      fits[[nm]], d, coh, nm, "main_cause_specific_cox"
    )
    model_status[[paste(coh, nm)]] <- tibble::tibble(
      cohort = coh, analysis = "main_cause_specific_cox", model = nm,
      converged = fs$converged,
      warning_count = fs$warning_count,
      warning_message = fs$warning_message,
      people = dplyr::n_distinct(d$candidate_id),
      intervals = nrow(d), events = sum(d$dementia_event),
      log_likelihood = fits[[nm]]$loglik[2L]
    )
  }
  ct <- contrast_test(fits$M2_history_deviation, coh, "main_cause_specific_cox")
  main_contrasts[[coh]] <- ct$individual
  main_joint[[coh]] <- ct$joint
  model_comparisons[[coh]] <- partial_likelihood_comparison(
    fits$M1_current_scores, fits$M2_history_deviation,
    coh, "main_cause_specific_cox"
  )
}

readr::write_csv(dplyr::bind_rows(main_coefficients), file.path(out_dir, "02_main_cox_coefficients.csv"), na = "")
readr::write_csv(dplyr::bind_rows(main_contrasts), file.path(out_dir, "03_history_deviation_contrasts.csv"), na = "")
readr::write_csv(dplyr::bind_rows(main_joint), file.path(out_dir, "04_joint_constraint_tests.csv"), na = "")
readr::write_csv(dplyr::bind_rows(model_comparisons), file.path(out_dir, "05_M1_M2_model_comparisons.csv"), na = "")

## Prespecified sensitivity analyses using M2 only.
make_sensitivity <- function(d, name) {
  if (name == "exclude_first_followup_wave_outcomes") {
    excluded_ids <- unique(d$candidate_id[d$dementia_event == 1L & d$positive_at_first_followup_wave])
    return(d |> dplyr::filter(!.data$candidate_id %in% excluded_ids))
  }
  if (name == "cognition_adjusted_complete_case") {
    return(d |> dplyr::filter(is.finite(.data$cognition_z)))
  }
  if (name == "exclude_same_date_outcome_death") {
    excluded_ids <- unique(d$candidate_id[d$same_date_event])
    return(d |> dplyr::filter(!.data$candidate_id %in% excluded_ids))
  }
  if (name == "direct_interview_endpoint") {
    proxy_event <- d$dementia_event == 1L & d$first_positive_proxy %in% TRUE
    d$dementia_event[proxy_event] <- 0L
    return(d)
  }
  if (name == "outcome_stable_no_later_negative") {
    excluded_ids <- unique(d$candidate_id[
      d$dementia_event == 1L & d$negative_after_first_positive %in% TRUE
    ])
    return(d |> dplyr::filter(!.data$candidate_id %in% excluded_ids))
  }
  stop("Unknown sensitivity analysis: ", name)
}

sensitivity_names <- c(
  "exclude_first_followup_wave_outcomes",
  "cognition_adjusted_complete_case",
  "exclude_same_date_outcome_death",
  "direct_interview_endpoint",
  "same_or_previous_wave_baseline_negative",
  "outcome_stable_no_later_negative",
  "strict_last_assessment_censoring",
  "strict_last_assessment_censoring_same_or_previous_wave"
)
sensitivity_coefficients <- list()
sensitivity_contrasts <- list()
sensitivity_joint <- list()

for (coh in c("ELSA", "SHARE")) {
  d0 <- dat |> dplyr::filter(.data$cohort == coh)
  for (sn in sensitivity_names) {
    ds <- if (sn == "same_or_previous_wave_baseline_negative") {
      dat_baseline_sensitivity |> dplyr::filter(.data$cohort == coh)
    } else if (sn == "strict_last_assessment_censoring") {
      dat_strict_observation |> dplyr::filter(.data$cohort == coh)
    } else if (sn == "strict_last_assessment_censoring_same_or_previous_wave") {
      dat_strict_observation_baseline_sensitivity |> dplyr::filter(.data$cohort == coh)
    } else {
      make_sensitivity(d0, sn)
    }
    terms <- exposure_m2
    if (sn == "cognition_adjusted_complete_case") terms <- c(terms, "cognition_z")
    fit <- fit_cox(ds, coh, terms)
    fs <- cox_fit_status(fit)
    sensitivity_coefficients[[paste(coh, sn)]] <- extract_cox(fit, ds, coh, "M2_history_deviation", sn)
    ct <- contrast_test(fit, coh, sn)
    sensitivity_contrasts[[paste(coh, sn)]] <- ct$individual
    sensitivity_joint[[paste(coh, sn)]] <- ct$joint
    model_status[[paste(coh, sn)]] <- tibble::tibble(
      cohort = coh, analysis = sn, model = "M2_history_deviation",
      converged = fs$converged,
      warning_count = fs$warning_count,
      warning_message = fs$warning_message,
      people = dplyr::n_distinct(ds$candidate_id),
      intervals = nrow(ds), events = sum(ds$dementia_event),
      log_likelihood = fit$loglik[2L]
    )
  }
}

readr::write_csv(dplyr::bind_rows(sensitivity_coefficients), file.path(out_dir, "06_sensitivity_cox_coefficients.csv"), na = "")
readr::write_csv(dplyr::bind_rows(sensitivity_contrasts), file.path(out_dir, "07_sensitivity_contrasts.csv"), na = "")
readr::write_csv(dplyr::bind_rows(sensitivity_joint), file.path(out_dir, "08_sensitivity_joint_tests.csv"), na = "")

observation_rule_coefficients <- dplyr::bind_rows(
  dplyr::bind_rows(main_coefficients) |>
    dplyr::filter(.data$model == "M2_history_deviation", .data$term %in% exposure_m2) |>
    dplyr::mutate(observation_rule = "retain_known_death"),
  dplyr::bind_rows(sensitivity_coefficients) |>
    dplyr::filter(
      .data$analysis == "strict_last_assessment_censoring",
      .data$term %in% exposure_m2
    ) |>
    dplyr::mutate(observation_rule = "strict_last_assessment_censoring")
) |>
  dplyr::select(
    .data$cohort, .data$observation_rule, .data$term,
    .data$hazard_ratio, .data$conf_low_95, .data$conf_high_95,
    .data$p_value, .data$p_holm_exposure_family,
    .data$people, .data$intervals, .data$events
  )
readr::write_csv(
  observation_rule_coefficients,
  file.path(out_dir, "09_observation_rule_coefficient_comparison.csv"),
  na = ""
)

observation_rule_joint <- dplyr::bind_rows(
  dplyr::bind_rows(main_joint) |>
    dplyr::mutate(observation_rule = "retain_known_death"),
  dplyr::bind_rows(sensitivity_joint) |>
    dplyr::filter(.data$analysis == "strict_last_assessment_censoring") |>
    dplyr::mutate(observation_rule = "strict_last_assessment_censoring")
) |>
  dplyr::select(
    .data$cohort, .data$observation_rule,
    .data$statistic, .data$df, .data$p_value
  )
readr::write_csv(
  observation_rule_joint,
  file.path(out_dir, "10_observation_rule_joint_test_comparison.csv"),
  na = ""
)

## Interval-sensitive complementary log-log model with clustered covariance.
fit_cloglog <- function(d, cohort, terms) {
  rhs <- c("factor(age_band)", covariates, terms)
  if (cohort == "SHARE") rhs <- c(rhs, "factor(country)")
  f <- stats::as.formula(paste(
    "dementia_event ~", paste(rhs, collapse = " + "),
    "+ offset(log(interval_years))"
  ))
  stats::glm(f, data = d, family = stats::binomial(link = "cloglog"))
}

extract_cloglog <- function(fit, d, cohort) {
  V <- sandwich::vcovCL(fit, cluster = d$candidate_id, type = "HC0")
  b <- stats::coef(fit)
  se <- sqrt(diag(V))
  keep <- names(b) %in% exposure_m2
  z <- b[keep] / se[keep]
  tibble::tibble(
    cohort = cohort,
    analysis = "wave_interval_cloglog_sensitivity",
    term = names(b)[keep], beta = b[keep], robust_se = se[keep],
    interval_hazard_ratio = exp(b[keep]),
    conf_low_95 = exp(b[keep] - 1.96 * se[keep]),
    conf_high_95 = exp(b[keep] + 1.96 * se[keep]),
    z_value = z,
    p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    people = dplyr::n_distinct(d$candidate_id),
    intervals = nrow(d), events = sum(d$dementia_event)
  ) |>
    dplyr::mutate(p_holm_exposure_family = stats::p.adjust(.data$p_value, method = "holm"))
}

cloglog_results <- list()
for (coh in c("ELSA", "SHARE")) {
  d <- dat |> dplyr::filter(.data$cohort == coh, .data$interval_years > 0)
  fit <- fit_cloglog(d, coh, exposure_m2)
  cloglog_results[[coh]] <- extract_cloglog(fit, d, coh)
  model_status[[paste(coh, "cloglog")]] <- tibble::tibble(
    cohort = coh, analysis = "wave_interval_cloglog_sensitivity",
    model = "M2_history_deviation", converged = isTRUE(fit$converged),
    warning_count = 0L,
    warning_message = "",
    people = dplyr::n_distinct(d$candidate_id),
    intervals = nrow(d), events = sum(d$dementia_event),
    log_likelihood = as.numeric(stats::logLik(fit))
  )
}
readr::write_csv(dplyr::bind_rows(cloglog_results), file.path(out_dir, "11_cloglog_sensitivity.csv"), na = "")

status <- dplyr::bind_rows(model_status)
readr::write_csv(status, file.path(out_dir, "12_model_status.csv"), na = "")
if (any(!status$converged)) stop("At least one model failed to converge; no result lock written.")

manifest <- tibble::tibble(
  role = c("main_input", "baseline_sensitivity_input", "strict_observation_input", "strict_observation_baseline_input", "analysis_specification"),
  path = c(input_path, baseline_sensitivity_path, strict_observation_path, strict_observation_baseline_sensitivity_path, protocol_path),
  md5 = unname(tools::md5sum(c(input_path, baseline_sensitivity_path, strict_observation_path, strict_observation_baseline_sensitivity_path, protocol_path)))
)
readr::write_csv(manifest, file.path(out_dir, "13_input_manifest.csv"), na = "")

writeLines(c(
  "Dynamic dementia-related outcome model family completed.",
  "ELSA and SHARE were analysed and reported separately because outcome wording differs.",
  "All prespecified results are retained irrespective of direction or statistical significance.",
  paste0("Completed at: ", format(Sys.time(), tz = "Asia/Shanghai"))
), lock_path)
message("Dementia-related outcome models completed: ", out_dir)
