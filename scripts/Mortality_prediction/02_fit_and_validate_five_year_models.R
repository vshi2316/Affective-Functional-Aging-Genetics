#!/usr/bin/env Rscript
# Fit three nested HRS models and evaluate frozen five-year mortality risk across cohorts.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
output_root <- Sys.getenv(
  "MORTALITY_PREDICTION_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_prediction")
)
landmark_dir <- file.path(output_root, "landmark")
landmark_path <- file.path(landmark_dir, "age50_landmark_dataset.rds")
landmark_lock <- file.path(landmark_dir, "landmark_lock.txt")
out_dir <- file.path(output_root, "models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lock_path <- file.path(out_dir, "mortality_prediction_models_lock.txt")
if (file.exists(lock_path)) stop("Five-year prediction model stage is already locked: ", lock_path)
if (!file.exists(landmark_path) || !file.exists(landmark_lock)) {
  stop("Run 01_build_age50_landmark_dataset.R before model fitting.")
}

pkgs <- c("survival", "dplyr", "tidyr", "readr", "tibble", "purrr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install before running: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
})

horizon_years <- 5
bootstrap_seed <- 20260828
bootstrap_replicates <- as.integer(
  Sys.getenv("MORTALITY_PREDICTION_BOOTSTRAP_B", unset = "1000")
)
if (!is.finite(bootstrap_replicates) || bootstrap_replicates < 100L) {
  stop("MORTALITY_PREDICTION_BOOTSTRAP_B must be at least 100.")
}

protocol <- tibble(
  field = c(
    "outcome", "horizon_years", "model_0", "model_1", "model_2",
    "primary_comparison", "development", "external_evaluation",
    "external_recalibration", "bootstrap_replicates", "bootstrap_seed"
  ),
  value = c(
    "all-cause mortality",
    as.character(horizon_years),
    "sex + standardized education + calendar year; attained age represented by baseline hazard",
    "model 0 + current affective score + current functional score",
    "model 0 + affective history + affective deviation + functional history + functional deviation",
    "model 2 minus model 1",
    "HRS participant-level 70% development split",
    "HRS 30% internal validation; MHAS and SHARE external; ELSA sensitivity; CHARLS support-only",
    "none; HRS coefficients and cumulative baseline hazards frozen",
    as.character(bootstrap_replicates),
    as.character(bootstrap_seed)
  )
)
write_csv(protocol, file.path(out_dir, "model_protocol.csv"), na = "")

landmark <- readRDS(landmark_path)
needed <- c(
  "cohort", "subject_id", "split", "start_age", "end_age", "event",
  "followup", "time5", "event5", "observed5", "current_A", "current_F",
  "H_A", "D_A", "H_F", "D_F", "sex", "education", "calendar_year"
)
if (!all(needed %in% names(landmark))) {
  stop("Landmark data missing required columns: ", paste(setdiff(needed, names(landmark)), collapse = ", "))
}
if (anyDuplicated(landmark %>% select(cohort, subject_id))) {
  stop("Model input contains repeated participants.")
}

rhs <- c(
  model_0 = "sex + education + calendar_year",
  model_1 = "sex + education + calendar_year + current_A + current_F",
  model_2 = "sex + education + calendar_year + H_A + D_A + H_F + D_F"
)
train <- landmark %>% filter(cohort == "HRS", split == "development")
if (nrow(train) < 1000 || sum(train$event5) < 100) {
  stop("Insufficient HRS development support.")
}

fit_one <- function(rhs_string) {
  coxph(
    as.formula(paste("Surv(start_age, end_age, event) ~", rhs_string)),
    data = train,
    ties = "efron",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
}

fits <- lapply(rhs, function(x) tryCatch(fit_one(x), error = function(e) e))
model_status <- imap_dfr(fits, function(fit, model) {
  ok <- !inherits(fit, "error") && length(coef(fit)) > 0L &&
    all(is.finite(coef(fit))) && all(is.finite(diag(vcov(fit))))
  tibble(
    model = model,
    success = ok,
    development_subjects = nrow(train),
    development_events_5y = sum(train$event5),
    error = if (inherits(fit, "error")) conditionMessage(fit) else if (!ok) "nonfinite_model" else ""
  )
})
write_csv(model_status, file.path(out_dir, "model_status.csv"), na = "")
if (!all(model_status$success)) stop("At least one required prediction model failed.")

coef_table <- imap_dfr(fits, function(fit, model) {
  beta <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  z <- beta / se
  tibble(
    model = model,
    term = names(beta),
    estimate = unname(beta),
    std_error = unname(se),
    z = unname(z),
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
    hazard_ratio = exp(estimate),
    ci_low = exp(estimate - 1.96 * std_error),
    ci_high = exp(estimate + 1.96 * std_error)
  )
})
write_csv(coef_table, file.path(out_dir, "development_coefficients.csv"), na = "")

baseline_hazards <- lapply(fits, function(fit) basehaz(fit, centered = FALSE))
hazard_max_age <- min(vapply(
  baseline_hazards,
  function(x) max(x$time, na.rm = TRUE),
  numeric(1)
))
landmark <- landmark %>%
  mutate(prediction_supported = start_age + horizon_years <= hazard_max_age + 1e-8)
support_audit <- landmark %>% count(cohort, split, prediction_supported, name = "subjects")
write_csv(support_audit, file.path(out_dir, "baseline_hazard_support_audit.csv"), na = "")

## Breslow cumulative baseline hazard is a right-continuous step function.
h0_at <- function(baseline_hazard, age) {
  out <- rep(NA_real_, length(age))
  inside <- is.finite(age) & age <= max(baseline_hazard$time)
  index <- findInterval(age[inside], baseline_hazard$time)
  value <- numeric(length(index))
  positive <- index > 0L
  value[positive] <- baseline_hazard$hazard[index[positive]]
  out[inside] <- value
  out
}

predict_risk <- function(fit, baseline_hazard, data) {
  linear_predictor <- as.numeric(
    predict(fit, newdata = data, type = "lp", reference = "zero")
  )
  cumulative_increment <-
    h0_at(baseline_hazard, data$start_age + horizon_years) -
    h0_at(baseline_hazard, data$start_age)
  pmin(pmax(1 - exp(-pmax(cumulative_increment, 0) * exp(linear_predictor)), 0), 1)
}

ipcw <- function(data) {
  censor_event <- as.integer(data$event5 == 0L & data$followup < horizon_years)
  censor_fit <- survfit(Surv(time5, censor_event) ~ 1, data = data)
  G <- function(time, left_limit = FALSE) {
    requested <- pmax(as.numeric(time), 0)
    index <- findInterval(requested, censor_fit$time, left.open = left_limit)
    value <- rep(1, length(requested))
    positive <- index > 0L
    value[positive] <- censor_fit$surv[index[positive]]
    value
  }
  G_event <- G(data$time5, left_limit = TRUE)
  G_horizon <- G(rep(horizon_years, nrow(data)))
  ifelse(
    data$event5 == 1L,
    1 / pmax(G_event, 1e-6),
    ifelse(data$followup >= horizon_years, 1 / pmax(G_horizon, 1e-6), 0)
  )
}

km_risk_5y <- function(data) {
  fit <- survfit(Surv(time5, event5) ~ 1, data = data)
  value <- summary(fit, times = horizon_years, extend = TRUE)$surv
  1 - as.numeric(value[[1]])
}

auc_ipcw <- function(data, predicted_risk, weight) {
  case <- data$event5 == 1L & weight > 0
  control <- data$followup >= horizon_years & data$event5 == 0L & weight > 0
  if (!any(case) || !any(control)) return(NA_real_)

  case_risk <- predicted_risk[case]
  control_risk <- predicted_risk[control]
  case_weight <- weight[case]
  control_weight <- weight[control]

  order_control <- order(control_risk)
  control_risk <- control_risk[order_control]
  control_weight <- control_weight[order_control]
  group <- cumsum(c(TRUE, diff(control_risk) != 0))
  unique_risk <- control_risk[c(TRUE, diff(control_risk) != 0)]
  unique_weight <- as.numeric(rowsum(control_weight, group, reorder = FALSE))
  cumulative_weight <- cumsum(unique_weight)

  index <- findInterval(case_risk, unique_risk)
  below <- numeric(length(case_risk))
  has_below <- index > 0L
  below[has_below] <- cumulative_weight[index[has_below]]
  equal <- numeric(length(case_risk))
  candidate <- index > 0L & index <= length(unique_risk)
  equal[candidate] <- ifelse(
    unique_risk[index[candidate]] == case_risk[candidate],
    unique_weight[index[candidate]],
    0
  )
  sum(case_weight * (below - 0.5 * equal)) /
    (sum(case_weight) * sum(control_weight))
}

set.seed(bootstrap_seed)
unit_data <- tibble(event5 = c(1L, 1L, 0L, 0L, 0L), followup = c(1, 2, 5, 5, 5))
unit_prediction <- c(0.9, 0.1, 0.8, 0.2, 0.1)
unit_weight <- c(1.2, 0.7, 1.0, 0.8, 1.1)
unit_case <- unit_data$event5 == 1L
unit_control <- unit_data$followup >= 5 & unit_data$event5 == 0L
unit_brute <- sum(
  outer(unit_prediction[unit_case], unit_prediction[unit_control], function(a, b) {
    as.numeric(a > b) + 0.5 * as.numeric(a == b)
  }) * outer(unit_weight[unit_case], unit_weight[unit_control])
) / (sum(unit_weight[unit_case]) * sum(unit_weight[unit_control]))
unit_rank <- auc_ipcw(unit_data, unit_prediction, unit_weight)
auc_test <- tibble(
  test = "weighted_auc_rank_vs_bruteforce",
  rank = unit_rank,
  brute_force = unit_brute,
  absolute_difference = abs(unit_rank - unit_brute),
  pass = abs(unit_rank - unit_brute) < 1e-12
)
write_csv(auc_test, file.path(out_dir, "auc_implementation_test.csv"), na = "")
if (!auc_test$pass) stop("Weighted AUC implementation test failed.")

calibration_summary <- function(data, predicted_risk, weight) {
  keep <- weight > 0 & is.finite(predicted_risk) & predicted_risk > 0 & predicted_risk < 1
  if (sum(keep) < 100) {
    return(tibble(
      mean_predicted_risk = NA_real_, observed_risk = NA_real_,
      calibration_slope = NA_real_, calibration_intercept = NA_real_
    ))
  }
  logit_prediction <- qlogis(pmin(pmax(predicted_risk[keep], 1e-6), 1 - 1e-6))
  outcome <- data$event5[keep]
  analysis_weight <- weight[keep]
  slope_fit <- suppressWarnings(
    glm(outcome ~ logit_prediction, weights = analysis_weight, family = binomial())
  )
  intercept_fit <- suppressWarnings(
    glm(outcome ~ offset(logit_prediction), weights = analysis_weight, family = binomial())
  )
  tibble(
    mean_predicted_risk = weighted.mean(predicted_risk[keep], analysis_weight),
    observed_risk = km_risk_5y(data),
    calibration_slope = unname(coef(slope_fit)[2]),
    calibration_intercept = unname(coef(intercept_fit)[1])
  )
}

calibration_groups <- function(data, predicted_risk, weight, groups = 10L) {
  working <- tibble(
    predicted_risk = predicted_risk,
    event5 = data$event5,
    time5 = data$time5,
    followup = data$followup,
    weight = weight
  ) %>%
    filter(is.finite(predicted_risk)) %>%
    mutate(risk_group = ntile(predicted_risk, groups))
  working %>%
    group_by(risk_group) %>%
    group_modify(~tibble(
      participants = nrow(.x),
      participants_with_observed_status = sum(.x$weight > 0),
      deaths_5y = sum(.x$event5),
      mean_predicted_risk = mean(.x$predicted_risk),
      observed_risk = km_risk_5y(.x),
      min_predicted_risk = min(.x$predicted_risk),
      max_predicted_risk = max(.x$predicted_risk)
    )) %>%
    ungroup()
}

evaluation_sets <- split(
  landmark,
  interaction(landmark$cohort, landmark$split, drop = TRUE)
)
metrics_parts <- list()
calibration_parts <- list()
prediction_cache <- list()

for (set_name in names(evaluation_sets)) {
  data <- evaluation_sets[[set_name]] %>% filter(prediction_supported)
  if (!nrow(data)) next
  evaluation_label <- paste(unique(data$cohort), unique(data$split), sep = "_")
  weight <- ipcw(data)
  predictions <- map2(fits, baseline_hazards, ~predict_risk(.x, .y, data))
  prediction_cache[[evaluation_label]] <- list(data = data, predictions = predictions)

  metrics_parts[[evaluation_label]] <- imap_dfr(predictions, function(prediction, model) {
    keep <- weight > 0 & is.finite(prediction)
    calibration <- calibration_summary(data, prediction, weight)
    tibble(
      cohort = unique(data$cohort),
      split = unique(data$split),
      model = model,
      subjects = nrow(data),
      events_5y = sum(data$event5),
      brier_5y = sum(weight[keep] * (data$event5[keep] - prediction[keep])^2) / sum(weight[keep]),
      auc_5y = auc_ipcw(data, prediction, weight)
    ) %>% bind_cols(calibration)
  })

  calibration_parts[[evaluation_label]] <- imap_dfr(predictions, function(prediction, model) {
    calibration_groups(data, prediction, weight) %>%
      mutate(
        cohort = unique(data$cohort), split = unique(data$split), model = model,
        .before = 1
      )
  })
}

metrics <- bind_rows(metrics_parts)
calibration_by_decile <- bind_rows(calibration_parts)
write_csv(metrics, file.path(out_dir, "three_model_validation_metrics.csv"), na = "")
write_csv(calibration_by_decile, file.path(out_dir, "calibration_by_decile.csv"), na = "")

fit_model_2 <- fits$model_2
coefficient_names <- names(coef(fit_model_2))
contrast_matrix <- matrix(
  0, nrow = 2, ncol = length(coefficient_names),
  dimnames = list(c("HA_minus_DA", "HF_minus_DF"), coefficient_names)
)
contrast_matrix[1, "H_A"] <- 1
contrast_matrix[1, "D_A"] <- -1
contrast_matrix[2, "H_F"] <- 1
contrast_matrix[2, "D_F"] <- -1
contrast_estimate <- as.numeric(contrast_matrix %*% coef(fit_model_2))
contrast_variance <- contrast_matrix %*% vcov(fit_model_2) %*% t(contrast_matrix)
wald_statistic <- as.numeric(
  t(contrast_estimate) %*% solve(contrast_variance, contrast_estimate)
)
wald <- tibble(
  test = "joint_model2_equals_model1",
  df = 2,
  statistic = wald_statistic,
  p_value = pchisq(wald_statistic, df = 2, lower.tail = FALSE),
  contrast_HA_DA = contrast_estimate[1],
  contrast_HF_DF = contrast_estimate[2]
)
write_csv(wald, file.path(out_dir, "nested_model2_vs_model1_wald.csv"), na = "")

likelihood_ratio_statistic <- 2 *
  (as.numeric(logLik(fits$model_2)) - as.numeric(logLik(fits$model_1)))
likelihood_ratio <- tibble(
  test = "likelihood_ratio_model2_vs_model1",
  df = 2,
  statistic = likelihood_ratio_statistic,
  p_value = pchisq(likelihood_ratio_statistic, df = 2, lower.tail = FALSE)
)
write_csv(
  likelihood_ratio,
  file.path(out_dir, "nested_model2_vs_model1_likelihood_ratio.csv"),
  na = ""
)

bootstrap_one <- function(data, predictions, replicates, seed) {
  if (nrow(data) != n_distinct(data$subject_id)) {
    stop("Bootstrap input contains repeated participants.")
  }
  set.seed(seed)
  output <- vector("list", replicates)
  for (iteration in seq_len(replicates)) {
    index <- sample.int(nrow(data), nrow(data), replace = TRUE)
    sampled_data <- data[index, , drop = FALSE]
    sampled_predictions <- lapply(predictions, `[`, index)
    sampled_weight <- ipcw(sampled_data)
    auc <- vapply(
      sampled_predictions,
      function(prediction) auc_ipcw(sampled_data, prediction, sampled_weight),
      numeric(1)
    )
    output[[iteration]] <- tibble(
      iteration = iteration,
      delta_1_minus_0 = auc["model_1"] - auc["model_0"],
      delta_2_minus_1 = auc["model_2"] - auc["model_1"]
    )
    if (iteration %% 100L == 0L) message("  ", iteration, "/", replicates)
  }
  bind_rows(output)
}

bootstrap_parts <- list()
bootstrap_summary_parts <- list()
partial_bootstrap_path <- file.path(out_dir, "paired_bootstrap_auc_replicates_partial.csv")
set_index <- 0L
for (evaluation_label in names(prediction_cache)) {
  set_index <- set_index + 1L
  cached <- prediction_cache[[evaluation_label]]
  data <- cached$data
  predictions <- cached$predictions
  message("Bootstrap AUC differences: ", evaluation_label)
  replicate_result <- bootstrap_one(
    data, predictions, bootstrap_replicates, bootstrap_seed + set_index
  ) %>%
    mutate(evaluation_set = evaluation_label, .before = 1)
  bootstrap_parts[[evaluation_label]] <- replicate_result

  weight <- ipcw(data)
  observed_auc <- vapply(
    predictions,
    function(prediction) auc_ipcw(data, prediction, weight),
    numeric(1)
  )
  observed_difference <- c(
    delta_1_minus_0 = unname(observed_auc["model_1"] - observed_auc["model_0"]),
    delta_2_minus_1 = unname(observed_auc["model_2"] - observed_auc["model_1"])
  )
  bootstrap_summary_parts[[evaluation_label]] <- map_dfr(
    names(observed_difference),
    function(comparison) {
      values <- replicate_result[[comparison]]
      tibble(
        evaluation_set = evaluation_label,
        comparison = comparison,
        estimate = unname(observed_difference[comparison]),
        ci_low = unname(quantile(values, 0.025, na.rm = TRUE)),
        ci_high = unname(quantile(values, 0.975, na.rm = TRUE)),
        bootstrap_replicates = sum(is.finite(values))
      )
    }
  )

  write_csv(
    bind_rows(bootstrap_parts),
    partial_bootstrap_path,
    na = ""
  )
}

bootstrap_replicates_data <- bind_rows(bootstrap_parts)
bootstrap_summary <- bind_rows(bootstrap_summary_parts)
write_csv(
  bootstrap_replicates_data,
  file.path(out_dir, "paired_bootstrap_auc_replicates.csv"),
  na = ""
)
write_csv(
  bootstrap_summary,
  file.path(out_dir, "paired_bootstrap_auc_differences.csv"),
  na = ""
)
if (file.exists(partial_bootstrap_path)) unlink(partial_bootstrap_path)

analysis_gates <- tibble(
  gate = c(
    "three_models_successful",
    "auc_implementation_test_passed",
    "six_evaluation_sets_present",
    "eighteen_model_metric_rows",
    "twelve_bootstrap_summary_rows",
    "all_bootstrap_replicates_complete",
    "nested_tests_finite",
    "external_coefficients_frozen"
  ),
  value = c(
    sum(model_status$success),
    as.numeric(auc_test$pass),
    n_distinct(paste(metrics$cohort, metrics$split)),
    nrow(metrics),
    nrow(bootstrap_summary),
    min(bootstrap_summary$bootstrap_replicates),
    sum(is.finite(c(wald$statistic, likelihood_ratio$statistic))),
    1
  ),
  pass = c(
    nrow(model_status) == 3L && all(model_status$success),
    auc_test$pass,
    n_distinct(paste(metrics$cohort, metrics$split)) == 6L,
    nrow(metrics) == 18L,
    nrow(bootstrap_summary) == 12L,
    all(bootstrap_summary$bootstrap_replicates == bootstrap_replicates),
    all(is.finite(c(wald$statistic, likelihood_ratio$statistic))),
    TRUE
  )
)
write_csv(analysis_gates, file.path(out_dir, "analysis_gates.csv"), na = "")
if (!all(analysis_gates$pass %in% TRUE)) {
  stop("Five-year prediction analysis gate failed; no model lock written.")
}

saveRDS(
  list(
    protocol = protocol,
    model_status = model_status,
    coefficients = coef_table,
    baseline_hazards = baseline_hazards,
    support_audit = support_audit,
    metrics = metrics,
    calibration_by_decile = calibration_by_decile,
    wald = wald,
    likelihood_ratio = likelihood_ratio,
    bootstrap_summary = bootstrap_summary,
    analysis_gates = analysis_gates
  ),
  file.path(out_dir, "results_bundle.rds")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
writeLines(
  c(
    "FIVE-YEAR MORTALITY PREDICTION COMPLETED",
    "Current affective-functional scores were compared with demographic covariates.",
    "History-deviation components were compared with current scores.",
    "HRS coefficients and cumulative baseline hazards were frozen for every evaluation set.",
    "Participant-level bootstrap quantified paired differences in five-year AUC.",
    "Submission outputs are generated in the next stage without model refitting."
  ),
  file.path(out_dir, "ANALYSIS_DECISION.txt"),
  useBytes = TRUE
)

script_path <- Sys.getenv(
  "MORTALITY_PREDICTION_MODEL_SCRIPT",
  unset = file.path(project_dir, "scripts", "Mortality_prediction", "02_fit_and_validate_five_year_models.R")
)
writeLines(
  c(
    paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
    paste0("Landmark MD5: ", unname(tools::md5sum(landmark_path))),
    paste0("Landmark lock MD5: ", unname(tools::md5sum(landmark_lock))),
    paste0("Model script MD5: ", if (file.exists(script_path)) unname(tools::md5sum(script_path)) else NA_character_),
    paste0("Bootstrap replicates: ", bootstrap_replicates)
  ),
  lock_path,
  useBytes = TRUE
)

message("Five-year mortality prediction completed: ", out_dir)
