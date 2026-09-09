#!/usr/bin/env Rscript
# Prespecified clinical-translation analysis:
# adults aged 50-64 years with at least one prior functional measurement.
# The primary model uses continuous history/deviation terms. Binary burden
# strata are descriptive only and are not used as the primary model.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

project_dir <- normalizePath(
  Sys.getenv("DCV_PROJECT_DIR", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
input_root <- Sys.getenv(
  "MORTALITY_PREDICTION_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_prediction")
)
landmark_path <- Sys.getenv(
  "CLINICAL_LANDMARK_INPUT",
  unset = Sys.getenv(
    "MORTALITY_PREDICTION_AGE50_LANDMARK_INPUT",
    unset = file.path(input_root, "landmark", "age50_landmark_dataset.rds")
  )
)
output_root <- Sys.getenv(
  "CLINICAL_PRIOR1_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "clinical_risk_prior1_continuous")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("survival", "dplyr", "tidyr", "readr", "tibble", "purrr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install before running: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

if (!file.exists(landmark_path)) stop("Missing landmark input: ", landmark_path)

min_age <- 50
max_age <- 64
horizon_years <- 5
min_prior_n <- 1L
bootstrap_replicates <- as.integer(Sys.getenv("CLINICAL_BOOTSTRAP_B", unset = "1000"))
bootstrap_seed <- as.integer(Sys.getenv("CLINICAL_BOOTSTRAP_SEED", unset = "20260909"))

protocol <- tibble(
  field = c(
    "target_population", "age_window", "minimum_prior_measurements",
    "primary_endpoint", "development", "internal_validation", "external_validation",
    "primary_model", "secondary_model", "threshold_rule", "dca_role"
  ),
  value = c(
    "Adults aged 50-64 years with at least one prior functional measurement",
    "50-64 years at landmark", as.character(min_prior_n),
    "five-year all-cause mortality", "HRS development split",
    "HRS internal validation", "MHAS and SHARE when event support is estimable",
    "continuous history and deviation terms",
    "current-score model without history/deviation terms",
    "cohort-specific H_F P75, descriptive only",
    "reported only when the evaluated set has sufficient event support"
  )
)
write_csv(protocol, file.path(output_root, "prior1_continuous_protocol.csv"), na = "")

dat <- readRDS(landmark_path)
required <- c(
  "cohort", "subject_id", "split", "prior_n", "start_age", "end_age", "event",
  "followup", "time5", "event5", "observed5", "current_A", "current_F", "H_A",
  "D_A", "H_F", "D_F", "sex", "education", "calendar_year"
)
missing_fields <- setdiff(required, names(dat))
if (length(missing_fields)) stop("Landmark data missing: ", paste(missing_fields, collapse = ", "))
if (anyDuplicated(dat %>% select(cohort, subject_id))) {
  stop("Landmark input contains duplicate participants.")
}

if (!"prediction_supported" %in% names(dat)) dat$prediction_supported <- TRUE

dat <- dat %>%
  mutate(
    cohort = as.character(cohort),
    split = as.character(split),
    prior_n = as.integer(prior_n),
    prediction_supported = as.logical(prediction_supported)
  ) %>%
  filter(
    start_age >= min_age,
    start_age <= max_age,
    prior_n >= min_prior_n
  ) %>%
  filter(if_all(
    c(start_age, end_age, event, followup, time5, event5, observed5,
      current_A, current_F, H_A, D_A, H_F, D_F, education, calendar_year),
    is.finite
  ))
if (!nrow(dat)) stop("No observations remain after age and prior-measurement filtering.")

preserve_cutpoints <- dat %>%
  group_by(cohort) %>%
  summarise(
    current_F_p25 = quantile(current_F, 0.25, na.rm = TRUE, names = FALSE),
    current_F_p75 = quantile(current_F, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

history_cutpoints <- dat %>%
  group_by(cohort) %>%
  summarise(
    H_F_p25 = quantile(H_F, 0.25, na.rm = TRUE, names = FALSE),
    H_F_p75 = quantile(H_F, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

dat <- dat %>%
  left_join(preserve_cutpoints, by = "cohort") %>%
  left_join(history_cutpoints, by = "cohort") %>%
  mutate(
    current_preserved = current_F >= current_F_p25 & current_F <= current_F_p75,
    high_HF = H_F >= H_F_p75
  )

event_audit <- dat %>%
  mutate(
    population = if_else(current_preserved, "all_prior1_preserved_current", "all_prior1_current_outside_P25_P75")
  ) %>%
  group_by(population, cohort, split) %>%
  summarise(
    participants = n(),
    prior_n_min = min(prior_n),
    prior_n_median = median(prior_n),
    deaths_5y = sum(event5 == 1L),
    observed_5y = sum(observed5 == 1L),
    high_HF_participants = sum(high_HF, na.rm = TRUE),
    low_HF_deaths = sum(event5 == 1L & !high_HF),
    high_HF_deaths = sum(event5 == 1L & high_HF),
    .groups = "drop"
  )
write_csv(event_audit, file.path(output_root, "prior1_event_audit.csv"), na = "")

analysis_data <- dat %>% filter(prediction_supported)
if (!nrow(analysis_data)) stop("No prediction-supported observations remain.")

event_gates <- analysis_data %>%
  group_by(cohort, split) %>%
  summarise(
    participants = n(),
    deaths_5y = sum(event5 == 1L),
    low_HF_participants = sum(!high_HF),
    high_HF_participants = sum(high_HF),
    low_HF_deaths = sum(event5 == 1L & !high_HF),
    high_HF_deaths = sum(event5 == 1L & high_HF),
    both_HF_groups_nonempty = low_HF_participants > 0 & high_HF_participants > 0,
    both_HF_groups_have_events = low_HF_deaths > 0 & high_HF_deaths > 0,
    .groups = "drop"
  )
write_csv(event_gates, file.path(output_root, "prior1_event_gates.csv"), na = "")

train <- analysis_data %>% filter(cohort == "HRS", split == "development")
if (!nrow(train) || sum(train$event5 == 1L) < 5) {
  stop("HRS development has fewer than five five-year deaths after prior_n >= 1 filtering.")
}
write_csv(
  tibble(
    cohort = "HRS", split = "development", participants = nrow(train),
    deaths_5y = sum(train$event5 == 1L),
    events_per_primary_parameter = sum(train$event5 == 1L) / 7
  ),
  file.path(output_root, "prior1_model_support.csv"), na = ""
)

rhs <- c(
  model_1_current = "sex + education + calendar_year + current_A + current_F",
  model_2_history_deviation = "sex + education + calendar_year + H_A + D_A + H_F + D_F"
)

fit_one <- function(rhs_string) {
  coxph(
    as.formula(paste("Surv(start_age, end_age, event) ~", rhs_string)),
    data = train, ties = "efron", model = TRUE, x = TRUE, y = TRUE
  )
}
fits <- lapply(rhs, function(x) tryCatch(fit_one(x), error = function(e) e))
if (any(vapply(fits, inherits, logical(1), what = "error"))) {
  bad <- which(vapply(fits, inherits, logical(1), what = "error"))
  stop("Clinical model failed: ", paste(names(rhs)[bad], collapse = ", "))
}

baseline_hazards <- lapply(fits, function(fit) basehaz(fit, centered = FALSE))
h0_at <- function(bh, age) {
  out <- rep(0, length(age))
  inside <- is.finite(age) & age <= max(bh$time)
  idx <- findInterval(age[inside], bh$time)
  ok <- idx > 0L
  out_inside <- numeric(length(idx))
  out_inside[ok] <- bh$hazard[idx[ok]]
  out[inside] <- out_inside
  out
}
predict_risk <- function(fit, bh, d) {
  lp <- as.numeric(predict(fit, newdata = d, type = "lp", reference = "zero"))
  inc <- h0_at(bh, d$start_age + horizon_years) - h0_at(bh, d$start_age)
  pmin(pmax(1 - exp(-pmax(inc, 0) * exp(lp)), 0), 1)
}

ipcw <- function(d) {
  censor <- as.integer(d$event5 == 0L & d$followup < horizon_years)
  sf <- survfit(Surv(time5, censor) ~ 1, data = d)
  g <- function(t, left = FALSE) {
    idx <- findInterval(pmax(as.numeric(t), 0), sf$time, left.open = left)
    out <- rep(1, length(idx))
    ok <- idx > 0L
    out[ok] <- sf$surv[idx[ok]]
    out
  }
  ge <- g(d$time5, left = TRUE)
  gh <- g(rep(horizon_years, nrow(d)))
  ifelse(
    d$event5 == 1L, 1 / pmax(ge, 1e-6),
    ifelse(d$followup >= horizon_years, 1 / pmax(gh, 1e-6), 0)
  )
}

auc_ipcw <- function(d, pred, w) {
  case <- d$event5 == 1L & w > 0
  ctrl <- d$followup >= horizon_years & d$event5 == 0L & w > 0
  if (!any(case) || !any(ctrl)) return(NA_real_)
  a <- pred[case]; b <- pred[ctrl]; wa <- w[case]; wb <- w[ctrl]
  sum(outer(a, b, function(x, y) as.numeric(x > y) + 0.5 * as.numeric(x == y)) * outer(wa, wb)) /
    (sum(wa) * sum(wb))
}

km_risk <- function(d) {
  if (!sum(d$event5 == 1L)) return(0)
  fit <- survfit(Surv(time5, event5) ~ 1, data = d)
  1 - as.numeric(summary(fit, times = horizon_years, extend = TRUE)$surv[[1]])
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, prob = prob, na.rm = TRUE, names = FALSE))
}

calibration <- function(d, pred, w) {
  keep <- w > 0 & is.finite(pred) & pred > 0 & pred < 1
  if (sum(keep) < 100 || length(unique(d$event5[keep])) < 2) {
    return(tibble(mean_predicted_risk = NA_real_, observed_risk = km_risk(d), calibration_slope = NA_real_, calibration_intercept = NA_real_))
  }
  z <- qlogis(pmin(pmax(pred[keep], 1e-6), 1 - 1e-6))
  y <- d$event5[keep]
  ww <- w[keep]
  slope_fit <- suppressWarnings(glm(y ~ z, weights = ww, family = binomial()))
  intercept_fit <- suppressWarnings(glm(y ~ offset(z), weights = ww, family = binomial()))
  tibble(
    mean_predicted_risk = weighted.mean(pred[keep], ww),
    observed_risk = km_risk(d),
    calibration_slope = unname(coef(slope_fit)[2]),
    calibration_intercept = unname(coef(intercept_fit)[1])
  )
}

evaluation_sets <- split(analysis_data, interaction(analysis_data$cohort, analysis_data$split, drop = TRUE))
metric_rows <- list()
bootstrap_rows <- list()
risk_gradient_rows <- list()
dca_rows <- list()
decision_thresholds <- seq(0.05, 0.15, by = 0.01)

net_benefit <- function(d, pred, w, threshold) {
  keep <- w > 0 & is.finite(pred)
  if (!any(keep)) return(NA_real_)
  y <- d$event5[keep]
  ww <- w[keep]
  positive <- pred[keep] >= threshold
  n <- sum(keep)
  tp <- sum(ww * positive * (y == 1L)) / n
  fp <- sum(ww * positive * (y == 0L)) / n
  tp - fp * threshold / (1 - threshold)
}

for (nm in names(evaluation_sets)) {
  d <- evaluation_sets[[nm]]
  if (!nrow(d)) next
  w <- ipcw(d)
  preds <- map2(fits, baseline_hazards, ~predict_risk(.x, .y, d))
  label <- paste(unique(d$cohort), unique(d$split), sep = "_")

  metric_rows[[label]] <- imap_dfr(preds, function(pred, model_name) {
    cal <- calibration(d, pred, w)
    tibble(
      cohort = unique(d$cohort), split = unique(d$split), model = model_name,
      subjects = nrow(d), deaths_5y = sum(d$event5 == 1L),
      auc_5y = auc_ipcw(d, pred, w),
      brier_5y = if (sum(w) > 0) sum(w * (d$event5 - pred)^2) / sum(w) else NA_real_
    ) %>% bind_cols(cal)
  })

  event_rate <- if (sum(w) > 0) sum(w * d$event5) / sum(w) else NA_real_
  dca_rows[[label]] <- map_dfr(names(preds), function(model_name) {
    map_dfr(decision_thresholds, function(threshold) tibble(
      cohort = unique(d$cohort), split = unique(d$split), model = model_name,
      threshold = threshold,
      net_benefit_per_1000 = 1000 * net_benefit(d, preds[[model_name]], w, threshold),
      estimable = is.finite(event_rate) && sum(d$event5 == 1L) > 0,
      role = "model"
    ))
  }) %>% bind_rows(
    map_dfr(decision_thresholds, function(threshold) tibble(
      cohort = unique(d$cohort), split = unique(d$split), model = "treat_all",
      threshold = threshold,
      net_benefit_per_1000 = if (is.finite(event_rate)) 1000 * (event_rate - (1 - event_rate) * threshold / (1 - threshold)) else NA_real_,
      estimable = is.finite(event_rate) && sum(d$event5 == 1L) > 0,
      role = "reference"
    )),
    map_dfr(decision_thresholds, function(threshold) tibble(
      cohort = unique(d$cohort), split = unique(d$split), model = "treat_none",
      threshold = threshold, net_benefit_per_1000 = 0,
      estimable = is.finite(event_rate) && sum(d$event5 == 1L) > 0,
      role = "reference"
    ))
  )

  set.seed(bootstrap_seed + length(bootstrap_rows))
  delta <- replicate(bootstrap_replicates, {
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
    d_b <- d[idx, , drop = FALSE]
    w_b <- ipcw(d_b)
    p1 <- preds$model_1_current[idx]
    p2 <- preds$model_2_history_deviation[idx]
    auc_ipcw(d_b, p2, w_b) - auc_ipcw(d_b, p1, w_b)
  })
  bootstrap_rows[[label]] <- tibble(
    cohort = unique(d$cohort), split = unique(d$split),
    comparison = "model_2_history_deviation_minus_model_1_current",
    estimate = auc_ipcw(d, preds$model_2_history_deviation, w) - auc_ipcw(d, preds$model_1_current, w),
    ci_low = safe_quantile(delta, 0.025),
    ci_high = safe_quantile(delta, 0.975),
    bootstrap_replicates = sum(is.finite(delta))
  )

  # Descriptive absolute-risk gradient in the current-function preserved subgroup.
  preserved <- d %>% filter(current_preserved)
  if (nrow(preserved)) {
    preserved$H_F_quartile <- factor(
      paste0("Q", dplyr::ntile(preserved$H_F, 4L)),
      levels = paste0("Q", 1:4)
    )
    risk_gradient_rows[[label]] <- preserved %>%
      group_by(cohort, split, H_F_quartile) %>%
      group_modify(~tibble(
        participants = nrow(.x),
        deaths_5y = sum(.x$event5 == 1L),
        observed_risk_5y = km_risk(.x),
        mean_H_F = mean(.x$H_F)
      )) %>%
      ungroup()
  }
}

metrics <- bind_rows(metric_rows)
auc_bootstrap <- bind_rows(bootstrap_rows)
risk_gradient <- bind_rows(risk_gradient_rows)
dca <- bind_rows(dca_rows)
write_csv(metrics, file.path(output_root, "prior1_continuous_model_metrics.csv"), na = "")
write_csv(auc_bootstrap, file.path(output_root, "prior1_continuous_auc_bootstrap.csv"), na = "")
write_csv(risk_gradient, file.path(output_root, "prior1_preserved_current_HF_gradient.csv"), na = "")
write_csv(dca, file.path(output_root, "prior1_continuous_decision_curve.csv"), na = "")

coefficients <- imap_dfr(fits, function(fit, model_name) {
  tibble(
    model = model_name,
    term = names(coef(fit)),
    estimate = unname(coef(fit)),
    hazard_ratio = exp(unname(coef(fit)))
  )
})
write_csv(coefficients, file.path(output_root, "prior1_continuous_model_coefficients.csv"), na = "")

if (nrow(risk_gradient)) {
  p_gradient <- risk_gradient %>%
    ggplot(aes(x = H_F_quartile, y = observed_risk_5y, fill = H_F_quartile)) +
    geom_col(width = 0.7) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Historical functional burden quartile", y = "Observed 5-year mortality", title = "Continuous historical burden gradient among participants with preserved current function") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none")
  ggsave(file.path(output_root, "Prior1_Continuous_HF_Gradient.pdf"), p_gradient, width = 180, height = 105, units = "mm", device = cairo_pdf)
  ggsave(file.path(output_root, "Prior1_Continuous_HF_Gradient.tiff"), p_gradient, width = 180, height = 105, units = "mm", dpi = 600, compression = "lzw", device = "tiff")
}

writeLines(c(
  paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  paste0("Input: ", landmark_path),
  "Population: age 50-64 years with prior_n >= 1.",
  "Primary evidence: continuous history/deviation Cox model.",
  "Binary high historical burden is descriptive and does not define the primary model.",
  "DCA was not used as the primary evidence because clinical utility requires adequate event support and external calibration.",
  "Any external set with zero events is reported as non-estimable rather than treated as validation."
), file.path(output_root, "prior1_continuous_analysis_lock.txt"), useBytes = TRUE)

message("Completed prior_n >= 1 continuous clinical-translation analysis: ", output_root)
