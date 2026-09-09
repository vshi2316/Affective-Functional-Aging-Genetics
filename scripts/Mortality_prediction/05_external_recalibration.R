#!/usr/bin/env Rscript
# External transportability analysis for the continuous history/deviation model.
# Coefficients are frozen in HRS development. Only the 5-year calibration
# intercept is updated within each evaluation cohort. No predictor coefficient
# is re-estimated in the recalibrated predictions.

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
  "CLINICAL_RECALIBRATION_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "clinical_risk_recalibration")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("survival", "dplyr", "tidyr", "readr", "tibble", "purrr", "ggplot2")
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
})

if (!file.exists(landmark_path)) stop("Missing landmark input: ", landmark_path)

min_age <- 50
max_age <- 64
horizon_years <- 5
min_prior_n <- 1L

protocol <- tibble(
  field = c(
    "target_population", "age_window", "minimum_prior_measurements",
    "development_cohort", "primary_model", "transportability_operation",
    "coefficient_handling", "recalibration_parameter", "primary_outputs"
  ),
  value = c(
    "Adults aged 50-64 years with at least one prior functional measurement",
    "50-64 years at landmark", as.character(min_prior_n),
    "HRS development", "continuous history and deviation model",
    "external-cohort transportability analysis",
    "HRS development coefficients remain fixed",
    "one cohort-specific 5-year calibration intercept; calibration slope fixed at 1",
    "calibration-in-the-large, calibration slope, observed and predicted baseline risk, Brier score, calibration plots"
  )
)
write_csv(protocol, file.path(output_root, "recalibration_protocol.csv"), na = "")

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
    prior_n >= min_prior_n,
    prediction_supported
  ) %>%
  filter(if_all(
    c(start_age, end_age, event, followup, time5, event5, observed5,
      current_A, current_F, H_A, D_A, H_F, D_F, education, calendar_year),
    is.finite
  ))
if (!nrow(dat)) stop("No observations remain after filtering.")

train <- dat %>% filter(cohort == "HRS", split == "development")
if (!nrow(train) || sum(train$event5 == 1L) < 5) {
  stop("HRS development has fewer than five five-year deaths.")
}

rhs <- c(
  model_1_current = "sex + education + calendar_year + current_A + current_F",
  model_2_history_deviation = "sex + education + calendar_year + H_A + D_A + H_F + D_F"
)
fits <- lapply(rhs, function(rhs_string) {
  coxph(
    as.formula(paste("Surv(start_age, end_age, event) ~", rhs_string)),
    data = train, ties = "efron", model = TRUE, x = TRUE, y = TRUE
  )
})
names(fits) <- names(rhs)

baseline_hazards <- lapply(fits, function(fit) basehaz(fit, centered = FALSE))
h0_at <- function(bh, age) {
  out <- rep(0, length(age))
  inside <- is.finite(age) & age <= max(bh$time)
  idx <- findInterval(age[inside], bh$time)
  val <- numeric(length(idx))
  ok <- idx > 0L
  val[ok] <- bh$hazard[idx[ok]]
  out[inside] <- val
  out
}
predict_risk <- function(fit, bh, d) {
  lp <- as.numeric(predict(fit, newdata = d, type = "lp", reference = "zero"))
  inc <- h0_at(bh, d$start_age + horizon_years) - h0_at(bh, d$start_age)
  pmin(pmax(1 - exp(-pmax(inc, 0) * exp(lp)), 1e-7), 1 - 1e-7)
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

km_risk <- function(d) {
  if (!sum(d$event5 == 1L)) return(0)
  fit <- survfit(Surv(time5, event5) ~ 1, data = d)
  1 - as.numeric(summary(fit, times = horizon_years, extend = TRUE)$surv[[1]])
}

weighted_observed_risk <- function(d, w) {
  keep <- w > 0 & is.finite(w)
  if (!any(keep)) return(NA_real_)
  sum(w[keep] * d$event5[keep]) / sum(w[keep])
}

safe_glm <- function(y, pred, w, offset_only = FALSE) {
  keep <- is.finite(y) & is.finite(pred) & is.finite(w) & w > 0
  if (sum(keep) < 50 || length(unique(y[keep])) < 2) return(NULL)
  z <- qlogis(pmin(pmax(pred[keep], 1e-7), 1 - 1e-7))
  if (offset_only) {
    suppressWarnings(glm(y[keep] ~ offset(z), weights = w[keep], family = binomial()))
  } else {
    suppressWarnings(glm(y[keep] ~ z, weights = w[keep], family = binomial()))
  }
}

recalibrate_intercept <- function(d, pred, w) {
  fit <- safe_glm(d$event5, pred, w, offset_only = TRUE)
  if (is.null(fit)) return(NA_real_)
  unname(coef(fit)[1])
}

calibration_slope <- function(d, pred, w) {
  fit <- safe_glm(d$event5, pred, w, offset_only = FALSE)
  if (is.null(fit)) return(NA_real_)
  unname(coef(fit)[2])
}

calibration_intercept <- function(d, pred, w) {
  fit <- safe_glm(d$event5, pred, w, offset_only = FALSE)
  if (is.null(fit)) return(NA_real_)
  unname(coef(fit)[1])
}

brier_ipcw <- function(d, pred, w) {
  keep <- w > 0 & is.finite(pred) & is.finite(w)
  if (!any(keep)) return(NA_real_)
  sum(w[keep] * (d$event5[keep] - pred[keep])^2) / sum(w[keep])
}

safe_mean <- function(x) {
  if (!length(x) || !any(is.finite(x))) return(NA_real_)
  mean(x[is.finite(x)])
}

evaluation_sets <- split(dat, interaction(dat$cohort, dat$split, drop = TRUE))
summary_rows <- list()
plot_rows <- list()

for (nm in names(evaluation_sets)) {
  d <- evaluation_sets[[nm]]
  if (!nrow(d)) next
  w <- ipcw(d)
  label <- paste(unique(d$cohort), unique(d$split), sep = "_")

  for (model_name in names(fits)) {
    pred <- predict_risk(fits[[model_name]], baseline_hazards[[model_name]], d)
    alpha <- recalibrate_intercept(d, pred, w)
    pred_recal <- plogis(qlogis(pred) + alpha)
    cal_slope <- calibration_slope(d, pred, w)
    cal_int <- calibration_intercept(d, pred, w)
    summary_rows[[paste(label, model_name, sep = "_")]] <- tibble(
      cohort = unique(d$cohort),
      split = unique(d$split),
      model = model_name,
      subjects = nrow(d),
      deaths_5y = sum(d$event5 == 1L),
      observed_km_risk_5y = km_risk(d),
      observed_ipcw_risk_5y = weighted_observed_risk(d, w),
      mean_predicted_risk_before = safe_mean(pred),
      mean_predicted_risk_after = safe_mean(pred_recal),
      calibration_in_the_large_before = cal_int,
      recalibration_intercept = alpha,
      calibration_slope = cal_slope,
      brier_before = brier_ipcw(d, pred, w),
      brier_after = brier_ipcw(d, pred_recal, w),
      coefficient_update = "none",
      recalibration_scope = "cohort-specific 5-year intercept only"
    )

    d_plot <- tibble(
      cohort = unique(d$cohort), split = unique(d$split), model = model_name,
      observed = d$event5, weight = w, predicted_before = pred,
      predicted_after = pred_recal
    ) %>%
      filter(weight > 0, is.finite(predicted_before), is.finite(predicted_after)) %>%
      mutate(decile = ntile(predicted_before, 10L)) %>%
      group_by(cohort, split, model, decile) %>%
      summarise(
        observed_risk = sum(weight * observed) / sum(weight),
        predicted_before = mean(predicted_before),
        predicted_after = mean(predicted_after),
        participants = n(),
        .groups = "drop"
      )
    plot_rows[[paste(label, model_name, sep = "_")]] <- d_plot
  }
}

summary_out <- bind_rows(summary_rows)
plot_out <- bind_rows(plot_rows)
write_csv(summary_out, file.path(output_root, "external_recalibration_summary.csv"), na = "")
write_csv(plot_out, file.path(output_root, "external_recalibration_decile_data.csv"), na = "")

calibration_plot <- plot_out %>%
  pivot_longer(
    cols = c(predicted_before, predicted_after),
    names_to = "prediction_version", values_to = "predicted_risk"
  ) %>%
  ggplot(aes(x = predicted_risk, y = observed_risk, colour = prediction_version)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 1.8) +
  facet_grid(model ~ cohort, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Predicted 5-year mortality risk",
    y = "Observed 5-year mortality risk",
    colour = NULL,
    title = "External calibration before and after cohort-specific intercept updating"
  ) +
  theme_bw(base_size = 9)

ggsave(file.path(output_root, "External_Recalibration_Calibration.pdf"), calibration_plot, width = 210, height = 145, units = "mm", device = cairo_pdf)
ggsave(file.path(output_root, "External_Recalibration_Calibration.tiff"), calibration_plot, width = 210, height = 145, units = "mm", dpi = 600, compression = "lzw", device = "tiff")

writeLines(c(
  paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  paste0("Input: ", landmark_path),
  "Population: age 50-64 years with prior_n >= 1.",
  "HRS development coefficients were frozen.",
  "Recalibration updated only the cohort-specific 5-year calibration intercept.",
  "Calibration slope was estimated for assessment and was not used to update predictions.",
  "This output is a transportability analysis, not an independent model validation or clinical deployment assessment."
), file.path(output_root, "external_recalibration_analysis_lock.txt"), useBytes = TRUE)

message("Completed external recalibration transportability analysis: ", output_root)
