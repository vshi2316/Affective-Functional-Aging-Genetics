#!/usr/bin/env Rscript
# Fit the frozen history-anchored dynamic mortality models once.
# Primary: HRS. Replication: MHAS and SHARE.
# Sensitivity: ELSA. Support-only: CHARLS.
# The attained-age time scale, age bands, four primary global tests,
# Holm correction, two reverse-causation analyses and one interaction
# are fixed before association results are read.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

output_root <- Sys.getenv(
  "MORTALITY_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_bridge")
)
candidate_path <- Sys.getenv(
  "MORTALITY_DYNAMIC_CANDIDATE",
  unset = file.path(output_root, "candidate", "dynamic_candidate_exact.rds")
)
candidate_gate_path <- Sys.getenv(
  "MORTALITY_DYNAMIC_CANDIDATE_GATES",
  unset = file.path(output_root, "candidate", "dynamic_candidate_gates.csv")
)
candidate_lock_path <- Sys.getenv(
  "MORTALITY_DYNAMIC_CANDIDATE_LOCK",
  unset = file.path(output_root, "candidate", "dynamic_candidate_lock.txt")
)
out_dir <- Sys.getenv(
  "MORTALITY_MODEL_OUT",
  unset = file.path(output_root, "models")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lock <- file.path(out_dir, "mortality_models_lock.txt")
if (file.exists(lock)) {
  stop("The frozen mortality models are already locked; do not rerun.")
}
required_inputs <- c(candidate_path, candidate_gate_path, candidate_lock_path)
if (any(!file.exists(required_inputs))) {
  stop("Missing candidate input, feasibility gates or candidate lock.")
}

pkgs <- c("survival", "dplyr", "tidyr", "purrr", "readr", "tibble", "metafor")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install before running: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(survival); library(dplyr); library(tidyr); library(purrr)
  library(readr); library(tibble); library(metafor)
})

main_cohorts <- c("HRS", "MHAS", "SHARE")
fit_cohorts <- c(main_cohorts, "ELSA")
band_levels <- c("lt65", "65_74", "75_84", "ge85")
band_labels <- c("<65", "65-74", "75-84", ">=85")
components <- c("H_A", "D_A", "H_F", "D_F")
band_terms <- as.vector(outer(components, band_levels, paste, sep = "__"))

protocol <- tribble(
  ~field, ~frozen_value,
  "analysis_name", "dynamic mortality analysis dynamic estimand-mortality bridge",
  "time_scale", "attained age; exact month-level visit dates",
  "risk_origin", "second explicit-direct joint measurement",
  "main_model", "H_A + D_A + H_F + D_F simultaneously, each varying over four attained-age bands",
  "covariates", "sex; education z-score; calendar year at dynamic entry",
  "variance", "subject-cluster robust sandwich variance",
  "ties", "Efron",
  "primary_testing_cohort", "HRS",
  "replication_cohorts", "MHAS;SHARE",
  "sensitivity_cohort", "ELSA",
  "support_only_cohort", "CHARLS; no age-varying model because 83 deaths and sparse bands",
  "primary_family_1", "global D_A=0 across four age bands",
  "primary_family_2", "global D_F=0 across four age bands",
  "primary_family_3", "global H_A-D_A=0 across four age bands",
  "primary_family_4", "global H_F-D_F=0 across four age bands",
  "multiplicity", "Holm correction across four HRS global tests",
  "age_variation", "four secondary 3-df Wald tests; BH correction within cohort",
  "interaction", "single constant D_A*D_F term; secondary; no age-specific interaction search",
  "reverse_causality_1", "12-month delayed risk entry after each current visit",
  "reverse_causality_2", "one-visit-lagged history and deviation components",
  "history_stability_sensitivity", "current-visit model restricted to intervals with at least two prior joint measurements (third-or-later measurement)",
  "meta_analysis", "HRS/MHAS/SHARE age-band coefficients and contrasts; REML Hartung-Knapp; prediction interval",
  "causal_claim", "prohibited",
  "rerun", "prohibited after result inspection"
)
write_csv(protocol, file.path(out_dir, "frozen_protocol.csv"), na = "")

claim_rules <- tribble(
  ~status, ~claim,
  "allowed_if_supported", "mortality tracked prior history burden and/or current within-person deviation",
  "allowed_if_supported", "history burden and current deviation had different prognostic associations",
  "allowed_if_supported", "associations varied across attained-age bands",
  "required", "report all cohorts, all age bands, heterogeneity and prediction intervals",
  "required", "call H components history-anchored burden, not pure trait-level between-person effects",
  "prohibited", "within-person deterioration caused death",
  "prohibited", "an association constituted dynamic clinical prediction or validated a causal mechanism",
  "prohibited", "call D_A*D_F synergy unless the frozen interaction is supported with replication-consistent direction"
)
write_csv(claim_rules, file.path(out_dir, "claim_rules.csv"), na = "")

g <- read_csv(candidate_gate_path, show_col_types = FALSE)
if (!all(g$pass %in% TRUE)) stop("The dynamic-candidate feasibility gates are not all passed.")
dat0 <- readRDS(candidate_path)
needed <- c(
  "cohort", "cohort_role", "subject_id", "interval_start_age", "interval_end_age",
  "interval_event", "H_A", "D_A", "H_F", "D_F",
  "lag1_H_A", "lag1_D_A", "lag1_H_F", "lag1_D_F", "lag1_available",
  "prior_n", "sex_baseline", "education_z", "calendar_year_c"
)
if (!all(needed %in% names(dat0))) stop("Candidate is missing: ", paste(setdiff(needed, names(dat0)), collapse = ", "))

dat <- dat0 %>% mutate(
  cohort = as.character(cohort), subject_id = as.character(subject_id),
  interval_event = as.logical(interval_event), lag1_available = as.logical(lag1_available)
) %>% group_by(cohort, subject_id) %>% arrange(interval_start_age, .by_group = TRUE) %>%
  mutate(entry_calendar_year_c = first(calendar_year_c)) %>% ungroup()

if (anyDuplicated(dat %>% select(cohort, subject_id, interval_start_age))) stop("Duplicate candidate interval starts.")
if (!all(main_cohorts %in% dat$cohort)) stop("A main cohort is absent.")

## Freeze support-only CHARLS information before excluding it from modelling.
charls_support <- dat %>% filter(cohort == "CHARLS") %>% summarise(
  subjects = n_distinct(subject_id), intervals = n(), deaths = sum(interval_event),
  lag1_subjects = n_distinct(subject_id[lag1_available]), lag1_deaths = sum(interval_event & lag1_available)
)
write_csv(charls_support, file.path(out_dir, "CHARLS_support_only.csv"), na = "")

prepare_variant <- function(d, variant) {
  x <- d
  if (variant == "lag_one_visit") {
    x <- x %>% filter(lag1_available) %>% mutate(
      H_A = lag1_H_A, D_A = lag1_D_A, H_F = lag1_H_F, D_F = lag1_D_F
    )
  }
  if (variant == "ge3_measurements") {
    x <- x %>% filter(prior_n >= 2L)
  }
  if (variant == "delay_12m") {
    x <- x %>% mutate(interval_start_age = interval_start_age + 1) %>%
      filter(interval_end_age > interval_start_age)
  }
  x <- x %>% select(
    cohort, subject_id, interval_start_age, interval_end_age, interval_event,
    H_A, D_A, H_F, D_F, sex_baseline, education_z, entry_calendar_year_c
  )
  if (!nrow(x)) return(x)
  sp <- survSplit(
    Surv(interval_start_age, interval_end_age, interval_event) ~ .,
    data = x, cut = c(65, 75, 85), episode = "age_episode"
  ) %>% mutate(
    age_band = factor(case_when(
      interval_start_age < 65 ~ "lt65", interval_start_age < 75 ~ "65_74",
      interval_start_age < 85 ~ "75_84", TRUE ~ "ge85"
    ), levels = band_levels),
    sex_f = factor(sex_baseline),
    DAxDF = D_A * D_F
  )
  for (cmp in components) for (b in band_levels) {
    nm <- paste(cmp, b, sep = "__")
    sp[[nm]] <- ifelse(sp$age_band == b, sp[[cmp]], 0)
  }
  sp
}

rhs_primary <- paste(c(band_terms, "sex_f", "education_z", "entry_calendar_year_c"), collapse = " + ")
rhs_interaction <- paste(rhs_primary, "+ DAxDF")

safe_fit <- function(cohort_name, variant, interaction = FALSE) {
  d <- dat %>% filter(cohort == cohort_name)
  sp <- tryCatch(prepare_variant(d, variant), error = function(e) e)
  if (inherits(sp, "error")) return(list(
    status = tibble(cohort = cohort_name, variant = variant, interaction = interaction,
                    success = FALSE, subjects = NA_integer_, rows = NA_integer_, deaths = NA_integer_,
                    error = conditionMessage(sp)), fit = NULL, data = NULL
  ))
  rhs <- if (interaction) rhs_interaction else rhs_primary
  form <- as.formula(paste("Surv(interval_start_age, interval_end_age, interval_event) ~", rhs))
  fit <- tryCatch(
    coxph(form, data = sp, ties = "efron", id = subject_id,
          robust = TRUE, model = FALSE, x = FALSE, y = FALSE),
    error = function(e) e
  )
  ok <- !inherits(fit, "error") && !is.null(fit$coefficients) &&
    all(is.finite(coef(fit))) && all(is.finite(diag(vcov(fit))))
  list(
    status = tibble(
      cohort = cohort_name, variant = variant, interaction = interaction,
      success = ok, subjects = n_distinct(sp$subject_id), rows = nrow(sp),
      deaths = sum(sp$interval_event),
      error = if (inherits(fit, "error")) conditionMessage(fit) else if (!ok) "nonfinite_or_nonconverged" else ""
    ),
    fit = if (ok) fit else NULL,
    data = NULL
  )
}

specs <- bind_rows(
  crossing(cohort = fit_cohorts,
           variant = c("primary", "delay_12m", "lag_one_visit", "ge3_measurements")) %>%
    mutate(interaction = FALSE),
  tibble(cohort = fit_cohorts, variant = "primary", interaction = TRUE)
)
fits <- pmap(specs, function(cohort, variant, interaction) safe_fit(cohort, variant, interaction))
status <- bind_rows(map(fits, "status"))
write_csv(status, file.path(out_dir, "model_status.csv"), na = "")

required_ok <- status %>% filter(variant == "primary", !interaction, cohort %in% fit_cohorts)
if (nrow(required_ok) != 4 || !all(required_ok$success)) {
  stop("At least one required exact-date primary model failed; inspect status. No lock written.")
}

get_fit <- function(cohort, variant = "primary", interaction = FALSE) {
  i <- which(specs$cohort == cohort & specs$variant == variant & specs$interaction == interaction)
  if (length(i) != 1) return(NULL)
  fits[[i]]$fit
}

extract_terms <- function(fit, cohort, variant, interaction) {
  if (is.null(fit)) return(tibble())
  b <- coef(fit); V <- vcov(fit); se <- sqrt(diag(V))
  tibble(
    cohort = cohort, variant = variant, interaction = interaction,
    term = names(b), beta = unname(b), se = unname(se),
    hr = exp(beta), ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
    z = beta / se, p = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}
estimates <- pmap_dfr(specs, function(cohort, variant, interaction) {
  extract_terms(get_fit(cohort, variant, interaction), cohort, variant, interaction)
})
write_csv(estimates, file.path(out_dir, "all_cohort_estimates.csv"), na = "")

wald_L <- function(fit, L, label) {
  b <- coef(fit); V <- vcov(fit)
  if (is.null(dim(L))) L <- matrix(L, nrow = 1)
  q <- as.numeric(L %*% b); S <- L %*% V %*% t(L)
  stat <- tryCatch(as.numeric(t(q) %*% solve(S, q)), error = function(e) NA_real_)
  tibble(test = label, df = nrow(L), statistic = stat,
         p = ifelse(is.finite(stat), pchisq(stat, df = nrow(L), lower.tail = FALSE), NA_real_))
}
selector <- function(fit, terms) {
  nm <- names(coef(fit)); L <- matrix(0, nrow = length(terms), ncol = length(nm), dimnames = list(NULL, nm))
  for (i in seq_along(terms)) {
    if (!terms[i] %in% nm) stop("Missing model term: ", terms[i])
    L[i, terms[i]] <- 1
  }
  L
}

global_tests_one <- function(cohort) {
  fit <- get_fit(cohort, "primary", FALSE); nm <- names(coef(fit))
  DA <- paste("D_A", band_levels, sep = "__"); DF <- paste("D_F", band_levels, sep = "__")
  HA <- paste("H_A", band_levels, sep = "__"); HF <- paste("H_F", band_levels, sep = "__")
  Lda <- selector(fit, DA); Ldf <- selector(fit, DF)
  Lha <- selector(fit, HA); Lhf <- selector(fit, HF)
  bind_rows(
    wald_L(fit, Lda, "global_D_A_zero"),
    wald_L(fit, Ldf, "global_D_F_zero"),
    wald_L(fit, Lha - Lda, "global_H_A_equals_D_A"),
    wald_L(fit, Lhf - Ldf, "global_H_F_equals_D_F")
  ) %>% mutate(cohort = cohort, p_holm_within_cohort = p.adjust(p, method = "holm"), .before = 1)
}
global_tests <- map_dfr(fit_cohorts, global_tests_one)
write_csv(global_tests, file.path(out_dir, "four_global_Wald_tests.csv"), na = "")

age_test_one <- function(cohort) {
  fit <- get_fit(cohort, "primary", FALSE)
  map_dfr(components, function(cmp) {
    tr <- paste(cmp, band_levels, sep = "__"); S <- selector(fit, tr)
    C <- rbind(c(1, -1, 0, 0), c(0, 1, -1, 0), c(0, 0, 1, -1))
    wald_L(fit, C %*% S, paste0(cmp, "_constant_across_age_bands"))
  }) %>% mutate(cohort = cohort, p_bh_within_cohort = p.adjust(p, method = "BH"), .before = 1)
}
age_tests <- map_dfr(fit_cohorts, age_test_one)
write_csv(age_tests, file.path(out_dir, "age_variation_Wald_tests.csv"), na = "")

## Per-band history-minus-deviation contrasts, retaining robust covariance.
contrast_one <- function(cohort, domain, age_band) {
  fit <- get_fit(cohort, "primary", FALSE); nm <- names(coef(fit))
  h <- paste0("H_", domain, "__", age_band); d <- paste0("D_", domain, "__", age_band)
  L <- rep(0, length(nm)); names(L) <- nm; L[h] <- 1; L[d] <- -1
  beta <- sum(L * coef(fit)); se <- sqrt(as.numeric(t(L) %*% vcov(fit) %*% L))
  tibble(cohort = cohort, domain = domain, age_band = age_band, contrast = "history_minus_deviation",
         beta = beta, se = se, ratio_of_HRs = exp(beta), ci_low = exp(beta - 1.96 * se),
         ci_high = exp(beta + 1.96 * se), p = 2 * pnorm(abs(beta / se), lower.tail = FALSE))
}
contrasts <- crossing(cohort = fit_cohorts, domain = c("A", "F"), age_band = band_levels) %>%
  pmap_dfr(contrast_one)
write_csv(contrasts, file.path(out_dir, "history_vs_deviation_contrasts.csv"), na = "")

run_meta <- function(d) {
  if (nrow(d) != 3 || any(!is.finite(d$beta) | !is.finite(d$se) | d$se <= 0)) return(tibble())
  m <- tryCatch(rma.uni(yi = d$beta, sei = d$se, method = "REML", test = "knha"), error = function(e) NULL)
  if (is.null(m)) return(tibble())
  pr <- tryCatch(predict(m), error = function(e) NULL)
  tibble(
    k = m$k, beta = as.numeric(m$b[1]), se = m$se,
    hr_or_ratio = exp(as.numeric(m$b[1])), ci_low = exp(m$ci.lb), ci_high = exp(m$ci.ub),
    p_hk = m$pval, tau2 = m$tau2, i2 = m$I2, q = m$QE, q_p = m$QEp,
    prediction_low = if (is.null(pr)) NA_real_ else exp(pr$pi.lb),
    prediction_high = if (is.null(pr)) NA_real_ else exp(pr$pi.ub)
  )
}

coef_main <- estimates %>% filter(
  cohort %in% main_cohorts, variant == "primary", !interaction, term %in% band_terms
) %>% separate(term, into = c("component", "age_band"), sep = "__", remove = FALSE)
meta_coefficients <- coef_main %>% group_by(component, age_band) %>%
  group_modify(~run_meta(.x)) %>% ungroup()
write_csv(meta_coefficients, file.path(out_dir, "main_age_band_meta.csv"), na = "")

meta_contrasts <- contrasts %>% filter(cohort %in% main_cohorts) %>%
  group_by(domain, age_band) %>% group_modify(~run_meta(.x)) %>% ungroup()
write_csv(meta_contrasts, file.path(out_dir, "main_history_deviation_contrast_meta.csv"), na = "")

interaction_results <- estimates %>% filter(variant == "primary", interaction, term == "DAxDF")
interaction_meta <- interaction_results %>% filter(cohort %in% main_cohorts) %>% run_meta()
write_csv(interaction_results, file.path(out_dir, "DA_DF_interaction_by_cohort.csv"), na = "")
write_csv(interaction_meta, file.path(out_dir, "DA_DF_interaction_meta.csv"), na = "")

## Main sensitivity summaries: all 16 exposure coefficients, no selective filter.
sensitivity_estimates <- estimates %>% filter(
  !interaction, variant %in% c("delay_12m", "lag_one_visit", "ge3_measurements"), term %in% band_terms
)
write_csv(sensitivity_estimates, file.path(out_dir, "reverse_causality_sensitivities.csv"), na = "")

gates <- tribble(
  ~gate, ~value, ~pass,
  "required_primary_models_4_of_4", sum(required_ok$success), nrow(required_ok) == 4 && all(required_ok$success),
  "HRS_four_global_tests_present", sum(global_tests$cohort == "HRS"), sum(global_tests$cohort == "HRS") == 4,
  "HRS_Holm_values_complete", sum(is.finite(global_tests$p_holm_within_cohort[global_tests$cohort == "HRS"])),
    all(is.finite(global_tests$p_holm_within_cohort[global_tests$cohort == "HRS"])),
  "main_meta_16_rows_k3", nrow(meta_coefficients), nrow(meta_coefficients) == 16 && all(meta_coefficients$k == 3),
  "main_contrast_meta_8_rows_k3", nrow(meta_contrasts), nrow(meta_contrasts) == 8 && all(meta_contrasts$k == 3),
  "CHARLS_not_modelled", sum(status$cohort == "CHARLS"), sum(status$cohort == "CHARLS") == 0
)
write_csv(gates, file.path(out_dir, "analysis_gates.csv"), na = "")
if (!all(gates$pass %in% TRUE)) stop("Final analysis-gate failure; inspect outputs. No lock written.")

hrs_tests <- global_tests %>% filter(cohort == "HRS")
decision <- c(
  "DYNAMIC ESTIMAND-MORTALITY BRIDGE COMPLETED AND FROZEN",
  paste0("HRS Holm-significant global tests: ", sum(hrs_tests$p_holm_within_cohort < 0.05), "/4."),
  paste0("Required primary models successful: ", sum(required_ok$success), "/4."),
  paste0("Reverse-causality sensitivity models successful: ",
         sum(status$variant %in% c("delay_12m", "lag_one_visit") & status$success), "/8."),
  paste0("At-least-three-measurement sensitivity models successful: ",
         sum(status$variant == "ge3_measurements" & status$success), "/4."),
  paste0("Interaction models successful: ", sum(status$interaction & status$success), "/4."),
  "Interpretation must use HRS primary tests, replication-cohort patterns, random-effects heterogeneity and prediction intervals together.",
  "No causal, genetic-mechanistic, or validated-prediction claim is permitted.",
  "Do not rerun or add outcomes, cut-points, subgroups, interactions or alternative lags."
)
writeLines(decision, file.path(out_dir, "DECISION.txt"), useBytes = TRUE)
saveRDS(list(protocol = protocol, status = status, estimates = estimates,
             global_tests = global_tests, age_tests = age_tests, contrasts = contrasts,
             meta_coefficients = meta_coefficients, meta_contrasts = meta_contrasts,
             interaction_results = interaction_results, interaction_meta = interaction_meta,
             gates = gates), file.path(out_dir, "results_bundle.rds"))
writeLines(c(capture.output(sessionInfo())), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
writeLines(c(
  paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  paste0("Candidate MD5: ", unname(tools::md5sum(candidate_path))),
  paste0("Candidate lock MD5: ", unname(tools::md5sum(candidate_lock_path))),
  "Do not rerun after result inspection."
), lock, useBytes = TRUE)

