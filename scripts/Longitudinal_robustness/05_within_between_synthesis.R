## Within-person/history-level synthesis and evidence freeze.
## Final current-data analysis before any external-data decision.
## 1) verifies frozen sample branches; 2) runs one explicit within-between
## decomposition; 3) consolidates evidence and decision gates.

options(stringsAsFactors = FALSE)
base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir)) stop("Set DCV_BASE_DIR before running this script.")
if (!dir.exists(base_dir)) stop("DCV_BASE_DIR does not exist: ", base_dir)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
robustness_dir <- file.path(base_dir, "robustness_analysis")
out_dir <- file.path(robustness_dir, "04_within_between_manuscript_freeze")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("dplyr", "readr", "tibble", "purrr", "nlme", "metafor")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install before running: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(nlme); library(metafor)
})

input_rds <- file.path(base_dir, "analysis_ready_core", "longitudinal_models_v1",
                       "affective_functional_longitudinal_input_v1.rds")
sample_file <- file.path(robustness_dir, "03_lagged_fallback_sample_reconciliation",
                         "sample_branch_reconciliation.csv")
qualification_file <- file.path(robustness_dir, "01_longitudinal_robustness",
                                 "model_qualification_summary.csv")
qualification_effects_file <- file.path(robustness_dir, "01_longitudinal_robustness",
                                         "model_qualification_effects.csv")
lag_file <- file.path(robustness_dir, "03_lagged_fallback_sample_reconciliation",
                      "lagged_target_terms.csv")
genomic_fit_file <- file.path(robustness_dir, "02_genomicsem_boundary",
                              "genomicsem_free_vs_rho1_fit.csv")
genomic_rho_file <- file.path(robustness_dir, "02_genomicsem_boundary",
                              "free_factor_correlation_parameter.csv")
required_files <- c(input_rds, sample_file, qualification_file,
                    qualification_effects_file, lag_file, genomic_fit_file, genomic_rho_file)
if (any(!file.exists(required_files))) stop("Missing required files: ", paste(required_files[!file.exists(required_files)], collapse = "; "))

## -------------------------------------------------------------------------
## 1. Freeze and verify sample definitions
## -------------------------------------------------------------------------
sample_branches <- read_csv(sample_file, show_col_types = FALSE)
expected <- tribble(
  ~branch, ~persons_expected, ~observations_expected,
  "affective_repeated_plus_baseline_covariates", 189732, 803080,
  "functional_repeated_plus_baseline_covariates", 208751, 953301,
  "joint_repeated_plus_baseline_covariates", 189718, 802984,
  "harmonization_eligible_repeated", 187973, 788393
)
sample_check <- sample_branches %>% filter(cohort == "TOTAL") %>%
  select(branch, persons_observed = persons, observations_observed = observations) %>%
  right_join(expected, by = "branch") %>%
  mutate(pass = persons_observed == persons_expected & observations_observed == observations_expected)
write_csv(sample_check, file.path(out_dir, "sample_freeze_check.csv"), na = "")
if (!all(sample_check$pass)) stop("Sample freeze failed. Do not continue or revise the manuscript.")

## -------------------------------------------------------------------------
## 2. Explicit within-between decomposition (the only new model family)
## -------------------------------------------------------------------------
dat <- readRDS(input_rds)
numeric_vars <- c("year", "wave", "age_baseline", "sex_baseline", "education_baseline",
                  "affective_harmonized_v1", "functional_harmonized_v1")
dat <- dat %>% mutate(across(all_of(numeric_vars), ~suppressWarnings(as.numeric(.x)))) %>%
  filter(if_all(all_of(c("year", "age_baseline", "sex_baseline", "education_baseline",
                         "affective_harmonized_v1", "functional_harmonized_v1")), ~!is.na(.x))) %>%
  arrange(cohort, subject_id, year, wave)
eligible <- dat %>% count(cohort, subject_id, name = "n_joint") %>% filter(n_joint >= 2)
dat <- dat %>% inner_join(eligible, by = c("cohort", "subject_id")) %>%
  group_by(cohort, subject_id) %>%
  mutate(
    baseline_year = min(year),
    time_since_baseline_years = year - baseline_year,
    affective_between = mean(affective_harmonized_v1),
    affective_within = affective_harmonized_v1 - affective_between,
    functional_between = mean(functional_harmonized_v1),
    functional_within = functional_harmonized_v1 - functional_between
  ) %>% ungroup()

specs <- tribble(
  ~outcome, ~model_id, ~formula_text, ~between_term, ~within_term,
  "functional_harmonized_v1", "functional_on_affective_within_between",
  "functional_harmonized_v1 ~ affective_between + affective_within + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "affective_between", "affective_within",
  "affective_harmonized_v1", "affective_on_functional_within_between",
  "affective_harmonized_v1 ~ functional_between + functional_within + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "functional_between", "functional_within"
)

fit_one <- function(d, cohort_name, outcome, model_id, formula_text, between_term, within_term) {
  x <- d %>% filter(cohort == cohort_name)
  n_obs <- nrow(x); n_subjects <- n_distinct(x$subject_id)
  fit <- tryCatch(
    nlme::lme(fixed = as.formula(formula_text), random = ~1 | subject_id,
      data = x, method = "ML", na.action = na.omit,
      control = nlme::lmeControl(msMaxIter = 200, opt = "optim", returnObject = FALSE)),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(
    effects = tibble(cohort = cohort_name, outcome, model_id, term = NA_character_, estimate = NA_real_,
      std_error = NA_real_, p_value = NA_real_, n_obs, n_subjects, status = "failed",
      error_message = conditionMessage(fit)),
    contrast = tibble(cohort = cohort_name, outcome, model_id, contrast = "between_minus_within",
      estimate = NA_real_, std_error = NA_real_, p_value = NA_real_, n_obs, n_subjects,
      status = "failed", error_message = conditionMessage(fit))
  ))
  tt <- as.data.frame(summary(fit)$tTable); tt$term <- rownames(tt); rownames(tt) <- NULL
  effects <- tibble(cohort = cohort_name, outcome, model_id, term = tt$term,
    estimate = tt$Value, std_error = tt$Std.Error, p_value = tt$`p-value`,
    n_obs, n_subjects, status = "ok", error_message = NA_character_)
  b <- nlme::fixef(fit); V <- vcov(fit)
  delta <- unname(b[between_term] - b[within_term])
  delta_se <- sqrt(V[between_term, between_term] + V[within_term, within_term] - 2 * V[between_term, within_term])
  contrast <- tibble(cohort = cohort_name, outcome, model_id, contrast = "between_minus_within",
    estimate = delta, std_error = delta_se, p_value = 2 * pnorm(abs(delta / delta_se), lower.tail = FALSE),
    n_obs, n_subjects, status = "ok", error_message = NA_character_)
  list(effects = effects, contrast = contrast)
}

cohorts <- sort(unique(dat$cohort))
fit_results <- pmap(specs, function(outcome, model_id, formula_text, between_term, within_term) {
  map(cohorts, ~fit_one(dat, .x, outcome, model_id, formula_text, between_term, within_term))
}) %>% flatten()
effects <- map_dfr(fit_results, "effects")
contrasts <- map_dfr(fit_results, "contrast")
target_terms <- c("affective_between", "affective_within", "functional_between", "functional_within")
effects <- effects %>% mutate(p_fdr_sensitivity_family = if_else(term %in% target_terms,
  p.adjust(if_else(term %in% target_terms, p_value, NA_real_), method = "BH"), NA_real_))
contrasts <- contrasts %>% mutate(p_fdr_sensitivity_family = p.adjust(p_value, method = "BH"))
write_csv(effects, file.path(out_dir, "within_between_effects.csv"), na = "")
write_csv(contrasts, file.path(out_dir, "between_minus_within_contrasts.csv"), na = "")

## Random-effects summaries are descriptive because k=5 and heterogeneity is expected.
meta_input <- effects %>% filter(status == "ok", term %in% target_terms, is.finite(estimate),
                                 is.finite(std_error), std_error > 0)
meta_summary <- meta_input %>% group_by(model_id, outcome, term) %>% group_modify(~{
  fit <- metafor::rma.uni(yi = .x$estimate, sei = .x$std_error, method = "REML")
  tibble(k = nrow(.x), estimate = as.numeric(fit$b), std_error = fit$se,
         ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, p_value = fit$pval,
         I2 = fit$I2, tau2 = fit$tau2)
}) %>% ungroup()
write_csv(meta_input, file.path(out_dir, "within_between_meta_input.csv"), na = "")
write_csv(meta_summary, file.path(out_dir, "within_between_meta_summary.csv"), na = "")

## -------------------------------------------------------------------------
## 3. Consolidate manuscript-ready evidence and frozen decision gates
## -------------------------------------------------------------------------
qualification <- read_csv(qualification_file, show_col_types = FALSE)
qualification_effects <- read_csv(qualification_effects_file, show_col_types = FALSE)
lag_targets <- read_csv(lag_file, show_col_types = FALSE)
genomic_fit <- read_csv(genomic_fit_file, show_col_types = FALSE)
genomic_rho <- read_csv(genomic_rho_file, show_col_types = FALSE)

model_status <- qualification %>% mutate(family = case_when(
  grepl("_RI$", model_id) ~ "random_intercept",
  grepl("_RS$", model_id) ~ "random_slope",
  grepl("_CAR1$", model_id) ~ "CAR1",
  TRUE ~ "other")) %>% count(family, status, name = "n_models")
write_csv(model_status, file.path(out_dir, "model_status_audit.csv"), na = "")

primary_terms <- qualification_effects %>% filter(
  (model_id == "functional_on_affective_RI" & term == "affective_harmonized_v1") |
  (model_id == "affective_on_functional_RI" & term == "functional_harmonized_v1") |
  (model_id == "functional_on_affective_RS" & term == "affective_harmonized_v1") |
  (model_id == "affective_on_functional_RS" & term == "functional_harmonized_v1"))
write_csv(primary_terms, file.path(out_dir, "RI_RS_target_effects.csv"), na = "")
write_csv(lag_targets, file.path(out_dir, "lag_delta_and_FE_target_effects.csv"), na = "")
write_csv(genomic_fit, file.path(out_dir, "genomic_model_fit.csv"), na = "")
write_csv(genomic_rho, file.path(out_dir, "genomic_factor_correlation.csv"), na = "")

car_rows <- qualification %>% filter(grepl("_CAR1$", model_id))
pa <- lag_targets %>% filter(estimand == "population_average_cluster_robust")
fe <- lag_targets %>% filter(estimand == "within_person_subject_and_wave_FE")
gfree <- genomic_fit %>% filter(model == "correlated_two_factor_free")
gbound <- genomic_fit %>% filter(model == "rho_equal_1_boundary")
decision_gates <- tribble(
  ~gate, ~value, ~pass, ~interpretation,
  "sample_denominators_reproduced", as.character(sum(sample_check$pass)), all(sample_check$pass), "All four branches match the frozen counts.",
  "RI_and_RS_models_successful", as.character(sum(qualification$status == "ok" & !grepl("CAR1", qualification$model_id))), all(qualification$status[!grepl("CAR1", qualification$model_id)] == "ok"), "RI and RS define the main longitudinal robustness framework.",
  "CAR1_success_count", paste0(sum(car_rows$status == "ok"), "/", nrow(car_rows)), sum(car_rows$status == "ok") == 8, "CAR1 is conditional sensitivity only; failures remain reported.",
  "population_average_lag_positive", paste0(sum(pa$estimate > 0), "/", nrow(pa)), all(pa$estimate > 0), "Supports recurrent population-average lagged associations.",
  "within_person_lag_positive", paste0(sum(fe$estimate > 0), "/", nrow(fe)), all(fe$estimate > 0), "Failure is expected and documents heterogeneity; not a stopping error.",
  "rho1_model_fit_worse", paste0("delta_CFI=", round(gbound$CFI - gfree$CFI, 4), "; delta_SRMR=", round(gbound$SRMR - gfree$SRMR, 4)),
  (gbound$CFI < gfree$CFI & gbound$SRMR > gfree$SRMR), "Retain two correlated domains within the specified covariance model; no biological-independence claim."
)
write_csv(decision_gates, file.path(out_dir, "frozen_decision_gates.csv"), na = "")

claim_map <- tribble(
  ~claim_id, ~claim, ~evidence_file, ~status,
  "C1", "Four parallel model-specific sample branches are reproducible.", "sample_freeze_check.csv", "supported",
  "C2", "Population-average affective-functional associations recur across five cohorts.", "RI_RS_target_effects.csv; lag_delta_and_FE_target_effects.csv", "supported_with_high_heterogeneity",
  "C3", "Within-person estimates vary across cohorts and directions.", "lag_delta_and_FE_target_effects.csv; within_between_effects.csv", "supported",
  "C4", "The free correlated two-factor model fits better than rho=1 constraint.", "genomic_model_fit.csv", "supported_as_model_comparison",
  "C5", "A universal reciprocal causal mechanism exists.", "none", "unsupported_do_not_claim",
  "C6", "The genetic structure independently replicates or predicts longitudinal change.", "none", "unsupported_requires_new_data"
)
write_csv(claim_map, file.path(out_dir, "claim_evidence_map.csv"), na = "")

writeLines(c(
  "FROZEN ONE-SENTENCE ARGUMENT:",
  "Across five ageing cohorts, affective-functional associations recurred at the population-average level but varied across cohorts and within-person estimands; in parallel, conceptually aligned GWAS traits were better represented by two strongly correlated genetic domains than by a unit-correlation model, without establishing a unified causal mechanism or independent genetic replication.",
  "",
  "MODEL ROLES:",
  "Primary longitudinal evidence: cohort-specific population-average random-intercept models.",
  "Robustness: random-slope models in all five cohorts.",
  "Conditional sensitivity: CAR(1), retaining all success/failure states.",
  "Minimal time sensitivity already completed: adjacent lagged observations with lag-by-delta-year interaction.",
  "Within-person sensitivity: subject-plus-wave fixed effects and the within component of the frozen decomposition.",
  "Genetic reference: common one-factor model.",
  "Selected genetic representation: free correlated two-factor model.",
  "Genetic sensitivity: residual-covariance model and rho=1 constrained comparison.",
  "",
  "STOP RULE:",
  "The prespecified current-data longitudinal model family is complete after this script. Any subsequent analysis requires a separately documented protocol."
), file.path(out_dir, "MANUSCRIPT_FREEZE.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Within-person/history-level synthesis completed: ", out_dir)
