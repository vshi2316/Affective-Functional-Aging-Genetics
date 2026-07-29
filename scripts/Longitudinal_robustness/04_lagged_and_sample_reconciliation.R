## Lagged estimands and sample reconciliation.
## Resolves the five singular lagged mixed models with two explicitly labelled
## estimands and reconciles the three manuscript sample denominators.

options(stringsAsFactors = FALSE)
base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir)) stop("Set DCV_BASE_DIR before running this script.")
if (!dir.exists(base_dir)) stop("DCV_BASE_DIR does not exist: ", base_dir)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
input_file <- file.path(base_dir, "analysis_ready_core", "longitudinal_models_v1", "affective_functional_longitudinal_input_v1.rds")
out_dir <- file.path(base_dir, "robustness_analysis", "03_lagged_fallback_sample_reconciliation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pkgs <- c("dplyr", "readr", "tibble", "purrr", "fixest")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install before running: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble); library(purrr); library(fixest)})

dat <- readRDS(input_file)
num <- c("year", "wave", "age_baseline", "sex_baseline", "education_baseline",
         "affective_harmonized_v1", "functional_harmonized_v1", "harmonization_eligibility_v1")
dat <- dat %>% mutate(across(all_of(num), ~suppressWarnings(as.numeric(.x)))) %>% arrange(cohort, subject_id, year, wave)

count_branch <- function(d, branch, complete_vars = character(), eligibility = NULL) {
  x <- d
  if (length(complete_vars)) x <- x %>% filter(if_all(all_of(complete_vars), ~!is.na(.x)))
  if (!is.null(eligibility)) x <- x %>% filter(.data[[eligibility]] == 1)
  eligible <- x %>% count(cohort, subject_id, name = "n_observations") %>% filter(n_observations >= 2)
  x %>% inner_join(eligible, by = c("cohort", "subject_id")) %>% group_by(cohort) %>%
    summarise(branch = branch, persons = n_distinct(subject_id), observations = n(), .groups = "drop")
}

branches <- bind_rows(
  count_branch(dat, "affective_repeated_plus_baseline_covariates",
               c("affective_harmonized_v1", "age_baseline", "sex_baseline", "education_baseline")),
  count_branch(dat, "functional_repeated_plus_baseline_covariates",
               c("functional_harmonized_v1", "age_baseline", "sex_baseline", "education_baseline")),
  count_branch(dat, "joint_repeated_plus_baseline_covariates",
               c("affective_harmonized_v1", "functional_harmonized_v1", "age_baseline", "sex_baseline", "education_baseline")),
  count_branch(dat, "harmonization_eligible_repeated", eligibility = "harmonization_eligibility_v1")
)
branches <- bind_rows(branches, branches %>% group_by(branch) %>% summarise(cohort = "TOTAL", persons = sum(persons), observations = sum(observations), .groups = "drop"))
write_csv(branches, file.path(out_dir, "sample_branch_reconciliation.csv"), na = "")

lagdat <- dat %>% group_by(cohort, subject_id) %>% arrange(year, wave, .by_group = TRUE) %>%
  mutate(lag_affective = lag(affective_harmonized_v1), lag_functional = lag(functional_harmonized_v1),
         delta_year = year - lag(year)) %>% ungroup() %>% filter(delta_year > 0, delta_year <= 10) %>%
  group_by(cohort) %>% mutate(delta_year_c = delta_year - median(delta_year, na.rm = TRUE)) %>% ungroup()

models <- tribble(
  ~outcome, ~estimand, ~formula_text, ~target_term,
  "functional", "population_average_cluster_robust",
  "functional_harmonized_v1 ~ lag_affective + lag_functional + delta_year_c + lag_affective:delta_year_c + age_baseline + sex_baseline + education_baseline",
  "lag_affective",
  "affective", "population_average_cluster_robust",
  "affective_harmonized_v1 ~ lag_functional + lag_affective + delta_year_c + lag_functional:delta_year_c + age_baseline + sex_baseline + education_baseline",
  "lag_functional",
  "functional", "within_person_subject_and_wave_FE",
  "functional_harmonized_v1 ~ lag_affective + lag_functional + delta_year_c + lag_affective:delta_year_c | subject_id + wave",
  "lag_affective",
  "affective", "within_person_subject_and_wave_FE",
  "affective_harmonized_v1 ~ lag_functional + lag_affective + delta_year_c + lag_functional:delta_year_c | subject_id + wave",
  "lag_functional"
)

extract_fixest <- function(fit, cohort_name, outcome, estimand, target_term, n_rows, n_persons) {
  if (inherits(fit, "error")) return(tibble(cohort = cohort_name, outcome, estimand, target_term,
    term = NA_character_, estimate = NA_real_, std_error = NA_real_, p_value = NA_real_,
    n_rows, n_persons, status = "failed", error_message = conditionMessage(fit)))
  terms <- names(coef(fit)); ses <- se(fit); ps <- pvalue(fit)
  tibble(cohort = cohort_name, outcome, estimand, target_term, term = terms,
    estimate = unname(coef(fit)), std_error = unname(ses[terms]), p_value = unname(ps[terms]),
    n_rows = nobs(fit), n_persons, status = "ok", error_message = NA_character_)
}

results <- pmap_dfr(models, function(outcome, estimand, formula_text, target_term) {
  map_dfr(sort(unique(lagdat$cohort)), function(cc) {
    d <- lagdat %>% filter(cohort == cc)
    fml <- as.formula(formula_text)
    model_vars <- intersect(all.vars(fml), names(d))
    analysis_rows <- complete.cases(d[, model_vars, drop = FALSE])
    n_persons <- n_distinct(d$subject_id[analysis_rows])
    fit <- tryCatch(fixest::feols(fml, data = d, vcov = ~subject_id,
      notes = TRUE, warn = TRUE), error = function(e) e)
    extract_fixest(fit, cc, outcome, estimand, target_term, nrow(d), n_persons)
  })
})
write_csv(results, file.path(out_dir, "lagged_cluster_and_FE_results.csv"), na = "")

targets <- results %>% filter(status == "ok", term == target_term) %>% mutate(direction = case_when(estimate > 0 ~ "positive", estimate < 0 ~ "negative", TRUE ~ "zero"))
write_csv(targets, file.path(out_dir, "lagged_target_terms.csv"), na = "")
writeLines(c(
  "population_average_cluster_robust retains baseline covariates and clusters SE by person.",
  "within_person_subject_and_wave_FE removes all time-invariant person-level confounding; baseline covariates are absorbed.",
  "The fixed-effect dynamic model is a sensitivity analysis and is not interpreted causally because lagged outcomes can induce dynamic-panel bias.",
  "External data must not be added until these results are reviewed against the frozen direction/heterogeneity rules."
), file.path(out_dir, "README.txt"))
message("Lagged and sample-reconciliation analysis complete: ", out_dir)
