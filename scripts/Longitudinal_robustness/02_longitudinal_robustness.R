## Longitudinal robustness models.
## Existing data only: model qualification, interval-aware reciprocal lagged models,
## and leave-one-cohort-out summaries. Failed models are recorded, never forced.

options(stringsAsFactors = FALSE)
base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir)) stop("Set DCV_BASE_DIR before running this script.")
if (!dir.exists(base_dir)) stop("DCV_BASE_DIR does not exist: ", base_dir)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
input_file <- file.path(base_dir, "analysis_ready_core", "longitudinal_models_v1", "affective_functional_longitudinal_input_v1.rds")
out_dir <- file.path(base_dir, "robustness_analysis", "01_longitudinal_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pkgs <- c("dplyr", "readr", "tibble", "purrr", "nlme", "metafor")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install before running: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble); library(purrr); library(nlme); library(metafor)})

dat <- readRDS(input_file)
na_to_na <- function(x) {x <- suppressWarnings(as.numeric(x)); x[x %in% c(-1,-3,-4,-9,88,99,999,9999)] <- NA_real_; x}
dat <- dat %>% mutate(across(c(year, wave, age_baseline, sex_baseline, education_baseline,
  affective_harmonized_v1, functional_harmonized_v1), na_to_na)) %>% arrange(subject_id, year, wave)

extract_fixed <- function(fit, cohort_name, model_id, n_obs, n_subjects, status = "ok", error_message = NA_character_) {
  if (status != "ok") return(tibble(cohort = cohort_name, model_id, term = NA_character_, estimate = NA_real_, std_error = NA_real_, p_value = NA_real_, n_obs, n_subjects, status, error_message))
  tt <- as.data.frame(summary(fit)$tTable); tt$term <- rownames(tt); rownames(tt) <- NULL
  tibble(cohort = cohort_name, model_id, term = tt$term, estimate = tt$Value, std_error = tt$Std.Error,
         p_value = tt$`p-value`, n_obs, n_subjects, status, error_message)
}

fit_one <- function(d, cohort_name, model_id, formula_text, structure = "RI") {
  d <- d %>% filter(cohort == cohort_name) %>% filter(complete.cases(.))
  minimum_visits <- if (structure == "RS") 3L else 2L
  eligible <- d %>% count(subject_id, name = "n") %>% filter(n >= minimum_visits)
  d <- inner_join(d, eligible, by = "subject_id")
  n_obs <- nrow(d); n_subjects <- n_distinct(d$subject_id)
  if (n_obs < 200 || n_subjects < 100) return(extract_fixed(NULL, cohort_name, model_id, n_obs, n_subjects, "skipped_low_sample", "pre-specified minimum not met"))
  random_form <- if (structure == "RS") ~time_since_baseline_years|subject_id else ~1|subject_id
  cor_form <- if (structure == "CAR1") nlme::corCAR1(form = ~time_since_baseline_years|subject_id) else NULL
  fit <- tryCatch(nlme::lme(fixed = as.formula(formula_text), random = random_form, correlation = cor_form,
    data = d, method = "ML", na.action = na.omit,
    control = nlme::lmeControl(msMaxIter = 200, opt = "optim", returnObject = FALSE)), error = function(e) e)
  if (inherits(fit, "error")) return(extract_fixed(NULL, cohort_name, model_id, n_obs, n_subjects, "failed", conditionMessage(fit)))
  extract_fixed(fit, cohort_name, model_id, n_obs, n_subjects)
}

cohorts <- sort(unique(dat$cohort))
base <- dat %>% select(cohort, subject_id, year, wave, age_baseline, sex_baseline, education_baseline,
  affective_harmonized_v1, functional_harmonized_v1) %>% group_by(subject_id) %>%
  mutate(time_since_baseline_years = year - min(year, na.rm = TRUE)) %>% ungroup()

specs <- tribble(
  ~model_id, ~formula_text, ~structure,
  "functional_on_affective_RI", "functional_harmonized_v1 ~ affective_harmonized_v1 + time_since_baseline_years + age_baseline + sex_baseline + education_baseline", "RI",
  "affective_on_functional_RI", "affective_harmonized_v1 ~ functional_harmonized_v1 + time_since_baseline_years + age_baseline + sex_baseline + education_baseline", "RI",
  "functional_on_affective_RS", "functional_harmonized_v1 ~ affective_harmonized_v1 + time_since_baseline_years + age_baseline + sex_baseline + education_baseline", "RS",
  "affective_on_functional_RS", "affective_harmonized_v1 ~ functional_harmonized_v1 + time_since_baseline_years + age_baseline + sex_baseline + education_baseline", "RS",
  "functional_on_affective_CAR1", "functional_harmonized_v1 ~ affective_harmonized_v1 + time_since_baseline_years + age_baseline + sex_baseline + education_baseline", "CAR1",
  "affective_on_functional_CAR1", "affective_harmonized_v1 ~ functional_harmonized_v1 + time_since_baseline_years + age_baseline + sex_baseline + education_baseline", "CAR1"
)
qual <- pmap_dfr(specs, function(model_id, formula_text, structure) {
  map_dfr(cohorts, ~fit_one(base, .x, model_id, formula_text, structure))
})
write_csv(qual, file.path(out_dir, "model_qualification_effects.csv"), na = "")
write_csv(qual %>% group_by(cohort, model_id, status) %>% summarise(n_terms = sum(!is.na(term)), n_obs = first(n_obs),
  n_subjects = first(n_subjects), error = first(error_message), .groups = "drop"),
  file.path(out_dir, "model_qualification_summary.csv"), na = "")

## Reciprocal lag models. delta_year is modeled explicitly, not inserted as a cosmetic covariate.
lagdat <- base %>% group_by(subject_id) %>% arrange(year, wave, .by_group = TRUE) %>%
  mutate(lag_affective = lag(affective_harmonized_v1), lag_functional = lag(functional_harmonized_v1),
         delta_year = year - lag(year)) %>% ungroup() %>% filter(delta_year > 0, delta_year <= 10) %>%
  group_by(cohort) %>% mutate(delta_year_c = delta_year - median(delta_year, na.rm = TRUE)) %>% ungroup()
lag_specs <- tribble(
  ~model_id, ~formula_text,
  "functional_lag", "functional_harmonized_v1 ~ lag_affective + lag_functional + delta_year_c + lag_affective:delta_year_c + age_baseline + sex_baseline + education_baseline",
  "affective_lag", "affective_harmonized_v1 ~ lag_functional + lag_affective + delta_year_c + lag_functional:delta_year_c + age_baseline + sex_baseline + education_baseline"
)
lag_res <- pmap_dfr(lag_specs, function(model_id, formula_text) map_dfr(cohorts, ~fit_one(lagdat, .x, model_id, formula_text, "RI")))
write_csv(lag_res, file.path(out_dir, "interval_aware_lagged_effects.csv"), na = "")

primary <- qual %>% filter(status == "ok", model_id %in% c("functional_on_affective_RI", "affective_on_functional_RI"),
  term %in% c("affective_harmonized_v1", "functional_harmonized_v1"), is.finite(estimate), is.finite(std_error), std_error > 0)
loo <- map_dfr(unique(primary$model_id), function(mid) {
  d <- primary %>% filter(model_id == mid)
  map_dfr(seq_len(nrow(d)), function(i) {
    x <- d[-i, , drop = FALSE]
    if (nrow(x) < 2) return(tibble(model_id = mid, omitted_cohort = d$cohort[i], k = nrow(x), estimate = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_, p_value = NA_real_, I2 = NA_real_))
    fit <- metafor::rma.uni(x$estimate, sei = x$std_error, method = "REML")
    tibble(model_id = mid, omitted_cohort = d$cohort[i], k = nrow(x), estimate = as.numeric(fit$b), ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, p_value = fit$pval, I2 = fit$I2)
  })
})
write_csv(loo, file.path(out_dir, "leave_one_cohort_out.csv"), na = "")
writeLines(c(
  "RS/CAR1 models are accepted only when status=ok and effect direction/scale are coherent with RI.",
  "Death/attrition is not addressed because the harmonized input has no verified mortality variable.",
  "A converged sensitivity model is not automatically scientifically preferred."
), file.path(out_dir, "longitudinal_readme.txt"))
message("Longitudinal robustness complete: ", out_dir)
