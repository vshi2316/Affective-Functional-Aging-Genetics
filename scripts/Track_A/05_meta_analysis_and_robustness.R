# ==============================================================================
# Track A compact analysis script 5
# Final longitudinal analysis components.
# ==============================================================================

# ------------------------------------------------------------------------------
# Component: run_meta_analysis_affective_functional.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it before running this script.")
  }
}

check_pkg("metafor")

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "longitudinal_models")
out_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "affective_functional_longitudinal_effects.csv")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

effects <- read_csv(input_file, show_col_types = FALSE)
effects <- effects %>%
  mutate(
    model_id = as.character(model_id),
    term = as.character(term),
    status = as.character(status)
  )

target_terms <- tribble(
  ~model_id, ~term,
  "affective_trajectory", "time_since_baseline_years",
  "functional_trajectory", "time_since_baseline_years",
  "functional_on_affective", "affective_harmonized",
  "affective_on_functional", "functional_harmonized"
) %>%
  mutate(
    model_id = as.character(model_id),
    term = as.character(term)
  )

meta_input <- effects %>%
  inner_join(target_terms, by = c("model_id", "term")) %>%
  filter(status == "ok", !is.na(estimate), !is.na(std.error), std.error > 0)

run_meta <- function(df, model_id_value, term_value) {
  dat <- df %>% filter(model_id == model_id_value, term == term_value)
  if (nrow(dat) < 2) {
    return(tibble(
      model_id = model_id_value,
      term = term_value,
      k = nrow(dat),
      pooled_estimate = NA_real_,
      pooled_se = NA_real_,
      ci_lb = NA_real_,
      ci_ub = NA_real_,
      p_value = NA_real_,
      i2 = NA_real_,
      tau2 = NA_real_,
      status = "insufficient_cohorts"
    ))
  }

  fit <- metafor::rma.uni(yi = dat$estimate, sei = dat$std.error, method = "REML", test = "knha")

  tibble(
    model_id = model_id_value,
    term = term_value,
    k = nrow(dat),
    pooled_estimate = as.numeric(fit$b[[1]]),
    pooled_se = fit$se,
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    p_value = fit$pval,
    i2 = fit$I2,
    tau2 = fit$tau2,
    status = "ok"
  )
}

meta_summary <- bind_rows(lapply(seq_len(nrow(target_terms)), function(i) {
  run_meta(meta_input, target_terms$model_id[[i]], target_terms$term[[i]])
}))

write_csv(meta_input, file.path(out_dir, "meta_input_affective_functional.csv"), na = "")
write_csv(meta_summary, file.path(out_dir, "meta_summary_affective_functional.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: run_meta_analysis_basic_cognition.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it before running this script.")
  }
}

check_pkg("metafor")

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "longitudinal_models")
out_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "basic_cognition_longitudinal_effects.csv")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

effects <- read_csv(input_file, show_col_types = FALSE)
effects <- effects %>%
  mutate(
    model_id = as.character(model_id),
    term = as.character(term),
    status = as.character(status)
  )

target_terms <- tribble(
  ~model_id, ~term,
  "cognition_trajectory", "time_since_baseline_years"
) %>%
  mutate(
    model_id = as.character(model_id),
    term = as.character(term)
  )

meta_input <- effects %>%
  inner_join(target_terms, by = c("model_id", "term")) %>%
  filter(status == "ok", !is.na(estimate), !is.na(std.error), std.error > 0)

run_meta <- function(df, model_id_value, term_value) {
  dat <- df %>% filter(model_id == model_id_value, term == term_value)
  if (nrow(dat) < 2) {
    return(tibble(
      model_id = model_id_value,
      term = term_value,
      k = nrow(dat),
      pooled_estimate = NA_real_,
      pooled_se = NA_real_,
      ci_lb = NA_real_,
      ci_ub = NA_real_,
      p_value = NA_real_,
      i2 = NA_real_,
      tau2 = NA_real_,
      status = "insufficient_cohorts"
    ))
  }

  fit <- metafor::rma.uni(yi = dat$estimate, sei = dat$std.error, method = "REML", test = "knha")

  tibble(
    model_id = model_id_value,
    term = term_value,
    k = nrow(dat),
    pooled_estimate = as.numeric(fit$b[[1]]),
    pooled_se = fit$se,
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    p_value = fit$pval,
    i2 = fit$I2,
    tau2 = fit$tau2,
    status = "ok"
  )
}

meta_summary <- bind_rows(lapply(seq_len(nrow(target_terms)), function(i) {
  run_meta(meta_input, target_terms$model_id[[i]], target_terms$term[[i]])
}))

write_csv(meta_input, file.path(out_dir, "meta_input_basic_cognition.csv"), na = "")
write_csv(meta_summary, file.path(out_dir, "meta_summary_basic_cognition.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: run_leave_one_cohort_out_meta.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it before running this script.")
  }
}

check_pkg("metafor")

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
out_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis_closure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

aff_file <- file.path(input_dir, "meta_input_affective_functional.csv")
cog_file <- file.path(input_dir, "meta_input_basic_cognition.csv")

if (!file.exists(aff_file)) stop("Missing required input file: ", aff_file)
if (!file.exists(cog_file)) stop("Missing required input file: ", cog_file)

run_loo <- function(df) {
  df <- df %>%
    mutate(
      cohort = as.character(cohort),
      model_id = as.character(model_id),
      term = as.character(term)
    )

  bind_rows(lapply(split(df, interaction(df$model_id, df$term, drop = TRUE)), function(dat) {
    dat <- as_tibble(dat)
    model_id_value <- dat$model_id[[1]]
    term_value <- dat$term[[1]]
    cohorts <- sort(unique(dat$cohort))

    bind_rows(lapply(cohorts, function(excluded) {
      loo_dat <- dat %>% filter(cohort != excluded)
      if (nrow(loo_dat) < 2) {
        return(tibble(
          model_id = model_id_value,
          term = term_value,
          excluded_cohort = excluded,
          k_remaining = nrow(loo_dat),
          pooled_estimate = NA_real_,
          pooled_se = NA_real_,
          ci_lb = NA_real_,
          ci_ub = NA_real_,
          p_value = NA_real_,
          i2 = NA_real_,
          tau2 = NA_real_,
          status = "insufficient_cohorts"
        ))
      }

      fit <- metafor::rma.uni(yi = loo_dat$estimate, sei = loo_dat$std.error, method = "REML", test = "knha")
      tibble(
        model_id = model_id_value,
        term = term_value,
        excluded_cohort = excluded,
        k_remaining = nrow(loo_dat),
        pooled_estimate = as.numeric(fit$b[[1]]),
        pooled_se = fit$se,
        ci_lb = fit$ci.lb,
        ci_ub = fit$ci.ub,
        p_value = fit$pval,
        i2 = fit$I2,
        tau2 = fit$tau2,
        status = "ok"
      )
    }))
  }))
}

aff_input <- read_csv(aff_file, show_col_types = FALSE)
cog_input <- read_csv(cog_file, show_col_types = FALSE)

aff_loo <- run_loo(aff_input)
cog_loo <- run_loo(cog_input)

write_csv(aff_loo, file.path(out_dir, "leave_one_cohort_out_affective_functional.csv"), na = "")
write_csv(cog_loo, file.path(out_dir, "leave_one_cohort_out_basic_cognition.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: summarize_cohort_directionality.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
long_dir <- file.path(base_dir, "analysis_ready_core", "longitudinal_models")
meta_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
out_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis_closure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

aff_file <- file.path(long_dir, "affective_functional_longitudinal_effects.csv")
cog_file <- file.path(long_dir, "basic_cognition_longitudinal_effects.csv")
meta_aff_file <- file.path(meta_dir, "meta_summary_affective_functional.csv")
meta_cog_file <- file.path(meta_dir, "meta_summary_basic_cognition.csv")

required <- c(aff_file, cog_file, meta_aff_file, meta_cog_file)
missing_required <- required[!file.exists(required)]
if (length(missing_required) > 0) {
  stop("Missing required input file(s): ", paste(missing_required, collapse = "; "))
}

bind_effects <- bind_rows(
  read_csv(aff_file, show_col_types = FALSE),
  read_csv(cog_file, show_col_types = FALSE)
) %>%
  filter(status == "ok", !is.na(term), !is.na(estimate), !is.na(std.error), std.error > 0) %>%
  mutate(
    model_id = as.character(model_id),
    term = as.character(term),
    cohort = as.character(cohort),
    ci_lb = estimate - 1.96 * std.error,
    ci_ub = estimate + 1.96 * std.error,
    direction = case_when(
      estimate > 0 ~ "positive",
      estimate < 0 ~ "negative",
      TRUE ~ "null"
    )
  )

meta_summary <- bind_rows(
  read_csv(meta_aff_file, show_col_types = FALSE),
  read_csv(meta_cog_file, show_col_types = FALSE)
) %>%
  mutate(
    model_id = as.character(model_id),
    term = as.character(term),
    pooled_direction = case_when(
      pooled_estimate > 0 ~ "positive",
      pooled_estimate < 0 ~ "negative",
      TRUE ~ "null"
    )
  ) %>%
  select(model_id, term, pooled_estimate, pooled_direction, i2, tau2, p_value)

directionality <- bind_effects %>%
  left_join(meta_summary, by = c("model_id", "term")) %>%
  mutate(
    same_as_pooled_direction = case_when(
      is.na(pooled_direction) ~ NA,
      direction == pooled_direction ~ "yes",
      TRUE ~ "no"
    )
  ) %>%
  select(
    cohort, model_id, outcome, term,
    estimate, std.error, ci_lb, ci_ub, direction,
    pooled_estimate, pooled_direction, same_as_pooled_direction,
    p.value, i2, tau2, n_obs, n_subjects
  ) %>%
  arrange(model_id, cohort)

write_csv(directionality, file.path(out_dir, "cohort_directionality_summary.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: build_forest_plots.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it before running this script.")
  }
}

check_pkg("metafor")

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
out_dir <- file.path(base_dir, "analysis_ready_core", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

aff_file <- file.path(input_dir, "meta_input_affective_functional.csv")
cog_file <- file.path(input_dir, "meta_input_basic_cognition.csv")

if (!file.exists(aff_file)) stop("Missing required input file: ", aff_file)
if (!file.exists(cog_file)) stop("Missing required input file: ", cog_file)

label_map <- c(
  affective_trajectory = "Affective trajectory",
  functional_trajectory = "Functional trajectory",
  functional_on_affective = "Functional on affective",
  affective_on_functional = "Affective on functional",
  cognition_trajectory = "Cognition trajectory"
)

file_map <- c(
  affective_trajectory = "forest_affective_trajectory.png",
  functional_trajectory = "forest_functional_trajectory.png",
  functional_on_affective = "forest_functional_on_affective.png",
  affective_on_functional = "forest_affective_on_functional.png",
  cognition_trajectory = "forest_cognition_trajectory.png"
)

plot_one <- function(dat, model_id_value) {
  model_dat <- dat %>% filter(model_id == model_id_value) %>% arrange(cohort)
  if (nrow(model_dat) < 2) {
    return(invisible(NULL))
  }

  fit <- metafor::rma.uni(yi = model_dat$estimate, sei = model_dat$std.error, method = "REML", test = "knha")
  out_file <- file.path(out_dir, file_map[[model_id_value]])

  png(out_file, width = 1600, height = 1000, res = 150)
  metafor::forest(
    fit,
    slab = model_dat$cohort,
    xlab = "Estimate",
    main = label_map[[model_id_value]],
    addpred = TRUE,
    header = c("Cohort", "Estimate [95% CI]")
  )
  dev.off()
}

aff_input <- read_csv(aff_file, show_col_types = FALSE) %>%
  mutate(model_id = as.character(model_id), cohort = as.character(cohort))
cog_input <- read_csv(cog_file, show_col_types = FALSE) %>%
  mutate(model_id = as.character(model_id), cohort = as.character(cohort))

all_input <- bind_rows(aff_input, cog_input)

for (model_id_value in unique(all_input$model_id)) {
  plot_one(all_input, model_id_value)
}

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: audit_share_cognition_meta_exclusion.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
long_dir <- file.path(base_dir, "analysis_ready_core", "longitudinal_models")
meta_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
out_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis_closure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

effects_file <- file.path(long_dir, "basic_cognition_longitudinal_effects.csv")
input_file <- file.path(long_dir, "basic_cognition_longitudinal_input.rds")
meta_input_file <- file.path(meta_dir, "meta_input_basic_cognition.csv")

required <- c(effects_file, input_file, meta_input_file)
missing_required <- required[!file.exists(required)]
if (length(missing_required) > 0) {
  stop("Missing required input file(s): ", paste(missing_required, collapse = "; "))
}

effects <- read_csv(effects_file, show_col_types = FALSE) %>%
  mutate(cohort = as.character(cohort))
meta_input <- read_csv(meta_input_file, show_col_types = FALSE) %>%
  mutate(cohort = as.character(cohort))
input_dat <- readRDS(input_file) %>%
  mutate(cohort = as.character(cohort))

cohort_input_summary <- input_dat %>%
  group_by(cohort) %>%
  summarise(
    n_rows_input = n(),
    n_subjects_input = n_distinct(paste(cohort, id, sep = "::")),
    non_missing_cognition = sum(!is.na(cognition_harmonized)),
    .groups = "drop"
  )

audit_tbl <- effects %>%
  filter(model_id == "cognition_trajectory") %>%
  select(cohort, status, n_obs, n_subjects, error_message) %>%
  left_join(cohort_input_summary, by = "cohort") %>%
  mutate(
    included_in_meta = ifelse(cohort %in% meta_input$cohort, "yes", "no"),
    audit_reason = case_when(
      included_in_meta == "yes" ~ "Included in cognition meta-analysis",
      status == "failed" ~ paste0("Model failed: ", ifelse(is.na(error_message), "unknown error", error_message)),
      status == "skipped_low_sample" ~ "Excluded due to low repeated-measures sample after filtering",
      TRUE ~ "Excluded before meta input construction; inspect cohort-specific cognition availability and model status"
    )
  ) %>%
  arrange(cohort)

share_row <- audit_tbl %>% filter(cohort == "SHARE")
audit_note <- if (nrow(share_row) == 0) {
  c(
    "SHARE cognition audit",
    "SHARE was not present in the cohort-specific cognition effects table.",
    "This means exclusion occurred before or during longitudinal cognition modeling.",
    "Inspect pooled_basic_cognition_harmonization_input.rds construction and SHARE cognition harmonization coverage."
  )
} else {
  c(
    "SHARE cognition audit",
    paste0("Meta included: ", share_row$included_in_meta[[1]]),
    paste0("Model status: ", share_row$status[[1]]),
    paste0("Audit reason: ", share_row$audit_reason[[1]])
  )
}

write_csv(audit_tbl, file.path(out_dir, "share_cognition_meta_audit.csv"), na = "")
writeLines(audit_note, file.path(out_dir, "share_cognition_meta_audit.txt"))

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: run_temporal_sensitivity_hrs_elsa.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it before running this script.")
  }
}

check_pkg("nlme")

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "longitudinal_models")
out_dir <- file.path(base_dir, "analysis_ready_core", "temporal_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "affective_functional_longitudinal_input.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

dat <- readRDS(input_file)

required_cols <- c(
  "cohort", "subject_id", "year",
  "affective_harmonized", "functional_harmonized",
  "time_since_baseline_years", "age_baseline", "sex_baseline", "education_baseline"
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns in affective_functional_longitudinal_input.rds: ",
       paste(missing_cols, collapse = ", "))
}

make_lagged <- function(df) {
  df %>%
    arrange(subject_id, year, wave) %>%
    group_by(subject_id) %>%
    mutate(
      lag_affective = dplyr::lag(affective_harmonized),
      lag_functional = dplyr::lag(functional_harmonized),
      delta_year = year - dplyr::lag(year)
    ) %>%
    ungroup()
}

extract_fixed <- function(fit, cohort_name, model_id, outcome_name, n_obs, n_subjects) {
  coef_tbl <- as.data.frame(summary(fit)$tTable)
  coef_tbl$term <- rownames(coef_tbl)
  rownames(coef_tbl) <- NULL
  names(coef_tbl) <- sub(" ", "_", names(coef_tbl), fixed = TRUE)

  tibble(
    cohort = cohort_name,
    model_id = model_id,
    outcome = outcome_name,
    term = coef_tbl$term,
    estimate = coef_tbl$Value,
    std.error = coef_tbl$Std.Error,
    statistic = coef_tbl$t.value,
    p.value = coef_tbl$p.value,
    n_obs = n_obs,
    n_subjects = n_subjects,
    status = "ok",
    error_message = NA_character_
  )
}

fit_model <- function(data, formula_text) {
  tryCatch(
    nlme::lme(
      fixed = as.formula(formula_text),
      random = ~ 1 | subject_id,
      data = data,
      method = "ML",
      na.action = na.omit,
      control = nlme::lmeControl(
        returnObject = TRUE,
        msMaxIter = 200,
        opt = "optim"
      )
    ),
    error = function(e) e
  )
}

run_one <- function(df, cohort_name, model_id, outcome_name, formula_text, required_vars) {
  dat_cohort <- df %>%
    filter(cohort == cohort_name) %>%
    filter(if_all(all_of(required_vars), ~ !is.na(.x))) %>%
    filter(if_all(all_of(c("age_baseline", "sex_baseline", "education_baseline")), ~ !is.na(.x)))

  subject_counts <- dat_cohort %>% count(subject_id, name = "n_rows")
  dat_cohort <- dat_cohort %>%
    inner_join(subject_counts %>% filter(n_rows >= 2), by = "subject_id")

  n_obs <- nrow(dat_cohort)
  n_subjects <- dplyr::n_distinct(dat_cohort$subject_id)

  if (n_obs < 200 || n_subjects < 100) {
    return(tibble(
      cohort = cohort_name,
      model_id = model_id,
      outcome = outcome_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_obs = n_obs,
      n_subjects = n_subjects,
      status = "skipped_low_sample",
      error_message = "Insufficient repeated-measures sample after lagged filtering."
    ))
  }

  fit <- fit_model(dat_cohort, formula_text)
  if (inherits(fit, "error")) {
    return(tibble(
      cohort = cohort_name,
      model_id = model_id,
      outcome = outcome_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_obs = n_obs,
      n_subjects = n_subjects,
      status = "failed",
      error_message = conditionMessage(fit)
    ))
  }

  extract_fixed(fit, cohort_name, model_id, outcome_name, n_obs, n_subjects)
}

model_specs <- tribble(
  ~model_id, ~outcome, ~formula_text, ~required_vars,
  "functional_on_lagged_affective", "functional_harmonized",
  "functional_harmonized ~ lag_affective + lag_functional + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "functional_harmonized|lag_affective|lag_functional|time_since_baseline_years",
  "affective_on_lagged_functional", "affective_harmonized",
  "affective_harmonized ~ lag_functional + lag_affective + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "affective_harmonized|lag_functional|lag_affective|time_since_baseline_years"
)

target_cohorts <- c("HRS", "ELSA")
lagged_dat <- dat %>%
  filter(cohort %in% target_cohorts) %>%
  make_lagged()

results <- bind_rows(lapply(seq_len(nrow(model_specs)), function(i) {
  req <- strsplit(model_specs$required_vars[[i]], "\\|")[[1]]
  bind_rows(lapply(target_cohorts, function(cc) {
    run_one(
      lagged_dat,
      cc,
      model_specs$model_id[[i]],
      model_specs$outcome[[i]],
      model_specs$formula_text[[i]],
      req
    )
  }))
}))

summary_tbl <- results %>%
  group_by(cohort, model_id, outcome, status) %>%
  summarise(
    n_terms = sum(!is.na(term)),
    n_obs = dplyr::first(n_obs),
    n_subjects = dplyr::first(n_subjects),
    error_message = dplyr::first(error_message),
    .groups = "drop"
  )

notes <- c(
  "Temporal sensitivity analysis in HRS and ELSA",
  "These are lagged mixed models, not full RI-CLPM.",
  "They are intended as a pragmatic temporal sensitivity layer for affective-functional co-evolution."
)

write_csv(results, file.path(out_dir, "temporal_sensitivity_hrs_elsa.csv"), na = "")
write_csv(summary_tbl, file.path(out_dir, "temporal_sensitivity_hrs_elsa_summary.csv"), na = "")
writeLines(notes, file.path(out_dir, "temporal_sensitivity_hrs_elsa_notes.txt"))

message("Done.")
message("Outputs written to: ", out_dir)

})


