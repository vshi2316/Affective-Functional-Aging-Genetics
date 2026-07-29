# ==============================================================================
# Track A compact analysis script 6
# Final longitudinal analysis components.
# ==============================================================================

# ------------------------------------------------------------------------------
# Component: prepare_results_tables.R
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
closure_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis_closure")
harm_dir <- file.path(base_dir, "analysis_ready_core", "harmonized_affective_functional")
out_dir <- file.path(base_dir, "analysis_ready_core", "results_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(long_dir, "affective_functional_longitudinal_effects.csv"),
  file.path(long_dir, "basic_cognition_longitudinal_effects.csv"),
  file.path(meta_dir, "meta_summary_affective_functional.csv"),
  file.path(meta_dir, "meta_summary_basic_cognition.csv"),
  file.path(closure_dir, "cohort_directionality_summary.csv"),
  file.path(harm_dir, "pooled_affective_functional_harmonized.rds")
)
missing_required <- required[!file.exists(required)]
if (length(missing_required) > 0) {
  stop("Missing required input file(s): ", paste(missing_required, collapse = "; "))
}

aff_effects <- read_csv(file.path(long_dir, "affective_functional_longitudinal_effects.csv"), show_col_types = FALSE)
cog_effects <- read_csv(file.path(long_dir, "basic_cognition_longitudinal_effects.csv"), show_col_types = FALSE)
meta_aff <- read_csv(file.path(meta_dir, "meta_summary_affective_functional.csv"), show_col_types = FALSE)
meta_cog <- read_csv(file.path(meta_dir, "meta_summary_basic_cognition.csv"), show_col_types = FALSE)
directionality <- read_csv(file.path(closure_dir, "cohort_directionality_summary.csv"), show_col_types = FALSE)
harmonized <- readRDS(file.path(harm_dir, "pooled_affective_functional_harmonized.rds")) %>%
  as_tibble()

analytic_rows <- harmonized %>%
  filter(
    harmonization_eligibility == 1,
    !is.na(year), !is.na(age), !is.na(sex_female),
    !is.na(education_level_clean),
    !is.na(affective_harmonized), !is.na(functional_harmonized)
  )

analytic_people <- analytic_rows %>%
  arrange(cohort, id, year, wave) %>%
  group_by(cohort, id) %>%
  summarise(
    baseline_age = first(age),
    female = first(sex_female),
    first_year = min(year),
    last_year = max(year),
    repeated_observations = n(),
    .groups = "drop"
  ) %>%
  filter(repeated_observations >= 2) %>%
  mutate(follow_up_years = last_year - first_year)

analytic_rows_repeated <- analytic_rows %>%
  semi_join(analytic_people, by = c("cohort", "id"))

baseline_summary <- analytic_people %>%
  group_by(cohort) %>%
  summarise(
    joint_complete_case_participants = n(),
    baseline_age_mean = mean(baseline_age),
    baseline_age_sd = sd(baseline_age),
    women_pct = 100 * mean(female == 1),
    follow_up_median = median(follow_up_years),
    follow_up_q1 = as.numeric(quantile(follow_up_years, 0.25)),
    follow_up_q3 = as.numeric(quantile(follow_up_years, 0.75)),
    .groups = "drop"
  ) %>%
  left_join(
    analytic_rows_repeated %>%
      group_by(cohort) %>%
      summarise(
        joint_complete_case_observations = n(),
        first_calendar_year = min(year),
        last_calendar_year = max(year),
        .groups = "drop"
      ),
    by = "cohort"
  )

table2_main <- bind_rows(meta_aff, meta_cog) %>%
  mutate(
    result_label = case_when(
      model_id == "functional_trajectory" ~ "Functional trajectory",
      model_id == "affective_trajectory" ~ "Affective trajectory",
      model_id == "functional_on_affective" ~ "Functional on affective",
      model_id == "affective_on_functional" ~ "Affective on functional",
      model_id == "cognition_trajectory" ~ "Cognition trajectory",
      TRUE ~ model_id
    )
  ) %>%
  select(
    result_label, model_id, term, k,
    pooled_estimate, pooled_se, ci_lb, ci_ub, p_value, i2, tau2, status
  ) %>%
  arrange(result_label)

supplementary_terms <- c(
  "time_since_baseline_years",
  "affective_harmonized",
  "functional_harmonized"
)

supplementary_estimates <- bind_rows(aff_effects, cog_effects) %>%
  filter(status == "ok", term %in% supplementary_terms) %>%
  arrange(model_id, cohort)

table1_cohort_characteristics <- bind_rows(
  aff_effects %>%
    group_by(cohort) %>%
    summarise(
      max_affective_n_obs = max(n_obs, na.rm = TRUE),
      max_affective_n_subjects = max(n_subjects, na.rm = TRUE),
      .groups = "drop"
    ),
  cog_effects %>%
    group_by(cohort) %>%
    summarise(
      max_cognition_n_obs = max(n_obs, na.rm = TRUE),
      max_cognition_n_subjects = max(n_subjects, na.rm = TRUE),
      .groups = "drop"
    )
) %>%
  group_by(cohort) %>%
  summarise(across(everything(), ~ suppressWarnings(max(.x, na.rm = TRUE))), .groups = "drop") %>%
  left_join(
    directionality %>%
      group_by(cohort) %>%
      summarise(
        n_models_with_directionality = n_distinct(model_id),
        .groups = "drop"
      ),
    by = "cohort"
  ) %>%
  left_join(baseline_summary, by = "cohort") %>%
  mutate(
    country_or_region = recode(
      cohort,
      CHARLS = "China", ELSA = "England", HRS = "United States",
      MHAS = "Mexico", SHARE = "Europe and Israel"
    ),
      analysis_role = if_else(
        cohort %in% c("ELSA", "HRS"),
        "Cohort stratum 1",
        "Cohort stratum 2"
      )
  ) %>%
  arrange(cohort)

write_csv(table1_cohort_characteristics, file.path(out_dir, "table1_cohort_characteristics.csv"), na = "")
write_csv(table2_main, file.path(out_dir, "table2_main_meta_results.csv"), na = "")
write_csv(supplementary_estimates, file.path(out_dir, "supplementary_cohort_specific_estimates.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: build_closure_status_report.R
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
meta_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
closure_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis_closure")
fig_dir <- file.path(base_dir, "analysis_ready_core", "figures")
temporal_dir <- file.path(base_dir, "analysis_ready_core", "temporal_sensitivity")
weighted_dir <- file.path(base_dir, "analysis_ready_core", "weighted_sensitivity")
out_dir <- file.path(base_dir, "analysis_ready_core", "results_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

checks <- tibble(
  criterion = c(
    "Main meta-analysis completed",
    "Leave-one-cohort-out completed",
    "Forest plots completed",
    "Cohort directionality table completed",
    "SHARE cognition gap explained",
    "At least one temporal or robustness sensitivity completed"
  ),
  status = c(
    all(file.exists(
      file.path(meta_dir, "meta_summary_affective_functional.csv"),
      file.path(meta_dir, "meta_summary_basic_cognition.csv")
    )),
    all(file.exists(
      file.path(closure_dir, "leave_one_cohort_out_affective_functional.csv"),
      file.path(closure_dir, "leave_one_cohort_out_basic_cognition.csv")
    )),
    all(file.exists(
      file.path(fig_dir, "forest_affective_trajectory.png"),
      file.path(fig_dir, "forest_functional_trajectory.png"),
      file.path(fig_dir, "forest_functional_on_affective.png"),
      file.path(fig_dir, "forest_affective_on_functional.png"),
      file.path(fig_dir, "forest_cognition_trajectory.png")
    )),
    file.exists(file.path(closure_dir, "cohort_directionality_summary.csv")),
    all(file.exists(
      file.path(closure_dir, "share_cognition_meta_audit.csv"),
      file.path(closure_dir, "share_cognition_meta_audit.txt")
    )),
    any(file.exists(
      file.path(temporal_dir, "temporal_sensitivity_hrs_elsa_summary.csv"),
      file.path(weighted_dir, "weighted_sensitivity_summary.csv")
    ))
  )
) %>%
  mutate(status = ifelse(status, "done", "pending"))

all_done <- all(checks$status == "done")

report_lines <- c(
  "Five-cohort closure status report v1",
  paste0("Overall closure status: ", ifelse(all_done, "READY_FOR_FINAL_DESIGN_FREEZE", "NOT_YET_CLOSED")),
  ""
)

for (i in seq_len(nrow(checks))) {
  report_lines <- c(report_lines, paste0("- ", checks$criterion[[i]], ": ", checks$status[[i]]))
}

write_csv(checks, file.path(out_dir, "closure_status_report.csv"), na = "")
writeLines(report_lines, file.path(out_dir, "closure_status_report.txt"))

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: build_trackA_evidence_package.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
domain_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
mi_dir <- file.path(base_dir, "analysis_ready_core", "measurement_invariance")
meta_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis")
closure_dir <- file.path(base_dir, "analysis_ready_core", "meta_analysis_closure")
temporal_dir <- file.path(base_dir, "analysis_ready_core", "temporal_sensitivity")
out_dir <- file.path(base_dir, "analysis_ready_core", "trackA_evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(domain_dir, "cohort_readiness_summary.csv"),
  file.path(domain_dir, "domain_coverage_report.csv"),
  file.path(mi_dir, "affective_functional_mi_fit_summary.csv"),
  file.path(meta_dir, "meta_summary_affective_functional.csv"),
  file.path(meta_dir, "meta_summary_basic_cognition.csv"),
  file.path(closure_dir, "cohort_directionality_summary.csv"),
  file.path(closure_dir, "leave_one_cohort_out_affective_functional.csv"),
  file.path(closure_dir, "leave_one_cohort_out_basic_cognition.csv")
)
missing_required <- required[!file.exists(required)]
if (length(missing_required) > 0) {
  stop("Missing required input file(s): ", paste(missing_required, collapse = "; "))
}

readiness <- read_csv(file.path(domain_dir, "cohort_readiness_summary.csv"), show_col_types = FALSE)
coverage <- read_csv(file.path(domain_dir, "domain_coverage_report.csv"), show_col_types = FALSE)
mi_fit <- read_csv(file.path(mi_dir, "affective_functional_mi_fit_summary.csv"), show_col_types = FALSE)
meta_summary <- bind_rows(
  read_csv(file.path(meta_dir, "meta_summary_affective_functional.csv"), show_col_types = FALSE),
  read_csv(file.path(meta_dir, "meta_summary_basic_cognition.csv"), show_col_types = FALSE)
)
directionality <- read_csv(file.path(closure_dir, "cohort_directionality_summary.csv"), show_col_types = FALSE)
loo_aff <- read_csv(file.path(closure_dir, "leave_one_cohort_out_affective_functional.csv"), show_col_types = FALSE)
loo_cog <- read_csv(file.path(closure_dir, "leave_one_cohort_out_basic_cognition.csv"), show_col_types = FALSE)

temporal_summary_file <- file.path(temporal_dir, "temporal_sensitivity_hrs_elsa.csv")
temporal_summary <- if (file.exists(temporal_summary_file)) {
  read_csv(temporal_summary_file, show_col_types = FALSE)
} else {
  tibble()
}

write_csv(readiness, file.path(out_dir, "eTable_trackA_cohort_readiness.csv"), na = "")
write_csv(coverage, file.path(out_dir, "eTable_trackA_domain_coverage.csv"), na = "")
write_csv(mi_fit, file.path(out_dir, "eTable_trackA_mi_fit_summary.csv"), na = "")
write_csv(meta_summary, file.path(out_dir, "eTable_trackA_meta_summary.csv"), na = "")
write_csv(directionality, file.path(out_dir, "eTable_trackA_directionality.csv"), na = "")
write_csv(bind_rows(loo_aff, loo_cog), file.path(out_dir, "eTable_trackA_leave_one_out.csv"), na = "")
if (nrow(temporal_summary) > 0) {
  write_csv(temporal_summary, file.path(out_dir, "eTable_trackA_temporal_sensitivity.csv"), na = "")
}

meta_plot_data <- meta_summary %>%
  mutate(
    result_label = case_when(
      model_id == "functional_trajectory" ~ "Functional trajectory",
      model_id == "affective_trajectory" ~ "Affective trajectory",
      model_id == "functional_on_affective" ~ "Functional on affective",
      model_id == "affective_on_functional" ~ "Affective on functional",
      model_id == "cognition_trajectory" ~ "Cognition trajectory",
      TRUE ~ model_id
    )
  )

p1 <- ggplot(meta_plot_data, aes(x = reorder(result_label, pooled_estimate), y = pooled_estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lb, ymax = ci_ub), width = 0.15) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  coord_flip() +
  labs(
    title = "Pooled Track A Effects",
    x = NULL,
    y = "Pooled estimate (95% CI)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "figure_trackA_pooled_effects.png"),
  plot = p1,
  width = 9,
  height = 5.5,
  dpi = 300
)

time_terms <- directionality %>%
  filter(term == "time_since_baseline_years") %>%
  mutate(
    result_label = case_when(
      model_id == "functional_trajectory" ~ "Functional trajectory",
      model_id == "affective_trajectory" ~ "Affective trajectory",
      model_id == "cognition_trajectory" ~ "Cognition trajectory",
      TRUE ~ model_id
    )
  )

p2 <- ggplot(time_terms, aes(x = cohort, y = estimate, color = result_label)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.4)) +
  geom_errorbar(
    aes(ymin = ci_lb, ymax = ci_ub),
    width = 0.12,
    position = position_dodge(width = 0.4)
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  labs(
    title = "Cohort-specific Time Slopes",
    x = NULL,
    y = "Estimate (95% CI)",
    color = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "figure_trackA_cohort_time_slopes.png"),
  plot = p2,
  width = 9,
  height = 5.5,
  dpi = 300
)

loo_plot_data <- bind_rows(loo_aff, loo_cog) %>%
  mutate(
    result_label = case_when(
      model_id == "functional_trajectory" ~ "Functional trajectory",
      model_id == "affective_trajectory" ~ "Affective trajectory",
      model_id == "functional_on_affective" ~ "Functional on affective",
      model_id == "affective_on_functional" ~ "Affective on functional",
      model_id == "cognition_trajectory" ~ "Cognition trajectory",
      TRUE ~ model_id
    )
  )

p3 <- ggplot(loo_plot_data, aes(x = excluded_cohort, y = pooled_estimate, color = result_label, group = result_label)) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~ result_label, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  labs(
    title = "Leave-one-cohort-out Stability",
    x = "Excluded cohort",
    y = "Re-estimated pooled effect",
    color = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave(
  filename = file.path(out_dir, "figure_trackA_leave_one_out.png"),
  plot = p3,
  width = 11,
  height = 7,
  dpi = 300
)

if (nrow(temporal_summary) > 0) {
  temporal_focus <- temporal_summary %>%
    filter(term %in% c("lag_affective", "lag_functional"))

  p4 <- ggplot(temporal_focus, aes(x = cohort, y = estimate, fill = term)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
    labs(
      title = "Temporal Sensitivity: Lagged Cross-domain Terms",
      x = NULL,
      y = "Lagged estimate",
      fill = NULL
    ) +
    theme_minimal(base_size = 12)

  ggsave(
    filename = file.path(out_dir, "figure_trackA_temporal_lagged_terms.png"),
    plot = p4,
    width = 8,
    height = 5,
    dpi = 300
  )
}

message("Done.")
message("Outputs written to: ", out_dir)

})


