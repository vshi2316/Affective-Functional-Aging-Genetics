## Build longitudinal submission figures from frozen outputs.
## This script performs no additional model fitting.

options(stringsAsFactors = FALSE)

base_dir <- normalizePath(
  Sys.getenv("DCV_BASE_DIR", unset = ""),
  winslash = "/", mustWork = TRUE
)
submission_dir <- Sys.getenv(
  "LONGITUDINAL_FIGURE_DIR",
  unset = file.path(base_dir, "submission_figures")
)
main_dir <- file.path(submission_dir, "main_figures")
supp_dir <- file.path(submission_dir, "supplementary_figures")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

packages <- c("dplyr", "readr", "ggplot2", "patchwork", "scales", "tidyr")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install before running: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(patchwork); library(scales)
})

robustness_root <- file.path(base_dir, "robustness_analysis")
sample <- read_csv(file.path(robustness_root, "03_lagged_fallback_sample_reconciliation", "sample_branch_reconciliation.csv"), show_col_types = FALSE)
ri_rs <- read_csv(file.path(robustness_root, "04_within_between_manuscript_freeze", "RI_RS_target_effects.csv"), show_col_types = FALSE)
qualification <- read_csv(file.path(robustness_root, "01_longitudinal_robustness", "model_qualification_summary.csv"), show_col_types = FALSE)
lagged <- read_csv(file.path(robustness_root, "04_within_between_manuscript_freeze", "lag_delta_and_FE_target_effects.csv"), show_col_types = FALSE)
within_between <- read_csv(file.path(robustness_root, "04_within_between_manuscript_freeze", "within_between_effects.csv"), show_col_types = FALSE)
meta <- read_csv(file.path(robustness_root, "04_within_between_manuscript_freeze", "within_between_meta_summary.csv"), show_col_types = FALSE)
contrast <- read_csv(file.path(robustness_root, "04_within_between_manuscript_freeze", "between_minus_within_contrasts.csv"), show_col_types = FALSE)

colours <- c("Affective to functional" = "#3B6FB6", "Functional to affective" = "#D05A47")
estimand_colours <- c("Population-average" = "#3569A8", "Within-person" = "#D05A47", "Between-person" = "#2A9D8F")

theme_nc <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = base_size + 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

save_pair <- function(plot, stem, width_mm, height_mm, directory) {
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot,
         width = width_mm / 25.4, height = height_mm / 25.4,
         device = grDevices::cairo_pdf)
  ggsave(file.path(directory, paste0(stem, ".png")), plot,
         width = width_mm / 25.4, height = height_mm / 25.4,
         dpi = 600, bg = "white")
}

direction_label <- function(model_id) {
  ifelse(grepl("functional_on_affective", model_id),
         "Affective to functional", "Functional to affective")
}

sample_total <- sample %>%
  filter(cohort == "TOTAL") %>%
  mutate(
    branch_label = recode(
      branch,
      affective_repeated_plus_baseline_covariates = "Affective branch",
      functional_repeated_plus_baseline_covariates = "Functional branch",
      joint_repeated_plus_baseline_covariates = "Joint branch",
      harmonization_eligible_repeated = "Complete harmonisation branch"
    ),
    branch_label = factor(branch_label, levels = rev(c(
      "Affective branch", "Functional branch", "Joint branch", "Complete harmonisation branch"
    )))
  )

ri <- ri_rs %>%
  filter(grepl("_RI$", model_id), status == "ok") %>%
  mutate(direction = direction_label(model_id), cohort = factor(cohort, levels = rev(c("CHARLS", "ELSA", "HRS", "MHAS", "SHARE"))))

meta_plot <- meta %>%
  mutate(
    direction = ifelse(grepl("functional_on_affective", model_id), "Affective to functional", "Functional to affective"),
    component = ifelse(grepl("between$", term), "Between-person", "Within-person"),
    label = paste(direction, component, sep = "\n"),
    label = factor(label, levels = rev(c(
      "Affective to functional\nBetween-person", "Affective to functional\nWithin-person",
      "Functional to affective\nBetween-person", "Functional to affective\nWithin-person"
    )))
  )

qual_plot <- qualification %>%
  mutate(structure = sub(".*_(RI|RS|CAR1)$", "\\1", model_id)) %>%
  count(structure, status, name = "n_models") %>%
  tidyr::complete(structure = c("RI", "RS", "CAR1"), status = c("ok", "failed"), fill = list(n_models = 0)) %>%
  mutate(structure = factor(structure, levels = c("RI", "RS", "CAR1")), status = recode(status, ok = "Successful", failed = "Failed"))

p_a <- ggplot(sample_total, aes(persons, branch_label)) +
  geom_segment(aes(x = 0, xend = persons, yend = branch_label), linewidth = 0.6, colour = "#B8C2CC") +
  geom_point(size = 2.5, colour = "#3B6FB6") +
  geom_text(aes(label = comma(persons)), hjust = -0.12, size = 2.3) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k"), expand = expansion(mult = c(0, 0.20))) +
  labs(x = "Participants", y = NULL, title = "Parallel analytic branches") + theme_nc()

p_b <- ggplot(ri, aes(estimate, cohort, colour = direction)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std_error, xmax = estimate + 1.96 * std_error), height = 0.14, linewidth = 0.45, position = position_dodge(width = 0.35)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = colours) +
  labs(x = "Random-intercept coefficient (95% CI)", y = NULL, title = "Population-average associations") + theme_nc()

p_c <- ggplot(meta_plot, aes(estimate, label, colour = component)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.14, linewidth = 0.5) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = estimand_colours) +
  labs(x = "Descriptive random-effects estimate (95% CI)", y = NULL, title = "Within-between decomposition") + theme_nc()

p_d <- ggplot(qual_plot, aes(structure, n_models, fill = status)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(n_models == 0, "", n_models)), position = position_stack(vjust = 0.5), size = 2.4, colour = "white") +
  scale_fill_manual(values = c("Successful" = "#2A9D8F", "Failed" = "#B94B5F")) +
  scale_y_continuous(breaks = seq(0, 10, 2), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Correlation structure", y = "Models", title = "Frozen fit audit") + theme_nc()

figure1 <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9), legend.position = "bottom")
save_pair(figure1, "Figure1_longitudinal_estimand_audit", 183, 165, main_dir)

ri_all <- ri_rs %>%
  filter(status == "ok") %>%
  mutate(
    direction = direction_label(model_id),
    structure = sub(".*_(RI|RS)$", "\\1", model_id),
    cohort = factor(cohort, levels = rev(c("CHARLS", "ELSA", "HRS", "MHAS", "SHARE")))
  )
p_s24 <- ggplot(ri_all, aes(estimate, cohort, colour = structure)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std_error, xmax = estimate + 1.96 * std_error), height = 0.12, position = position_dodge(width = 0.35)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.7) +
  facet_wrap(~direction, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = c(RI = "#3B6FB6", RS = "#D05A47")) +
  labs(x = "Coefficient (95% CI)", y = NULL) + theme_nc()
save_pair(p_s24, "Supplementary_Figure_24_RI_RS_robustness", 183, 90, supp_dir)

lag_plot <- lagged %>%
  filter(status == "ok") %>%
  mutate(
    direction = ifelse(outcome == "functional", "Affective to functional", "Functional to affective"),
    estimand_label = recode(estimand,
      population_average_cluster_robust = "Population-average",
      within_person_subject_and_wave_FE = "Within-person"),
    cohort = factor(cohort, levels = rev(c("CHARLS", "ELSA", "HRS", "MHAS", "SHARE")))
  )
p_s25 <- ggplot(lag_plot, aes(estimate, cohort, colour = estimand_label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std_error, xmax = estimate + 1.96 * std_error), height = 0.12, position = position_dodge(width = 0.35)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.7) +
  facet_wrap(~direction, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = estimand_colours) +
  labs(x = "Lagged coefficient (95% CI)", y = NULL) + theme_nc()
save_pair(p_s25, "Supplementary_Figure_25_lagged_estimands", 183, 90, supp_dir)

wb_plot <- within_between %>%
  filter(status == "ok", term %in% c("affective_between", "affective_within", "functional_between", "functional_within")) %>%
  mutate(
    direction = ifelse(grepl("functional_on_affective", model_id), "Affective to functional", "Functional to affective"),
    component = ifelse(grepl("between$", term), "Between-person", "Within-person"),
    cohort = factor(cohort, levels = rev(c("CHARLS", "ELSA", "HRS", "MHAS", "SHARE")))
  )
p_s26 <- ggplot(wb_plot, aes(estimate, cohort, colour = component)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std_error, xmax = estimate + 1.96 * std_error), height = 0.12, position = position_dodge(width = 0.35)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.7) +
  facet_wrap(~direction, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = estimand_colours) +
  labs(x = "Within-between coefficient (95% CI)", y = NULL) + theme_nc()
save_pair(p_s26, "Supplementary_Figure_26_within_between", 183, 90, supp_dir)

main_source <- bind_rows(
  sample_total %>% transmute(panel = "a", branch = as.character(branch_label), persons, observations),
  ri %>% transmute(panel = "b", cohort = as.character(cohort), direction, estimate, std_error, p_value),
  meta_plot %>% transmute(panel = "c", direction, component, estimate, std_error, ci_lb, ci_ub, p_value, I2, tau2),
  qual_plot %>% transmute(panel = "d", structure = as.character(structure), status, n_models)
)
write_csv(main_source, file.path(main_dir, "Figure1_source_data.csv"), na = "")
write_csv(ri_all, file.path(supp_dir, "Supplementary_Figure_24_source_data.csv"), na = "")
write_csv(lag_plot, file.path(supp_dir, "Supplementary_Figure_25_source_data.csv"), na = "")
write_csv(wb_plot, file.path(supp_dir, "Supplementary_Figure_26_source_data.csv"), na = "")
write_csv(contrast, file.path(supp_dir, "between_minus_within_contrasts.csv"), na = "")

writeLines(capture.output(sessionInfo()), file.path(submission_dir, "figure_sessionInfo.txt"))
message("Submission figures written to: ", submission_dir)
