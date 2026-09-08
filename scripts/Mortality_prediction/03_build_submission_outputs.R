#!/usr/bin/env Rscript
# Build manuscript tables, Supplementary Figure 34 and Supplementary Tables S54-S57.
# This stage reads locked aggregate results and does not refit any model.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
output_root <- Sys.getenv(
  "MORTALITY_PREDICTION_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_prediction")
)
landmark_dir <- file.path(output_root, "landmark")
model_dir <- file.path(output_root, "models")
out_dir <- file.path(output_root, "submission")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lock_path <- file.path(out_dir, "submission_outputs_lock.txt")
if (file.exists(lock_path)) stop("Prediction submission outputs are already locked: ", lock_path)

required <- c(
  landmark_lock = file.path(landmark_dir, "landmark_lock.txt"),
  landmark_flow = file.path(landmark_dir, "landmark_flow.csv"),
  history_depth = file.path(landmark_dir, "landmark_history_depth.csv"),
  landmark_integrity = file.path(landmark_dir, "landmark_integrity_audit.csv"),
  endpoint_integrity = file.path(landmark_dir, "endpoint_integrity_audit.csv"),
  input_manifest = file.path(landmark_dir, "input_manifest.csv"),
  model_lock = file.path(model_dir, "mortality_prediction_models_lock.txt"),
  coefficients = file.path(model_dir, "development_coefficients.csv"),
  metrics = file.path(model_dir, "three_model_validation_metrics.csv"),
  bootstrap = file.path(model_dir, "paired_bootstrap_auc_differences.csv"),
  wald = file.path(model_dir, "nested_model2_vs_model1_wald.csv"),
  likelihood_ratio = file.path(model_dir, "nested_model2_vs_model1_likelihood_ratio.csv"),
  support = file.path(model_dir, "baseline_hazard_support_audit.csv"),
  auc_test = file.path(model_dir, "auc_implementation_test.csv"),
  analysis_gates = file.path(model_dir, "analysis_gates.csv")
)
if (any(!file.exists(required))) {
  stop("Missing locked prediction output(s): ", paste(names(required)[!file.exists(required)], collapse = ", "))
}

pkgs <- c(
  "dplyr", "tidyr", "readr", "tibble", "ggplot2", "patchwork", "scales",
  "openxlsx", "flextable", "officer", "svglite", "ragg"
)
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install before running: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
  library(flextable)
  library(officer)
  library(svglite)
  library(ragg)
})

gates <- read_csv(required[["analysis_gates"]], show_col_types = FALSE)
landmark_integrity <- read_csv(required[["landmark_integrity"]], show_col_types = FALSE)
if (!all(gates$pass %in% TRUE) || !all(landmark_integrity$pass %in% TRUE)) {
  stop("Locked analysis gates are not all passed.")
}

flow <- read_csv(required[["landmark_flow"]], show_col_types = FALSE)
history_depth <- read_csv(required[["history_depth"]], show_col_types = FALSE)
endpoint_integrity <- read_csv(required[["endpoint_integrity"]], show_col_types = FALSE)
input_manifest <- read_csv(required[["input_manifest"]], show_col_types = FALSE)
coefficients <- read_csv(required[["coefficients"]], show_col_types = FALSE)
metrics <- read_csv(required[["metrics"]], show_col_types = FALSE)
bootstrap <- read_csv(required[["bootstrap"]], show_col_types = FALSE)
wald <- read_csv(required[["wald"]], show_col_types = FALSE)
likelihood_ratio <- read_csv(required[["likelihood_ratio"]], show_col_types = FALSE)
support <- read_csv(required[["support"]], show_col_types = FALSE)
auc_test <- read_csv(required[["auc_test"]], show_col_types = FALSE)

evaluation_levels <- c(
  "HRS_development", "HRS_internal_validation", "MHAS_external",
  "SHARE_external", "ELSA_external", "CHARLS_external"
)
evaluation_labels <- c(
  HRS_development = "HRS development",
  HRS_internal_validation = "HRS internal validation",
  MHAS_external = "MHAS external",
  SHARE_external = "SHARE external",
  ELSA_external = "ELSA sensitivity",
  CHARLS_external = "CHARLS support-only"
)
evaluation_roles <- c(
  HRS_development = "Development",
  HRS_internal_validation = "Internal validation",
  MHAS_external = "External validation",
  SHARE_external = "External validation",
  ELSA_external = "Sensitivity",
  CHARLS_external = "Support-only"
)
model_labels <- c(
  model_0 = "Model 0: demographic",
  model_1 = "Model 1: current scores",
  model_2 = "Model 2: history-deviation"
)
comparison_labels <- c(
  delta_1_minus_0 = "Model 1 - Model 0",
  delta_2_minus_1 = "Model 2 - Model 1"
)

metrics <- metrics %>%
  mutate(
    evaluation_set = paste(cohort, split, sep = "_"),
    evaluation_label = unname(evaluation_labels[evaluation_set]),
    evaluation_role = unname(evaluation_roles[evaluation_set]),
    model_label = unname(model_labels[model])
  )
bootstrap <- bootstrap %>%
  mutate(
    evaluation_label = unname(evaluation_labels[evaluation_set]),
    evaluation_role = unname(evaluation_roles[evaluation_set]),
    comparison_label = unname(comparison_labels[comparison])
  )

format_ci <- function(estimate, low, high, digits = 4) {
  paste0(
    formatC(estimate, format = "f", digits = digits),
    " (", formatC(low, format = "f", digits = digits),
    " to ", formatC(high, format = "f", digits = digits), ")"
  )
}

## Main Table 2 ---------------------------------------------------------------
table2_sets <- c("HRS_internal_validation", "MHAS_external", "SHARE_external", "ELSA_external")
table2_metric <- metrics %>%
  filter(evaluation_set %in% table2_sets) %>%
  select(evaluation_set, evaluation_label, model, subjects, events_5y, auc_5y) %>%
  pivot_wider(names_from = model, values_from = auc_5y) %>%
  arrange(match(evaluation_set, table2_sets))
table2_delta <- bootstrap %>%
  filter(evaluation_set %in% table2_sets) %>%
  transmute(
    evaluation_set,
    comparison,
    formatted = format_ci(estimate, ci_low, ci_high, digits = 4)
  ) %>%
  pivot_wider(names_from = comparison, values_from = formatted)
table2 <- table2_metric %>%
  left_join(table2_delta, by = "evaluation_set") %>%
  transmute(
    Cohort = evaluation_label,
    Participants = subjects,
    `Five-year deaths` = events_5y,
    `AUC Model 0` = model_0,
    `AUC Model 1` = model_1,
    `AUC Model 2` = model_2,
    `Delta AUC Model 1 - Model 0 (95% CI)` = delta_1_minus_0,
    `Delta AUC Model 2 - Model 1 (95% CI)` = delta_2_minus_1
  )
write_csv(table2, file.path(out_dir, "Table2_FiveYear_Mortality_Discrimination.csv"), na = "")

table2_wb <- createWorkbook(creator = "Affective-Functional Aging Genetics")
addWorksheet(table2_wb, "Table 2", gridLines = FALSE)
writeData(table2_wb, "Table 2", "Table 2 | Five-year mortality discrimination across internal and external evaluation cohorts", startRow = 1)
mergeCells(table2_wb, "Table 2", cols = 1:ncol(table2), rows = 1)
writeData(table2_wb, "Table 2", table2, startRow = 3, withFilter = TRUE)
title_style <- createStyle(
  fontName = "Arial", fontSize = 11, fontColour = "#FFFFFF",
  fgFill = "#1F4E78", textDecoration = "bold", halign = "left", valign = "center"
)
header_style <- createStyle(
  fontName = "Arial", fontSize = 9, fgFill = "#D9EAF7",
  textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE,
  border = "bottom", borderColour = "#5B9BD5"
)
body_style <- createStyle(fontName = "Arial", fontSize = 9, valign = "center")
addStyle(table2_wb, "Table 2", title_style, rows = 1, cols = 1:ncol(table2), gridExpand = TRUE)
addStyle(table2_wb, "Table 2", header_style, rows = 3, cols = 1:ncol(table2), gridExpand = TRUE)
addStyle(table2_wb, "Table 2", body_style, rows = 4:(3 + nrow(table2)), cols = 1:ncol(table2), gridExpand = TRUE)
addStyle(
  table2_wb, "Table 2", createStyle(numFmt = "0.000"),
  rows = 4:(3 + nrow(table2)), cols = 4:6, gridExpand = TRUE, stack = TRUE
)
setColWidths(table2_wb, "Table 2", cols = 1:8, widths = c(23, 13, 14, 12, 12, 12, 27, 27))
setRowHeights(table2_wb, "Table 2", rows = 1, heights = 24)
freezePane(table2_wb, "Table 2", firstActiveRow = 4)
saveWorkbook(
  table2_wb,
  file.path(out_dir, "Table2_FiveYear_Mortality_Discrimination.xlsx"),
  overwrite = TRUE
)

table2_display <- table2 %>%
  mutate(across(starts_with("AUC"), ~formatC(.x, format = "f", digits = 3)))
ft <- flextable(table2_display) %>%
  theme_booktabs() %>%
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 8.5, part = "all") %>%
  bold(part = "header") %>%
  align(j = 2:6, align = "center", part = "all") %>%
  valign(valign = "center", part = "all") %>%
  autofit()
doc <- read_docx() %>%
  body_add_par(
    "Table 2. Five-year mortality discrimination across internal and external evaluation cohorts",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft) %>%
  body_add_par(
    paste(
      "Model 0 included sex, standardized education and calendar year, with attained age represented by the Cox baseline hazard.",
      "Model 1 added current affective and functional scores. Model 2 used affective history, affective deviation, functional history and functional deviation.",
      "Differences in AUC and 95% confidence intervals were estimated by 1,000 paired participant-level bootstrap repetitions."
    ),
    style = "Normal"
  )
print(doc, target = file.path(out_dir, "Table2_FiveYear_Mortality_Discrimination.docx"))

## Supplementary Figure 34 ----------------------------------------------------
figure_levels <- setdiff(evaluation_levels, "CHARLS_external")
figure_metrics <- metrics %>%
  filter(evaluation_set %in% figure_levels) %>%
  mutate(
    evaluation_set = factor(evaluation_set, levels = rev(figure_levels)),
    evaluation_label = factor(evaluation_label, levels = rev(unname(evaluation_labels[figure_levels]))),
    model = factor(model, levels = names(model_labels)),
    model_label = factor(model_label, levels = unname(model_labels))
  )
auc_ranges <- figure_metrics %>%
  group_by(evaluation_set, evaluation_label) %>%
  summarise(auc_min = min(auc_5y), auc_max = max(auc_5y), .groups = "drop")

model_palette <- c(
  "Model 0: demographic" = "#767676",
  "Model 1: current scores" = "#3775BA",
  "Model 2: history-deviation" = "#1B9E77"
)

p_a <- ggplot(figure_metrics, aes(x = auc_5y, y = evaluation_label)) +
  geom_segment(
    data = auc_ranges,
    aes(x = auc_min, xend = auc_max, y = evaluation_label, yend = evaluation_label),
    inherit.aes = FALSE, colour = "#D8D8D8", linewidth = 0.7
  ) +
  geom_point(aes(colour = model_label, shape = model_label), size = 2.4, stroke = 0.6) +
  scale_colour_manual(values = model_palette, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_x_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.04, 0.04))) +
  labs(x = "IPCW AUC at five years", y = NULL) +
  theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    plot.margin = margin(5, 5, 5, 5)
  )

figure_delta <- bootstrap %>%
  filter(comparison == "delta_2_minus_1", evaluation_set %in% figure_levels) %>%
  mutate(
    evaluation_set = factor(evaluation_set, levels = rev(figure_levels)),
    evaluation_label = factor(evaluation_label, levels = rev(unname(evaluation_labels[figure_levels])))
  )
p_b <- ggplot(
  figure_delta,
  aes(x = estimate, y = evaluation_label, xmin = ci_low, xmax = ci_high)
) +
  geom_vline(xintercept = 0, colour = "#767676", linetype = 2, linewidth = 0.45) +
  geom_errorbarh(height = 0, colour = "#1B9E77", linewidth = 0.7) +
  geom_point(shape = 21, size = 2.6, stroke = 0.6, fill = "#1B9E77", colour = "white") +
  scale_x_continuous(labels = label_number(accuracy = 0.001)) +
  labs(x = "Delta AUC: Model 2 - Model 1", y = NULL) +
  theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.margin = margin(5, 5, 5, 5)
  )

figure <- p_a + p_b +
  plot_layout(widths = c(1.15, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 9, face = "bold", family = "Arial"),
    legend.position = "bottom"
  )

figure_base <- file.path(out_dir, "Supplementary_Figure_34_FiveYear_Mortality_Prediction")
width_mm <- 183
height_mm <- 95
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite(paste0(figure_base, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()
cairo_pdf(paste0(figure_base, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(figure)
dev.off()
agg_tiff(
  paste0(figure_base, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure)
dev.off()
agg_png(
  paste0(figure_base, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300
)
print(figure)
dev.off()

figure_source <- bind_rows(
  figure_metrics %>% transmute(
    panel = "a", evaluation_set = as.character(evaluation_set),
    evaluation_label = as.character(evaluation_label), model = as.character(model),
    model_label = as.character(model_label), estimate = auc_5y,
    ci_low = NA_real_, ci_high = NA_real_, interval = "none"
  ),
  figure_delta %>% transmute(
    panel = "b", evaluation_set = as.character(evaluation_set),
    evaluation_label = as.character(evaluation_label), model = comparison,
    model_label = comparison_label, estimate, ci_low, ci_high,
    interval = "paired participant bootstrap 95% CI"
  )
)
write_csv(
  figure_source,
  file.path(out_dir, "Source_Data_Supplementary_Figure_34.csv"),
  na = ""
)

figure_legend <- paste(
  "Supplementary Fig. 34 | Five-year mortality discrimination across internal and external evaluation cohorts.",
  "a, Inverse-probability-of-censoring-weighted five-year AUC for the demographic, current-score and history-deviation models.",
  "b, Difference in five-year AUC between the history-deviation and current-score models; points show observed differences and lines show 95% confidence intervals from 1,000 paired participant-level bootstrap repetitions.",
  "All coefficients and cumulative baseline hazards were estimated in the HRS development sample and frozen for subsequent evaluation.",
  "ELSA supplied sensitivity evidence; CHARLS support-only results are reported in Supplementary Tables S56 and S57."
)
writeLines(
  figure_legend,
  file.path(out_dir, "Supplementary_Figure_34_Legend.txt"),
  useBytes = TRUE
)

## Supplementary Tables S54-S57 ----------------------------------------------
support_summary <- support %>%
  group_by(cohort, split) %>%
  summarise(
    prediction_evaluable = sum(subjects[prediction_supported %in% TRUE]),
    excluded_outside_hazard_support = sum(subjects[prediction_supported %in% FALSE]),
    .groups = "drop"
  )
s54a <- flow %>%
  left_join(support_summary, by = c("cohort", "split")) %>%
  mutate(
    evaluation_set = paste(cohort, split, sep = "_"),
    evaluation_label = unname(evaluation_labels[evaluation_set]),
    validation_role = unname(evaluation_roles[evaluation_set])
  ) %>%
  transmute(
    Cohort = evaluation_label,
    `Validation role` = validation_role,
    `Eligible landmark participants` = subjects,
    `Prediction-evaluable participants` = prediction_evaluable,
    `Excluded outside HRS hazard support` = excluded_outside_hazard_support,
    `Deaths within five years` = events_5y,
    `Observed five-year status` = observed_5y,
    `Censored before five years` = censored_before_5y,
    `Median follow-up, years` = median_followup,
    `Landmark age, minimum` = min_landmark_age,
    `Landmark age, median` = median_landmark_age,
    `Landmark age, maximum` = max_landmark_age
  )
s54b <- history_depth %>%
  mutate(
    evaluation_set = paste(cohort, split, sep = "_"),
    evaluation_label = unname(evaluation_labels[evaluation_set])
  ) %>%
  transmute(
    Cohort = evaluation_label,
    N = subjects,
    `Prior measurements, minimum` = min_prior_measurements,
    `Prior measurements, Q1` = q1_prior_measurements,
    `Prior measurements, median` = median_prior_measurements,
    `Prior measurements, Q3` = q3_prior_measurements,
    `Prior measurements, maximum` = max_prior_measurements,
    `One prior measurement, n` = one_prior_measurement,
    `One prior measurement, %` = one_prior_measurement_pct,
    `Two or more prior measurements, n` = two_or_more_prior_measurements,
    `Two or more prior measurements, %` = two_or_more_prior_measurements_pct
  )

s55a <- coefficients %>%
  mutate(model_label = unname(model_labels[model])) %>%
  transmute(
    Model = model_label, Term = term, Beta = estimate, SE = std_error,
    Z = z, `P value` = p_value, HR = hazard_ratio,
    `95% CI lower` = ci_low, `95% CI upper` = ci_high
  )
s55b <- bind_rows(
  wald %>% transmute(
    Comparison = "Model 2 versus Model 1", Test = "Joint Wald", df, Statistic = statistic,
    `P value` = p_value, `H_A - D_A` = contrast_HA_DA, `H_F - D_F` = contrast_HF_DF
  ),
  likelihood_ratio %>% transmute(
    Comparison = "Model 2 versus Model 1", Test = "Likelihood-ratio", df, Statistic = statistic,
    `P value` = p_value, `H_A - D_A` = NA_real_, `H_F - D_F` = NA_real_
  )
)

s56 <- metrics %>%
  arrange(match(evaluation_set, evaluation_levels), model) %>%
  transmute(
    Cohort = evaluation_label,
    `Validation role` = evaluation_role,
    Model = model_label,
    Participants = subjects,
    `Five-year deaths` = events_5y,
    `Five-year AUC` = auc_5y,
    `Five-year Brier score` = brier_5y,
    `Mean predicted risk` = mean_predicted_risk,
    `Observed five-year risk` = observed_risk,
    `Calibration slope` = calibration_slope,
    `Calibration intercept` = calibration_intercept
  )

s57a <- bootstrap %>%
  arrange(match(evaluation_set, evaluation_levels), comparison) %>%
  transmute(
    `Evaluation set` = evaluation_label,
    `Validation role` = evaluation_role,
    Comparison = comparison_label,
    `Delta AUC` = estimate,
    `95% CI lower` = ci_low,
    `95% CI upper` = ci_high,
    `Bootstrap repetitions` = bootstrap_replicates
  )
s57b <- bind_rows(
  landmark_integrity %>% transmute(
    `Audit domain` = "Landmark integrity", Scope = "All cohorts",
    Check = check, Value = as.character(value), Pass = pass, Detail = ""
  ),
  endpoint_integrity %>% transmute(
    `Audit domain` = "Endpoint integrity", Scope = cohort,
    Check = "multiple/nonterminal event subjects",
    Value = paste0(subjects_with_multiple_event_rows, "/", nonterminal_event_subjects),
    Pass = subjects_with_multiple_event_rows == 0 & nonterminal_event_subjects == 0,
    Detail = paste0("event subjects=", event_subjects)
  ),
  support %>% transmute(
    `Audit domain` = "Baseline-hazard support", Scope = paste(cohort, split),
    Check = ifelse(prediction_supported, "supported", "excluded outside support"),
    Value = as.character(subjects), Pass = TRUE, Detail = ""
  ),
  auc_test %>% transmute(
    `Audit domain` = "AUC implementation", Scope = "Unit test", Check = test,
    Value = as.character(absolute_difference), Pass = pass,
    Detail = paste0("rank=", rank, "; brute_force=", brute_force)
  ),
  gates %>% transmute(
    `Audit domain` = "Analysis gate", Scope = "Model stage", Check = gate,
    Value = as.character(value), Pass = pass, Detail = ""
  ),
  input_manifest %>% transmute(
    `Audit domain` = "Input manifest", Scope = input, Check = "MD5",
    Value = md5, Pass = exists, Detail = basename(path)
  )
)

write_csv(s54a, file.path(out_dir, "Supplementary_Table_S54_Panel_A.csv"), na = "")
write_csv(s54b, file.path(out_dir, "Supplementary_Table_S54_Panel_B.csv"), na = "")
write_csv(s55a, file.path(out_dir, "Supplementary_Table_S55_Panel_A.csv"), na = "")
write_csv(s55b, file.path(out_dir, "Supplementary_Table_S55_Panel_B.csv"), na = "")
write_csv(s56, file.path(out_dir, "Supplementary_Table_S56.csv"), na = "")
write_csv(s57a, file.path(out_dir, "Supplementary_Table_S57_Panel_A.csv"), na = "")
write_csv(s57b, file.path(out_dir, "Supplementary_Table_S57_Panel_B.csv"), na = "")

deep_blue <- "#1F4E78"
section_blue <- "#5B9BD5"
light_blue <- "#D9EAF7"
title_style <- createStyle(
  fontName = "Arial", fontSize = 11, fontColour = "#FFFFFF",
  fgFill = deep_blue, textDecoration = "bold", halign = "left", valign = "center"
)
section_style <- createStyle(
  fontName = "Arial", fontSize = 9, fontColour = "#FFFFFF",
  fgFill = section_blue, textDecoration = "bold", halign = "left", valign = "center"
)
header_style <- createStyle(
  fontName = "Arial", fontSize = 9, fgFill = light_blue,
  textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE,
  border = "bottom", borderColour = section_blue
)
body_style <- createStyle(fontName = "Arial", fontSize = 9, valign = "center")

write_panel <- function(wb, sheet, start_row, label, data) {
  writeData(wb, sheet, label, startRow = start_row, startCol = 1)
  mergeCells(wb, sheet, cols = 1:ncol(data), rows = start_row)
  addStyle(wb, sheet, section_style, rows = start_row, cols = 1:ncol(data), gridExpand = TRUE)
  writeData(wb, sheet, data, startRow = start_row + 1L, withFilter = FALSE)
  addStyle(wb, sheet, header_style, rows = start_row + 1L, cols = 1:ncol(data), gridExpand = TRUE)
  if (nrow(data)) {
    addStyle(
      wb, sheet, body_style,
      rows = (start_row + 2L):(start_row + 1L + nrow(data)),
      cols = 1:ncol(data), gridExpand = TRUE
    )
  }
  start_row + nrow(data) + 3L
}

write_supplementary_sheets <- function(wb) {
  sheets <- names(wb)
  if (any(c("S54", "S55", "S56", "S57") %in% sheets)) {
    stop("Target workbook already contains one or more of S54-S57.")
  }
  titles <- c(
    S54 = "Age-50 landmark prediction sample, follow-up, and history depth",
    S55 = "HRS development coefficients and nested model comparisons",
    S56 = "Five-year mortality discrimination, prediction error, and calibration across three nested models",
    S57 = "Paired participant-bootstrap differences in five-year AUC and analytic audit"
  )
  for (sheet in names(titles)) addWorksheet(wb, sheet, gridLines = FALSE)

  write_title <- function(sheet, title, columns) {
    writeData(wb, sheet, paste0("Supplementary Table ", sheet, " | ", title), startRow = 1)
    mergeCells(wb, sheet, cols = 1:columns, rows = 1)
    addStyle(wb, sheet, title_style, rows = 1, cols = 1:columns, gridExpand = TRUE)
    setRowHeights(wb, sheet, rows = 1, heights = 24)
  }

  write_title("S54", titles[["S54"]], max(ncol(s54a), ncol(s54b)))
  row <- write_panel(wb, "S54", 3, "A. Landmark sample and follow-up", s54a)
  write_panel(wb, "S54", row, "B. History depth", s54b)

  write_title("S55", titles[["S55"]], max(ncol(s55a), ncol(s55b)))
  row <- write_panel(wb, "S55", 3, "A. HRS development coefficients", s55a)
  write_panel(wb, "S55", row, "B. Nested model comparisons", s55b)

  write_title("S56", titles[["S56"]], ncol(s56))
  write_panel(wb, "S56", 3, "Three-model performance", s56)

  write_title("S57", titles[["S57"]], max(ncol(s57a), ncol(s57b)))
  row <- write_panel(wb, "S57", 3, "A. Paired participant-bootstrap differences in AUC", s57a)
  write_panel(wb, "S57", row, "B. Analytic and reproducibility audit", s57b)

  sheet_columns <- c(
    S54 = max(ncol(s54a), ncol(s54b)),
    S55 = max(ncol(s55a), ncol(s55b)),
    S56 = ncol(s56),
    S57 = max(ncol(s57a), ncol(s57b))
  )
  for (sheet in names(sheet_columns)) {
    setColWidths(wb, sheet, cols = seq_len(sheet_columns[[sheet]]), widths = "auto")
    widths <- sapply(seq_len(sheet_columns[[sheet]]), function(i) {
      values <- readWorkbook(wb, sheet, cols = i, rows = 1:200, colNames = FALSE)
      max_chars <- max(nchar(as.character(unlist(values))), na.rm = TRUE)
      min(max(28, max_chars + 2), 42)
    })
    setColWidths(wb, sheet, cols = seq_len(sheet_columns[[sheet]]), widths = widths)
    freezePane(wb, sheet, firstActiveRow = 4)
  }
  titles
}

standalone_wb <- createWorkbook(creator = "Affective-Functional Aging Genetics")
addWorksheet(standalone_wb, "Content", gridLines = FALSE)
titles <- write_supplementary_sheets(standalone_wb)
content_data <- tibble(Table = names(titles), Title = unname(titles))
writeData(standalone_wb, "Content", content_data, startRow = 1, withFilter = TRUE)
addStyle(standalone_wb, "Content", header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
addStyle(standalone_wb, "Content", body_style, rows = 2:(nrow(content_data) + 1), cols = 1:2, gridExpand = TRUE)
setColWidths(standalone_wb, "Content", cols = c(1, 2), widths = c(12, 72))
saveWorkbook(
  standalone_wb,
  file.path(out_dir, "Supplementary_Tables_S54-S57.xlsx"),
  overwrite = TRUE
)

base_supplementary <- Sys.getenv("MORTALITY_PREDICTION_SUPPLEMENTARY_TABLES_INPUT", unset = "")
if (nzchar(base_supplementary) && file.exists(base_supplementary)) {
  integrated_wb <- loadWorkbook(base_supplementary)
  titles <- write_supplementary_sheets(integrated_wb)
  if ("Content" %in% names(integrated_wb)) {
    existing_content <- readWorkbook(integrated_wb, "Content")
    start_row <- nrow(existing_content) + 2L
    append_content <- data.frame(Table = names(titles), Title = unname(titles))
    writeData(
      integrated_wb, "Content", append_content,
      startRow = start_row, startCol = 1, colNames = FALSE
    )
    addStyle(
      integrated_wb, "Content", body_style,
      rows = start_row:(start_row + nrow(append_content) - 1L), cols = 1:2, gridExpand = TRUE
    )
  }
  saveWorkbook(
    integrated_wb,
    file.path(out_dir, "Supplementary_Tables.xlsx"),
    overwrite = TRUE
  )
}

## QA, reporting text and manifest -------------------------------------------
reporting_text <- c(
  "METHODS HEADING: Secondary five-year mortality prediction analysis",
  "RESULTS HEADING: Five-year mortality discrimination",
  "MAIN TABLE: Table 2. Five-year mortality discrimination across internal and external evaluation cohorts",
  "SUPPLEMENTARY FIGURE: Supplementary Fig. 34. Five-year mortality discrimination across internal and external evaluation cohorts",
  "SUPPLEMENTARY TABLES: S54-S57",
  "PRIMARY RESULT: Current affective-functional scores improved five-year discrimination across HRS, MHAS, SHARE and ELSA.",
  "INCREMENTAL RESULT: History-deviation separation retained additional discrimination in HRS internal validation and SHARE external validation."
)
writeLines(reporting_text, file.path(out_dir, "Mortality_Prediction_Reporting_Summary.txt"), useBytes = TRUE)

qa_record <- c(
  paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  "Backend: R only (ggplot2 and patchwork)",
  paste0("Figure size: ", width_mm, " mm x ", height_mm, " mm"),
  "Figure exports: SVG, PDF, 600 dpi LZW TIFF and 300 dpi PNG preview",
  paste0("Table 2 rows: ", nrow(table2)),
  paste0("S54 panel rows: ", nrow(s54a), " and ", nrow(s54b)),
  paste0("S55 panel rows: ", nrow(s55a), " and ", nrow(s55b)),
  paste0("S56 rows: ", nrow(s56)),
  paste0("S57 panel rows: ", nrow(s57a), " and ", nrow(s57b)),
  "No model was fitted or selected in the submission-output stage."
)
writeLines(qa_record, file.path(out_dir, "Mortality_Prediction_QA_Record.txt"), useBytes = TRUE)

output_files <- list.files(out_dir, full.names = TRUE)
output_files <- output_files[basename(output_files) != basename(lock_path)]
manifest <- tibble(
  file = basename(output_files),
  bytes = file.info(output_files)$size,
  md5 = unname(tools::md5sum(output_files))
)
write_csv(manifest, file.path(out_dir, "Mortality_Prediction_Output_Manifest.csv"), na = "")
script_path <- Sys.getenv(
  "MORTALITY_PREDICTION_PRESENTATION_SCRIPT",
  unset = file.path(project_dir, "scripts", "Mortality_prediction", "03_build_submission_outputs.R")
)
writeLines(
  c(
    paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
    paste0("Landmark lock MD5: ", unname(tools::md5sum(required[["landmark_lock"]]))),
    paste0("Model lock MD5: ", unname(tools::md5sum(required[["model_lock"]]))),
    paste0("Presentation script MD5: ", if (file.exists(script_path)) unname(tools::md5sum(script_path)) else NA_character_),
    paste0("Output files before manifest: ", nrow(manifest)),
    "Submission outputs were generated from locked aggregate results without model refitting."
  ),
  lock_path,
  useBytes = TRUE
)

message("Mortality prediction submission outputs completed: ", out_dir)
