
#!/usr/bin/env Rscript
# Build the dementia-related clinical-event figure and aggregate source data.
# No model is fitted in this stage.

options(stringsAsFactors = FALSE)

required <- c(
  "dplyr", "tidyr", "readr", "tibble", "ggplot2", "patchwork",
  "scales", "svglite", "ragg", "openxlsx"
)
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install required packages before running: ", paste(missing_pkgs, collapse = ", "))
}

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
results_dir <- Sys.getenv(
  "DEMENTIA_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "dementia_bridge")
)
input_dir <- Sys.getenv(
  "DEMENTIA_EVENT_DATA_OUT",
  unset = file.path(results_dir, "event_dataset")
)
model_dir <- Sys.getenv(
  "DEMENTIA_MODEL_OUT",
  unset = file.path(results_dir, "models")
)
out_dir <- Sys.getenv(
  "DEMENTIA_FIGURE_OUT",
  unset = file.path(results_dir, "figure")
)
source_dir <- file.path(out_dir, "Source_Data_CSV")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  support = file.path(input_dir, "model_support.csv"),
  censoring = file.path(input_dir, "observation_censoring_summary.csv"),
  death_audit = file.path(input_dir, "death_after_last_assessment_summary.csv"),
  endpoint_quality = file.path(input_dir, "endpoint_quality.csv"),
  flow = file.path(model_dir, "01_analysis_flow.csv"),
  main_coef = file.path(model_dir, "02_main_cox_coefficients.csv"),
  contrasts = file.path(model_dir, "03_history_deviation_contrasts.csv"),
  joint = file.path(model_dir, "04_joint_constraint_tests.csv"),
  comparison = file.path(model_dir, "05_M1_M2_model_comparisons.csv"),
  sensitivity_coef = file.path(model_dir, "06_sensitivity_cox_coefficients.csv"),
  sensitivity_contrast = file.path(model_dir, "07_sensitivity_contrasts.csv"),
  sensitivity_joint = file.path(model_dir, "08_sensitivity_joint_tests.csv"),
  observation_coef = file.path(model_dir, "09_observation_rule_coefficient_comparison.csv"),
  observation_joint = file.path(model_dir, "10_observation_rule_joint_test_comparison.csv"),
  cloglog = file.path(model_dir, "11_cloglog_sensitivity.csv"),
  status = file.path(model_dir, "12_model_status.csv")
)
if (any(!file.exists(paths))) {
  stop("Missing result file(s): ", paste(names(paths)[!file.exists(paths)], collapse = ", "))
}

read_result <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
dat <- lapply(paths, read_result)

if (any(!dat$status$converged) || any(dat$status$warning_count > 0, na.rm = TRUE)) {
  stop("At least one model failed or retained a warning. Figure not generated.")
}

main_flow <- dat$flow |>
  dplyr::mutate(
    cohort = factor(.data$cohort, levels = c("ELSA", "SHARE")),
    outcome_label = dplyr::if_else(
      .data$cohort == "ELSA",
      "First survey-reported dementia",
      "First survey-reported Alzheimer disease or dementia"
    )
  )

main_coef <- dat$main_coef |>
  dplyr::mutate(cohort = factor(.data$cohort, levels = c("ELSA", "SHARE")))

panel_b_data <- main_coef |>
  dplyr::filter(
    .data$model == "M1_current_scores",
    .data$term %in% c("current_affective", "current_functional")
  ) |>
  dplyr::mutate(
    domain = dplyr::if_else(grepl("affective", .data$term), "Affective", "Functional"),
    display_term = factor(
      .data$term,
      levels = c("current_functional", "current_affective"),
      labels = c("Current functional score", "Current affective score")
    ),
    estimate_label = sprintf("%.2f (%.2f–%.2f)", .data$hazard_ratio, .data$conf_low_95, .data$conf_high_95)
  )

panel_c_data <- main_coef |>
  dplyr::filter(
    .data$model == "M2_history_deviation",
    .data$term %in% c("H_A", "D_A", "H_F", "D_F")
  ) |>
  dplyr::mutate(
    domain = dplyr::if_else(.data$term %in% c("H_A", "D_A"), "Affective", "Functional"),
    component = dplyr::if_else(.data$term %in% c("H_A", "H_F"), "History", "Current deviation"),
    display_term = factor(
      .data$term,
      levels = c("D_F", "H_F", "D_A", "H_A"),
      labels = c(
        "Functional deviation", "Functional history",
        "Affective deviation", "Affective history"
      )
    ),
    estimate_label = sprintf("%.2f (%.2f–%.2f)", .data$hazard_ratio, .data$conf_low_95, .data$conf_high_95)
  )

panel_d_data <- dat$contrasts |>
  dplyr::mutate(
    cohort = factor(.data$cohort, levels = c("ELSA", "SHARE")),
    domain = dplyr::if_else(grepl("affective", .data$contrast), "Affective", "Functional"),
    display_term = factor(.data$domain, levels = c("Functional", "Affective")),
    estimate_label = sprintf(
      "%.2f (%.2f–%.2f)",
      .data$ratio_of_hazard_ratios, .data$conf_low_95, .data$conf_high_95
    )
  ) |>
  dplyr::left_join(
    dat$joint |>
      dplyr::select(.data$cohort, joint_p = .data$p_value),
    by = "cohort"
  )

sup_sensitivity <- dat$sensitivity_coef |>
  dplyr::filter(.data$term %in% c("H_A", "D_A", "H_F", "D_F")) |>
  dplyr::mutate(
    cohort = factor(.data$cohort, levels = c("ELSA", "SHARE")),
    domain = dplyr::if_else(.data$term %in% c("H_A", "D_A"), "Affective", "Functional"),
    component = dplyr::if_else(.data$term %in% c("H_A", "H_F"), "History", "Current deviation"),
    term_label = factor(
      .data$term,
      levels = c("H_A", "D_A", "H_F", "D_F"),
      labels = c("Affective history", "Affective deviation", "Functional history", "Functional deviation")
    ),
    analysis_label = dplyr::recode(
      .data$analysis,
      exclude_first_followup_wave_outcomes = "Exclude first follow-up events",
      cognition_adjusted_complete_case = "Cognition adjusted",
      exclude_same_date_outcome_death = "Exclude same-date records",
      direct_interview_endpoint = "Direct interviews only",
      same_or_previous_wave_baseline_negative = "Same or previous-wave negative",
      outcome_stable_no_later_negative = "No later negative reversal",
      strict_last_assessment_censoring = "Strict observation censoring",
      strict_last_assessment_censoring_same_or_previous_wave = "Strict censoring + previous-wave negative"
    )
  )

format_p <- function(p) {
  superscript <- function(x) {
    chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(x))
  }
  vapply(p, function(z) {
    if (!is.finite(z)) return("P unavailable")
    if (z >= 0.001) return(paste0("P = ", formatC(z, format = "f", digits = 3)))
    exponent <- floor(log10(z))
    mantissa <- z / (10^exponent)
    paste0("P = ", formatC(mantissa, format = "f", digits = 2), " × 10", superscript(exponent))
  }, character(1))
}

palette <- c(
  Affective = "#D97745",
  Functional = "#287F8E",
  neutral_dark = "#30343B",
  neutral_mid = "#747B83",
  neutral_light = "#D9DEE3",
  strict = "#7A6F9B"
)

theme_pub <- function(base_size = 7.2) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#30343B"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "#30343B"),
      axis.text = ggplot2::element_text(colour = "#30343B"),
      axis.title = ggplot2::element_text(colour = "#30343B"),
      strip.background = ggplot2::element_rect(fill = "#F2F4F6", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.6, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.2, colour = "#565D65"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 5, 4, 5, unit = "pt")
    )
}

## Panel A: study-time and observation-rule schematic.
timeline <- tibble::tibble(
  x = c(1.0, 2.0, 3.0, 4.2),
  y = 1,
  label = c("Prior visit 1", "Prior visit 2", "Prior visit 3", "Current visit"),
  type = c("history", "history", "history", "current")
)
p_a <- ggplot2::ggplot() +
  ggplot2::annotate("segment", x = 0.8, xend = 7.4, y = 1, yend = 1,
                    linewidth = 0.45, colour = palette["neutral_dark"],
                    arrow = ggplot2::arrow(length = grid::unit(2.2, "mm"))) +
  ggplot2::geom_point(
    data = timeline,
    ggplot2::aes(.data$x, .data$y, fill = .data$type),
    shape = 21, size = 3.0, stroke = 0.5, colour = "white"
  ) +
  ggplot2::scale_fill_manual(values = c(history = palette["neutral_mid"], current = palette["Functional"])) +
  ggplot2::annotate("text", x = 2.0, y = 1.43, label = "Historical mean (H)",
                    family = "Arial", fontface = "bold", size = 2.55, colour = palette["neutral_dark"]) +
  ggplot2::annotate("segment", x = 0.95, xend = 3.05, y = 1.31, yend = 1.31,
                    linewidth = 0.4, colour = palette["neutral_mid"]) +
  ggplot2::annotate("text", x = 4.2, y = 1.43, label = "Current deviation (D)",
                    family = "Arial", fontface = "bold", size = 2.55, colour = palette["Functional"]) +
  ggplot2::annotate("text", x = timeline$x, y = 0.74, label = timeline$label,
                    family = "Arial", size = 2.15, colour = palette["neutral_dark"]) +
  ggplot2::annotate("segment", x = 5.25, xend = 5.25, y = 0.43, yend = 1.58,
                    linetype = 2, linewidth = 0.4, colour = palette["neutral_mid"]) +
  ggplot2::annotate("text", x = 5.25, y = 0.28,
                    label = "Last observed outcome\nassessment",
                    family = "Arial", size = 2.1, lineheight = 0.92,
                    colour = palette["neutral_mid"]) +
  ggplot2::annotate("point", x = 6.35, y = 1.18, shape = 22, size = 3.0,
                    fill = palette["Affective"], colour = "white", stroke = 0.5) +
  ggplot2::annotate("text", x = 6.35, y = 1.50,
                    label = "First detected dementia-related\noutcome",
                    family = "Arial", size = 2.1, lineheight = 0.92,
                    colour = palette["neutral_dark"]) +
  ggplot2::annotate("point", x = 6.35, y = 0.67, shape = 24, size = 3.1,
                    fill = palette["neutral_dark"], colour = "white", stroke = 0.5) +
  ggplot2::annotate("text", x = 6.35, y = 0.39,
                    label = "Known death within\ndynamic support",
                    family = "Arial", size = 2.1, lineheight = 0.92,
                    colour = palette["neutral_dark"]) +
  ggplot2::coord_cartesian(xlim = c(0.6, 7.55), ylim = c(0.12, 1.74), clip = "off") +
  ggplot2::labs(title = "Dynamic exposure and subsequent clinical-event follow-up") +
  ggplot2::theme_void(base_family = "Arial", base_size = 7.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 7.8, hjust = 0),
    plot.margin = ggplot2::margin(5, 7, 2, 7, unit = "pt")
  )

forest_base <- function(
  data, xvar, low, high, x_limits, x_breaks,
  label_x, label_column_start, title, subtitle
) {
  if (max(data[[high]], na.rm = TRUE) >= label_column_start) {
    stop("The estimate-label column overlaps a confidence interval in ", title, ".")
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data[[xvar]], y = .data$display_term, colour = .data$domain)
  ) +
    ggplot2::annotate(
      "rect", xmin = label_column_start, xmax = Inf,
      ymin = -Inf, ymax = Inf, fill = "#F8F9FA", colour = NA
    ) +
    ggplot2::geom_vline(xintercept = 1, linewidth = 0.4, linetype = 2, colour = palette["neutral_mid"]) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data[[low]], xmax = .data[[high]]),
      height = 0.16, linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_text(
      ggplot2::aes(x = label_x, label = .data$estimate_label),
      hjust = 1, family = "Arial", size = 2.25, colour = palette["neutral_dark"],
      show.legend = FALSE
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$cohort), scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_manual(values = palette[c("Affective", "Functional")]) +
    ggplot2::scale_x_log10(limits = x_limits, breaks = x_breaks) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Hazard ratio (95% CI)", y = NULL) +
    theme_pub() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 6.5),
      panel.spacing.y = grid::unit(7, "pt"),
      strip.text.y.right = ggplot2::element_text(
        angle = 270, face = "bold", margin = ggplot2::margin(2, 2, 2, 2)
      )
    )
}

p_b <- forest_base(
  panel_b_data, "hazard_ratio", "conf_low_95", "conf_high_95",
  x_limits = c(0.72, 2.70), x_breaks = c(0.75, 1, 1.25, 1.5, 2),
  label_x = 2.58, label_column_start = 1.72,
  title = "Current-score model",
  subtitle = "M1: current affective and functional scores"
)

p_c <- forest_base(
  panel_c_data, "hazard_ratio", "conf_low_95", "conf_high_95",
  x_limits = c(0.68, 2.85), x_breaks = c(0.7, 1, 1.25, 1.5, 2),
  label_x = 2.72, label_column_start = 1.78,
  title = "History–deviation model",
  subtitle = "M2: four components entered simultaneously"
)

joint_labels <- dat$joint |>
  dplyr::mutate(
    cohort = factor(.data$cohort, levels = c("ELSA", "SHARE")),
    label = paste0("Joint Wald: ", format_p(.data$p_value)),
    x = 2.05,
    display_term = factor("Affective", levels = c("Functional", "Affective"))
  )

if (max(panel_d_data$conf_high_95, na.rm = TRUE) >= 1.68) {
  stop("The Panel D estimate-label column overlaps a confidence interval.")
}

p_d <- ggplot2::ggplot(
  panel_d_data,
  ggplot2::aes(x = .data$ratio_of_hazard_ratios, y = .data$display_term, colour = .data$domain)
) +
  ggplot2::annotate(
    "rect", xmin = 1.68, xmax = Inf,
    ymin = -Inf, ymax = Inf, fill = "#F8F9FA", colour = NA
  ) +
  ggplot2::geom_vline(xintercept = 1, linewidth = 0.4, linetype = 2, colour = palette["neutral_mid"]) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = .data$conf_low_95, xmax = .data$conf_high_95),
    height = 0.16, linewidth = 0.55
  ) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_text(
    ggplot2::aes(x = 2.05, label = .data$estimate_label),
    hjust = 1, family = "Arial", size = 2.2, colour = palette["neutral_dark"],
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = joint_labels,
    ggplot2::aes(x = .data$x, y = .data$display_term, label = .data$label),
    inherit.aes = FALSE, hjust = 1, nudge_y = 0.34,
    family = "Arial", fontface = "bold",
    size = 2.25, colour = palette["neutral_dark"]
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(.data$cohort), scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = palette[c("Affective", "Functional")]) +
  ggplot2::scale_x_log10(limits = c(0.88, 2.15), breaks = c(0.9, 1, 1.2, 1.5, 2)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.42, 0.72))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "History versus deviation",
    subtitle = "Ratio of HRs; values >1 favour stronger history association",
    x = "Ratio of hazard ratios (95% CI)", y = NULL
  ) +
  theme_pub() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 6.5),
    panel.spacing.y = grid::unit(7, "pt"),
    strip.text.y.right = ggplot2::element_text(
      angle = 270, face = "bold", margin = ggplot2::margin(2, 2, 2, 2)
    ),
    plot.margin = ggplot2::margin(8, 5, 4, 5, unit = "pt")
  )

layout_design <- "
AAAA
BBCC
DDCC
"
figure_3 <- p_a + p_b + p_c + p_d +
  patchwork::plot_layout(design = layout_design, heights = c(0.76, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 9),
    plot.tag.position = c(0.005, 0.995)
  )

save_pub <- function(plot, stem, width_mm = 183, height_mm = 166, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw", background = "white"
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    paste0(stem, "_preview.png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

save_pub(
  figure_3,
  file.path(out_dir, "Figure3_Dementia_Related_Clinical_Events")
)

## Supplementary sensitivity forest plot.
analysis_order <- c(
  "Exclude first follow-up events", "Cognition adjusted", "Direct interviews only",
  "Exclude same-date records", "Same or previous-wave negative",
  "No later negative reversal", "Strict observation censoring",
  "Strict censoring + previous-wave negative"
)
sup_sensitivity <- sup_sensitivity |>
  dplyr::mutate(
    analysis_label = factor(.data$analysis_label, levels = rev(analysis_order))
  )

supp_figure <- ggplot2::ggplot(
  sup_sensitivity,
  ggplot2::aes(
    x = .data$hazard_ratio, y = .data$analysis_label,
    colour = .data$domain, shape = .data$component
  )
) +
  ggplot2::geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2, colour = palette["neutral_mid"]) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = .data$conf_low_95, xmax = .data$conf_high_95),
    height = 0.14, linewidth = 0.45,
    position = ggplot2::position_dodge(width = 0.55)
  ) +
  ggplot2::geom_point(size = 1.8, position = ggplot2::position_dodge(width = 0.55)) +
  ggplot2::facet_grid(rows = ggplot2::vars(.data$cohort), cols = ggplot2::vars(.data$term_label)) +
  ggplot2::scale_colour_manual(values = palette[c("Affective", "Functional")]) +
  ggplot2::scale_shape_manual(values = c(History = 16, `Current deviation` = 17)) +
  ggplot2::scale_x_log10(limits = c(0.62, 1.90), breaks = c(0.7, 1, 1.25, 1.5)) +
  ggplot2::labs(x = "Hazard ratio (95% CI)", y = NULL) +
  theme_pub(base_size = 7) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 6.2),
    strip.text = ggplot2::element_text(size = 6.3, face = "bold"),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.box = "horizontal"
  )

save_pub(
  supp_figure,
  file.path(out_dir, "Supplementary_Figure_Dementia_Sensitivity_Analyses"),
  width_mm = 183, height_mm = 175, dpi = 600
)

## Figure source data.
panel_a_source <- tibble::tribble(
  ~element, ~definition,
  "Historical mean (H)", "Mean of all eligible measurements preceding the current visit",
  "Current deviation (D)", "Current score minus the preceding historical mean",
  "Primary observation rule", "Known death retained as a competing event when supported by a dynamic interval",
  "Strict observation sensitivity", "Follow-up censored at the last observed dementia-related assessment"
)
panel_d_source <- panel_d_data |>
  dplyr::select(
    .data$cohort, .data$domain, .data$beta_difference, .data$robust_se,
    .data$ratio_of_hazard_ratios, .data$conf_low_95, .data$conf_high_95,
    .data$p_value, .data$p_holm_contrast_family, .data$joint_p
  )

readr::write_csv(panel_a_source, file.path(source_dir, "Figure3A_timeline_definitions.csv"), na = "")
readr::write_csv(panel_b_data, file.path(source_dir, "Figure3B_current_score_model.csv"), na = "")
readr::write_csv(panel_c_data, file.path(source_dir, "Figure3C_history_deviation_model.csv"), na = "")
readr::write_csv(panel_d_source, file.path(source_dir, "Figure3D_history_deviation_contrasts.csv"), na = "")
readr::write_csv(sup_sensitivity, file.path(source_dir, "Supplementary_Figure_sensitivity_coefficients.csv"), na = "")

wb_source <- openxlsx::createWorkbook()
for (sheet in c("Figure3A", "Figure3B", "Figure3C", "Figure3D", "SupplementaryFigure")) {
  openxlsx::addWorksheet(wb_source, sheet)
}
openxlsx::writeData(wb_source, "Figure3A", panel_a_source)
openxlsx::writeData(wb_source, "Figure3B", panel_b_data)
openxlsx::writeData(wb_source, "Figure3C", panel_c_data)
openxlsx::writeData(wb_source, "Figure3D", panel_d_source)
openxlsx::writeData(wb_source, "SupplementaryFigure", sup_sensitivity)
openxlsx::freezePane(wb_source, "Figure3B", firstRow = TRUE)
openxlsx::freezePane(wb_source, "Figure3C", firstRow = TRUE)
openxlsx::freezePane(wb_source, "Figure3D", firstRow = TRUE)
openxlsx::freezePane(wb_source, "SupplementaryFigure", firstRow = TRUE)
openxlsx::saveWorkbook(
  wb_source,
  file.path(out_dir, "Figure3_Source_Data.xlsx"),
  overwrite = TRUE
)

## Submission tables assembled from the model outputs.
wb_sup <- openxlsx::createWorkbook()
submission_sheets <- c("Sample_flow", "Models", "Contrasts", "Sensitivity", "Ascertainment")
for (sheet in submission_sheets) openxlsx::addWorksheet(wb_sup, sheet)

header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1")
write_block <- function(wb, sheet, title, data, start_row) {
  openxlsx::writeData(wb, sheet, title, startRow = start_row, startCol = 1)
  openxlsx::addStyle(wb, sheet, header_style, rows = start_row, cols = 1, gridExpand = TRUE)
  openxlsx::writeData(wb, sheet, data, startRow = start_row + 1, startCol = 1, withFilter = FALSE)
  start_row + nrow(data) + 4L
}

r <- 1L
r <- write_block(wb_sup, "Sample_flow", "Main analysis flow", dat$flow, r)
r <- write_block(wb_sup, "Sample_flow", "Observation-rule sample support", dat$support, r)
r <- write_block(wb_sup, "Sample_flow", "Observation censoring audit", dat$censoring, r)

write_block(wb_sup, "Models", "M0, M1 and M2 cause-specific Cox models", dat$main_coef, 1L)

r <- 1L
r <- write_block(wb_sup, "Contrasts", "History-deviation contrasts", dat$contrasts, r)
r <- write_block(wb_sup, "Contrasts", "Joint robust Wald tests", dat$joint, r)
r <- write_block(wb_sup, "Contrasts", "M1 versus M2 model comparisons", dat$comparison, r)

r <- 1L
r <- write_block(wb_sup, "Sensitivity", "Sensitivity coefficients", dat$sensitivity_coef, r)
r <- write_block(wb_sup, "Sensitivity", "Sensitivity history-deviation contrasts", dat$sensitivity_contrast, r)
r <- write_block(wb_sup, "Sensitivity", "Sensitivity joint tests", dat$sensitivity_joint, r)
r <- write_block(wb_sup, "Sensitivity", "Complementary log-log sensitivity", dat$cloglog, r)

r <- 1L
r <- write_block(wb_sup, "Ascertainment", "Endpoint quality audit", dat$endpoint_quality, r)
r <- write_block(wb_sup, "Ascertainment", "Deaths after last observed assessment", dat$death_audit, r)
r <- write_block(wb_sup, "Ascertainment", "Observation-rule coefficient comparison", dat$observation_coef, r)
r <- write_block(wb_sup, "Ascertainment", "Observation-rule joint-test comparison", dat$observation_joint, r)

for (sheet in submission_sheets) {
  openxlsx::setColWidths(wb_sup, sheet, cols = 1:40, widths = "auto")
  openxlsx::freezePane(wb_sup, sheet, firstRow = TRUE)
}
openxlsx::saveWorkbook(
  wb_sup,
  file.path(out_dir, "Dementia_Related_Outcomes_Submission_Tables.xlsx"),
  overwrite = TRUE
)

## Caption and QA record.
caption <- c(
  "Figure 3 | Historical burden, current deviation, and subsequent dementia-related outcomes.",
  "A, Historical burden was calculated from eligible measurements preceding the current visit, and current deviation represented the current score relative to that history. The primary observation rule retained a known death as a competing event when supported by a dynamic interval; strict censoring at the last observed dementia-related assessment was evaluated separately.",
  "B, Cause-specific hazard ratios and 95% confidence intervals for current affective and functional scores (M1).",
  "C, Mutually adjusted cause-specific hazard ratios and 95% confidence intervals for affective history, affective deviation, functional history, and functional deviation (M2).",
  "D, Ratios of the history and deviation hazard ratios within each domain. Ratios above one indicate a stronger association for historical burden. Labels report the two-degree-of-freedom robust Wald test of the affective and functional equal-coefficient constraints.",
  "Models used attained age as the time scale, participant-clustered robust standard errors, and adjustment for sex, education, and calendar year; SHARE models used country-stratified baseline hazards. ELSA and SHARE were analysed separately because their dementia-related outcome definitions differed. Hazard ratios correspond to a one-standard-deviation higher score or component."
)
writeLines(caption, file.path(out_dir, "Figure3_Legend.txt"), useBytes = TRUE)

qa <- c(
  "Core conclusion: Historical burden and current deviation contain prognostic information not fully represented by the corresponding current scores.",
  "Archetype: schematic-led quantitative composite.",
  "Backend: R only (ggplot2 and patchwork).",
  "Final size: 183 mm x 166 mm.",
  "Exports: editable SVG, vector PDF, 600-dpi TIFF, 300-dpi PNG preview.",
  "Panel labels: uppercase A-D to match manuscript citations.",
  "Statistics: hazard ratios with participant-clustered robust 95% confidence intervals; Holm-adjusted exposure-family results retained in source data.",
  "Review boundary: cause-specific prognostic associations, not absolute risk prediction, clinical utility, diagnosis, or causality.",
  "Required visual check after local rendering: no overlap, readable text at final size, selectable PDF/SVG text, and intact en dashes and minus signs."
)
writeLines(qa, file.path(out_dir, "Figure3_QA_Record.txt"), useBytes = TRUE)

manifest_paths <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
manifest_paths <- manifest_paths[file.info(manifest_paths)$isdir %in% FALSE]
manifest <- tibble::tibble(
  file = normalizePath(manifest_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(manifest_paths)$size,
  md5 = unname(tools::md5sum(manifest_paths))
)
readr::write_csv(manifest, file.path(out_dir, "Figure3_Output_Manifest.csv"), na = "")

message("Figure 3 and supporting files generated in: ", out_dir)
