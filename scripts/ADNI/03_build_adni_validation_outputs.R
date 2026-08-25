#!/usr/bin/env Rscript
# Confirmatory diagnosis-by-time sensitivity model, Figure 6 assets, and submission tables.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

required <- c("dplyr", "readr", "tibble", "lme4", "lmerTest", "broom.mixed", "ggplot2", "openxlsx", "svglite", "ragg")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))

project_dir <- normalizePath(Sys.getenv("DCV_PROJECT_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
output_root <- Sys.getenv("ADNI_OUTPUT_ROOT", unset = file.path(project_dir, "results", "adni_validation"))
input_file <- file.path(output_root, "dataset", "03_analysis_dataset_locked.csv")
out_dir <- file.path(output_root, "validation_outputs")
figure_dir <- file.path(out_dir, "figure")
source_dir <- file.path(figure_dir, "Source_Data_CSV")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
runtime_tmp <- tempdir()
dir.create(runtime_tmp, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(runtime_tmp)) stop("Unable to recreate R temporary directory for figure and Excel output.")
if (!file.exists(input_file)) stop("Run 01_build_locked_cdrsb_dataset.R first.")

d <- readr::read_csv(input_file, show_col_types = FALSE) |>
  dplyr::filter(!is.na(.data$baseline_APOE4)) |>
  dplyr::mutate(PTID = as.character(.data$PTID), sex = factor(.data$sex), baseline_diagnosis = factor(.data$baseline_diagnosis)) |>
  dplyr::filter(dplyr::if_all(c("PGS_affective_z", "PGS_functional_z", "time_years",
                                "baseline_age_c", "education_c", dplyr::starts_with("PC")), is.finite))

fit <- lmerTest::lmer(
  CDRSB ~ PGS_affective_z * time_years + PGS_functional_z * time_years + baseline_diagnosis * time_years +
    baseline_age_c + sex + education_c + baseline_APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    (1 | PTID),
  data = d, REML = FALSE,
  control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)
fixed <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) |>
  dplyr::mutate(n_rows = nrow(d), n_people = dplyr::n_distinct(d$PTID),
                n_people_2plus_visits = sum(table(d$PTID) >= 2L), singular = lme4::isSingular(fit, tol = 1e-4))
focal <- fixed |>
  dplyr::filter(.data$term %in% c("PGS_affective_z:time_years", "time_years:PGS_functional_z")) |>
  dplyr::mutate(p_fdr_two_primary_terms = stats::p.adjust(.data$p.value, method = "BH"))

participant_followup <- d |>
  dplyr::group_by(.data$PTID, .data$baseline_diagnosis) |>
  dplyr::summarise(n_observations = dplyr::n(), followup_years = max(.data$time_years),
                   baseline_CDRSB = .data$CDRSB[which.min(.data$time_years)][1], .groups = "drop")
sample_summary <- participant_followup |>
  dplyr::group_by(.data$baseline_diagnosis) |>
  dplyr::summarise(participants = dplyr::n(), observations = sum(.data$n_observations),
                   participants_2plus_visits = sum(.data$n_observations >= 2L),
                   median_followup_years = median(.data$followup_years),
                   q1_followup_years = stats::quantile(.data$followup_years, 0.25, names = FALSE),
                   q3_followup_years = stats::quantile(.data$followup_years, 0.75, names = FALSE),
                   baseline_CDRSB_median = median(.data$baseline_CDRSB), .groups = "drop")
diagnosis_slopes <- fixed |>
  dplyr::filter(grepl("^time_years:baseline_diagnosis", .data$term)) |>
  dplyr::mutate(diagnosis = sub("^time_years:baseline_diagnosis", "", .data$term),
                reference_diagnosis = levels(d$baseline_diagnosis)[1])

readr::write_csv(fixed, file.path(out_dir, "01_diagnosis_time_model_coefficients.csv"), na = "")
readr::write_csv(focal, file.path(out_dir, "02_primary_interaction_summary.csv"), na = "")
readr::write_csv(sample_summary, file.path(out_dir, "03_diagnosis_sample_summary.csv"), na = "")
readr::write_csv(diagnosis_slopes, file.path(out_dir, "04_diagnosis_time_interactions.csv"), na = "")

forest_data <- focal |>
  dplyr::transmute(
    score = dplyr::case_when(
      .data$term == "PGS_affective_z:time_years" ~ "Affective-factor PGS",
      .data$term == "time_years:PGS_functional_z" ~ "Functional-factor PGS",
      TRUE ~ .data$term
    ),
    estimate = .data$estimate, conf.low = .data$conf.low, conf.high = .data$conf.high,
    p.value = .data$p.value, p_fdr_two_primary_terms = .data$p_fdr_two_primary_terms
  ) |>
  dplyr::arrange(.data$estimate) |>
  dplyr::mutate(score = factor(.data$score, levels = .data$score))

theme_pub <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#30343B"),
                   axis.text = ggplot2::element_text(colour = "#30343B"),
                   axis.title = ggplot2::element_text(colour = "#30343B"),
                   plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
                   plot.margin = ggplot2::margin(4, 5, 4, 5, unit = "pt"))
}
figure_6 <- ggplot2::ggplot(forest_data, ggplot2::aes(x = .data$estimate, y = .data$score)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "#747B83", linewidth = 0.4) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high), height = 0.16, linewidth = 0.55, colour = "#287F8E") +
  ggplot2::geom_point(size = 2.3, colour = "#287F8E") +
  ggplot2::scale_x_continuous(limits = c(-0.060, 0.060), breaks = c(-0.05, -0.025, 0, 0.025, 0.05)) +
  ggplot2::labs(title = "Independent ADNI clinical validation", subtitle = "Joint model with baseline-diagnosis-specific CDR-SB time slopes",
                x = "Difference in annual CDR-SB change per 1-SD higher PGS (95% CI)", y = NULL) +
  theme_pub()

save_pub <- function(plot, stem, width_mm = 183, height_mm = 70, dpi = 600) {
  width_in <- width_mm / 25.4; height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in); print(plot); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial"); print(plot); grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw", background = "white"); print(plot); grDevices::dev.off()
  ragg::agg_png(paste0(stem, "_preview.png"), width = width_in, height = height_in, units = "in", res = 300, background = "white"); print(plot); grDevices::dev.off()
}
save_pub(figure_6, file.path(figure_dir, "Figure6_ADNI_Factor_PGS_Validation"))

diagnostic_data <- tibble::tibble(fitted = stats::fitted(fit), residual = stats::resid(fit))
residual_plot <- ggplot2::ggplot(diagnostic_data, ggplot2::aes(.data$fitted, .data$residual)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "#747B83") + ggplot2::geom_point(alpha = 0.25, size = 0.7, colour = "#287F8E") +
  ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "#D97745", linewidth = 0.6) +
  ggplot2::labs(x = "Fitted CDR-SB", y = "Conditional residual") + theme_pub()
ggplot2::ggsave(file.path(figure_dir, "Supplementary_Figure_ADNI_Residuals.png"), residual_plot, width = 6.5, height = 4.5, dpi = 600)

readr::write_csv(forest_data, file.path(source_dir, "Figure6_factor_PGS_time_effects.csv"), na = "")
readr::write_csv(diagnostic_data, file.path(source_dir, "Supplementary_Figure_residual_diagnostics.csv"), na = "")

wb <- openxlsx::createWorkbook()
for (sheet in c("Table_S54_sample", "Table_S55_model", "Table_S55_diagnosis_slopes")) openxlsx::addWorksheet(wb, sheet)
openxlsx::writeData(wb, "Table_S54_sample", sample_summary)
openxlsx::writeData(wb, "Table_S55_model", fixed)
openxlsx::writeData(wb, "Table_S55_diagnosis_slopes", diagnosis_slopes)
header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1")
for (sheet in names(wb)) {
  openxlsx::addStyle(wb, sheet, header_style, rows = 1, cols = 1:40, gridExpand = TRUE)
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet, cols = 1:40, widths = "auto")
}
openxlsx::saveWorkbook(wb, file.path(out_dir, "ADNI_CDRSB_Submission_Tables.xlsx"), overwrite = TRUE)

legend <- c(
  "Figure 6 | Independent ADNI clinical validation of affective- and functional-factor polygenic scores.",
  "Points show the difference in annual CDR-SB change per one-standard-deviation higher factor-PGS; horizontal lines show 95% confidence intervals.",
  "The linear mixed model estimated both factor-PGS-by-time terms jointly and adjusted for baseline age, sex, education, APOE4, genetic principal components, baseline diagnosis and diagnosis-specific time slopes, with a participant-specific random intercept.",
  "Higher CDR-SB indicates greater clinical severity. The model estimates longitudinal associations and does not establish causal effects or clinical prediction utility."
)
writeLines(legend, file.path(figure_dir, "Figure6_Legend.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
message("ADNI sensitivity model, Figure 6 assets, and submission tables written to: ", out_dir)
