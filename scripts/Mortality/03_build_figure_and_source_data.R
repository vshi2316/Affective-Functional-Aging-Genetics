#!/usr/bin/env Rscript
# Build the prespecified mortality-bridge figure and its tabular source data.
# This script is read-only with respect to analysis results.

options(stringsAsFactors = FALSE, warn = 1)

project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
output_root <- Sys.getenv(
  "MORTALITY_OUTPUT_ROOT",
  unset = file.path(project_dir, "results", "mortality_bridge")
)
model_dir <- Sys.getenv(
  "MORTALITY_MODEL_OUT",
  unset = file.path(output_root, "models")
)
figure_dir <- Sys.getenv(
  "MORTALITY_FIGURE_OUT",
  unset = file.path(output_root, "figures")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  estimates = file.path(model_dir, "all_cohort_estimates.csv"),
  meta = file.path(model_dir, "main_age_band_meta.csv"),
  contrasts = file.path(model_dir, "main_history_deviation_contrast_meta.csv"),
  global_tests = file.path(model_dir, "four_global_Wald_tests.csv"),
  status = file.path(model_dir, "model_status.csv"),
  gates = file.path(model_dir, "analysis_gates.csv"),
  lock = file.path(model_dir, "mortality_models_lock.txt")
)
if (any(!file.exists(required))) {
  stop("Missing frozen model output(s): ",
       paste(names(required)[!file.exists(required)], collapse = ", "))
}

pkgs <- c("dplyr", "tidyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install before running: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

gates <- read_csv(required[["gates"]], show_col_types = FALSE)
if (!all(gates$pass %in% TRUE)) stop("Frozen model gates are not all passed.")

est <- read_csv(required[["estimates"]], show_col_types = FALSE)
meta <- read_csv(required[["meta"]], show_col_types = FALSE)
contrast <- read_csv(required[["contrasts"]], show_col_types = FALSE)
tests <- read_csv(required[["global_tests"]], show_col_types = FALSE)

age_labels <- c(
  lt65 = "<65",
  `65_74` = "65–74",
  `75_84` = "75–84",
  ge85 = "≥85"
)
component_labels <- c(
  H_A = "Affective history",
  D_A = "Current affective deviation",
  H_F = "Functional history",
  D_F = "Current functional deviation"
)
cohort_order <- c("HRS", "MHAS", "SHARE", "ELSA")
palette <- c(
  H_A = "#0072B2",
  D_A = "#56B4E9",
  H_F = "#D55E00",
  D_F = "#E69F00"
)

cohort_plot_data <- est %>%
  filter(
    cohort %in% cohort_order,
    variant == "primary",
    !interaction,
    grepl("^(H_A|D_A|H_F|D_F)__", term)
  ) %>%
  separate(term, into = c("component", "age_band"), sep = "__", remove = FALSE) %>%
  mutate(
    cohort = factor(cohort, levels = rev(cohort_order)),
    component = factor(component, levels = names(component_labels)),
    age_band = factor(age_band, levels = names(age_labels)),
    component_label = factor(component_labels[as.character(component)],
                             levels = unname(component_labels)),
    age_label = factor(age_labels[as.character(age_band)], levels = unname(age_labels))
  )

meta_plot_data <- meta %>%
  mutate(
    component = factor(component, levels = names(component_labels)),
    age_band = factor(age_band, levels = names(age_labels)),
    component_label = factor(component_labels[as.character(component)],
                             levels = unname(component_labels)),
    age_label = factor(age_labels[as.character(age_band)], levels = unname(age_labels))
  )

contrast_plot_data <- contrast %>%
  mutate(
    domain = recode(domain, A = "Affective", F = "Functional"),
    age_band = factor(age_band, levels = names(age_labels)),
    age_label = factor(age_labels[as.character(age_band)], levels = unname(age_labels))
  )

p_a <- ggplot(
  cohort_plot_data,
  aes(x = hr, y = cohort, xmin = ci_low, xmax = ci_high, colour = component)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2, colour = "grey45") +
  geom_errorbarh(height = 0, linewidth = 0.45, position = position_dodge(width = 0.55)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.55)) +
  facet_grid(component_label ~ age_label, scales = "free_y") +
  scale_x_log10() +
  scale_colour_manual(values = palette, guide = "none") +
  labs(x = "Hazard ratio per baseline SD", y = NULL, title = "a  Cohort-specific associations") +
  theme_bw(base_size = 8.5) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(size = 7.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 9.2)
  )

p_b <- ggplot(
  meta_plot_data,
  aes(x = hr_or_ratio, y = component_label, xmin = ci_low, xmax = ci_high, colour = component)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2, colour = "grey45") +
  geom_errorbarh(height = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~age_label, nrow = 1) +
  scale_x_log10() +
  scale_colour_manual(values = palette, guide = "none") +
  labs(x = "Random-effects hazard ratio", y = NULL,
       title = "b  Main-cohort meta-analysis") +
  theme_bw(base_size = 8.5) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 9.2)
  )

p_c <- ggplot(
  contrast_plot_data,
  aes(x = hr_or_ratio, y = domain, xmin = ci_low, xmax = ci_high, colour = domain)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2, colour = "grey45") +
  geom_errorbarh(height = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~age_label, nrow = 1) +
  scale_x_log10() +
  scale_colour_manual(values = c(Affective = "#0072B2", Functional = "#D55E00"),
                      guide = "none") +
  labs(x = "Ratio of hazard ratios: history / current deviation", y = NULL,
       title = "c  Formal history–deviation contrasts") +
  theme_bw(base_size = 8.5) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 9.2)
  )

figure <- p_a / p_b / p_c + plot_layout(heights = c(2.1, 1, 1))

pdf_path <- file.path(figure_dir, "Figure_mortality_bridge.pdf")
svg_path <- file.path(figure_dir, "Figure_mortality_bridge.svg")
tiff_path <- file.path(figure_dir, "Figure_mortality_bridge.tiff")
ggsave(pdf_path, figure, width = 180, height = 250, units = "mm", device = cairo_pdf)
ggsave(svg_path, figure, width = 180, height = 250, units = "mm", device = "svg")
ggsave(
  tiff_path, figure, width = 180, height = 250, units = "mm",
  dpi = 600, compression = "lzw", device = "tiff"
)

source_data <- bind_rows(
  cohort_plot_data %>% transmute(
    panel = "a", cohort = as.character(cohort), age_band = as.character(age_band),
    component = as.character(component), estimate = hr, ci_low, ci_high,
    p_value = p, tau2 = NA_real_, i2 = NA_real_, prediction_low = NA_real_,
    prediction_high = NA_real_
  ),
  meta_plot_data %>% transmute(
    panel = "b", cohort = "HRS_MHAS_SHARE_meta", age_band = as.character(age_band),
    component = as.character(component), estimate = hr_or_ratio, ci_low, ci_high,
    p_value = p_hk, tau2, i2, prediction_low, prediction_high
  ),
  contrast_plot_data %>% transmute(
    panel = "c", cohort = "HRS_MHAS_SHARE_meta", age_band = as.character(age_band),
    component = paste0(domain, "_history_minus_deviation"),
    estimate = hr_or_ratio, ci_low, ci_high, p_value = p_hk, tau2, i2,
    prediction_low, prediction_high
  )
)
write_csv(
  source_data,
  file.path(figure_dir, "Source_Data_Figure_mortality_bridge.csv"),
  na = ""
)
write_csv(
  tests,
  file.path(figure_dir, "Source_Data_primary_global_tests.csv"),
  na = ""
)

manifest <- tibble(
  file = basename(c(pdf_path, svg_path, tiff_path,
                    file.path(figure_dir, "Source_Data_Figure_mortality_bridge.csv"),
                    file.path(figure_dir, "Source_Data_primary_global_tests.csv"))),
  md5 = unname(tools::md5sum(file.path(figure_dir, file)))
)
write_csv(manifest, file.path(figure_dir, "figure_manifest_md5.csv"), na = "")
writeLines(
  c(
    paste0("Completed: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
    paste0("Frozen-model lock MD5: ", unname(tools::md5sum(required[["lock"]]))),
    "No model was refitted and no result was selected for display."
  ),
  file.path(figure_dir, "figure_build_lock.txt"),
  useBytes = TRUE
)
