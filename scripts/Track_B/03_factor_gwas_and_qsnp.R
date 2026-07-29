options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(data.table)
  library(GenomicSEM)
})

resolve_base_dir <- function() {
  base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(base_dir) || !dir.exists(base_dir)) {
    stop("Set DCV_BASE_DIR to the project directory containing analysis_ready_core.")
  }
  normalizePath(base_dir, winslash = "/", mustWork = TRUE)
}

find_output_column <- function(column_names, exact_candidates, regex_candidates = character()) {
  exact_hit <- exact_candidates[exact_candidates %in% column_names]
  if (length(exact_hit) > 0) return(exact_hit[[1]])
  for (pattern in regex_candidates) {
    regex_hit <- grep(pattern, column_names, value = TRUE, ignore.case = TRUE)
    if (length(regex_hit) > 0) return(regex_hit[[1]])
  }
  NA_character_
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

make_qq_data <- function(p_values, maximum_points = 100000L) {
  p_values <- suppressWarnings(as.numeric(p_values))
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values <= 1]
  p_values <- sort(p_values)
  n <- length(p_values)
  if (n < 100) return(tibble())
  index <- unique(round(seq(1, n, length.out = min(n, maximum_points))))
  tibble(
    expected = -log10(ppoints(n)[index]),
    observed = -log10(p_values[index])
  )
}

make_qq_plot <- function(p_values, title, output_file) {
  plot_data <- make_qq_data(p_values)
  if (nrow(plot_data) == 0) return(invisible(NULL))
  plot <- ggplot(plot_data, aes(expected, observed)) +
    geom_point(size = 0.55, alpha = 0.60, colour = "#2C6E9E") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#B23A48") +
    labs(title = title, x = "Expected -log10(P)", y = "Observed -log10(P)") +
    theme_classic(base_size = 11)
  ggsave(output_file, plot, width = 5.8, height = 5.4, dpi = 600)
}

make_manhattan_plot <- function(results, p_column, title, output_file) {
  plot_data <- results[
    is.finite(CHR) & is.finite(BP) & is.finite(get(p_column)) &
      CHR >= 1 & CHR <= 22 & BP > 0 & get(p_column) > 0 & get(p_column) <= 1,
    .(CHR, BP, SNP, p_value = get(p_column))
  ]
  if (nrow(plot_data) == 0) return(invisible(NULL))

  plot_data[, neglog10_p := -log10(p_value)]
  significant <- plot_data[p_value < 5e-8]
  nonsignificant <- plot_data[p_value >= 5e-8]
  set.seed(20260713)
  if (nrow(nonsignificant) > 400000) {
    nonsignificant <- nonsignificant[sample.int(nrow(nonsignificant), 400000)]
  }
  plot_data <- rbindlist(list(significant, nonsignificant), use.names = TRUE)
  setorder(plot_data, CHR, BP)

  chromosome_lengths <- plot_data[, .(max_bp = max(BP)), by = CHR][order(CHR)]
  chromosome_lengths[, offset := shift(cumsum(max_bp), fill = 0)]
  chromosome_lengths[, centre := offset + max_bp / 2]
  plot_data <- merge(plot_data, chromosome_lengths[, .(CHR, offset)], by = "CHR")
  plot_data[, cumulative_bp := BP + offset]

  plot <- ggplot(plot_data, aes(cumulative_bp, neglog10_p, colour = factor(CHR %% 2))) +
    geom_point(size = 0.40, alpha = 0.65) +
    geom_hline(yintercept = -log10(5e-8), linetype = 2, colour = "#B23A48") +
    scale_x_continuous(breaks = chromosome_lengths$centre, labels = chromosome_lengths$CHR) +
    scale_colour_manual(values = c("#2C6E9E", "#6AAED6")) +
    labs(title = title, x = "Chromosome", y = "-log10(P)") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", panel.grid = element_blank())
  ggsave(output_file, plot, width = 13.5, height = 5.2, dpi = 600)
}

base_dir <- resolve_base_dir()
trackb_dir <- file.path(base_dir, "analysis_ready_core", "trackB")
source_dir <- file.path(trackb_dir, "factor_gwas_qsnp")
comparison_dir <- file.path(trackb_dir, "genomicsem_competing_models")
output_dir <- file.path(trackb_dir, "twofactor_gwas_qsnp")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

covstruc_file <- file.path(source_dir, "trackB_afffunc_ldsc_covstruc.rds")
snp_file <- file.path(source_dir, "trackB_afffunc_sumstats.rds")
fit_file <- file.path(comparison_dir, "trackB_competing_model_fit.csv")

required_files <- c(covstruc_file, snp_file, fit_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = "; "))
}

fit_audit <- read_csv(fit_file, show_col_types = FALSE)
two_factor_audit <- fit_audit %>%
  filter(model_name == "correlated_affective_functional")

if (
  nrow(two_factor_audit) != 1 ||
    two_factor_audit$status[[1]] != "estimated" ||
    !isTRUE(two_factor_audit$proper_solution[[1]])
) {
  stop("The correlated two-factor model did not pass the completed model audit.")
}

covstruc <- readRDS(covstruc_file)
snp_sumstats <- readRDS(snp_file)

required_snp_columns <- c(
  "CHR", "BP", "SNP", "A1", "A2", "MAF",
  "beta.mdd", "se.mdd", "beta.lonely", "se.lonely",
  "beta.frailty", "se.frailty", "beta.srh", "se.srh",
  "beta.walkpace", "se.walkpace"
)
missing_columns <- setdiff(required_snp_columns, names(snp_sumstats))
if (length(missing_columns) > 0) {
  stop("The SNP summary object is missing columns: ", paste(missing_columns, collapse = ", "))
}

two_factor_snp_model <- "
  affective_domain =~ mdd + lonely
  functional_domain =~ frailty + srh + walkpace
  affective_domain ~~ functional_domain
  affective_domain ~ SNP
  functional_domain ~ SNP
"

message("Running the correlated two-factor userGWAS on ", nrow(snp_sumstats), " variants.")
message("The common projection GWAS outputs are treated as read-only inputs.")

two_factor_raw <- GenomicSEM::userGWAS(
  covstruc = covstruc,
  SNPs = snp_sumstats,
  estimation = "DWLS",
  model = two_factor_snp_model,
  printwarn = TRUE,
  sub = FALSE,
  parallel = FALSE,
  std.lv = TRUE,
  fix_measurement = TRUE,
  Q_SNP = TRUE,
  analytic = TRUE,
  batch_size = 50000
)

results <- as.data.table(two_factor_raw)
writeLines(names(results), file.path(output_dir, "trackB_twofactor_output_columns.txt"))

column_map <- list(
  beta_affective = find_output_column(
    names(results),
    c("beta_affective_domain"),
    c("^beta.*affective_domain$")
  ),
  se_affective = find_output_column(
    names(results),
    c("SE_affective_domain", "se_affective_domain"),
    c("^se.*affective_domain$")
  ),
  z_affective = find_output_column(
    names(results),
    c("Z_beta_affective_domain", "z_beta_affective_domain"),
    c("^z.*affective_domain$")
  ),
  p_affective = find_output_column(
    names(results),
    c("p_val_affective_domain", "p_affective_domain"),
    c("^p.*affective_domain$")
  ),
  beta_functional = find_output_column(
    names(results),
    c("beta_functional_domain"),
    c("^beta.*functional_domain$")
  ),
  se_functional = find_output_column(
    names(results),
    c("SE_functional_domain", "se_functional_domain"),
    c("^se.*functional_domain$")
  ),
  z_functional = find_output_column(
    names(results),
    c("Z_beta_functional_domain", "z_beta_functional_domain"),
    c("^z.*functional_domain$")
  ),
  p_functional = find_output_column(
    names(results),
    c("p_val_functional_domain", "p_functional_domain"),
    c("^p.*functional_domain$")
  ),
  q_snp = find_output_column(
    names(results),
    c("Q_omnibus"),
    c("^Q_omnibus$")
  ),
  q_snp_df = find_output_column(
    names(results),
    c("Q_omnibus_df"),
    c("^Q_omnibus_df$")
  ),
  q_snp_p = find_output_column(
    names(results),
    c("Q_omnibus_pval"),
    c("^Q_omnibus_p")
  )
)

required_mapped <- c(
  "beta_affective", "se_affective", "p_affective",
  "beta_functional", "se_functional", "p_functional",
  "q_snp", "q_snp_p"
)
missing_mapped <- required_mapped[is.na(unlist(column_map[required_mapped]))]
if (length(missing_mapped) > 0) {
  stop(
    "GenomicSEM output columns could not be mapped for: ",
    paste(missing_mapped, collapse = ", "),
    ". Inspect trackB_twofactor_output_columns.txt."
  )
}

for (target in names(column_map)) {
  source <- column_map[[target]]
  if (!is.na(source) && source %in% names(results) && source != target) {
    setnames(results, source, target)
  }
}

numeric_columns <- intersect(
  c(
    "CHR", "BP", "MAF", "beta_affective", "se_affective", "z_affective",
    "p_affective", "beta_functional", "se_functional", "z_functional",
    "p_functional", "q_snp", "q_snp_df", "q_snp_p"
  ),
  names(results)
)
results[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]

if (!"z_affective" %in% names(results)) {
  results[, z_affective := beta_affective / se_affective]
}
if (!"z_functional" %in% names(results)) {
  results[, z_functional := beta_functional / se_functional]
}

results[, affective_gws := is.finite(p_affective) & p_affective < 5e-8]
results[, functional_gws := is.finite(p_functional) & p_functional < 5e-8]
results[, q_snp_gws := is.finite(q_snp_p) & q_snp_p < 5e-8]
results[, gws_class := fifelse(
  affective_gws & functional_gws, "both_factors",
  fifelse(affective_gws, "affective_only", fifelse(functional_gws, "functional_only", "not_factor_gws"))
)]

summary_table <- tibble(
  n_variants = nrow(results),
  n_affective_gws = sum(results$affective_gws, na.rm = TRUE),
  n_functional_gws = sum(results$functional_gws, na.rm = TRUE),
  n_both_factor_gws = sum(results$affective_gws & results$functional_gws, na.rm = TRUE),
  n_affective_only_gws = sum(results$affective_gws & !results$functional_gws, na.rm = TRUE),
  n_functional_only_gws = sum(!results$affective_gws & results$functional_gws, na.rm = TRUE),
  n_q_snp_gws = sum(results$q_snp_gws, na.rm = TRUE),
  min_p_affective = safe_min(results$p_affective),
  min_p_functional = safe_min(results$p_functional),
  min_p_q_snp = safe_min(results$q_snp_p),
  z_correlation_all = suppressWarnings(cor(
    results$z_affective,
    results$z_functional,
    use = "complete.obs"
  )),
  beta_correlation_all = suppressWarnings(cor(
    results$beta_affective,
    results$beta_functional,
    use = "complete.obs"
  ))
)

significant_results <- results[
  affective_gws | functional_gws | q_snp_gws
]

fwrite(
  results,
  file.path(output_dir, "trackB_twofactor_gwas_qsnp.csv"),
  na = ""
)
fwrite(
  significant_results,
  file.path(output_dir, "trackB_twofactor_significant_results.csv"),
  na = ""
)
write_csv(
  summary_table,
  file.path(output_dir, "trackB_twofactor_gwas_summary.csv"),
  na = ""
)
saveRDS(
  results,
  file.path(output_dir, "trackB_twofactor_gwas_qsnp.rds"),
  compress = FALSE
)

make_manhattan_plot(
  results,
  "p_affective",
  "Affective factor GWAS",
  file.path(figure_dir, "figure_trackB_affective_factor_manhattan.png")
)
make_manhattan_plot(
  results,
  "p_functional",
  "Functional factor GWAS",
  file.path(figure_dir, "figure_trackB_functional_factor_manhattan.png")
)
make_manhattan_plot(
  results,
  "q_snp_p",
  "Two-factor residual heterogeneity",
  file.path(figure_dir, "figure_trackB_twofactor_qsnp_manhattan.png")
)
make_qq_plot(
  results$p_affective,
  "Affective factor GWAS",
  file.path(figure_dir, "figure_trackB_affective_factor_qq.png")
)
make_qq_plot(
  results$p_functional,
  "Functional factor GWAS",
  file.path(figure_dir, "figure_trackB_functional_factor_qq.png")
)
make_qq_plot(
  results$q_snp_p,
  "Two-factor residual heterogeneity",
  file.path(figure_dir, "figure_trackB_twofactor_qsnp_qq.png")
)

writeLines(
  c(
    "Track B correlated two-factor GWAS and Q SNP v1",
    "The analysis reuses the completed five-trait SNP summary object and LDSC covariance structure.",
    "The affective factor is defined by major depression and loneliness.",
    "The functional factor is defined by frailty, self-rated health and walking pace.",
    "Both factors are regressed on SNP simultaneously and their covariance is freely estimated.",
    "Measurement parameters are fixed to the model fitted without SNP effects.",
    "The common projection GWAS outputs are treated as read-only inputs.",
    "Do not update downstream candidate genes until the two-factor results have been compared with the one-factor results."
  ),
  file.path(output_dir, "trackB_twofactor_gwas_note.txt")
)

message("Two-factor GWAS completed.")
message("Outputs written to: ", output_dir)
print(summary_table)


