options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(GenomicSEM)
})

resolve_base_dir <- function() {
  base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(base_dir) || !dir.exists(base_dir)) {
    stop("Set DCV_BASE_DIR to the project directory containing analysis_ready_core.")
  }
  normalizePath(base_dir, winslash = "/", mustWork = TRUE)
}

safe_usermodel <- function(covstruc, model_syntax) {
  tryCatch(
    GenomicSEM::usermodel(
      covstruc = covstruc,
      estimation = "DWLS",
      model = model_syntax,
      CFIcalc = TRUE,
      std.lv = TRUE
    ),
    error = function(error) error
  )
}

extract_fit <- function(fit, model_name, model_role) {
  if (inherits(fit, "error") || !("modelfit" %in% names(fit))) {
    return(tibble(
      model_name = model_name,
      model_role = model_role,
      chisq = NA_real_, df = NA_real_, p_chisq = NA_real_,
      AIC = NA_real_, CFI = NA_real_, SRMR = NA_real_,
      status = "failed",
      error_message = if (inherits(fit, "error")) conditionMessage(fit) else "modelfit unavailable"
    ))
  }

  model_fit <- as.data.frame(fit$modelfit)
  get_value <- function(column) {
    if (!column %in% names(model_fit)) return(NA_real_)
    suppressWarnings(as.numeric(model_fit[[column]][1]))
  }

  tibble(
    model_name = model_name,
    model_role = model_role,
    chisq = get_value("chisq"),
    df = get_value("df"),
    p_chisq = get_value("p_chisq"),
    AIC = get_value("AIC"),
    CFI = get_value("CFI"),
    SRMR = get_value("SRMR"),
    status = "estimated",
    error_message = NA_character_
  )
}

extract_parameters <- function(fit, model_name) {
  if (inherits(fit, "error") || !("results" %in% names(fit))) return(tibble())
  parameters <- as_tibble(fit$results)
  transmute(
    parameters,
    model_name = model_name,
    lhs = as.character(lhs),
    op = as.character(op),
    rhs = as.character(rhs),
    unstand_est = suppressWarnings(as.numeric(Unstand_Est)),
    unstand_se = suppressWarnings(as.numeric(Unstand_SE)),
    std_genotype = suppressWarnings(as.numeric(STD_Genotype)),
    std_genotype_se = suppressWarnings(as.numeric(STD_Genotype_SE)),
    std_all = suppressWarnings(as.numeric(STD_All)),
    p_value = suppressWarnings(as.numeric(p_value))
  )
}

base_dir <- resolve_base_dir()
trackb_dir <- file.path(base_dir, "analysis_ready_core", "trackB")
input_file <- file.path(
  trackb_dir,
  "factor_gwas_qsnp",
  "trackB_afffunc_ldsc_covstruc.rds"
)
output_dir <- file.path(trackb_dir, "genomicsem_competing_models")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Missing LDSC covariance object: ", input_file)
covstruc <- readRDS(input_file)

# Use `functional_domain` to avoid parser ambiguity with the reserved R word
# `function` in lavaan and GenomicSEM model syntax.
model_syntax <- list(
  one_factor_primary = "
    aff_func =~ mdd + lonely + frailty + srh + walkpace
  ",
  correlated_affective_functional = "
    affective_domain =~ mdd + lonely
    functional_domain =~ frailty + srh + walkpace
    affective_domain ~~ functional_domain
  ",
  one_factor_functional_residual = "
    aff_func =~ mdd + lonely + frailty + srh + walkpace
    frailty ~~ srh
  "
)

model_roles <- c(
  one_factor_primary = "primary completed factor-GWAS model",
  correlated_affective_functional = "prespecified dimensional sensitivity",
  one_factor_functional_residual = "substantive residual-covariance sensitivity"
)

fits <- lapply(model_syntax, function(syntax) safe_usermodel(covstruc, syntax))
fit_summary <- bind_rows(Map(
  function(fit, name) extract_fit(fit, name, model_roles[[name]]),
  fits,
  names(fits)
))
parameter_summary <- bind_rows(Map(extract_parameters, fits, names(fits)))

diagnostics <- parameter_summary %>%
  group_by(model_name) %>%
  summarise(
    n_loading_abs_gt_1p05 = sum(op == "=~" & is.finite(std_all) & abs(std_all) > 1.05),
    n_negative_residual_variance = sum(op == "~~" & lhs == rhs & is.finite(std_all) & std_all < 0),
    n_nonfinite_standardised = sum(!is.finite(std_all)),
    n_nonsignificant_loadings = sum(op == "=~" & is.finite(p_value) & p_value >= 0.05),
    .groups = "drop"
  ) %>%
  mutate(
    proper_solution = n_loading_abs_gt_1p05 == 0 &
      n_negative_residual_variance == 0 &
      n_nonfinite_standardised == 0
  )

fit_summary <- fit_summary %>%
  left_join(diagnostics, by = "model_name") %>%
  mutate(
    delta_CFI_vs_primary = CFI - CFI[model_name == "one_factor_primary"][1],
    delta_SRMR_vs_primary = SRMR - SRMR[model_name == "one_factor_primary"][1]
  )

if (!"S" %in% names(covstruc)) {
  stop("The LDSC object does not contain the genetic covariance matrix `S`.")
}
genetic_covariance <- as.matrix(covstruc[["S"]])
traits <- c("mdd", "lonely", "frailty", "srh", "walkpace")
if (is.null(colnames(genetic_covariance)) || !all(traits %in% colnames(genetic_covariance))) {
  stop("The five expected traits are not present in the LDSC genetic covariance matrix `S`.")
}
if (is.null(rownames(genetic_covariance))) {
  rownames(genetic_covariance) <- colnames(genetic_covariance)
}
if (!all(traits %in% rownames(genetic_covariance))) {
  stop("The row order of the genetic covariance matrix cannot be matched to the five traits.")
}
genetic_covariance <- genetic_covariance[traits, traits, drop = FALSE]
if (any(!is.finite(genetic_covariance)) || any(diag(genetic_covariance) <= 0)) {
  stop("The genetic covariance matrix contains nonfinite values or nonpositive variances.")
}
genetic_sd <- sqrt(diag(genetic_covariance))
genetic_correlation <- genetic_covariance / outer(genetic_sd, genetic_sd)
eigenvalues <- eigen(genetic_correlation, symmetric = TRUE)$values

eigen_summary <- tibble(
  component = seq_along(eigenvalues),
  eigenvalue = eigenvalues,
  variance_proportion = eigenvalues / sum(eigenvalues),
  cumulative_variance = cumsum(eigenvalues / sum(eigenvalues))
)

factor_correlation <- parameter_summary %>%
  filter(
    model_name == "correlated_affective_functional",
    op == "~~",
    (lhs == "affective_domain" & rhs == "functional_domain") |
      (lhs == "functional_domain" & rhs == "affective_domain")
  ) %>%
  transmute(
    model_name,
    factor_correlation = std_all,
    standard_error = std_genotype_se,
    p_value
  )

primary_fit <- fit_summary %>% filter(model_name == "one_factor_primary")
two_factor_fit <- fit_summary %>% filter(model_name == "correlated_affective_functional")
residual_fit <- fit_summary %>% filter(model_name == "one_factor_functional_residual")

material_improvement <- function(candidate) {
  nrow(candidate) == 1 &&
    candidate$status == "estimated" &&
    isTRUE(candidate$proper_solution) &&
    (
      (!is.na(candidate$delta_CFI_vs_primary) && candidate$delta_CFI_vs_primary >= 0.01) ||
        (!is.na(candidate$delta_SRMR_vs_primary) && candidate$delta_SRMR_vs_primary <= -0.01)
    )
}

decision <- case_when(
  nrow(two_factor_fit) != 1 || two_factor_fit$status != "estimated" ~
    "The dimensional sensitivity model did not complete. Do not change the primary factor interpretation.",
  !isTRUE(two_factor_fit$proper_solution) ~
    "The dimensional sensitivity model was inadmissible. Retain the one-factor model and report the failed sensitivity.",
  material_improvement(two_factor_fit) ~
    "The two-factor model materially improved fit. Treat affective and functional dimensions as distinguishable and rerun factor GWAS before changing the primary manuscript.",
  material_improvement(residual_fit) ~
    "The two-factor model did not materially improve fit, but the prespecified functional residual covariance did. Retain the shared factor as the working model and run a residual-covariance factor-GWAS sensitivity analysis before final interpretation.",
  TRUE ~
    "The two-factor model did not materially improve fit. The one-factor model remains a parsimonious working representation, with residual heterogeneity reported through Q SNP."
)

write_csv(fit_summary, file.path(output_dir, "trackB_competing_model_fit.csv"), na = "")
write_csv(parameter_summary, file.path(output_dir, "trackB_competing_model_parameters.csv"), na = "")
write_csv(eigen_summary, file.path(output_dir, "trackB_five_trait_eigenstructure.csv"), na = "")
write_csv(
  as.data.frame(genetic_correlation) %>% rownames_to_column("trait"),
  file.path(output_dir, "trackB_five_trait_genetic_correlation_matrix.csv"),
  na = ""
)
write_csv(factor_correlation, file.path(output_dir, "trackB_two_factor_correlation.csv"), na = "")
writeLines(unname(decision), file.path(output_dir, "trackB_competing_model_decision.txt"))

plot_data <- fit_summary %>%
  filter(status == "estimated") %>%
  select(model_name, CFI, SRMR) %>%
  tidyr::pivot_longer(c(CFI, SRMR), names_to = "fit_index", values_to = "value") %>%
  mutate(
    model_label = recode(
      model_name,
      correlated_affective_functional = "Correlated two-factor",
      one_factor_functional_residual = "One-factor plus residual",
      one_factor_primary = "One-factor primary"
    ),
    reference = if_else(fit_index == "CFI", 0.95, 0.08),
    fit_index = factor(fit_index, levels = c("CFI", "SRMR"))
  )

fit_plot <- ggplot(plot_data, aes(model_label, value, fill = fit_index)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(
    data = distinct(plot_data, fit_index, reference),
    aes(yintercept = reference),
    linetype = 2,
    colour = "grey45"
  ) +
  facet_wrap(~ fit_index, scales = "free_y") +
  labs(x = NULL, y = "Fit index", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(
  file.path(figure_dir, "figure_trackB_competing_model_fit.png"),
  fit_plot,
  width = 8,
  height = 5,
  dpi = 600
)

message("Done. Outputs written to: ", output_dir)
message("Decision: ", decision)


