## ============================================================================
## Track C Phase 1b: covariate-adjusted endpoint reassessment
## User-run only. This script does not alter Phase-1 outputs.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(logistf)
  library(jsonlite)
  library(digest)
})

options(stringsAsFactors = FALSE)

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
trackc_root <- Sys.getenv(
  "TRACKC_ROOT",
  unset = file.path(base_dir, "analysis_ready_core", "trackC")
)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
trackc_root <- normalizePath(trackc_root, winslash = "/", mustWork = TRUE)
trackb <- file.path(base_dir, "analysis_ready_core", "trackB")
frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
spec_dir <- file.path(trackc_root, "01_analysis_specification")
phase1_dir <- file.path(trackc_root, "02_endpoint_reassessment")
out_dir <- Sys.getenv(
  "TRACKC_PHASE1B_OUT_DIR",
  unset = file.path(trackc_root, "03_endpoint_covariate_adjusted")
)

required_markers <- c(
  file.path(frozen_dir, "_STAGE_COMPLETE.json"),
  file.path(spec_dir, "_STAGE_COMPLETE.json"),
  file.path(phase1_dir, "_STAGE_COMPLETE.json")
)
if (!all(file.exists(required_markers))) {
  stop("Phase 1 is incomplete. Required COMPLETE markers are missing.")
}
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Overwrite protection: Phase 1b output directory already contains files: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_json(
  list(stage = "phase1b_adjusted", status = "INCOMPLETE_OR_RUNNING", started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  file.path(out_dir, "_STAGE_INCOMPLETE.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

write_csv_utf8 <- function(x, path) fwrite(x, path, bom = TRUE, na = "")
bh_adjust <- function(p) p.adjust(as.numeric(p), method = "BH")
sha256_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)

scale_robust <- function(x) {
  x <- as.numeric(x)
  finite <- is.finite(x)
  fill <- median(x[finite], na.rm = TRUE)
  x[!finite] <- fill
  med <- median(x, na.rm = TRUE)
  s <- 1.4826 * median(abs(x - med), na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- 1
  (x - med) / s
}

read_frozen_set <- function(set_name, version) {
  path <- file.path(frozen_dir, paste0(set_name, "_", version, ".csv"))
  unique(toupper(trimws(as.character(fread(path)$gene_symbol))))
}

prepare_source_data <- function() {
  paths <- c(
    D2_TPMI = file.path(trackb, "trackD2_tpmi_transancestry_validation", "tables", "trackD2_tpmi_gene_window_summary_all_phenotypes.csv"),
    D5_FinnGen_positive = file.path(trackb, "trackD5_finngen_clinical_endpoint_replication", "tables", "trackD5_finngen_gene_window_summary_all_endpoints.csv"),
    D6_FinnGen_controls = file.path(trackb, "trackD6_finngen_endpoint_specificity", "tables", "trackD6_negative_control_gene_window_summary_all_endpoints.csv")
  )
  data <- lapply(paths, fread)
  for (nm in names(data)) {
    data[[nm]][, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
    data[[nm]] <- data[[nm]][tolower(as.character(in_background_for_test)) %in% c("true", "1")]
    data[[nm]][, phenocode := as.character(phenocode)]
  }
  list(paths = paths, data = data)
}

build_model_data <- function(endpoint, version, set_genes, covariates, source_data) {
  resource <- if (endpoint$source == "D2_TPMI") "TPMI" else "FinnGen"
  flag <- if (resource == "TPMI") "in_tpmi_gene_window_detectable" else "in_finngen_gene_window_detectable"
  snp_col <- if (resource == "TPMI") "tpmi_median_window_snps" else "finngen_median_window_snps"
  background <- copy(covariates[get(flag) == TRUE])
  if (version == "non_mhc") background <- background[is_extended_mhc == FALSE]

  endpoint_rows <- copy(source_data[[endpoint$source]])[
    as.character(phenocode) == as.character(endpoint$phenocode)
  ]
  endpoint_rows[, sidak_gene_p_num := as.numeric(sidak_gene_p)]
  setorder(endpoint_rows, sidak_gene_p_num)
  endpoint_rows <- unique(endpoint_rows, by = "gene_symbol")
  hit_map <- endpoint_rows[, .(gene_symbol, hit = sidak_gene_p_num < 0.05)]
  model <- merge(background, hit_map, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  model[is.na(hit), hit := FALSE]
  ## Keep the exposure numeric so logistf uses the stable coefficient name
  ## "in_set" (a logical exposure is expanded to "in_setTRUE").
  model[, in_set := as.integer(gene_symbol %in% set_genes)]
  model[, chr_factor := factor(chr)]
  model[, z_log_gene_length := scale_robust(log1p(as.numeric(gene_length_bp)))]
  model[, z_log_gene_density := scale_robust(log1p(as.numeric(gene_density_plusminus_1mb)))]
  model[, z_log_magma_nsnps := scale_robust(log1p(as.numeric(magma_nsnps_median)))]
  model[, z_effective_parameter_ratio := scale_robust(as.numeric(effective_parameter_ratio))]
  model[, z_log_resource_window_snps := scale_robust(log1p(as.numeric(get(snp_col))))]
  list(resource = resource, data = model)
}

fit_firth_model <- function(model_data) {
  formula <- hit ~ in_set + z_log_gene_length + z_log_gene_density +
    z_log_magma_nsnps + z_effective_parameter_ratio +
    z_log_resource_window_snps + chr_factor
  fit <- logistf(
    formula,
    data = model_data,
    pl = TRUE,
    plconf = 2,
    control = logistf.control(maxit = 50, maxstep = 5),
    plcontrol = logistpl.control(maxit = 100, maxstep = 5)
  )
  fit_conv <- max(abs(as.numeric(fit$conv)), na.rm = TRUE)
  profile_conv <- max(abs(as.numeric(fit$pl.conv)), na.rm = TRUE)
  coefficient_name <- if ("in_set" %in% names(fit$coefficients)) {
    "in_set"
  } else {
    candidates <- grep("^in_set", names(fit$coefficients), value = TRUE)
    if (length(candidates) != 1L) {
      stop(
        "Could not identify the frozen-set coefficient. Available coefficients: ",
        paste(names(fit$coefficients), collapse = ", ")
      )
    }
    candidates
  }
  beta <- unname(fit$coefficients[coefficient_name])
  p_value <- unname(fit$prob[coefficient_name])
  ci_low <- exp(unname(fit$ci.lower[coefficient_name]))
  ci_high <- exp(unname(fit$ci.upper[coefficient_name]))
  if (anyNA(c(beta, p_value, ci_low, ci_high))) {
    stop("Firth model completed but the frozen-set estimate could not be extracted.")
  }
  data.table(
    adjusted_log_or = beta,
    adjusted_or = exp(beta),
    adjusted_ci95_low = ci_low,
    adjusted_ci95_high = ci_high,
    adjusted_ci_ratio = ci_high / pmax(ci_low, .Machine$double.xmin),
    adjusted_log_ci_width = log(ci_high) - log(pmax(ci_low, .Machine$double.xmin)),
    adjusted_ci_ratio_gt_100 = ci_high / pmax(ci_low, .Machine$double.xmin) > 100,
    adjusted_p = p_value,
    firth_iterations = paste(fit$iter, collapse = ";"),
    profile_iterations = paste(fit$pl.iter, collapse = ";"),
    fit_convergence_max = fit_conv,
    profile_convergence_max = profile_conv,
    fit_and_profile_converged_1e_5 = fit_conv <= 1e-5 & profile_conv <= 1e-5
  )
}

calculate_design_diagnostics <- function(model_data) {
  focal <- c(
    "in_set", "z_log_gene_length", "z_log_gene_density",
    "z_log_magma_nsnps", "z_effective_parameter_ratio",
    "z_log_resource_window_snps"
  )
  vif_rows <- rbindlist(lapply(focal, function(variable) {
    others <- setdiff(focal, variable)
    formula <- as.formula(paste(variable, "~", paste(c(others, "chr_factor"), collapse = " + ")))
    fit <- lm(formula, data = model_data)
    r2 <- summary(fit)$r.squared
    data.table(
      predictor = variable,
      auxiliary_r_squared = r2,
      vif = 1 / pmax(1 - r2, .Machine$double.eps),
      vif_gt_5 = 1 / pmax(1 - r2, .Machine$double.eps) > 5,
      vif_gt_10 = 1 / pmax(1 - r2, .Machine$double.eps) > 10
    )
  }))
  numeric_matrix <- as.matrix(model_data[, ..focal])
  correlation <- cor(numeric_matrix, use = "pairwise.complete.obs")
  correlation_long <- as.data.table(as.table(correlation))
  setnames(correlation_long, c("predictor_1", "predictor_2", "correlation"))
  design <- model.matrix(
    ~ in_set + z_log_gene_length + z_log_gene_density + z_log_magma_nsnps +
      z_effective_parameter_ratio + z_log_resource_window_snps + chr_factor,
    data = model_data
  )
  condition <- data.table(
    design_rows = nrow(design),
    design_columns = ncol(design),
    design_matrix_rank = qr(design)$rank,
    design_condition_number = kappa(design, exact = FALSE),
    rank_deficient = qr(design)$rank < ncol(design)
  )
  list(vif = vif_rows, correlation = correlation_long, condition = condition)
}

covariates_path <- file.path(spec_dir, "magma_symbol_universe_and_matching_covariates.csv")
endpoints_path <- file.path(spec_dir, "frozen_endpoint_manifest.csv")
phase1_results_path <- file.path(phase1_dir, "trackC_endpoint_reassessment_results.csv")
covariates <- fread(covariates_path)
covariates[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
endpoints <- fread(endpoints_path, colClasses = list(character = "phenocode"))
phase1_results <- fread(phase1_results_path, colClasses = list(character = "phenocode"))
source_bundle <- prepare_source_data()

set_names <- c("cross_model_core", "affective_prioritized", "functional_prioritized")
versions <- c("full", "non_mhc")
results <- list()
failures <- list()
design_vif <- list()
design_correlations <- list()
design_conditions <- list()
diagnosed_design_keys <- character()
counter <- 0L

message("Phase 1b: fitting 138 Firth logistic models. This may take several minutes.")
for (i in seq_len(nrow(endpoints))) {
  endpoint <- endpoints[i]
  for (set_name in set_names) {
    for (version in versions) {
      counter <- counter + 1L
      message(sprintf("[%d/138] %s | %s | %s", counter, endpoint$phenocode, set_name, version))
      set_genes <- read_frozen_set(set_name, version)
      built <- build_model_data(endpoint, version, set_genes, covariates, source_bundle$data)
      design_key <- paste(built$resource, set_name, version, sep = "::")
      if (!design_key %in% diagnosed_design_keys) {
        design_diag <- calculate_design_diagnostics(built$data)
        id <- data.table(resource = built$resource, set_name = set_name, set_version = version)
        design_vif[[length(design_vif) + 1L]] <- cbind(id, design_diag$vif)
        design_correlations[[length(design_correlations) + 1L]] <- cbind(id, design_diag$correlation)
        design_conditions[[length(design_conditions) + 1L]] <- cbind(id, design_diag$condition)
        diagnosed_design_keys <- c(diagnosed_design_keys, design_key)
      }
      base_row <- data.table(
        test_family = endpoint$test_family,
        source = endpoint$source,
        phenocode = endpoint$phenocode,
        phenotype = endpoint$phenotype,
        validation_domain = endpoint$validation_domain,
        endpoint_family = endpoint$endpoint_family,
        endpoint_specificity_class = endpoint$endpoint_specificity_class,
        primary_role = endpoint$primary_role,
        resource = built$resource,
        set_name = set_name,
        set_version = version,
        background_mhc_excluded = version == "non_mhc",
        background_n = nrow(built$data),
        set_n_in_background = sum(built$data$in_set),
        set_hit_n = sum(built$data$hit & built$data$in_set),
        total_hit_n = sum(built$data$hit)
      )
      fitted <- tryCatch(
        fit_firth_model(built$data),
        error = function(e) e
      )
      if (inherits(fitted, "error")) {
        failures[[length(failures) + 1L]] <- cbind(
          base_row,
          data.table(error_message = conditionMessage(fitted))
        )
      } else {
        results[[length(results) + 1L]] <- cbind(base_row, fitted)
      }
    }
  }
}

result <- rbindlist(results, fill = TRUE)
failure_table <- rbindlist(failures, fill = TRUE)
if (!ncol(failure_table)) {
  failure_table <- data.table(
    source = character(), phenocode = character(), phenotype = character(),
    set_name = character(), set_version = character(), error_message = character()
  )
}
if (nrow(result)) {
  result[, adjusted_bh_within_family_version := bh_adjust(adjusted_p), by = .(test_family, set_version)]
  result[, adjusted_global_bh_primary := NA_real_]
  result[set_version == "full", adjusted_global_bh_primary := bh_adjust(adjusted_p)]

  compare_columns <- c(
    "source", "phenocode", "set_name", "set_version",
    "fisher_exact_conditional_mle_odds_ratio", "fisher_exact_ci95_low",
    "fisher_exact_ci95_high", "fisher_p_greater",
    "fisher_bh_within_family_version", "fisher_global_bh_primary_endpoints",
    "matched_permutation_p", "permutation_bh_within_family_version",
    "permutation_global_bh_primary_endpoints"
  )
  comparison <- merge(
    result,
    phase1_results[, ..compare_columns],
    by = c("source", "phenocode", "set_name", "set_version"),
    all.x = TRUE,
    sort = FALSE
  )
  comparison[, adjusted_family_significant := adjusted_bh_within_family_version < 0.05]
  comparison[, fisher_family_significant := fisher_bh_within_family_version < 0.05]
  comparison[, permutation_family_significant := permutation_bh_within_family_version < 0.05]
  comparison[, three_method_concordance :=
    adjusted_family_significant == fisher_family_significant &
    adjusted_family_significant == permutation_family_significant]
  comparison[, fisher_firth_or_direction_concordant := fifelse(
    is.finite(fisher_exact_conditional_mle_odds_ratio) & fisher_exact_conditional_mle_odds_ratio > 0,
    sign(log(fisher_exact_conditional_mle_odds_ratio)) == sign(adjusted_log_or),
    NA
  )]
  comparison[, interpretation_change := fifelse(
    adjusted_family_significant & !fisher_family_significant,
    "emerges_after_covariate_adjustment",
    fifelse(
      !adjusted_family_significant & fisher_family_significant,
      "attenuates_after_covariate_adjustment",
      "same_family_significance_state"
    )
  )]
  setorder(comparison, test_family, set_version, set_name, adjusted_p)
  write_csv_utf8(result, file.path(out_dir, "trackC_firth_adjusted_results.csv"))
  write_csv_utf8(comparison, file.path(out_dir, "trackC_firth_vs_phase1_comparison.csv"))

  specificity <- comparison[, .(
    tests = .N,
    adjusted_family_significant_n = sum(adjusted_family_significant, na.rm = TRUE),
    fisher_family_significant_n = sum(fisher_family_significant, na.rm = TRUE),
    permutation_family_significant_n = sum(permutation_family_significant, na.rm = TRUE),
    minimum_adjusted_p = min(adjusted_p, na.rm = TRUE),
    minimum_adjusted_bh = min(adjusted_bh_within_family_version, na.rm = TRUE)
  ), by = .(source, set_name, set_version)]
  write_csv_utf8(specificity, file.path(out_dir, "trackC_adjusted_specificity_summary.csv"))

  comparator_summary <- comparison[, .(
    positive_endpoint_significant_n = sum(source == "D5_FinnGen_positive" & adjusted_family_significant, na.rm = TRUE),
    cross_domain_comparator_significant_n = sum(source == "D6_FinnGen_controls" & adjusted_family_significant, na.rm = TRUE),
    tpmi_endpoint_significant_n = sum(source == "D2_TPMI" & adjusted_family_significant, na.rm = TRUE),
    endpoint_selectivity_interpretation = fifelse(
      sum(source == "D6_FinnGen_controls" & adjusted_family_significant, na.rm = TRUE) > 0,
      "cross_domain_comparator_enrichment_detected",
      "no_adjusted_comparator_endpoint_enrichment_detected"
    )
  ), by = .(set_name, set_version)]
  write_csv_utf8(comparator_summary, file.path(out_dir, "trackC_adjusted_comparator_endpoint_summary.csv"))
}

write_csv_utf8(failure_table, file.path(out_dir, "trackC_firth_model_failures.csv"))
write_csv_utf8(rbindlist(design_vif, fill = TRUE), file.path(out_dir, "phase1b_predictor_vif.csv"))
write_csv_utf8(rbindlist(design_correlations, fill = TRUE), file.path(out_dir, "phase1b_predictor_correlations.csv"))
write_csv_utf8(rbindlist(design_conditions, fill = TRUE), file.path(out_dir, "phase1b_design_condition_numbers.csv"))
input_paths <- c(covariates_path, endpoints_path, phase1_results_path, source_bundle$paths)
input_audit <- data.table(
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv_utf8(input_audit, file.path(out_dir, "phase1b_input_provenance_sha256.csv"))

audit <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  planned_models = 138L,
  completed_models = nrow(result),
  completed_models_with_nonmissing_estimates = sum(
    complete.cases(result[, .(adjusted_log_or, adjusted_or, adjusted_ci95_low, adjusted_ci95_high, adjusted_p)])
  ),
  failed_models = nrow(failure_table),
  method = "Firth penalized logistic regression with profile penalized-likelihood CI for frozen-set membership",
  outcome = "gene-window Sidak P < 0.05",
  covariates = c(
    "chromosome factor", "log gene length", "log local gene density",
    "log MAGMA NSNPS", "MAGMA NPARAM/NSNPS", "log resource median window SNP count"
  ),
  multiple_testing = "BH within scientific family and set version; global BH for full-set primary analyses",
  non_mhc_rule = "MHC removed from target membership and the full model background",
  diagnostics = "model/profile convergence, profile CI width, Fisher-Firth OR direction, predictor correlations, VIF, rank and condition number",
  comparator_endpoint_label = "D6 endpoints are cross-domain/non-neuroimmune comparators, not strict biological negative controls",
  role = "covariate-adjusted robustness analysis; does not repair or replace the Phase-1 matched permutation"
)
write_json(audit, file.path(out_dir, "phase1b_adjusted_audit.json"), pretty = TRUE, auto_unbox = TRUE)

valid_estimate_rows <- sum(
  complete.cases(result[, .(adjusted_log_or, adjusted_or, adjusted_ci95_low, adjusted_ci95_high, adjusted_p)])
)
if (nrow(failure_table) > 0L || nrow(result) != 138L || valid_estimate_rows != 138L) {
  stop("Phase 1b produced failed or missing models. See trackC_firth_model_failures.csv. INCOMPLETE marker retained.")
}

write_json(
  list(
    stage = "phase1b_adjusted",
    status = "COMPLETE",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    completed_models = nrow(result),
    completed_models_with_nonmissing_estimates = valid_estimate_rows,
    failed_models = 0L
  ),
  file.path(out_dir, "_STAGE_COMPLETE.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
unlink(file.path(out_dir, "_STAGE_INCOMPLETE.json"))
message("Track C Phase 1b completed: ", out_dir)


