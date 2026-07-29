## ============================================================================
## Track C C3b: uniform covariate adjustment for every imbalanced C3 combo
## User-run only. Adjustment scope is fixed by C3 balance triggers, not results.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(logistf)
  library(jsonlite)
  library(digest)
})

options(stringsAsFactors = FALSE)

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
trackc_root <- Sys.getenv("TRACKC_ROOT",
                      unset = file.path(base_dir, "analysis_ready_core", "trackC"))
c3_dir <- Sys.getenv("TRACKC_C3_DIR", unset = file.path(trackc_root, "06_C3_external_omics"))
out_dir <- Sys.getenv("TRACKC_C3B_OUT_DIR", unset = file.path(trackc_root, "07_C3b_uniform_adjustment"))

base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
trackc_root <- normalizePath(trackc_root, winslash = "/", mustWork = TRUE)
c3_dir <- normalizePath(c3_dir, winslash = "/", mustWork = TRUE)
frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
spec_dir <- file.path(trackc_root, "01_analysis_specification")
c3_dir <- file.path(base_dir, "analysis_ready_core", "trackB",
                       "trackC3_triangular_validation_immune_spatial", "tables")
spatial_m7_path <- file.path(base_dir, "外部验证", "空间蛋白组资源",
                             "41586_2026_10660_MOESM7_ESM.xlsx")

if (!file.exists(file.path(c3_dir, "_STAGE_COMPLETE.json")) ||
    !file.exists(file.path(c3_dir, "_C3B_REQUIRED.json"))) {
  stop("C3 COMPLETE/C3B_REQUIRED markers are missing.")
}
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Overwrite protection: output directory already contains files: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_json(
  list(stage = "C3b_uniform_adjustment", status = "INCOMPLETE_OR_RUNNING",
       started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  file.path(out_dir, "_STAGE_INCOMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)

write_csv_utf8 <- function(x, path) fwrite(x, path, bom = TRUE, na = "")
sha256_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)
bh_adjust <- function(p) p.adjust(as.numeric(p), method = "BH")

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

covariates_path <- file.path(spec_dir, "magma_symbol_universe_and_matching_covariates.csv")
covariates <- fread(covariates_path)
covariates[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]

background_paths <- c(
  immune_age_reported_DEG = file.path(spec_dir, "background_immune_age_reported_gene_universe.csv"),
  immune_marker = file.path(spec_dir, "background_immune_marker_reported_gene_universe.csv"),
  spatial_proteome = file.path(spec_dir, "background_spatial_proteome_detectable.csv")
)
backgrounds <- lapply(background_paths, function(path) {
  unique(toupper(trimws(as.character(fread(path)$gene_symbol))))
})

c3_results_path <- file.path(c3_dir, "C3_external_omics_enrichment_results.csv")
c3_trigger_path <- file.path(c3_dir, "C3_combination_balance_and_C3b_triggers.csv")
c3_results <- fread(c3_results_path)
triggers <- fread(c3_trigger_path)[requires_uniform_C3b_adjustment == TRUE]
if (nrow(triggers) != 12L) stop("Expected exactly 12 frozen C3b trigger combinations.")
if (any(triggers$coverage_n < 10L)) stop("A C3b trigger has insufficient inferential coverage.")
write_csv_utf8(triggers, file.path(out_dir, "C3b_frozen_trigger_combinations.csv"))

read_frozen_set <- function(set_name, set_version) {
  path <- file.path(frozen_dir, paste0(set_name, "_", set_version, ".csv"))
  unique(toupper(trimws(as.character(fread(path)$gene_symbol))))
}

## Reconstruct the exact C3 signatures without selecting on C3 P values.
age_path <- file.path(c3_dir, "trackC3_immune_age_deg_all_linear.csv")
marker_path <- file.path(c3_dir, "trackC3_immune_marker_all.csv")
age <- fread(age_path)
marker <- fread(marker_path)
age[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
marker[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
age_rows <- unique(age[sig_age_deg == TRUE, .(
  resource = "immune_age_reported_DEG",
  signature_id = paste("age", sex, celltype, direction, sep = "|"), gene_symbol
)])
marker_rows <- unique(marker[sig_marker == TRUE, .(
  resource = "immune_marker",
  signature_id = paste("marker", celltype, direction_nhood, sep = "|"), gene_symbol
)])

spatial_raw <- as.data.table(read_excel(spatial_m7_path, sheet = "A. All tissue enriched proteins"))
setnames(spatial_raw, "Gene", "gene_symbol")
spatial_raw[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
id_columns <- c("Tissue enriched protein", "gene_symbol")
tissue_columns <- setdiff(names(spatial_raw), id_columns)
for (v in tissue_columns) spatial_raw[, (v) := as.numeric(get(v))]
spatial_long <- melt(spatial_raw, id.vars = id_columns, measure.vars = tissue_columns,
                     variable.name = "tissue", value.name = "enrichment_score")
spatial_tissue_rows <- unique(spatial_long[is.finite(enrichment_score) & enrichment_score > 0,
  .(resource = "spatial_proteome", signature_id = paste0("tissue|", tissue), gene_symbol)])
spatial_categories <- list(
  brain_neural = c("Brain", "Nerve", "Spinal cord"),
  immune_lymphoid = c("Bone marrow", "Lymph node", "Thymus", "Spleen", "Tonsil"),
  musculoskeletal = c("Bone", "Cartilage", "Skeletal muscle", "Smooth muscle", "Tendon")
)
spatial_category_rows <- rbindlist(lapply(names(spatial_categories), function(category) {
  unique(spatial_long[tissue %in% spatial_categories[[category]] &
                        is.finite(enrichment_score) & enrichment_score > 0,
                      .(resource = "spatial_proteome",
                        signature_id = paste0("category|", category), gene_symbol)])
}))
signature_rows <- unique(rbindlist(list(age_rows, marker_rows,
                                        spatial_tissue_rows, spatial_category_rows)))

signature_manifest <- unique(c3_results[, .(
  resource, signature_id, signature_label, signature_group,
  analysis_layer, correction_family, interpretation_boundary
)])
reconstructed_counts <- signature_rows[, .(reconstructed_genes = uniqueN(gene_symbol)),
                                       by = .(resource, signature_id)]
manifest_check <- merge(signature_manifest, reconstructed_counts,
                        by = c("resource", "signature_id"), all.x = TRUE)
if (nrow(manifest_check) != 145L || anyNA(manifest_check$reconstructed_genes)) {
  stop("Reconstructed signatures do not reproduce the 145-signature C3 manifest.")
}
write_csv_utf8(manifest_check, file.path(out_dir, "C3b_signature_reconstruction_check.csv"))

prepare_universe <- function(resource_name, set_name, set_version) {
  background_genes <- backgrounds[[resource_name]]
  if (set_version == "non_mhc") {
    background_genes <- setdiff(background_genes,
                                covariates[is_extended_mhc == TRUE, gene_symbol])
  }
  d <- copy(covariates[gene_symbol %in% background_genes])
  set_genes <- read_frozen_set(set_name, set_version)
  d[, in_set := as.integer(gene_symbol %in% set_genes)]
  d[, chr_factor := factor(chr)]
  d[, z_log_gene_length := scale_robust(log1p(as.numeric(gene_length_bp)))]
  d[, z_log_gene_density := scale_robust(log1p(as.numeric(gene_density_plusminus_1mb)))]
  d[, z_log_magma_nsnps := scale_robust(log1p(as.numeric(magma_nsnps_median)))]
  d[, z_effective_parameter_ratio := scale_robust(as.numeric(effective_parameter_ratio))]
  d
}

fit_firth <- function(d) {
  fit <- logistf(
    hit ~ in_set + z_log_gene_length + z_log_gene_density +
      z_log_magma_nsnps + z_effective_parameter_ratio + chr_factor,
    data = d, pl = TRUE, plconf = 2,
    control = logistf.control(maxit = 50, maxstep = 5),
    plcontrol = logistpl.control(maxit = 100, maxstep = 5)
  )
  nm <- if ("in_set" %in% names(fit$coefficients)) "in_set" else
    grep("^in_set", names(fit$coefficients), value = TRUE)
  if (length(nm) != 1L) stop("Frozen-set coefficient not uniquely identified.")
  beta <- unname(fit$coefficients[nm])
  low <- exp(unname(fit$ci.lower[nm]))
  high <- exp(unname(fit$ci.upper[nm]))
  p <- unname(fit$prob[nm])
  fit_conv <- max(abs(as.numeric(fit$conv)), na.rm = TRUE)
  profile_conv <- max(abs(as.numeric(fit$pl.conv)), na.rm = TRUE)
  if (anyNA(c(beta, low, high, p))) stop("Missing focal Firth estimate.")
  data.table(
    adjusted_log_or = beta, adjusted_or = exp(beta),
    adjusted_ci95_low = low, adjusted_ci95_high = high,
    adjusted_ci_ratio = high / pmax(low, .Machine$double.xmin),
    adjusted_ci_ratio_gt_100 = high / pmax(low, .Machine$double.xmin) > 100,
    adjusted_p = p,
    fit_convergence_max = fit_conv, profile_convergence_max = profile_conv,
    fit_and_profile_converged_1e_5 = fit_conv <= 1e-5 & profile_conv <= 1e-5
  )
}

design_diagnostics <- function(d) {
  focal <- c("in_set", "z_log_gene_length", "z_log_gene_density",
             "z_log_magma_nsnps", "z_effective_parameter_ratio")
  vif <- rbindlist(lapply(focal, function(variable) {
    others <- setdiff(focal, variable)
    fit <- lm(as.formula(paste(variable, "~", paste(c(others, "chr_factor"), collapse = "+"))),
              data = d)
    r2 <- summary(fit)$r.squared
    data.table(predictor = variable, auxiliary_r_squared = r2,
               vif = 1 / pmax(1 - r2, .Machine$double.eps))
  }))
  design <- model.matrix(~ in_set + z_log_gene_length + z_log_gene_density +
                           z_log_magma_nsnps + z_effective_parameter_ratio + chr_factor,
                         data = d)
  condition <- data.table(
    design_rows = nrow(design), design_columns = ncol(design),
    design_matrix_rank = qr(design)$rank,
    design_condition_number = kappa(design, exact = FALSE),
    rank_deficient = qr(design)$rank < ncol(design)
  )
  list(vif = vif, condition = condition)
}

results <- list()
failures <- list()
vif_rows <- list()
condition_rows <- list()
counter <- 0L
expected_models <- sum(vapply(seq_len(nrow(triggers)), function(i) {
  sum(signature_manifest$resource == triggers$resource[i])
}, integer(1)))

for (i in seq_len(nrow(triggers))) {
  trigger <- triggers[i]
  d <- prepare_universe(trigger$resource, trigger$set_name, trigger$set_version)
  if (sum(d$in_set) != trigger$coverage_n) {
    stop("C3b set coverage does not reproduce C3 for trigger row ", i)
  }
  diag <- design_diagnostics(d)
  id <- trigger[, .(resource, set_name, set_version, coverage_n)]
  vif_rows[[length(vif_rows) + 1L]] <- cbind(id, diag$vif)
  condition_rows[[length(condition_rows) + 1L]] <- cbind(id, diag$condition)
  signatures <- signature_manifest[resource == trigger$resource]
  for (j in seq_len(nrow(signatures))) {
    counter <- counter + 1L
    signature <- signatures[j]
    genes <- signature_rows[resource == trigger$resource &
                              signature_id == signature$signature_id, gene_symbol]
    model <- copy(d)
    model[, hit := gene_symbol %in% genes]
    base_row <- data.table(
      resource = trigger$resource, signature_id = signature$signature_id,
      signature_label = signature$signature_label,
      signature_group = signature$signature_group,
      analysis_layer = signature$analysis_layer,
      correction_family = signature$correction_family,
      interpretation_boundary = signature$interpretation_boundary,
      set_name = trigger$set_name, set_version = trigger$set_version,
      background_mhc_excluded = trigger$set_version == "non_mhc",
      background_n = nrow(model), set_n_in_background = sum(model$in_set),
      signature_n_in_background = sum(model$hit),
      overlap_n = sum(model$hit & model$in_set == 1L)
    )
    ## A Firth coefficient can be returned for extremely rare signatures, but
    ## fewer than 10 outcome genes provides no stable adjusted enrichment test.
    ## Retain the row and overlap counts while withholding inferential estimates.
    if (sum(model$hit) < 10L || sum(model$hit) == nrow(model)) {
      fitted <- data.table(
        adjusted_log_or = NA_real_, adjusted_or = NA_real_,
        adjusted_ci95_low = NA_real_, adjusted_ci95_high = NA_real_,
        adjusted_ci_ratio = NA_real_, adjusted_ci_ratio_gt_100 = NA,
        adjusted_p = NA_real_, fit_convergence_max = NA_real_,
        profile_convergence_max = NA_real_, fit_and_profile_converged_1e_5 = NA,
        model_status = "noninferential_sparse_signature_lt10"
      )
    } else {
      fitted <- tryCatch(fit_firth(model), error = function(e) e)
      if (!inherits(fitted, "error")) fitted[, model_status := "estimated"]
    }
    if (inherits(fitted, "error")) {
      failures[[length(failures) + 1L]] <- cbind(base_row,
                                                 error_message = conditionMessage(fitted))
    } else {
      results[[length(results) + 1L]] <- cbind(base_row, fitted)
    }
    message(sprintf("[%d/%d] %s | %s | %s | %s", counter, expected_models,
                    trigger$resource, trigger$set_name, trigger$set_version,
                    signature$signature_id))
  }
}

result <- rbindlist(results, fill = TRUE)
failure_table <- rbindlist(failures, fill = TRUE)
if (!ncol(failure_table)) {
  failure_table <- data.table(resource = character(), signature_id = character(),
                              set_name = character(), set_version = character(),
                              error_message = character())
}

result[, adjusted_bh_within_family_version := bh_adjust(adjusted_p),
       by = .(correction_family, set_version)]
result[, adjusted_global_bh_primary_full := NA_real_]
primary_rows <- result$set_version == "full" &
  result$analysis_layer %in% c("primary_single_cell", "primary_spatial_category")
result[which(primary_rows), adjusted_global_bh_primary_full := bh_adjust(adjusted_p)]

c3_compare <- c3_results[, .(
  resource, signature_id, set_name, set_version,
  fisher_odds_ratio, fisher_ci95_low, fisher_ci95_high,
  fisher_p_greater, fisher_bh_within_family_version,
  matched_permutation_p, permutation_bh_within_family_version,
  fisher_global_bh_full, permutation_global_bh_full
)]
comparison <- merge(result, c3_compare,
                    by = c("resource", "signature_id", "set_name", "set_version"),
                    all.x = TRUE, sort = FALSE)
comparison[, fisher_firth_direction_concordant := fifelse(
  is.finite(fisher_odds_ratio) & fisher_odds_ratio > 0,
  sign(log(fisher_odds_ratio)) == sign(adjusted_log_or), NA
)]
comparison[, `:=`(
  adjusted_family_significant = adjusted_bh_within_family_version < 0.05,
  fisher_family_significant = fisher_bh_within_family_version < 0.05,
  permutation_family_significant = permutation_bh_within_family_version < 0.05
)]
comparison[, three_method_family_concordance :=
             adjusted_family_significant == fisher_family_significant &
             adjusted_family_significant == permutation_family_significant]

summary_table <- comparison[, .(
  tests = .N,
  noninferential_sparse_or_structural_n = sum(model_status != "estimated"),
  adjusted_family_significant_n = sum(adjusted_family_significant, na.rm = TRUE),
  fisher_family_significant_n = sum(fisher_family_significant, na.rm = TRUE),
  permutation_family_significant_n = sum(permutation_family_significant, na.rm = TRUE),
  three_method_concordant_n = sum(three_method_family_concordance, na.rm = TRUE),
  minimum_adjusted_p = if (all(is.na(adjusted_p))) NA_real_ else min(adjusted_p, na.rm = TRUE),
  minimum_adjusted_bh = if (all(is.na(adjusted_bh_within_family_version))) NA_real_ else
    min(adjusted_bh_within_family_version, na.rm = TRUE)
), by = .(resource, analysis_layer, set_name, set_version)]

setorder(comparison, correction_family, set_version, set_name, adjusted_p)
write_csv_utf8(result, file.path(out_dir, "C3b_uniform_adjusted_results.csv"))
write_csv_utf8(comparison, file.path(out_dir, "C3b_adjusted_vs_C3_comparison.csv"))
write_csv_utf8(summary_table, file.path(out_dir, "C3b_uniform_adjusted_summary.csv"))
write_csv_utf8(rbindlist(vif_rows, fill = TRUE), file.path(out_dir, "C3b_predictor_vif.csv"))
write_csv_utf8(rbindlist(condition_rows, fill = TRUE), file.path(out_dir, "C3b_design_condition_numbers.csv"))
write_csv_utf8(failure_table, file.path(out_dir, "C3b_model_failures.csv"))

input_paths <- c(covariates_path, background_paths, c3_results_path, c3_trigger_path,
                 age_path, marker_path, spatial_m7_path,
                 file.path(frozen_dir,
                           paste0(rep(c("cross_model_core", "affective_prioritized",
                                        "functional_prioritized"), each = 2), "_",
                                  rep(c("full", "non_mhc"), 3), ".csv")))
input_audit <- data.table(
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv_utf8(input_audit, file.path(out_dir, "C3b_input_provenance_sha256.csv"))

valid_models <- sum(complete.cases(result[, .(adjusted_log_or, adjusted_or,
                                               adjusted_ci95_low, adjusted_ci95_high,
                                               adjusted_p)]))
noninferential_models <- sum(result$model_status != "estimated")
audit <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  trigger_rule = "all signatures in each of the 12 imbalanced inferential C3 resource-set-version combinations",
  expected_models = expected_models, completed_models = nrow(result),
  models_with_nonmissing_estimates = valid_models,
  noninferential_sparse_or_structural_models = noninferential_models,
  adjusted_outcome_coverage_rule = "external signatures with fewer than 10 background genes retain descriptive overlaps but no adjusted inferential estimate",
  failed_models = nrow(failure_table),
  model = "Firth penalized logistic regression",
  outcome = "membership in the fixed external signature",
  covariates = c("chromosome factor", "log gene length", "log local gene density",
                 "log MAGMA NSNPS", "MAGMA effective-parameter ratio"),
  multiple_testing = "BH within original correction family and set version; global BH across full-set primary C3b models",
  role = "uniform covariate-adjusted robustness analysis; does not repair matched permutations",
  frozen_gene_membership_changed = FALSE,
  d3_rule = "D3 remains blocked until C3b results and QA are interpreted"
)
write_json(audit, file.path(out_dir, "C3b_audit.json"), pretty = TRUE, auto_unbox = TRUE)

if (nrow(failure_table) > 0L || nrow(result) != expected_models ||
    valid_models + noninferential_models != expected_models ||
    !all(result[model_status == "estimated"]$fit_and_profile_converged_1e_5)) {
  stop("C3b failed QA. INCOMPLETE marker retained; inspect diagnostics.")
}
write_json(
  list(stage = "C3b_uniform_adjustment", status = "COMPLETE",
       completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
       triggered_combinations = nrow(triggers), models = nrow(result),
       estimated_models = valid_models,
       noninferential_sparse_or_structural_models = noninferential_models, failures = 0L),
  file.path(out_dir, "_STAGE_COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)
unlink(file.path(out_dir, "_STAGE_INCOMPLETE.json"))
message("Track C C3b completed: ", out_dir)
