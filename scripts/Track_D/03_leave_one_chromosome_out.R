## ============================================================================
## Track C Phase 1c-A: targeted leave-one-chromosome-out influence analysis
## User-run only. Reads frozen non-MHC sets and valid Phase-1b final outputs.
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
phase1b_dir <- Sys.getenv(
  "TRACKC_PHASE1B_DIR",
  unset = file.path(trackc_root, "03_endpoint_covariate_adjusted")
)
out_dir <- Sys.getenv(
  "TRACKC_PHASE1C_OUT_DIR",
  unset = file.path(trackc_root, "04_targeted_loco")
)

base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
trackc_root <- normalizePath(trackc_root, winslash = "/", mustWork = TRUE)
phase1b_dir <- normalizePath(phase1b_dir, winslash = "/", mustWork = TRUE)
trackb <- file.path(base_dir, "analysis_ready_core", "trackB")
frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
spec_dir <- file.path(trackc_root, "01_analysis_specification")

required_markers <- c(
  file.path(frozen_dir, "_STAGE_COMPLETE.json"),
  file.path(spec_dir, "_STAGE_COMPLETE.json"),
  file.path(phase1b_dir, "_STAGE_COMPLETE.json")
)
if (!all(file.exists(required_markers))) {
  stop("Required freeze/spec/Phase-1b COMPLETE markers are missing.")
}
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Overwrite protection: output directory already contains files: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_json(
  list(stage = "phase1c_targeted_loco", status = "INCOMPLETE_OR_RUNNING",
       started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  file.path(out_dir, "_STAGE_INCOMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)

write_csv_utf8 <- function(x, path) fwrite(x, path, bom = TRUE, na = "")
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

## Ten representative, non-MHC-supported targets were frozen after Phase 1b.
## Overlapping FinnGen dementia definitions are represented by F5_DEMENTIA.
target_manifest <- data.table(
  target_id = sprintf("T%02d", 1:10),
  source = c(
    "D2_TPMI", "D6_FinnGen_controls",
    "D5_FinnGen_positive", "D5_FinnGen_positive",
    "D2_TPMI", "D2_TPMI", "D5_FinnGen_positive",
    "D5_FinnGen_positive", "D6_FinnGen_controls", "D6_FinnGen_controls"
  ),
  phenocode = c(
    "290.1", "T2D", "F5_DEMENTIA", "M13_RHEUMA",
    "290.1", "290.11", "F5_DEMENTIA",
    "M13_SPONDYLINFLAMNAS", "T2D", "M13_ARTHROSIS_KNEE"
  ),
  set_name = c(
    "cross_model_core", "cross_model_core",
    "affective_prioritized", "affective_prioritized",
    "functional_prioritized", "functional_prioritized",
    "functional_prioritized", "functional_prioritized",
    "functional_prioritized", "functional_prioritized"
  ),
  rationale = c(
    "cross-resource dementia anchor",
    "cross-domain comparator retained after non-MHC exclusion",
    "primary FinnGen dementia endpoint; overlapping dementia definitions not repeated",
    "non-MHC immune endpoint",
    "TPMI dementia anchor",
    "TPMI Alzheimer disease sensitivity endpoint",
    "primary FinnGen dementia endpoint",
    "non-MHC inflammatory endpoint",
    "cross-domain metabolic comparator",
    "cross-domain musculoskeletal comparator"
  ),
  set_version = "non_mhc"
)
write_csv_utf8(target_manifest, file.path(out_dir, "phase1c_frozen_target_manifest.csv"))

covariates_path <- file.path(spec_dir, "magma_symbol_universe_and_matching_covariates.csv")
endpoints_path <- file.path(spec_dir, "frozen_endpoint_manifest.csv")
phase1b_path <- file.path(phase1b_dir, "trackC_firth_vs_phase1_comparison.csv")
covariates <- fread(covariates_path)
covariates[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
endpoints <- fread(endpoints_path, colClasses = list(character = "phenocode"))
phase1b <- fread(phase1b_path, colClasses = list(character = "phenocode"))

if (nrow(phase1b) != 138L || anyNA(phase1b$adjusted_log_or)) {
  stop("Phase-1b input is not the valid 138-estimate r2 result.")
}

source_paths <- c(
  D2_TPMI = file.path(trackb, "trackD2_tpmi_transancestry_validation", "tables",
                      "trackD2_tpmi_gene_window_summary_all_phenotypes.csv"),
  D5_FinnGen_positive = file.path(trackb, "trackD5_finngen_clinical_endpoint_replication", "tables",
                                  "trackD5_finngen_gene_window_summary_all_endpoints.csv"),
  D6_FinnGen_controls = file.path(trackb, "trackD6_finngen_endpoint_specificity", "tables",
                                  "trackD6_negative_control_gene_window_summary_all_endpoints.csv")
)
source_data <- lapply(source_paths, fread)
for (nm in names(source_data)) {
  source_data[[nm]][, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
  source_data[[nm]] <- source_data[[nm]][
    tolower(as.character(in_background_for_test)) %in% c("true", "1")
  ]
  source_data[[nm]][, phenocode := as.character(phenocode)]
}

read_set <- function(set_name) {
  path <- file.path(frozen_dir, paste0(set_name, "_non_mhc.csv"))
  unique(toupper(trimws(as.character(fread(path)$gene_symbol))))
}

build_model_data <- function(source, phenocode_value, set_name) {
  resource <- if (source == "D2_TPMI") "TPMI" else "FinnGen"
  flag <- if (resource == "TPMI") "in_tpmi_gene_window_detectable" else "in_finngen_gene_window_detectable"
  snp_col <- if (resource == "TPMI") "tpmi_median_window_snps" else "finngen_median_window_snps"
  background <- copy(covariates[get(flag) == TRUE & is_extended_mhc == FALSE])
  endpoint_rows <- copy(source_data[[source]])[
    get("phenocode") == as.character(phenocode_value)
  ]
  endpoint_rows[, sidak_gene_p_num := as.numeric(sidak_gene_p)]
  setorder(endpoint_rows, sidak_gene_p_num)
  endpoint_rows <- unique(endpoint_rows, by = "gene_symbol")
  hit_map <- endpoint_rows[, .(gene_symbol, hit = sidak_gene_p_num < 0.05)]
  model <- merge(background, hit_map, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  model[is.na(hit), hit := FALSE]
  set_genes <- read_set(set_name)
  model[, in_set := as.integer(gene_symbol %in% set_genes)]
  model[, chr_factor := factor(chr)]
  model[, z_log_gene_length := scale_robust(log1p(as.numeric(gene_length_bp)))]
  model[, z_log_gene_density := scale_robust(log1p(as.numeric(gene_density_plusminus_1mb)))]
  model[, z_log_magma_nsnps := scale_robust(log1p(as.numeric(magma_nsnps_median)))]
  model[, z_effective_parameter_ratio := scale_robust(as.numeric(effective_parameter_ratio))]
  model[, z_log_resource_window_snps := scale_robust(log1p(as.numeric(get(snp_col))))]
  list(resource = resource, data = model)
}

fit_firth <- function(d) {
  d <- copy(d)
  d[, chr_factor := droplevels(factor(chr))]
  fit <- logistf(
    hit ~ in_set + z_log_gene_length + z_log_gene_density +
      z_log_magma_nsnps + z_effective_parameter_ratio +
      z_log_resource_window_snps + chr_factor,
    data = d, pl = TRUE, plconf = 2,
    control = logistf.control(maxit = 50, maxstep = 5),
    plcontrol = logistpl.control(maxit = 100, maxstep = 5)
  )
  nm <- if ("in_set" %in% names(fit$coefficients)) "in_set" else grep("^in_set", names(fit$coefficients), value = TRUE)
  if (length(nm) != 1L) stop("Frozen-set coefficient not uniquely identified.")
  beta <- unname(fit$coefficients[nm])
  low <- exp(unname(fit$ci.lower[nm]))
  high <- exp(unname(fit$ci.upper[nm]))
  p <- unname(fit$prob[nm])
  fit_conv <- max(abs(as.numeric(fit$conv)), na.rm = TRUE)
  profile_conv <- max(abs(as.numeric(fit$pl.conv)), na.rm = TRUE)
  if (anyNA(c(beta, low, high, p))) stop("Missing focal estimate after fitting.")
  data.table(
    log_or = beta, odds_ratio = exp(beta), ci95_low = low, ci95_high = high,
    p_value = p, fit_convergence_max = fit_conv,
    profile_convergence_max = profile_conv,
    converged_1e_5 = fit_conv <= 1e-5 & profile_conv <= 1e-5
  )
}

results <- list()
failures <- list()
counter <- 0L

for (i in seq_len(nrow(target_manifest))) {
  target <- target_manifest[i]
  endpoint_meta <- endpoints[
    source == target$source & as.character(phenocode) == as.character(target$phenocode)
  ]
  if (nrow(endpoint_meta) != 1L) stop("Endpoint manifest lookup was not unique for ", target$target_id)
  built <- build_model_data(target$source, target$phenocode, target$set_name)
  d <- built$data
  chromosomes <- sort(unique(as.character(d[in_set == 1L, chr])))
  deletion_units <- c("BASELINE", chromosomes)

  for (unit in deletion_units) {
    counter <- counter + 1L
    analysis_data <- if (unit == "BASELINE") d else d[as.character(chr) != unit]
    base_row <- data.table(
      target_id = target$target_id,
      source = target$source,
      phenocode = target$phenocode,
      phenotype = endpoint_meta$phenotype,
      set_name = target$set_name,
      set_version = "non_mhc",
      resource = built$resource,
      deletion_type = if (unit == "BASELINE") "baseline" else "leave_one_chromosome_out",
      omitted_chromosome = if (unit == "BASELINE") NA_character_ else unit,
      background_n = nrow(analysis_data),
      set_n = sum(analysis_data$in_set),
      set_hit_n = sum(analysis_data$in_set == 1L & analysis_data$hit),
      total_hit_n = sum(analysis_data$hit),
      removed_background_n = nrow(d) - nrow(analysis_data),
      removed_set_n = sum(d$in_set) - sum(analysis_data$in_set),
      removed_set_hit_n = sum(d$in_set == 1L & d$hit) - sum(analysis_data$in_set == 1L & analysis_data$hit)
    )
    fitted <- tryCatch(fit_firth(analysis_data), error = function(e) e)
    if (inherits(fitted, "error")) {
      failures[[length(failures) + 1L]] <- cbind(base_row, error_message = conditionMessage(fitted))
    } else {
      results[[length(results) + 1L]] <- cbind(base_row, fitted)
    }
    message(sprintf("[%d] %s | %s | omit %s", counter, target$target_id, target$set_name, unit))
  }
}

result <- rbindlist(results, fill = TRUE)
failure_table <- rbindlist(failures, fill = TRUE)
if (!ncol(failure_table)) {
  failure_table <- data.table(target_id = character(), source = character(), phenocode = character(),
                              set_name = character(), omitted_chromosome = character(), error_message = character())
}

baseline <- result[deletion_type == "baseline", .(
  target_id, baseline_log_or = log_or, baseline_or = odds_ratio,
  baseline_ci95_low = ci95_low, baseline_ci95_high = ci95_high, baseline_p = p_value
)]
result <- merge(result, baseline, by = "target_id", all.x = TRUE, sort = FALSE)
result[, relative_log_or_attenuation := (baseline_log_or - log_or) / pmax(abs(baseline_log_or), .Machine$double.eps)]
result[, direction_reversal_vs_baseline := sign(log_or) != sign(baseline_log_or)]

phase1b_reference <- phase1b[
  set_version == "non_mhc",
  .(source, phenocode = as.character(phenocode), set_name, phase1b_log_or = adjusted_log_or)
]
baseline_check <- merge(
  result[deletion_type == "baseline", .(target_id, source, phenocode = as.character(phenocode), set_name, baseline_log_or = log_or)],
  phase1b_reference,
  by = c("source", "phenocode", "set_name"), all.x = TRUE, sort = FALSE
)
baseline_check[, absolute_log_or_difference := abs(baseline_log_or - phase1b_log_or)]
baseline_check[, reproduces_phase1b_within_1e_5 := absolute_log_or_difference <= 1e-5]

loco <- result[deletion_type == "leave_one_chromosome_out"]
summary_table <- loco[, {
  worst <- which.max(relative_log_or_attenuation)
  .(
    chromosomes_tested = .N,
    baseline_or = baseline_or[1],
    minimum_loco_or = min(odds_ratio),
    maximum_loco_or = max(odds_ratio),
    maximum_relative_log_or_attenuation = max(relative_log_or_attenuation),
    most_influential_chromosome = omitted_chromosome[worst],
    any_direction_reversal = any(direction_reversal_vs_baseline),
    any_attenuation_ge_50_percent = any(relative_log_or_attenuation >= 0.5),
    all_loco_or_above_one = all(odds_ratio > 1),
    minimum_loco_ci95_low = min(ci95_low),
    interpretation = fifelse(
      any(direction_reversal_vs_baseline), "chromosome_sensitive_direction_reversal",
      fifelse(any(relative_log_or_attenuation >= 0.5),
              "chromosome_sensitive_ge50pct_logOR_attenuation",
              "no_single_chromosome_major_influence_detected")
    )
  )
}, by = .(target_id, source, phenocode, phenotype, set_name, set_version, resource)]

setorder(result, target_id, deletion_type, omitted_chromosome)
write_csv_utf8(result, file.path(out_dir, "phase1c_targeted_loco_all_models.csv"))
write_csv_utf8(summary_table, file.path(out_dir, "phase1c_targeted_loco_summary.csv"))
write_csv_utf8(baseline_check, file.path(out_dir, "phase1c_baseline_reproduction_check.csv"))
write_csv_utf8(failure_table, file.path(out_dir, "phase1c_model_failures.csv"))

input_paths <- c(covariates_path, endpoints_path, phase1b_path, source_paths,
                 file.path(frozen_dir, paste0(unique(target_manifest$set_name), "_non_mhc.csv")))
input_audit <- data.table(
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv_utf8(input_audit, file.path(out_dir, "phase1c_input_provenance_sha256.csv"))

valid_models <- sum(complete.cases(result[, .(log_or, odds_ratio, ci95_low, ci95_high, p_value)]))
all_reproduced <- nrow(baseline_check) == 10L && all(baseline_check$reproduces_phase1b_within_1e_5)
audit <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  role = "descriptive targeted chromosome influence analysis; not a new confirmatory test family",
  target_count = 10L,
  fitted_models = nrow(result),
  models_with_nonmissing_estimates = valid_models,
  failed_models = nrow(failure_table),
  baseline_models_reproducing_phase1b = sum(baseline_check$reproduces_phase1b_within_1e_5),
  selection_rule = "representative non-MHC-supported anchors; overlapping FinnGen dementia definitions not repeated",
  influence_flags = c("direction reversal", "at least 50% attenuation of baseline log OR"),
  multiple_testing = "none; deletion results are influence diagnostics and P values are descriptive",
  next_stage_rule = "leave-one-locus-out only within chromosomes flagged by direction reversal or >=50% log-OR attenuation"
)
write_json(audit, file.path(out_dir, "phase1c_targeted_loco_audit.json"), pretty = TRUE, auto_unbox = TRUE)

if (nrow(failure_table) > 0L || valid_models != nrow(result) || !all_reproduced ||
    !all(result$converged_1e_5)) {
  stop("Phase 1c-A failed QA. INCOMPLETE marker retained; inspect failure and baseline-check files.")
}

write_json(
  list(stage = "phase1c_targeted_loco", status = "COMPLETE",
       completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
       target_count = 10L, fitted_models = nrow(result), failed_models = 0L),
  file.path(out_dir, "_STAGE_COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)
unlink(file.path(out_dir, "_STAGE_INCOMPLETE.json"))
message("Track C Phase 1c-A completed: ", out_dir)


