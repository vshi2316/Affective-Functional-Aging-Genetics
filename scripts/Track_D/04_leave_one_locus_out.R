## ============================================================================
## Track C Phase 1c-B: targeted leave-one-locus-out influence analysis
## User-run only. Triggered strictly by Phase 1c-A LOCO results.
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
phase1c_loco_dir <- Sys.getenv(
  "TRACKC_PHASE1C_LOCO_DIR",
  unset = file.path(trackc_root, "04_targeted_loco")
)
out_dir <- Sys.getenv(
  "TRACKC_PHASE1C_LOLO_OUT_DIR",
  unset = file.path(trackc_root, "05_targeted_lolo")
)

base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
trackc_root <- normalizePath(trackc_root, winslash = "/", mustWork = TRUE)
phase1b_dir <- normalizePath(phase1b_dir, winslash = "/", mustWork = TRUE)
phase1c_loco_dir <- normalizePath(phase1c_loco_dir, winslash = "/", mustWork = TRUE)
trackb <- file.path(base_dir, "analysis_ready_core", "trackB")
frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
spec_dir <- file.path(trackc_root, "01_analysis_specification")

required_markers <- c(
  file.path(frozen_dir, "_STAGE_COMPLETE.json"),
  file.path(spec_dir, "_STAGE_COMPLETE.json"),
  file.path(phase1b_dir, "_STAGE_COMPLETE.json"),
  file.path(phase1c_loco_dir, "_STAGE_COMPLETE.json")
)
if (!all(file.exists(required_markers))) stop("Required upstream COMPLETE markers are missing.")
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Overwrite protection: output directory already contains files: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_json(
  list(stage = "phase1c_targeted_lolo", status = "INCOMPLETE_OR_RUNNING",
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

locus_gap_bp <- 1000000L
locus_deletion_flank_bp <- 1000000L

## These five target-chromosome pairs are exactly those triggered by Phase 1c-A.
trigger_manifest <- data.table(
  target_id = c("T01", "T04", "T05", "T06", "T08"),
  source = c("D2_TPMI", "D5_FinnGen_positive", "D2_TPMI", "D2_TPMI", "D5_FinnGen_positive"),
  phenocode = c("290.1", "M13_RHEUMA", "290.1", "290.11", "M13_SPONDYLINFLAMNAS"),
  set_name = c("cross_model_core", "affective_prioritized", "functional_prioritized",
               "functional_prioritized", "functional_prioritized"),
  trigger_chromosome = c("3", "11", "3", "3", "16"),
  trigger_reason = c("logOR attenuation >=50%", "logOR attenuation >=50%",
                     "direction reversal", "direction reversal", "logOR attenuation >=50%"),
  set_version = "non_mhc"
)

loco_summary_path <- file.path(phase1c_loco_dir, "phase1c_targeted_loco_summary.csv")
loco_summary <- fread(loco_summary_path, colClasses = list(character = c("phenocode", "most_influential_chromosome")))
derived_triggers <- loco_summary[
  any_direction_reversal == TRUE | any_attenuation_ge_50_percent == TRUE,
  .(target_id, source, phenocode = as.character(phenocode), set_name,
    trigger_chromosome = as.character(most_influential_chromosome))
]
expected_triggers <- trigger_manifest[, .(target_id, source, phenocode, set_name, trigger_chromosome)]
setorder(derived_triggers, target_id)
setorder(expected_triggers, target_id)
if (!fsetequal(derived_triggers, expected_triggers)) {
  stop("Frozen Phase 1c-B trigger manifest does not exactly reproduce Phase 1c-A trigger results.")
}
write_csv_utf8(trigger_manifest, file.path(out_dir, "phase1c_lolo_frozen_trigger_manifest.csv"))

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

build_model_data <- function(source_name, phenocode_value, set_name_value) {
  resource <- if (source_name == "D2_TPMI") "TPMI" else "FinnGen"
  flag <- if (resource == "TPMI") "in_tpmi_gene_window_detectable" else "in_finngen_gene_window_detectable"
  snp_col <- if (resource == "TPMI") "tpmi_median_window_snps" else "finngen_median_window_snps"
  background <- copy(covariates[get(flag) == TRUE & is_extended_mhc == FALSE])
  endpoint_rows <- copy(source_data[[source_name]])[
    get("phenocode") == as.character(phenocode_value)
  ]
  endpoint_rows[, sidak_gene_p_num := as.numeric(sidak_gene_p)]
  setorder(endpoint_rows, sidak_gene_p_num)
  endpoint_rows <- unique(endpoint_rows, by = "gene_symbol")
  hit_map <- endpoint_rows[, .(gene_symbol, sidak_gene_p = sidak_gene_p_num,
                               hit = sidak_gene_p_num < 0.05)]
  model <- merge(background, hit_map, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  model[is.na(hit), hit := FALSE]
  set_genes <- read_set(set_name_value)
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

make_loci <- function(d, chromosome_value, target_id_value) {
  genes <- copy(d[in_set == 1L & as.character(chr) == as.character(chromosome_value),
                  .(gene_symbol, chr = as.character(chr), start_grch37, stop_grch37,
                    center, hit, sidak_gene_p)])
  setorder(genes, start_grch37, stop_grch37, gene_symbol)
  if (!nrow(genes)) stop("No target-set genes on triggered chromosome for ", target_id_value)
  previous_running_end <- shift(cummax(genes$stop_grch37), fill = genes$start_grch37[1])
  genes[, new_locus := .I == 1L | (start_grch37 - previous_running_end) > locus_gap_bp]
  genes[, locus_number := cumsum(new_locus)]
  loci <- genes[, .(
    chromosome = first(chr),
    locus_start = min(start_grch37),
    locus_stop = max(stop_grch37),
    deletion_start = max(0, min(start_grch37) - locus_deletion_flank_bp),
    deletion_stop = max(stop_grch37) + locus_deletion_flank_bp,
    set_gene_n = .N,
    set_hit_gene_n = sum(hit),
    set_genes = paste(gene_symbol, collapse = ";"),
    set_hit_genes = paste(gene_symbol[hit], collapse = ";")
  ), by = locus_number]
  loci[, `:=`(
    target_id = target_id_value,
    locus_id = sprintf("%s_chr%s_L%02d", target_id_value, chromosome, locus_number)
  )]
  loci[]
}

results <- list()
failures <- list()
locus_manifests <- list()
counter <- 0L

for (i in seq_len(nrow(trigger_manifest))) {
  target <- trigger_manifest[i]
  endpoint_meta <- endpoints[
    source == target$source & as.character(phenocode) == as.character(target$phenocode)
  ]
  if (nrow(endpoint_meta) != 1L) stop("Endpoint manifest lookup was not unique for ", target$target_id)
  built <- build_model_data(target$source, target$phenocode, target$set_name)
  d <- built$data
  loci <- make_loci(d, target$trigger_chromosome, target$target_id)
  locus_output <- loci[, .(locus_id, locus_number, chromosome, locus_start, locus_stop,
                           deletion_start, deletion_stop, set_gene_n, set_hit_gene_n,
                           set_genes, set_hit_genes)]
  locus_output[, `:=`(
    target_id = target$target_id,
    source = target$source,
    phenocode = target$phenocode,
    set_name = target$set_name,
    trigger_chromosome = target$trigger_chromosome
  )]
  setcolorder(locus_output, c("target_id", "source", "phenocode", "set_name",
                             "trigger_chromosome", setdiff(names(locus_output),
                             c("target_id", "source", "phenocode", "set_name", "trigger_chromosome"))))
  locus_manifests[[length(locus_manifests) + 1L]] <- locus_output

  deletion_units <- c("BASELINE", loci$locus_id)
  for (unit in deletion_units) {
    counter <- counter + 1L
    if (unit == "BASELINE") {
      locus <- NULL
      analysis_data <- d
    } else {
      locus <- loci[locus_id == unit]
      remove_row <- as.character(d$chr) == locus$chromosome &
        d$center >= locus$deletion_start & d$center <= locus$deletion_stop
      analysis_data <- d[!remove_row]
    }
    base_row <- data.table(
      target_id = target$target_id, source = target$source,
      phenocode = target$phenocode, phenotype = endpoint_meta$phenotype,
      set_name = target$set_name, set_version = "non_mhc", resource = built$resource,
      trigger_chromosome = target$trigger_chromosome,
      deletion_type = if (unit == "BASELINE") "baseline" else "leave_one_locus_out",
      omitted_locus = if (unit == "BASELINE") NA_character_ else unit,
      deletion_start = if (unit == "BASELINE") NA_real_ else locus$deletion_start,
      deletion_stop = if (unit == "BASELINE") NA_real_ else locus$deletion_stop,
      background_n = nrow(analysis_data), set_n = sum(analysis_data$in_set),
      set_hit_n = sum(analysis_data$in_set == 1L & analysis_data$hit),
      total_hit_n = sum(analysis_data$hit),
      removed_background_n = nrow(d) - nrow(analysis_data),
      removed_set_n = sum(d$in_set) - sum(analysis_data$in_set),
      removed_set_hit_n = sum(d$in_set == 1L & d$hit) -
        sum(analysis_data$in_set == 1L & analysis_data$hit)
    )
    fitted <- tryCatch(fit_firth(analysis_data), error = function(e) e)
    if (inherits(fitted, "error")) {
      failures[[length(failures) + 1L]] <- cbind(base_row, error_message = conditionMessage(fitted))
    } else {
      results[[length(results) + 1L]] <- cbind(base_row, fitted)
    }
    message(sprintf("[%d] %s | omit %s", counter, target$target_id, unit))
  }
}

locus_manifest <- rbindlist(locus_manifests, fill = TRUE)
result <- rbindlist(results, fill = TRUE)
failure_table <- rbindlist(failures, fill = TRUE)
if (!ncol(failure_table)) {
  failure_table <- data.table(target_id = character(), source = character(), phenocode = character(),
                              set_name = character(), omitted_locus = character(), error_message = character())
}

baseline <- result[deletion_type == "baseline", .(
  target_id, baseline_log_or = log_or, baseline_or = odds_ratio,
  baseline_ci95_low = ci95_low, baseline_ci95_high = ci95_high, baseline_p = p_value
)]
result <- merge(result, baseline, by = "target_id", all.x = TRUE, sort = FALSE)
result[, relative_log_or_attenuation := (baseline_log_or - log_or) /
         pmax(abs(baseline_log_or), .Machine$double.eps)]
result[, direction_reversal_vs_baseline := sign(log_or) != sign(baseline_log_or)]

phase1b_reference <- phase1b[
  set_version == "non_mhc",
  .(source, phenocode = as.character(phenocode), set_name, phase1b_log_or = adjusted_log_or)
]
baseline_check <- merge(
  result[deletion_type == "baseline",
         .(target_id, source, phenocode = as.character(phenocode), set_name, baseline_log_or = log_or)],
  phase1b_reference, by = c("source", "phenocode", "set_name"), all.x = TRUE, sort = FALSE
)
baseline_check[, absolute_log_or_difference := abs(baseline_log_or - phase1b_log_or)]
baseline_check[, reproduces_phase1b_within_1e_5 := absolute_log_or_difference <= 1e-5]

lolo <- result[deletion_type == "leave_one_locus_out"]
summary_table <- lolo[, {
  worst <- which.max(relative_log_or_attenuation)
  .(
    loci_tested = .N,
    baseline_or = baseline_or[1],
    minimum_lolo_or = min(odds_ratio),
    maximum_lolo_or = max(odds_ratio),
    maximum_relative_log_or_attenuation = max(relative_log_or_attenuation),
    most_influential_locus = omitted_locus[worst],
    most_influential_locus_removed_set_n = removed_set_n[worst],
    most_influential_locus_removed_set_hit_n = removed_set_hit_n[worst],
    any_direction_reversal = any(direction_reversal_vs_baseline),
    any_attenuation_ge_50_percent = any(relative_log_or_attenuation >= 0.5),
    all_lolo_or_above_one = all(odds_ratio > 1),
    interpretation = fifelse(
      any(direction_reversal_vs_baseline), "single_locus_sensitive_direction_reversal",
      fifelse(any(relative_log_or_attenuation >= 0.5),
              "single_locus_sensitive_ge50pct_logOR_attenuation",
              "no_single_locus_major_influence_detected")
    )
  )
}, by = .(target_id, source, phenocode, phenotype, set_name, set_version,
          resource, trigger_chromosome)]

setorder(result, target_id, deletion_type, omitted_locus)
write_csv_utf8(locus_manifest, file.path(out_dir, "phase1c_lolo_locus_manifest.csv"))
write_csv_utf8(result, file.path(out_dir, "phase1c_targeted_lolo_all_models.csv"))
write_csv_utf8(summary_table, file.path(out_dir, "phase1c_targeted_lolo_summary.csv"))
write_csv_utf8(baseline_check, file.path(out_dir, "phase1c_lolo_baseline_reproduction_check.csv"))
write_csv_utf8(failure_table, file.path(out_dir, "phase1c_lolo_model_failures.csv"))

input_paths <- c(covariates_path, endpoints_path, phase1b_path, loco_summary_path, source_paths,
                 file.path(frozen_dir, paste0(unique(trigger_manifest$set_name), "_non_mhc.csv")))
input_audit <- data.table(
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv_utf8(input_audit, file.path(out_dir, "phase1c_lolo_input_provenance_sha256.csv"))

valid_models <- sum(complete.cases(result[, .(log_or, odds_ratio, ci95_low, ci95_high, p_value)]))
all_reproduced <- nrow(baseline_check) == 5L && all(baseline_check$reproduces_phase1b_within_1e_5)
audit <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  role = "descriptive targeted locus influence analysis; not a new confirmatory test family",
  triggered_target_count = 5L,
  locus_gap_bp = locus_gap_bp,
  deletion_flank_bp = locus_deletion_flank_bp,
  locus_definition = "target-set genes on the triggered chromosome clustered when successive intervals are <=1 Mb apart",
  deletion_rule = "remove all background genes with centers inside the target-locus span plus/minus 1 Mb",
  fitted_models = nrow(result), models_with_nonmissing_estimates = valid_models,
  failed_models = nrow(failure_table),
  baseline_models_reproducing_phase1b = sum(baseline_check$reproduces_phase1b_within_1e_5),
  influence_flags = c("direction reversal", "at least 50% attenuation of baseline log OR"),
  multiple_testing = "none; deletion results are descriptive influence diagnostics"
)
write_json(audit, file.path(out_dir, "phase1c_targeted_lolo_audit.json"), pretty = TRUE, auto_unbox = TRUE)

if (nrow(failure_table) > 0L || valid_models != nrow(result) || !all_reproduced ||
    !all(result$converged_1e_5)) {
  stop("Phase 1c-B failed QA. INCOMPLETE marker retained; inspect output diagnostics.")
}

write_json(
  list(stage = "phase1c_targeted_lolo", status = "COMPLETE",
       completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
       triggered_target_count = 5L, loci_tested = nrow(locus_manifest),
       fitted_models = nrow(result), failed_models = 0L),
  file.path(out_dir, "_STAGE_COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)
unlink(file.path(out_dir, "_STAGE_INCOMPLETE.json"))
message("Track C Phase 1c-B completed: ", out_dir)


