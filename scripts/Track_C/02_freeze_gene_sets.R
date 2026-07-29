## ============================================================================
## Freeze Track C gene sets before re-analysis of external resources
##
## Primary sets
##   1. Cross-model core
##   2. Affective-prioritized
##   3. Functional-prioritized
##
## Each set is exported as a full and non-MHC version. The script reads only
## FUMA/MAGMA annotation outputs and pre-existing input manifests. It does not
## read external validation results.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

resolve_project_base <- function() {
  candidates <- unique(c(
    Sys.getenv("DCV_BASE_DIR", unset = ""),
    Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
  ))
  candidates <- candidates[nzchar(candidates)]
  for (x in candidates) {
    if (dir.exists(file.path(x, "analysis_ready_core", "trackB"))) {
      return(normalizePath(x, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Could not resolve project base. Set DCV_BASE_DIR to the project root.")
}

as_flag <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

safe_fread <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Missing required input: ", path)
    return(data.table())
  }
  fread(path, data.table = TRUE, showProgress = FALSE)
}

file_hash <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(list(algorithm = "SHA256", value = digest::digest(file = path, algo = "sha256")))
  }
  list(algorithm = "MD5", value = unname(tools::md5sum(path)))
}

hash_manifest <- function(paths, role) {
  stopifnot(length(paths) == length(role))
  rbindlist(lapply(seq_along(paths), function(i) {
    p <- normalizePath(paths[[i]], winslash = "/", mustWork = TRUE)
    h <- file_hash(p)
    info <- file.info(p)
    data.table(
      role = role[[i]],
      path = p,
      bytes = as.numeric(info$size),
      file_mtime = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
      hash_algorithm = h$algorithm,
      hash = h$value
    )
  }))
}

clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x) | x %in% c("", "NA", ".")] <- NA_character_
  x
}

collapse_magma <- function(path, domain) {
  x <- safe_fread(path)
  required <- c("GENE", "SYMBOL", "CHR", "START", "STOP", "P")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(domain, " MAGMA table lacks columns: ", paste(missing, collapse = ", "))
  }

  x[, `:=`(
    GENE = as.character(GENE),
    symbol = clean_symbol(SYMBOL),
    CHR = suppressWarnings(as.integer(CHR)),
    START = suppressWarnings(as.numeric(START)),
    STOP = suppressWarnings(as.numeric(STOP)),
    P = suppressWarnings(as.numeric(P))
  )]
  x <- x[is.finite(P) & P > 0 & P <= 1]
  n_tested_records <- nrow(x)
  threshold <- 0.05 / n_tested_records
  x[, bonferroni_significant := P < threshold]

  missing_sig_symbol <- x[bonferroni_significant == TRUE & is.na(symbol)]
  if (nrow(missing_sig_symbol)) {
    stop(domain, " has Bonferroni-significant Ensembl records without a gene symbol.")
  }

  setorder(x, symbol, P, GENE)
  representative <- x[!is.na(symbol), .SD[1], by = symbol]
  aggregate <- x[!is.na(symbol), .(
    ensembl_ids = paste(sort(unique(GENE)), collapse = ";"),
    n_ensembl_records = uniqueN(GENE),
    n_magma_records = .N,
    min_p = min(P),
    bonferroni_significant = any(bonferroni_significant),
    significant_ensembl_ids = paste(sort(unique(GENE[bonferroni_significant])), collapse = ";")
  ), by = symbol]

  representative <- representative[, .(
    symbol,
    representative_ensembl_id = GENE,
    chr = CHR,
    start = START,
    stop = STOP,
    representative_p = P
  )]
  collapsed <- merge(aggregate, representative, by = "symbol", all.x = TRUE)

  duplicate_audit <- aggregate[n_ensembl_records > 1 | n_magma_records > 1]
  duplicate_audit[, domain := domain]

  list(
    raw = x,
    collapsed = collapsed,
    duplicate_audit = duplicate_audit,
    n_tested_records = n_tested_records,
    n_tested_symbols = uniqueN(x$symbol, na.rm = TRUE),
    threshold = threshold,
    n_significant_records = x[bonferroni_significant == TRUE, .N],
    n_significant_symbols = collapsed[bonferroni_significant == TRUE, .N]
  )
}

read_positional_symbols <- function(path, domain) {
  x <- safe_fread(path)
  if (!"symbol" %in% names(x)) stop(domain, " positional table lacks symbol column.")
  x[, symbol := clean_symbol(symbol)]
  unique(x[!is.na(symbol), .(symbol)])
}

base_dir <- resolve_project_base()
trackb_dir <- file.path(base_dir, "analysis_ready_core", "trackB")
fuma_summary_dir <- file.path(base_dir, "FUMA", "FUMA_JOB06_07_twofactor_visualization")
fuma_table_dir <- file.path(fuma_summary_dir, "tables")

out_dir <- file.path(trackb_dir, "trackC_twofactor_frozen_gene_sets")
table_dir <- file.path(out_dir, "tables")
protocol_dir <- file.path(out_dir, "validation_protocol")
note_dir <- file.path(out_dir, "notes")
lock_file <- file.path(out_dir, "FREEZE_LOCK.txt")

allow_rebuild <- as_flag(Sys.getenv("TRACKC_ALLOW_REBUILD", unset = "FALSE"))
if (file.exists(lock_file) && !allow_rebuild) {
  stop(
    "A Track C freeze lock already exists: ", lock_file,
    "\nThe frozen gene sets cannot be rebuilt after external results have been reviewed."
  )
}
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) && !allow_rebuild) {
  stop("Output directory already contains files but has no accepted rebuild flag: ", out_dir)
}

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(protocol_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(note_dir, recursive = TRUE, showWarnings = FALSE)

domains <- c("common", "affective", "functional")
magma_paths <- setNames(
  file.path(fuma_table_dir, paste0("fuma_", domains, "_magma_gene_results.csv")),
  domains
)
positional_paths <- setNames(
  file.path(fuma_table_dir, paste0("fuma_", domains, "_positionally_mapped_genes.csv")),
  domains
)
parameter_audit_path <- file.path(fuma_table_dir, "fuma_parameter_audit.csv")
parameter_consistency_path <- file.path(fuma_table_dir, "fuma_parameter_consistency.csv")

required_inputs <- c(magma_paths, positional_paths, parameter_audit_path, parameter_consistency_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing FUMA inputs:\n", paste(missing_inputs, collapse = "\n"))

parameter_consistency <- safe_fread(parameter_consistency_path)
if (!all(as_flag(parameter_consistency$consistent))) {
  stop("FUMA parameters are not consistent across common, affective, and functional tasks.")
}

parameter_audit <- safe_fread(parameter_audit_path)
expected_numeric_parameters <- c(
  genome_build_grch38 = 0,
  lead_p = 5e-8,
  candidate_p = 0.05,
  independent = 0.6,
  lead = 0.1,
  merge_distance_kb = 250,
  positional_mapping = 1,
  magma = 1
)
for (nm in names(expected_numeric_parameters)) {
  observed <- unique(suppressWarnings(as.numeric(parameter_audit[[nm]])))
  expected <- unname(expected_numeric_parameters[[nm]])
  if (length(observed) != 1L || !is.finite(observed) ||
      !isTRUE(all.equal(observed, expected, tolerance = 1e-12))) {
    stop("Unexpected FUMA parameter ", nm, ": ", paste(observed, collapse = " | "))
  }
}

expected_character_parameters <- c(reference_panel = "1KG/Phase3", ancestry = "EUR")
for (nm in names(expected_character_parameters)) {
  observed <- unique(as.character(parameter_audit[[nm]]))
  expected <- unname(expected_character_parameters[[nm]])
  if (length(observed) != 1L || observed != expected) {
    stop("Unexpected FUMA parameter ", nm, ": ", paste(observed, collapse = " | "))
  }
}

input_manifest <- hash_manifest(
  paths = required_inputs,
  role = c(
    paste0(domains, "_magma_gene_results"),
    paste0(domains, "_positionally_mapped_genes"),
    "fuma_parameter_audit",
    "fuma_parameter_consistency"
  )
)
fwrite(input_manifest, file.path(table_dir, "trackC_freeze_input_hash_manifest.csv"))
fwrite(parameter_audit, file.path(table_dir, "trackC_fuma_parameter_snapshot.csv"))

magma <- lapply(domains, function(d) collapse_magma(magma_paths[[d]], d))
names(magma) <- domains

n_tested <- vapply(magma, `[[`, integer(1), "n_tested_records")
thresholds <- vapply(magma, `[[`, numeric(1), "threshold")
if (length(unique(n_tested)) != 1L || unique(n_tested) != 18516L) {
  stop("Expected 18,516 MAGMA-tested records in each model; observed: ", paste(n_tested, collapse = ", "))
}
if (length(unique(signif(thresholds, 12))) != 1L) {
  stop("MAGMA Bonferroni thresholds differ across models.")
}

ensembl_sets <- lapply(magma, function(z) sort(unique(z$raw$GENE)))
if (!identical(ensembl_sets$common, ensembl_sets$affective) ||
    !identical(ensembl_sets$common, ensembl_sets$functional)) {
  stop("The three MAGMA analyses did not test identical Ensembl gene records.")
}

rename_domain <- function(x, domain) {
  x <- copy(x)
  keep <- setdiff(names(x), "symbol")
  setnames(x, keep, paste0(keep, "_", domain))
  x
}

master <- Reduce(
  function(x, y) merge(x, y, by = "symbol", all = TRUE),
  lapply(domains, function(d) rename_domain(magma[[d]]$collapsed, d))
)

for (d in domains) {
  sig_col <- paste0("bonferroni_significant_", d)
  master[is.na(get(sig_col)), (sig_col) := FALSE]
}

master[, `:=`(
  chr = fcoalesce(chr_common, chr_affective, chr_functional),
  start = fcoalesce(start_common, start_affective, start_functional),
  stop = fcoalesce(stop_common, stop_affective, stop_functional)
)]
master[, is_mhc_grch37 := chr == 6L & start <= 34000000 & stop >= 25000000]
master[is.na(is_mhc_grch37), is_mhc_grch37 := FALSE]

positional <- lapply(domains, function(d) read_positional_symbols(positional_paths[[d]], d))
names(positional) <- domains
for (d in domains) {
  master[, (paste0("positionally_mapped_", d)) := symbol %in% positional[[d]]$symbol]
}

master[, `:=`(
  cross_model_core = bonferroni_significant_common &
    bonferroni_significant_affective & bonferroni_significant_functional,
  affective_prioritized = bonferroni_significant_affective &
    !bonferroni_significant_functional,
  functional_prioritized = bonferroni_significant_functional &
    !bonferroni_significant_affective,
  common_projection_significant = bonferroni_significant_common
)]

master[, primary_set := fifelse(
  cross_model_core, "cross_model_core",
  fifelse(
    affective_prioritized, "affective_prioritized",
    fifelse(functional_prioritized, "functional_prioritized", "not_in_primary_sets")
  )
)]
master[, `:=`(
  freeze_version = "TrackC",
  genome_build = "GRCh37",
  mhc_rule = "chr6:25000000-34000000; overlap with gene interval",
  membership_interpretation = fifelse(
    primary_set == "cross_model_core",
    "Within-study convergence across correlated model parameterizations; not independent replication",
    fifelse(
      primary_set %in% c("affective_prioritized", "functional_prioritized"),
      "Prioritized by the pattern of MAGMA significance; not a formal cross-factor effect difference",
      "Reference universe gene"
    )
  )
)]

setorder(master, primary_set, symbol)
fwrite(master, file.path(table_dir, "trackC_master_gene_membership_matrix.csv"))

record_summary <- rbindlist(lapply(domains, function(d) {
  z <- magma[[d]]
  data.table(
    model = d,
    tested_magma_records = z$n_tested_records,
    tested_unique_symbols = z$n_tested_symbols,
    bonferroni_threshold = z$threshold,
    significant_magma_records = z$n_significant_records,
    significant_unique_symbols = z$n_significant_symbols
  )
}))
fwrite(record_summary, file.path(table_dir, "trackC_magma_record_symbol_reconciliation.csv"))

duplicate_audit <- rbindlist(lapply(magma, `[[`, "duplicate_audit"), fill = TRUE)
setcolorder(duplicate_audit, c("domain", setdiff(names(duplicate_audit), "domain")))
fwrite(duplicate_audit, file.path(table_dir, "trackC_ensembl_duplicate_symbol_audit.csv"))

significant_record_audit <- rbindlist(lapply(domains, function(d) {
  x <- copy(magma[[d]]$raw[bonferroni_significant == TRUE])
  x[, domain := d]
  x
}), fill = TRUE)
fwrite(significant_record_audit, file.path(table_dir, "trackC_all_significant_ensembl_records.csv"))

set_definitions <- data.table(
  set_id = c("cross_model_core", "affective_prioritized", "functional_prioritized"),
  operational_definition = c(
    "Bonferroni significant in common, affective, and functional MAGMA analyses",
    "Bonferroni significant in affective MAGMA and not Bonferroni significant in functional MAGMA",
    "Bonferroni significant in functional MAGMA and not Bonferroni significant in affective MAGMA"
  ),
  prohibited_interpretation = c(
    "Independent replication",
    "Affective-specific effect or causal gene",
    "Functional-specific effect or causal gene"
  )
)
fwrite(set_definitions, file.path(table_dir, "trackC_operational_set_definitions.csv"))

set_specs <- list(
  cross_model_core = "cross_model_core",
  affective_prioritized = "affective_prioritized",
  functional_prioritized = "functional_prioritized"
)

export_columns <- c(
  "symbol", "primary_set", "chr", "start", "stop", "is_mhc_grch37",
  "representative_ensembl_id_common", "representative_ensembl_id_affective",
  "representative_ensembl_id_functional", "ensembl_ids_common", "ensembl_ids_affective",
  "ensembl_ids_functional", "min_p_common", "min_p_affective", "min_p_functional",
  "bonferroni_significant_common", "bonferroni_significant_affective",
  "bonferroni_significant_functional", "positionally_mapped_common",
  "positionally_mapped_affective", "positionally_mapped_functional",
  "freeze_version", "genome_build", "mhc_rule", "membership_interpretation"
)

set_manifest <- rbindlist(lapply(names(set_specs), function(set_name) {
  flag <- set_specs[[set_name]]
  full <- copy(master[get(flag) == TRUE, ..export_columns])
  non_mhc <- full[is_mhc_grch37 == FALSE]
  fwrite(full, file.path(table_dir, paste0("trackC_gene_set_", set_name, "_full.csv")))
  fwrite(non_mhc, file.path(table_dir, paste0("trackC_gene_set_", set_name, "_non_mhc.csv")))
  data.table(
    set_id = set_name,
    full_n = nrow(full),
    non_mhc_n = nrow(non_mhc),
    mhc_n = sum(full$is_mhc_grch37),
    full_members = paste(full$symbol, collapse = ";"),
    non_mhc_members = paste(non_mhc$symbol, collapse = ";")
  )
}))

fwrite(set_manifest, file.path(table_dir, "trackC_frozen_gene_set_manifest.csv"))

expected_counts <- c(
  cross_model_core = 17L,
  affective_prioritized = 41L,
  functional_prioritized = 259L
)
assert_expected <- as_flag(Sys.getenv("TRACKC_ASSERT_EXPECTED_COUNTS", unset = "TRUE"))
observed_counts <- setNames(set_manifest$full_n, set_manifest$set_id)
if (assert_expected && !identical(as.integer(observed_counts[names(expected_counts)]), as.integer(expected_counts))) {
  stop(
    "Frozen set counts do not match the audited expectations. Expected ",
    paste(names(expected_counts), expected_counts, sep = "=", collapse = ", "),
    "; observed ",
    paste(names(observed_counts), observed_counts, sep = "=", collapse = ", ")
  )
}

## Snapshot external endpoint/resource manifests without reading result tables.
protocol_sources <- c(
  c3_resource_manifest = file.path(
    trackb_dir, "trackC3_triangular_validation_immune_spatial", "tables",
    "trackC3_external_resource_manifest.csv"
  ),
  d2_tpmi_endpoint_manifest = file.path(
    trackb_dir, "trackD2_tpmi_transancestry_validation", "tables",
    "trackD2_tpmi_input_manifest.csv"
  ),
  d5_finngen_endpoint_manifest = file.path(
    trackb_dir, "trackD5_finngen_clinical_endpoint_replication", "tables",
    "trackD5_finngen_input_manifest.csv"
  ),
  d6_finngen_negative_control_manifest = file.path(
    trackb_dir, "trackD6_finngen_endpoint_specificity", "tables",
    "trackD6_finngen_negative_control_input_manifest.csv"
  )
)

protocol_snapshot <- rbindlist(lapply(names(protocol_sources), function(role) {
  src <- protocol_sources[[role]]
  if (!file.exists(src)) {
    return(data.table(
      role = role, source_path = src, snapshot_path = NA_character_, exists = FALSE,
      hash_algorithm = NA_character_, hash = NA_character_
    ))
  }
  dst <- file.path(protocol_dir, paste0(role, "_frozen.csv"))
  if (!file.copy(src, dst, overwrite = TRUE)) stop("Failed to snapshot protocol input: ", src)
  h <- file_hash(src)
  data.table(
    role = role,
    source_path = normalizePath(src, winslash = "/", mustWork = TRUE),
    snapshot_path = normalizePath(dst, winslash = "/", mustWork = TRUE),
    exists = TRUE,
    hash_algorithm = h$algorithm,
    hash = h$value
  )
}))
fwrite(protocol_snapshot, file.path(protocol_dir, "trackC_external_manifest_snapshot_audit.csv"))

test_families <- data.table(
  family_id = c(
    "C3_single_cell", "C3_spatial_proteome", "D2_TPMI_endpoints",
    "D5_D6_FinnGen_endpoints", "candidate_eQTL", "candidate_pQTL"
  ),
  primary_correction = "Benjamini-Hochberg within the pre-specified family",
  global_sensitivity = "Benjamini-Hochberg across all primary external-validation tests",
  resource_specific_background = c(
    "MAGMA-tested genes intersected with genes entering the single-cell expression test",
    "MAGMA-tested genes intersected with proteins detectable in the spatial resource",
    "MAGMA-tested genes with valid GRCh38 windows and adequate TPMI SNP coverage",
    "MAGMA-tested genes with valid GRCh38 windows and adequate FinnGen SNP coverage",
    "Tested candidate genes intersected with available genes in each tissue and eQTL resource",
    "Tested candidate genes intersected with available proteins in each pQTL resource"
  ),
  matched_permutation_variables = c(
    "gene-set size; mean expression; detection rate",
    "gene-set size; protein measurability",
    "gene-set size; gene length; effective window SNP count; chromosome or LD stratum",
    "gene-set size; gene length; effective window SNP count; chromosome or LD stratum",
    "not applicable to candidate-level association tests",
    "not applicable to candidate-level association tests"
  ),
  permitted_label = c(
    "external orthogonal support", "external orthogonal support",
    "cross-resource endpoint reassessment", "cross-resource endpoint reassessment",
    "candidate-driven molecular follow-up", "candidate-driven molecular follow-up"
  )
)
fwrite(test_families, file.path(protocol_dir, "trackC_prespecified_test_families.csv"))

qtl_candidates <- data.table(
  gene_symbol = c("PHF2", "GPX1", "TCF4"),
  candidate_role = c(
    "affective-prioritized candidate",
    "functional-prioritized candidate",
    "cross-model core candidate"
  ),
  primary_factor_gwas = c("affective", "functional", "affective and functional in parallel"),
  analysis_label = "candidate-driven molecular follow-up",
  multiplicity_scope = c(
    "all PHF2 tissue and QTL-resource combinations",
    "all GPX1 tissue and QTL-resource combinations",
    "all TCF4 factor, tissue, and QTL-resource combinations"
  ),
  causal_claim_permitted = FALSE
)
fwrite(qtl_candidates, file.path(protocol_dir, "trackC_candidate_qtl_followup_plan.csv"))

output_files_before_lock <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
output_files_before_lock <- output_files_before_lock[file.info(output_files_before_lock)$isdir == FALSE]
output_manifest <- hash_manifest(
  output_files_before_lock,
  paste0("frozen_output/", substring(normalizePath(output_files_before_lock, winslash = "/"), nchar(normalizePath(out_dir, winslash = "/")) + 2L))
)
fwrite(output_manifest, file.path(table_dir, "trackC_frozen_output_hash_manifest.csv"))

freeze_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
lock_lines <- c(
  "Track C frozen gene-set lock",
  paste0("Frozen at: ", freeze_time),
  paste0("Project base: ", base_dir),
  paste0("MAGMA records tested per model: ", unique(n_tested)),
  paste0("MAGMA Bonferroni threshold: ", format(unique(thresholds), scientific = TRUE, digits = 12)),
  "MHC rule: GRCh37 chr6:25,000,000-34,000,000, based on gene-interval overlap",
  paste0("Cross-model core full/non-MHC: ", set_manifest[set_id == "cross_model_core", full_n], "/",
         set_manifest[set_id == "cross_model_core", non_mhc_n]),
  paste0("Affective-prioritized full/non-MHC: ", set_manifest[set_id == "affective_prioritized", full_n], "/",
         set_manifest[set_id == "affective_prioritized", non_mhc_n]),
  paste0("Functional-prioritized full/non-MHC: ", set_manifest[set_id == "functional_prioritized", full_n], "/",
         set_manifest[set_id == "functional_prioritized", non_mhc_n]),
  "The prioritized labels encode MAGMA significance patterns and do not establish cross-factor effect differences.",
  "External results must not be used to add or remove genes from these sets."
)
writeLines(lock_lines, lock_file, useBytes = TRUE)

note <- c(
  "# Track C frozen gene sets",
  "",
  "The primary sets were defined before re-analysis of external validation resources.",
  "Cross-model core denotes within-study convergence across correlated model parameterizations.",
  "Affective-prioritized and functional-prioritized denote patterns of MAGMA significance.",
  "No set is labelled domain-specific, causal, or independently replicated.",
  "",
  "Duplicate-symbol rule:",
  "A unique symbol enters a set when any corresponding Ensembl MAGMA record passes the model-specific Bonferroni threshold.",
  "All Ensembl records are retained in the audit table; the minimum P value is used for symbol-level display.",
  "",
  "External validation:",
  "Resource-specific backgrounds are intersections between the MAGMA-tested universe and genes measurable in each resource.",
  "Primary multiplicity correction is pre-specified within scientific test families.",
  "A global BH correction across primary validation tests is retained as a sensitivity analysis."
)
writeLines(note, file.path(note_dir, "trackC_freeze_method_note.md"), useBytes = TRUE)

message("Track C freeze completed: ", out_dir)
message("Freeze lock: ", lock_file)




