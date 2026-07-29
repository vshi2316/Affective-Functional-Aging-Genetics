## ============================================================================
## Track C phase 1
## User-run R pipeline: JOB01/JOB06/JOB07 freeze, analysis lock,
## and TPMI/FinnGen endpoint reassessment.
##
## Run with:
## Sys.setenv(
##   DCV_BASE_DIR = "/path/to/project-data",
##   TRACKC_ROOT = "/path/to/project-data/analysis_ready_core/trackC",
##   TRACKC_STAGES = "freeze,spec,endpoints"
## )
## source("run_trackC_phase1.R", encoding = "UTF-8")
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(jsonlite)
  library(digest)
})

options(stringsAsFactors = FALSE)

resolve_project_base <- function() {
  candidates <- unique(c(
    Sys.getenv("DCV_BASE_DIR", unset = ""),
    Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
  ))
  candidates <- candidates[nzchar(candidates)]
  for (x in candidates) {
    if (dir.exists(file.path(x, "analysis_ready_core"))) {
      return(normalizePath(x, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot resolve project base. Set DCV_BASE_DIR before source().")
}

base_dir <- resolve_project_base()
core_dir <- file.path(base_dir, "analysis_ready_core")
trackb <- file.path(core_dir, "trackB")
trackc_root <- Sys.getenv(
  "TRACKC_ROOT",
  unset = file.path(core_dir, "trackC_formal_run")
)
trackc_root <- normalizePath(trackc_root, winslash = "/", mustWork = FALSE)

requested_stages <- strsplit(
  Sys.getenv("TRACKC_STAGES", unset = "freeze,spec,endpoints"),
  ",",
  fixed = TRUE
)[[1]]
requested_stages <- trimws(tolower(requested_stages))
allowed_stages <- c("freeze", "spec", "endpoints")
if (!all(requested_stages %in% allowed_stages)) {
  stop("TRACKC_STAGES may contain only: freeze,spec,endpoints")
}

seed <- 20260714L
n_permutations <- 10000L
mhc_chr <- 6L
mhc_start <- 25000000L
mhc_stop <- 34000000L

assert_new_directory <- function(path) {
  if (dir.exists(path)) {
    existing <- list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    if (length(existing)) {
      stop("Overwrite protection: output directory already contains files: ", path)
    }
  } else {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

start_stage <- function(path, stage) {
  marker <- list(
    stage = stage,
    status = "INCOMPLETE_OR_RUNNING",
    started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write_json(marker, file.path(path, "_STAGE_INCOMPLETE.json"), pretty = TRUE, auto_unbox = TRUE)
}

complete_stage <- function(path, stage, details = list()) {
  marker <- c(list(
    stage = stage,
    status = "COMPLETE",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), details)
  write_json(marker, file.path(path, "_STAGE_COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE)
  incomplete <- file.path(path, "_STAGE_INCOMPLETE.json")
  if (file.exists(incomplete)) unlink(incomplete)
}

write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, path, bom = TRUE, na = "")
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

file_provenance <- function(paths) {
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  info <- file.info(paths)
  data.table(
    path = paths,
    bytes = as.numeric(info$size),
    file_mtime = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    sha256 = vapply(paths, sha256_file, character(1))
  )
}

read_params <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  section <- ""
  out <- list()
  for (line in trimws(x)) {
    if (grepl("^\\[.*\\]$", line)) {
      section <- sub("^\\[(.*)\\]$", "\\1", line)
    } else if (nzchar(line) && !startsWith(line, "#") && grepl("=", line, fixed = TRUE)) {
      pieces <- strsplit(line, "=", fixed = TRUE)[[1]]
      key <- trimws(pieces[1])
      value <- trimws(paste(pieces[-1], collapse = "="))
      out[[paste(section, key, sep = ".")]] <- value
    }
  }
  out
}

bh_adjust <- function(p) p.adjust(as.numeric(p), method = "BH")

or_ci <- function(a, b, c, d) {
  cells <- as.numeric(c(a, b, c, d))
  corrected <- any(cells == 0)
  if (corrected) cells <- cells + 0.5
  odds <- cells[1] * cells[4] / (cells[2] * cells[3])
  se <- sqrt(sum(1 / cells))
  data.table(
    odds_ratio_corrected_if_zero = odds,
    odds_ratio_ci95_low = exp(log(odds) - 1.96 * se),
    odds_ratio_ci95_high = exp(log(odds) + 1.96 * se),
    haldane_anscombe_correction = corrected
  )
}

## ----------------------------------------------------------------------------
## Stage 00: freeze JOB01/JOB06/JOB07 MAGMA gene sets
## ----------------------------------------------------------------------------

run_freeze_stage <- function() {
  out_dir <- file.path(trackc_root, "00_frozen_gene_sets")
  assert_new_directory(out_dir)
  start_stage(out_dir, "freeze")

  job_dirs <- c(
    common = file.path(base_dir, "FUMA", "FUMA_JOB01_affective_functional_factor_GWAS"),
    affective = file.path(base_dir, "FUMA", "JOB06"),
    functional = file.path(base_dir, "FUMA", "JOB07")
  )
  magma_paths <- file.path(job_dirs, "magma.genes.out")
  params_paths <- file.path(job_dirs, "params.config")
  if (!all(file.exists(c(magma_paths, params_paths)))) {
    stop("One or more JOB01/JOB06/JOB07 source files are missing.")
  }

  magma <- lapply(magma_paths, fread)
  names(magma) <- names(job_dirs)
  tested_counts <- vapply(magma, nrow, integer(1))
  if (!all(tested_counts == 18516L)) {
    stop("Unexpected MAGMA tested-gene counts: ", paste(tested_counts, collapse = ", "))
  }
  gene_universes <- lapply(magma, function(x) sort(as.character(x$GENE)))
  if (!identical(gene_universes$common, gene_universes$affective) ||
      !identical(gene_universes$common, gene_universes$functional)) {
    stop("MAGMA Ensembl tested-gene universes are not identical.")
  }

  threshold <- 0.05 / 18516
  significant <- lapply(magma, function(x) {
    y <- copy(x)[as.numeric(P) < threshold]
    y[, gene_symbol := fifelse(
      !is.na(SYMBOL) & nzchar(trimws(as.character(SYMBOL))),
      toupper(trimws(as.character(SYMBOL))),
      as.character(GENE)
    )]
    y[]
  })
  symbols <- lapply(significant, function(x) unique(x$gene_symbol))

  frozen_sets <- list(
    cross_model_core = Reduce(intersect, symbols),
    affective_prioritized = setdiff(symbols$affective, symbols$functional),
    functional_prioritized = setdiff(symbols$functional, symbols$affective)
  )
  observed <- vapply(frozen_sets, length, integer(1))
  expected <- c(
    cross_model_core = 17L,
    affective_prioritized = 41L,
    functional_prioritized = 259L
  )
  if (!identical(observed[names(expected)], expected)) {
    stop("Frozen-set counts do not equal 17/41/259: ", paste(observed, collapse = ", "))
  }

  union_symbols <- sort(unique(unlist(symbols)))
  mapping <- rbindlist(lapply(union_symbols, function(symbol) {
    rows <- rbindlist(lapply(names(significant), function(domain) {
      significant[[domain]][gene_symbol == symbol][, domain := domain][]
    }), fill = TRUE)
    exemplar <- rows[1]
    data.table(
      gene_symbol = symbol,
      ensembl_ids = paste(sort(unique(rows$GENE)), collapse = ";"),
      n_ensembl_ids = uniqueN(rows$GENE),
      chr = as.integer(exemplar$CHR),
      start_grch37 = min(as.integer(rows$START), na.rm = TRUE),
      stop_grch37 = max(as.integer(rows$STOP), na.rm = TRUE),
      is_extended_mhc = as.integer(exemplar$CHR) == mhc_chr &
        min(as.integer(rows$START), na.rm = TRUE) <= mhc_stop &
        max(as.integer(rows$STOP), na.rm = TRUE) >= mhc_start,
      common_bonferroni = symbol %in% symbols$common,
      affective_bonferroni = symbol %in% symbols$affective,
      functional_bonferroni = symbol %in% symbols$functional,
      frozen_set = paste(names(frozen_sets)[vapply(frozen_sets, function(z) symbol %in% z, logical(1))], collapse = ";")
    )
  }))
  write_csv_utf8(mapping, file.path(out_dir, "ensembl_symbol_mapping_audit.csv"))

  non_mhc_counts <- integer()
  for (set_name in names(frozen_sets)) {
    full <- mapping[gene_symbol %in% frozen_sets[[set_name]]]
    non_mhc <- full[is_extended_mhc == FALSE]
    setorder(full, gene_symbol)
    setorder(non_mhc, gene_symbol)
    write_csv_utf8(full, file.path(out_dir, paste0(set_name, "_full.csv")))
    write_csv_utf8(non_mhc, file.path(out_dir, paste0(set_name, "_non_mhc.csv")))
    non_mhc_counts[set_name] <- nrow(non_mhc)
  }

  duplicate_symbols <- lapply(significant, function(x) {
    x[, .(ensembl_ids = paste(sort(unique(GENE)), collapse = ";"), n = uniqueN(GENE)), by = gene_symbol][n > 1]
  })
  duplicate_rows <- rbindlist(lapply(names(duplicate_symbols), function(domain) {
    duplicate_symbols[[domain]][, domain := domain][]
  }), fill = TRUE)
  write_csv_utf8(duplicate_rows, file.path(out_dir, "duplicate_symbol_audit.csv"))

  params <- lapply(params_paths, read_params)
  names(params) <- names(job_dirs)
  all_keys <- sort(unique(unlist(lapply(params, names))))
  parameter_audit <- rbindlist(lapply(all_keys, function(key) {
    vals <- vapply(params, function(x) if (is.null(x[[key]])) "<missing>" else x[[key]], character(1))
    data.table(
      parameter = key,
      common = vals["common"],
      affective = vals["affective"],
      functional = vals["functional"],
      consistent = uniqueN(vals) == 1L
    )
  }))
  write_csv_utf8(parameter_audit, file.path(out_dir, "fuma_parameter_audit_complete.csv"))

  provenance <- file_provenance(c(magma_paths, params_paths))
  provenance[, domain := rep(names(job_dirs), 2)]
  write_csv_utf8(provenance, file.path(out_dir, "input_provenance_sha256.csv"))

  manifest <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    status = "FROZEN_BEFORE_EXTERNAL_REASSESSMENT",
    magma_tested_genes = as.list(tested_counts),
    magma_bonferroni_threshold = threshold,
    full_counts = as.list(observed),
    non_mhc_counts = as.list(non_mhc_counts),
    mhc_definition = "GRCh37 chr6:25,000,000-34,000,000; gene interval overlap",
    mapping_rule = "FUMA Ensembl v102 SYMBOL; blank symbol falls back to Ensembl ID; exact symbol deduplication retains all IDs",
    interpretation_constraints = c(
      "Prioritized denotes a significance pattern, not a formally tested latent-factor effect difference.",
      "Cross-model core is within-study convergence across correlated models, not independent replication.",
      "The three frozen gene sets remain unchanged throughout endpoint analysis.",
      "No genes may be added or removed after viewing C3 or Track D results."
    )
  )
  write_json(manifest, file.path(out_dir, "freeze_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

  audit_lines <- c(
    "# Track C frozen gene-set audit",
    "",
    paste0("Generated: ", manifest$generated_at),
    "",
    sprintf("Identical MAGMA universe: 18,516 Ensembl IDs; Bonferroni threshold %.15g.", threshold),
    sprintf("Full sets: core %d; affective-prioritized %d; functional-prioritized %d.", observed[1], observed[2], observed[3]),
    sprintf("Non-MHC sets: core %d; affective-prioritized %d; functional-prioritized %d.", non_mhc_counts[1], non_mhc_counts[2], non_mhc_counts[3]),
    "",
    "JOB01 differs from JOB06/07 in genetype and optional expression/positional annotation resources. The MAGMA tested-gene universe is identical, so the MAGMA comparison is valid; all FUMA outputs must not be described as completely parameter-identical.",
    "",
    "These are significance-pattern sets, not domain-specific effects."
  )
  writeLines(audit_lines, file.path(out_dir, "FREEZE_AUDIT.md"), useBytes = TRUE)
  complete_stage(out_dir, "freeze", list(
    full_counts = as.list(observed),
    non_mhc_counts = as.list(non_mhc_counts)
  ))
  message("Completed freeze stage: ", out_dir)
  invisible(out_dir)
}

## ----------------------------------------------------------------------------
## Stage 01: lock backgrounds, endpoints, test families and matching variables
## ----------------------------------------------------------------------------

load_detectable_window_genes <- function(path) {
  x <- fread(path, select = c("gene_symbol", "n_snps_in_window", "in_background_for_test"))
  x[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
  x <- x[tolower(as.character(in_background_for_test)) %in% c("true", "1")]
  x[, .(median_window_snps = median(as.numeric(n_snps_in_window), na.rm = TRUE)), by = gene_symbol]
}

audit_window_backgrounds <- function(named_paths) {
  sets <- list()
  rows <- list()
  for (source in names(named_paths)) {
    x <- fread(named_paths[[source]], select = c("gene_symbol", "phenocode", "in_background_for_test"))
    x[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
    x <- x[tolower(as.character(in_background_for_test)) %in% c("true", "1")]
    for (code in unique(as.character(x$phenocode))) {
      genes <- sort(unique(x[as.character(phenocode) == code, gene_symbol]))
      key <- paste(source, code, sep = "::")
      sets[[key]] <- genes
      rows[[length(rows) + 1L]] <- data.table(
        source = source,
        phenocode = code,
        background_n = length(genes),
        gene_set_sha256 = digest(paste(genes, collapse = "\n"), algo = "sha256", serialize = FALSE)
      )
    }
  }
  if (length(sets) != 23L) stop("Expected 23 endpoint-specific window backgrounds; found ", length(sets))
  reference <- sets[[1]]
  identical_flags <- vapply(sets, identical, logical(1), y = reference)
  if (!all(identical_flags)) {
    bad <- names(identical_flags)[!identical_flags]
    stop("The 23 endpoint backgrounds are not identical: ", paste(bad, collapse = ", "))
  }
  list(
    shared_genes = reference,
    audit = rbindlist(rows)[, identical_to_shared_background := TRUE][]
  )
}

run_spec_stage <- function() {
  frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
  if (!file.exists(file.path(frozen_dir, "freeze_manifest.json"))) {
    stop("Freeze stage output is missing. Run freeze first.")
  }
  out_dir <- file.path(trackc_root, "01_analysis_specification")
  assert_new_directory(out_dir)
  start_stage(out_dir, "spec")

  magma_paths <- c(
    common = file.path(base_dir, "FUMA", "FUMA_JOB01_affective_functional_factor_GWAS", "magma.genes.out"),
    affective = file.path(base_dir, "FUMA", "JOB06", "magma.genes.out"),
    functional = file.path(base_dir, "FUMA", "JOB07", "magma.genes.out")
  )
  magma <- lapply(magma_paths, fread)
  names(magma) <- names(magma_paths)
  common <- copy(magma$common)
  common[, gene_symbol := fifelse(!is.na(SYMBOL) & nzchar(trimws(as.character(SYMBOL))), toupper(trimws(as.character(SYMBOL))), as.character(GENE))]

  universe <- common[, .(
    ensembl_ids = paste(sort(unique(GENE)), collapse = ";"),
    n_ensembl_ids = uniqueN(GENE),
    chr = as.integer(first(CHR)),
    start_grch37 = min(as.integer(START)),
    stop_grch37 = max(as.integer(STOP))
  ), by = gene_symbol]
  universe[, gene_length_bp := stop_grch37 - start_grch37 + 1]
  universe[, center := floor((start_grch37 + stop_grch37) / 2)]
  universe[, gene_density_plusminus_1mb := vapply(center, function(z) sum(abs(center - z) <= 1000000), integer(1)), by = chr]

  for (domain in names(magma)) {
    x <- copy(magma[[domain]])
    x[, gene_symbol := fifelse(!is.na(SYMBOL) & nzchar(trimws(as.character(SYMBOL))), toupper(trimws(as.character(SYMBOL))), as.character(GENE))]
    agg <- x[, .(
      nsnps = median(as.numeric(NSNPS), na.rm = TRUE),
      nparam = median(as.numeric(NPARAM), na.rm = TRUE)
    ), by = gene_symbol]
    setnames(agg, c("nsnps", "nparam"), paste0(domain, c("_nsnps", "_nparam")))
    universe <- merge(universe, agg, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  }
  universe[, magma_nsnps_median := apply(.SD, 1, median, na.rm = TRUE), .SDcols = patterns("_nsnps$")]
  universe[, magma_nparam_median := apply(.SD, 1, median, na.rm = TRUE), .SDcols = patterns("_nparam$")]
  universe[, effective_parameter_ratio := magma_nparam_median / pmax(magma_nsnps_median, 1)]
  universe[, is_extended_mhc := chr == mhc_chr & start_grch37 <= mhc_stop & stop_grch37 >= mhc_start]

  c3_tables <- file.path(trackb, "trackC3_triangular_validation_immune_spatial", "tables")
  immune_age_path <- file.path(c3_tables, "trackC3_immune_age_deg_all_linear.csv")
  immune_marker_path <- file.path(c3_tables, "trackC3_immune_marker_all.csv")
  immune_age <- fread(immune_age_path, select = "gene_symbol")
  immune_marker <- fread(immune_marker_path, select = "gene_symbol")
  immune_age_genes <- unique(toupper(trimws(immune_age$gene_symbol)))
  immune_marker_genes <- unique(toupper(trimws(immune_marker$gene_symbol)))

  spatial_dir <- file.path(base_dir, "外部验证", "空间蛋白组资源")
  moesm6 <- file.path(spatial_dir, "41586_2026_10660_MOESM6_ESM.xlsx")
  moesm7 <- file.path(spatial_dir, "41586_2026_10660_MOESM7_ESM.xlsx")
  spatial_median <- as.data.table(read_excel(moesm6, sheet = "E.specificity.median", .name_repair = "unique"))
  spatial_compare <- as.data.table(read_excel(moesm7, sheet = "B.All comparison dataset", .name_repair = "unique"))
  spatial_genes <- unique(toupper(trimws(c(as.character(spatial_median$gene), as.character(spatial_compare$Genes)))))
  spatial_genes <- spatial_genes[nzchar(spatial_genes) & !is.na(spatial_genes)]

  d2_tables <- file.path(trackb, "trackD2_tpmi_transancestry_validation", "tables")
  d5_tables <- file.path(trackb, "trackD5_finngen_clinical_endpoint_replication", "tables")
  d6_tables <- file.path(trackb, "trackD6_finngen_endpoint_specificity", "tables")
  tpmi_window_path <- file.path(d2_tables, "trackD2_tpmi_gene_window_summary_all_phenotypes.csv")
  fg_pos_window_path <- file.path(d5_tables, "trackD5_finngen_gene_window_summary_all_endpoints.csv")
  fg_neg_window_path <- file.path(d6_tables, "trackD6_negative_control_gene_window_summary_all_endpoints.csv")
  window_background_audit <- audit_window_backgrounds(c(
    D2_TPMI = tpmi_window_path,
    D5_FinnGen_positive = fg_pos_window_path,
    D6_FinnGen_controls = fg_neg_window_path
  ))
  if (length(window_background_audit$shared_genes) != 16354L) {
    stop("Shared endpoint window background is not 16,354 genes: ", length(window_background_audit$shared_genes))
  }
  write_csv_utf8(
    window_background_audit$audit,
    file.path(out_dir, "endpoint_window_background_consistency_audit.csv")
  )
  write_csv_utf8(
    data.table(gene_symbol = window_background_audit$shared_genes),
    file.path(out_dir, "shared_23_endpoint_window_background_raw.csv")
  )
  tpmi_detect <- load_detectable_window_genes(tpmi_window_path)
  fg_pos_detect <- load_detectable_window_genes(fg_pos_window_path)
  fg_neg_detect <- load_detectable_window_genes(fg_neg_window_path)
  fg_detect <- merge(fg_pos_detect, fg_neg_detect, by = "gene_symbol", suffixes = c("_pos", "_neg"))
  fg_detect[, median_window_snps := apply(.SD, 1, median, na.rm = TRUE), .SDcols = patterns("median_window_snps_")]
  if (!identical(sort(tpmi_detect$gene_symbol), sort(fg_detect$gene_symbol))) {
    stop("TPMI and FinnGen resource-level detectable gene backgrounds are not identical.")
  }

  universe[, in_immune_age_reported_universe := gene_symbol %in% immune_age_genes]
  universe[, in_immune_marker_reported_universe := gene_symbol %in% immune_marker_genes]
  universe[, in_spatial_proteome_detectable := gene_symbol %in% spatial_genes]
  universe[, in_tpmi_gene_window_detectable := gene_symbol %in% tpmi_detect$gene_symbol]
  universe[, in_finngen_gene_window_detectable := gene_symbol %in% fg_detect$gene_symbol]
  universe <- merge(universe, tpmi_detect[, .(gene_symbol, tpmi_median_window_snps = median_window_snps)], by = "gene_symbol", all.x = TRUE, sort = FALSE)
  universe <- merge(universe, fg_detect[, .(gene_symbol, finngen_median_window_snps = median_window_snps)], by = "gene_symbol", all.x = TRUE, sort = FALSE)
  write_csv_utf8(universe, file.path(out_dir, "magma_symbol_universe_and_matching_covariates.csv"))

  resource_rules <- list(
    immune_age_reported_gene_universe = universe[in_immune_age_reported_universe == TRUE, gene_symbol],
    immune_marker_reported_gene_universe = universe[in_immune_marker_reported_universe == TRUE, gene_symbol],
    spatial_proteome_detectable = universe[in_spatial_proteome_detectable == TRUE, gene_symbol],
    tpmi_gene_window_detectable = universe[in_tpmi_gene_window_detectable == TRUE, gene_symbol],
    finngen_gene_window_detectable = universe[in_finngen_gene_window_detectable == TRUE, gene_symbol]
  )
  background_summary <- rbindlist(lapply(names(resource_rules), function(resource) {
    genes <- sort(unique(resource_rules[[resource]]))
    write_csv_utf8(data.table(gene_symbol = genes), file.path(out_dir, paste0("background_", resource, ".csv")))
    data.table(
      resource = resource,
      magma_intersection_count = length(genes),
      background_rule = "intersection with the JOB01/JOB06/JOB07 identical MAGMA universe after unique-symbol mapping"
    )
  }))
  write_csv_utf8(background_summary, file.path(out_dir, "resource_background_summary.csv"))

  manifest_specs <- list(
    list(source = "D2_TPMI", family = "tpmi_endpoints", path = file.path(d2_tables, "trackD2_tpmi_input_manifest.csv")),
    list(source = "D5_FinnGen_positive", family = "finngen_endpoints", path = file.path(d5_tables, "trackD5_finngen_input_manifest.csv")),
    list(source = "D6_FinnGen_controls", family = "finngen_endpoints", path = file.path(d6_tables, "trackD6_finngen_negative_control_input_manifest.csv"))
  )
  endpoints <- rbindlist(lapply(manifest_specs, function(spec) {
    x <- fread(spec$path)
    data.table(
      source = spec$source,
      test_family = spec$family,
      phenocode = as.character(x$phenocode),
      phenotype = if ("phenotype_label" %in% names(x)) as.character(x$phenotype_label) else as.character(x$phenotype),
      validation_domain = as.character(x$validation_domain),
      endpoint_family = if ("endpoint_family" %in% names(x)) as.character(x$endpoint_family) else "",
      endpoint_specificity_class = if ("endpoint_specificity_class" %in% names(x)) as.character(x$endpoint_specificity_class) else "positive_or_primary",
      primary_role = as.character(x$primary_role),
      input_manifest = normalizePath(spec$path, winslash = "/", mustWork = TRUE)
    )
  }), fill = TRUE)
  if (nrow(endpoints[source == "D2_TPMI"]) != 8L ||
      nrow(endpoints[source == "D5_FinnGen_positive"]) != 10L ||
      nrow(endpoints[source == "D6_FinnGen_controls"]) != 5L) {
    stop("Endpoint counts are not 8 TPMI + 10 FinnGen positive + 5 FinnGen controls.")
  }
  write_csv_utf8(endpoints, file.path(out_dir, "frozen_endpoint_manifest.csv"))

  test_families <- data.table(
    family = c(
      "C3_single_cell_celltype_enrichment", "C3_spatial_proteome_enrichment",
      "D2_TPMI_endpoints", "D5_D6_FinnGen_endpoints",
      "candidate_eqtl_followup", "candidate_pqtl_followup"
    ),
    correction = "BH within scientific family; global BH sensitivity analysis",
    primary_version = "full frozen gene sets",
    sensitivity_version = "non-MHC frozen gene sets"
  )
  write_csv_utf8(test_families, file.path(out_dir, "frozen_test_family_manifest.csv"))

  source_paths <- c(
    magma_paths, file.path(frozen_dir, "freeze_manifest.json"),
    immune_age_path, immune_marker_path, moesm6, moesm7,
    tpmi_window_path, fg_pos_window_path, fg_neg_window_path,
    vapply(manifest_specs, `[[`, character(1), "path")
  )
  write_csv_utf8(file_provenance(source_paths), file.path(out_dir, "input_provenance_sha256.csv"))

  specification <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    status = "LOCKED_BEFORE_TESTING",
    gene_sets = c("cross_model_core", "affective_prioritized", "functional_prioritized"),
    primary_versions = "full",
    sensitivity_versions = "non-MHC",
    endpoint_count = nrow(endpoints),
    shared_raw_window_background_count = length(window_background_audit$shared_genes),
    shared_raw_window_background_expected_count = 16354L,
    formal_background_rule = "shared 23-endpoint window background intersected with the identical MAGMA unique-symbol universe; report the actual intersection count",
    matching_covariates = c(
      "chromosome", "gene length", "gene density within +/-1 Mb",
      "MAGMA median NSNPS", "MAGMA median NPARAM", "NPARAM/NSNPS",
      "resource detectability", "resource median gene-window SNP count"
    ),
    permutations = n_permutations,
    seed = seed,
    external_evidence_label = "cross-resource validation / external endpoint reassessment; not blinded independent replication",
    limitations = c(
      "The immune-age workbook exposes reported DEG genes rather than a complete detectable transcriptome.",
      "NPARAM/NSNPS is a matching proxy, not a direct LD-score estimate.",
      "The resources were viewed in v1 and therefore do not constitute blinded replication."
    )
  )
  write_json(specification, file.path(out_dir, "analysis_specification.json"), pretty = TRUE, auto_unbox = TRUE)
  decision_lines <- c(
    "# Track C method decision lock",
    "",
    paste0("Locked: ", specification$generated_at),
    "",
    "1. The 17/41/259 operational definitions are retained.",
    "2. All 23 endpoint-specific window backgrounds must be identical and contain 16,354 raw valid genes.",
    "3. Formal endpoint enrichment uses the shared raw window background intersected with the MAGMA unique-symbol universe; the actual intersection count is reported.",
    "4. A non-MHC analysis removes MHC genes from the target set, Fisher background, endpoint hit matrix and permutation control pool.",
    "5. Matched permutations report pre-match and post-match SMD values, maximum absolute SMD, chromosome counts, fallback and failure diagnostics.",
    "6. Absolute SMD below 0.1 is an ideal reference, not a post hoc exclusion rule.",
    "7. Fisher conditional exact 95% confidence intervals accompany Woolf/Haldane-Anscombe intervals.",
    "8. Only the three frozen gene sets enter the primary endpoint correction denominator.",
    "9. Immune age-DEG enrichment remains limited by the absence of a complete detectable-transcriptome background.",
    "10. Results are external endpoint reassessment / cross-resource validation, not blinded independent replication."
  )
  writeLines(decision_lines, file.path(out_dir, "METHOD_DECISION_LOCK.md"), useBytes = TRUE)
  complete_stage(out_dir, "spec", list(endpoint_count = nrow(endpoints)))
  message("Completed specification stage: ", out_dir)
  invisible(out_dir)
}

## ----------------------------------------------------------------------------
## Stage 02: D2/D5/D6 endpoint reassessment
## ----------------------------------------------------------------------------

make_balance_strata <- function(universe, target_genes, resource_snp_col) {
  x <- copy(universe)
  x[, is_target := gene_symbol %in% target_genes]
  x[, log_gene_length := log1p(as.numeric(gene_length_bp))]
  x[, log_magma_nsnps := log1p(as.numeric(magma_nsnps_median))]
  x[, resource_window_snps := log1p(as.numeric(get(resource_snp_col)))]
  covars <- c(
    "log_gene_length", "gene_density_plusminus_1mb", "log_magma_nsnps",
    "effective_parameter_ratio", "resource_window_snps"
  )
  for (v in covars) {
    z <- as.numeric(x[[v]])
    z[!is.finite(z)] <- median(z[is.finite(z)], na.rm = TRUE)
    med <- median(z, na.rm = TRUE)
    scale <- 1.4826 * median(abs(z - med), na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- sd(z, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    x[[paste0(v, "_z")]] <- (z - med) / scale
  }
  zcols <- paste0(covars, "_z")
  direction <- colMeans(x[is_target == TRUE, ..zcols]) - colMeans(x[is_target == FALSE, ..zcols])
  x[, balance_score := as.numeric(as.matrix(.SD) %*% direction), .SDcols = zcols]

  strata <- rbindlist(lapply(sort(unique(x$chr)), function(chrom) {
    y <- x[chr == chrom]
    if (!any(y$is_target)) return(NULL)
    y[, matching_fallback_to_whole_chromosome := FALSE]
    probs <- seq(0, 1, length.out = 6)
    breaks <- unique(as.numeric(quantile(y$balance_score, probs = probs, na.rm = TRUE, type = 8)))
    if (length(breaks) < 2) {
      y[, score_bin := 1L]
    } else {
      breaks[1] <- -Inf
      breaks[length(breaks)] <- Inf
      y[, score_bin := as.integer(cut(balance_score, breaks = breaks, include.lowest = TRUE, labels = FALSE))]
      counts <- y[, .(target_n = sum(is_target), control_n = sum(!is_target)), by = score_bin]
      if (any(counts[target_n > 0, control_n < target_n])) {
        y[, score_bin := 1L]
        y[, matching_fallback_to_whole_chromosome := TRUE]
      }
    }
    y[, stratum := paste(chr, score_bin, sep = ":")]
    y
  }), fill = TRUE)
  attr(strata, "matching_covariates") <- zcols
  strata
}

build_null_counts <- function(universe, target_genes, resource_snp_col, endpoint_hit_matrix) {
  matched <- make_balance_strata(universe, target_genes, resource_snp_col)
  balance_covariates <- attr(matched, "matching_covariates")
  target <- matched[is_target == TRUE]
  controls <- matched[is_target == FALSE]
  if (!nrow(target)) stop("No frozen-set genes intersect the resource background.")
  control_index <- match(controls$gene_symbol, rownames(endpoint_hit_matrix))
  if (anyNA(control_index)) stop("Endpoint matrix is missing matched control genes.")
  controls[, matrix_row := control_index]
  groups <- target[, .(target_n = .N), by = stratum]
  pools <- merge(groups, controls[, .(pool = list(matrix_row)), by = stratum], by = "stratum", all.x = TRUE)
  if (any(lengths(pools$pool) < pools$target_n)) stop("Insufficient controls in a matching stratum.")

  null_counts <- matrix(0L, nrow = n_permutations, ncol = ncol(endpoint_hit_matrix))
  colnames(null_counts) <- colnames(endpoint_hit_matrix)
  post_smd <- matrix(NA_real_, nrow = n_permutations, ncol = length(balance_covariates))
  colnames(post_smd) <- sub("_z$", "", balance_covariates)
  target_cov <- as.matrix(target[, ..balance_covariates])
  all_control_cov <- as.matrix(controls[, ..balance_covariates])
  target_mean <- colMeans(target_cov)
  target_var <- apply(target_cov, 2, var)
  control_mean_pre <- colMeans(all_control_cov)
  control_var_pre <- apply(all_control_cov, 2, var)
  pre_denom <- sqrt((target_var + control_var_pre) / 2)
  pre_smd <- (target_mean - control_mean_pre) / pre_denom
  pre_smd[!is.finite(pre_smd)] <- NA_real_
  control_lookup <- controls$matrix_row
  for (p in seq_len(n_permutations)) {
    selected <- unlist(Map(function(pool, k) sample(pool, size = k, replace = FALSE), pools$pool, pools$target_n), use.names = FALSE)
    null_counts[p, ] <- colSums(endpoint_hit_matrix[selected, , drop = FALSE])
    selected_cov <- all_control_cov[match(selected, control_lookup), , drop = FALSE]
    selected_mean <- colMeans(selected_cov)
    selected_var <- apply(selected_cov, 2, var)
    denom <- sqrt((target_var + selected_var) / 2)
    post_smd[p, ] <- (target_mean - selected_mean) / denom
  }
  post_smd[!is.finite(post_smd)] <- NA_real_
  balance <- data.table(
    covariate = colnames(post_smd),
    pre_match_smd = as.numeric(pre_smd),
    pre_match_abs_smd = abs(as.numeric(pre_smd)),
    post_match_mean_smd = colMeans(post_smd, na.rm = TRUE),
    post_match_mean_abs_smd = colMeans(abs(post_smd), na.rm = TRUE),
    post_match_max_abs_smd = apply(abs(post_smd), 2, max, na.rm = TRUE),
    post_match_prop_abs_smd_lt_0_1 = colMeans(abs(post_smd) < 0.1, na.rm = TRUE),
    ideal_abs_smd_reference = 0.1
  )
  chromosome_balance <- merge(
    target[, .(target_gene_n = .N), by = chr],
    controls[, .(available_control_gene_n = .N), by = chr],
    by = "chr",
    all.x = TRUE
  )
  chromosome_balance[, `:=`(
    matched_gene_n_each_permutation = target_gene_n,
    matching_preserves_chromosome_count = TRUE
  )]
  list(
    null_counts = null_counts,
    target_genes = target$gene_symbol,
    balance = balance,
    chromosome_balance = chromosome_balance,
    diagnostics = data.table(
      target_genes_in_background = nrow(target),
      control_genes = nrow(controls),
      matching_strata = nrow(pools),
      whole_chromosome_fallback_count = uniqueN(matched[matching_fallback_to_whole_chromosome == TRUE, chr]),
      insufficient_pool_count = 0L,
      matching_failure_count = 0L,
      unique_controls_within_permutation = TRUE,
      ideal_abs_smd_reference = 0.1,
      analyses_retained_regardless_of_balance_threshold = TRUE
    )
  )
}

run_endpoint_stage <- function() {
  frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
  spec_dir <- file.path(trackc_root, "01_analysis_specification")
  if (!file.exists(file.path(spec_dir, "analysis_specification.json"))) {
    stop("Specification stage output is missing. Run freeze and spec first.")
  }
  out_dir <- file.path(trackc_root, "02_endpoint_reassessment")
  assert_new_directory(out_dir)
  start_stage(out_dir, "endpoints")
  set.seed(seed)

  covariates <- fread(file.path(spec_dir, "magma_symbol_universe_and_matching_covariates.csv"))
  covariates[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
  set_names <- c("cross_model_core", "affective_prioritized", "functional_prioritized")
  versions <- c("full", "non_mhc")
  frozen_sets <- setNames(lapply(set_names, function(set_name) {
    setNames(lapply(versions, function(version) {
      x <- fread(file.path(frozen_dir, paste0(set_name, "_", version, ".csv")))
      unique(toupper(trimws(as.character(x$gene_symbol))))
    }), versions)
  }), set_names)

  d2_path <- file.path(trackb, "trackD2_tpmi_transancestry_validation", "tables", "trackD2_tpmi_gene_window_summary_all_phenotypes.csv")
  d5_path <- file.path(trackb, "trackD5_finngen_clinical_endpoint_replication", "tables", "trackD5_finngen_gene_window_summary_all_endpoints.csv")
  d6_path <- file.path(trackb, "trackD6_finngen_endpoint_specificity", "tables", "trackD6_negative_control_gene_window_summary_all_endpoints.csv")
  source_data <- list(
    D2_TPMI = fread(d2_path),
    D5_FinnGen_positive = fread(d5_path),
    D6_FinnGen_controls = fread(d6_path)
  )
  for (nm in names(source_data)) {
    source_data[[nm]][, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
    source_data[[nm]] <- source_data[[nm]][tolower(as.character(in_background_for_test)) %in% c("true", "1")]
  }
  endpoints <- fread(file.path(spec_dir, "frozen_endpoint_manifest.csv"), colClasses = list(character = "phenocode"))

  build_resource_matrix <- function(resource) {
    sources <- if (resource == "TPMI") "D2_TPMI" else c("D5_FinnGen_positive", "D6_FinnGen_controls")
    endpoint_sub <- endpoints[source %in% sources]
    background_flag <- if (resource == "TPMI") "in_tpmi_gene_window_detectable" else "in_finngen_gene_window_detectable"
    background <- covariates[get(background_flag) %in% TRUE]
    genes <- background$gene_symbol
    hit <- matrix(FALSE, nrow = length(genes), ncol = nrow(endpoint_sub), dimnames = list(genes, paste(endpoint_sub$source, endpoint_sub$phenocode, sep = "::")))
    tables <- list()
    for (i in seq_len(nrow(endpoint_sub))) {
      e <- endpoint_sub[i]
      x <- copy(source_data[[e$source]])[as.character(phenocode) == as.character(e$phenocode)]
      x[, sidak_gene_p_num := as.numeric(sidak_gene_p)]
      setorder(x, sidak_gene_p_num)
      x <- unique(x, by = "gene_symbol")
      x <- x[gene_symbol %in% genes]
      hit[x$gene_symbol, i] <- as.numeric(x$sidak_gene_p) < 0.05
      tables[[colnames(hit)[i]]] <- x
    }
    list(background = background, endpoints = endpoint_sub, hit = hit, tables = tables)
  }

  resources <- list(TPMI = build_resource_matrix("TPMI"), FinnGen = build_resource_matrix("FinnGen"))
  results <- list()
  diagnostics <- list()
  balance_audit <- list()
  chromosome_audit <- list()
  counter <- 0L
  for (resource_name in names(resources)) {
    resource <- resources[[resource_name]]
    snp_col <- if (resource_name == "TPMI") "tpmi_median_window_snps" else "finngen_median_window_snps"
    for (set_name in set_names) {
      for (version in versions) {
        analysis_background <- if (version == "non_mhc") {
          resource$background[is_extended_mhc %in% FALSE]
        } else {
          copy(resource$background)
        }
        analysis_hit <- resource$hit[analysis_background$gene_symbol, , drop = FALSE]
        if (version == "non_mhc" && any(analysis_background$is_extended_mhc %in% TRUE)) {
          stop("Internal error: non-MHC background still contains MHC genes.")
        }
        null <- build_null_counts(
          analysis_background,
          frozen_sets[[set_name]][[version]],
          snp_col,
          analysis_hit
        )
        diagnostics[[length(diagnostics) + 1L]] <- cbind(
          data.table(
            resource = resource_name,
            set_name = set_name,
            set_version = version,
            background_mhc_excluded = version == "non_mhc",
            analysis_background_n = nrow(analysis_background)
          ),
          null$diagnostics
        )
        balance_audit[[length(balance_audit) + 1L]] <- cbind(
          data.table(resource = resource_name, set_name = set_name, set_version = version),
          null$balance
        )
        chromosome_audit[[length(chromosome_audit) + 1L]] <- cbind(
          data.table(resource = resource_name, set_name = set_name, set_version = version),
          null$chromosome_balance
        )
        selected <- intersect(null$target_genes, rownames(analysis_hit))
        for (j in seq_len(nrow(resource$endpoints))) {
          e <- resource$endpoints[j]
          observed_hits <- sum(analysis_hit[selected, j])
          background_nonset <- setdiff(rownames(analysis_hit), selected)
          a <- observed_hits
          b <- length(selected) - a
          c <- sum(analysis_hit[background_nonset, j])
          d <- length(background_nonset) - c
          contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
          ft <- fisher.test(contingency, alternative = "greater")
          ft_exact_ci <- fisher.test(contingency, alternative = "two.sided", conf.int = TRUE)
          ci <- or_ci(a, b, c, d)
          null_hits <- null$null_counts[, j]
          counter <- counter + 1L
          results[[counter]] <- cbind(
            data.table(
              test_family = e$test_family,
              source = e$source,
              phenocode = e$phenocode,
              phenotype = e$phenotype,
              validation_domain = e$validation_domain,
              endpoint_family = e$endpoint_family,
              endpoint_specificity_class = e$endpoint_specificity_class,
              primary_role = e$primary_role,
              set_name = set_name,
              set_version = version,
              background_mhc_excluded = version == "non_mhc",
              background_n = nrow(analysis_hit),
              set_n_in_background = length(selected),
              set_hit_n_sidak_lt_0_05 = a,
              background_nonset_hit_n = c,
              fisher_odds_ratio = unname(ft$estimate),
              fisher_exact_conditional_mle_odds_ratio = unname(ft_exact_ci$estimate),
              fisher_exact_ci95_low = unname(ft_exact_ci$conf.int[1]),
              fisher_exact_ci95_high = unname(ft_exact_ci$conf.int[2]),
              fisher_p_greater = ft$p.value,
              matched_permutation_p = (1 + sum(null_hits >= observed_hits)) / (n_permutations + 1),
              matched_null_mean_hits = mean(null_hits),
              matched_null_sd_hits = sd(null_hits),
              n_permutations = n_permutations
            ),
            ci
          )
        }
      }
    }
  }
  result <- rbindlist(results, fill = TRUE)
  result[, fisher_bh_within_family_version := bh_adjust(fisher_p_greater), by = .(test_family, set_version)]
  result[, permutation_bh_within_family_version := bh_adjust(matched_permutation_p), by = .(test_family, set_version)]
  result[, fisher_global_bh_primary_endpoints := NA_real_]
  result[set_version == "full", fisher_global_bh_primary_endpoints := bh_adjust(fisher_p_greater)]
  result[, permutation_global_bh_primary_endpoints := NA_real_]
  result[set_version == "full", permutation_global_bh_primary_endpoints := bh_adjust(matched_permutation_p)]
  setorder(result, test_family, set_version, set_name, matched_permutation_p, fisher_p_greater)
  write_csv_utf8(result, file.path(out_dir, "trackC_endpoint_reassessment_results.csv"))
  write_csv_utf8(rbindlist(diagnostics, fill = TRUE), file.path(out_dir, "matched_permutation_diagnostics.csv"))
  write_csv_utf8(rbindlist(balance_audit, fill = TRUE), file.path(out_dir, "matched_covariate_balance.csv"))
  write_csv_utf8(rbindlist(chromosome_audit, fill = TRUE), file.path(out_dir, "matched_chromosome_balance.csv"))

  summary <- result[, .(
    tests = .N,
    fisher_bh_lt_0_05 = sum(fisher_bh_within_family_version < 0.05, na.rm = TRUE),
    permutation_bh_lt_0_05 = sum(permutation_bh_within_family_version < 0.05, na.rm = TRUE),
    min_fisher_p = min(fisher_p_greater, na.rm = TRUE),
    min_permutation_p = min(matched_permutation_p, na.rm = TRUE)
  ), by = .(test_family, set_version, set_name)]
  write_csv_utf8(summary, file.path(out_dir, "trackC_endpoint_reassessment_summary.csv"))

  audit <- list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed = seed,
    permutations = n_permutations,
    hit_definition = "gene-window Sidak P < 0.05",
    fisher_alternative = "greater",
    confidence_interval = "Woolf log-OR 95% CI; Haldane-Anscombe 0.5 correction if any cell is zero",
    exact_confidence_interval = "Fisher conditional maximum-likelihood odds ratio and two-sided exact 95% confidence interval",
    matching = "chromosome-exact adaptive multivariate-balance-score strata; without-replacement controls",
    non_mhc_rule = "MHC genes removed from the target set, Fisher background, endpoint hit matrix and permutation control pool",
    balance_reporting = "pre-match SMD, post-match mean signed/absolute SMD, maximum absolute SMD, proportion below 0.1, chromosome counts and failure diagnostics",
    results_rows = nrow(result),
    primary_full_results_rows = nrow(result[set_version == "full"])
  )
  write_json(audit, file.path(out_dir, "endpoint_reassessment_audit.json"), pretty = TRUE, auto_unbox = TRUE)
  complete_stage(out_dir, "endpoints", list(
    results_rows = nrow(result),
    permutations = n_permutations,
    non_mhc_background_synchronized = TRUE,
    balance_audit_written = TRUE
  ))
  message("Completed endpoint stage: ", out_dir)
  invisible(out_dir)
}

## ----------------------------------------------------------------------------
## Execute requested stages when source() is called
## ----------------------------------------------------------------------------

message("Track C output root: ", trackc_root)
message("Requested stages: ", paste(requested_stages, collapse = ", "))
if ("freeze" %in% requested_stages) run_freeze_stage()
if ("spec" %in% requested_stages) run_spec_stage()
if ("endpoints" %in% requested_stages) run_endpoint_stage()
message("Track C phase 1 finished.")
