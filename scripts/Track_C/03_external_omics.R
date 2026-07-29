## ============================================================================
## Track C C3: frozen-set external single-cell/spatial-omics reassessment
## User-run only. Does not modify v1 C3/D3 outputs or re-rank frozen genes.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(jsonlite)
  library(digest)
})

options(stringsAsFactors = FALSE)

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
trackc_root <- Sys.getenv(
  "TRACKC_ROOT",
  unset = file.path(base_dir, "analysis_ready_core", "trackC")
)
phase1c_lolo_dir <- Sys.getenv(
  "TRACKC_PHASE1C_LOLO_DIR",
  unset = file.path(trackc_root, "05_targeted_lolo")
)
out_dir <- Sys.getenv(
  "TRACKC_C3_OUT_DIR",
  unset = file.path(trackc_root, "06_C3_external_omics")
)

base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
trackc_root <- normalizePath(trackc_root, winslash = "/", mustWork = TRUE)
phase1c_lolo_dir <- normalizePath(phase1c_lolo_dir, winslash = "/", mustWork = TRUE)
frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
spec_dir <- file.path(trackc_root, "01_analysis_specification")
c3_dir <- file.path(base_dir, "analysis_ready_core", "trackB",
                       "trackC3_triangular_validation_immune_spatial", "tables")
spatial_m7_path <- file.path(base_dir, "外部验证", "空间蛋白组资源",
                             "41586_2026_10660_MOESM7_ESM.xlsx")

required_markers <- c(
  file.path(frozen_dir, "_STAGE_COMPLETE.json"),
  file.path(spec_dir, "_STAGE_COMPLETE.json"),
  file.path(phase1c_lolo_dir, "_STAGE_COMPLETE.json")
)
if (!all(file.exists(required_markers))) stop("Required upstream COMPLETE markers are missing.")
if (!file.exists(spatial_m7_path)) stop("Spatial proteome MOESM7 workbook is missing.")
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Overwrite protection: output directory already contains files: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_json(
  list(stage = "C3_external_omics", status = "INCOMPLETE_OR_RUNNING",
       started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  file.path(out_dir, "_STAGE_INCOMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)

n_permutations <- as.integer(Sys.getenv("TRACKC_C3_PERMUTATIONS", unset = "10000"))
random_seed <- 20260714L
set.seed(random_seed)
write_csv_utf8 <- function(x, path) fwrite(x, path, bom = TRUE, na = "")
sha256_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)
bh_adjust <- function(p) p.adjust(as.numeric(p), method = "BH")

set_names <- c("cross_model_core", "affective_prioritized", "functional_prioritized")
set_versions <- c("full", "non_mhc")

read_frozen_set <- function(set_name, set_version) {
  path <- file.path(frozen_dir, paste0(set_name, "_", set_version, ".csv"))
  unique(toupper(trimws(as.character(fread(path)$gene_symbol))))
}

covariates_path <- file.path(spec_dir, "magma_symbol_universe_and_matching_covariates.csv")
covariates <- fread(covariates_path)
covariates[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]

background_paths <- c(
  immune_age = file.path(spec_dir, "background_immune_age_reported_gene_universe.csv"),
  immune_marker = file.path(spec_dir, "background_immune_marker_reported_gene_universe.csv"),
  spatial_proteome = file.path(spec_dir, "background_spatial_proteome_detectable.csv")
)
backgrounds <- lapply(background_paths, function(path) {
  unique(toupper(trimws(as.character(fread(path)$gene_symbol))))
})
expected_background_n <- c(immune_age = 2266L, immune_marker = 5787L, spatial_proteome = 9031L)
actual_background_n <- vapply(backgrounds, length, integer(1))
if (!identical(actual_background_n[names(expected_background_n)], expected_background_n)) {
  stop("Frozen C3 full-background counts do not match the locked specification.")
}

age_path <- file.path(c3_dir, "trackC3_immune_age_deg_all_linear.csv")
marker_path <- file.path(c3_dir, "trackC3_immune_marker_all.csv")
age <- fread(age_path)
marker <- fread(marker_path)
age[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
marker[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]

## Immune-age sheets expose reported DEGs, not a complete detectable transcriptome.
## These tests therefore assess concentration across reported age-DEG signatures.
age_signatures <- unique(age[sig_age_deg == TRUE, .(
  test_family = "C3_single_cell_celltype_enrichment",
  resource = "immune_age_reported_DEG",
  signature_id = paste("age", sex, celltype, direction, sep = "|"),
  signature_label = paste(sex, celltype, direction, sep = " | "),
  signature_group = "reported_age_DEG_conditional",
  analysis_layer = "conditional_secondary",
  correction_family = "C3_reported_age_DEG_conditional",
  gene_symbol
)])

## Marker table contains reported tested marker rows and a prespecified sig_marker flag.
## Sex-stratified rows are collapsed to cell type x direction for the primary tests;
## this avoids treating identical male/female marker lists as independent evidence.
marker_signatures <- unique(marker[sig_marker == TRUE, .(
  test_family = "C3_single_cell_celltype_enrichment",
  resource = "immune_marker",
  signature_id = paste("marker", celltype, direction_nhood, sep = "|"),
  signature_label = paste(celltype, direction_nhood, sep = " | "),
  signature_group = "sex_collapsed_celltype_marker",
  analysis_layer = "primary_single_cell",
  correction_family = "C3_primary_single_cell_marker",
  gene_symbol
)])

## Spatial signatures include every named tissue plus three categories fixed a priori.
spatial_raw <- as.data.table(read_excel(spatial_m7_path, sheet = "A. All tissue enriched proteins"))
setnames(spatial_raw, "Gene", "gene_symbol")
spatial_raw[, gene_symbol := toupper(trimws(as.character(gene_symbol)))]
id_columns <- c("Tissue enriched protein", "gene_symbol")
tissue_columns <- setdiff(names(spatial_raw), id_columns)
for (v in tissue_columns) spatial_raw[, (v) := as.numeric(get(v))]
spatial_long <- melt(
  spatial_raw, id.vars = id_columns, measure.vars = tissue_columns,
  variable.name = "tissue", value.name = "enrichment_score"
)
spatial_tissue_signatures <- unique(spatial_long[
  is.finite(enrichment_score) & enrichment_score > 0,
  .(test_family = "C3_spatial_proteome_enrichment",
    resource = "spatial_proteome", signature_id = paste0("tissue|", tissue),
    signature_label = as.character(tissue), signature_group = "individual_tissue",
    analysis_layer = "exploratory_spatial_tissue",
    correction_family = "C3_exploratory_individual_tissues",
    gene_symbol)
])

spatial_categories <- list(
  brain_neural = c("Brain", "Nerve", "Spinal cord"),
  immune_lymphoid = c("Bone marrow", "Lymph node", "Thymus", "Spleen", "Tonsil"),
  musculoskeletal = c("Bone", "Cartilage", "Skeletal muscle", "Smooth muscle", "Tendon")
)
if (!all(unique(unlist(spatial_categories)) %in% tissue_columns)) {
  stop("One or more prespecified spatial category tissues are absent from MOESM7.")
}
spatial_category_signatures <- rbindlist(lapply(names(spatial_categories), function(category) {
  tissues <- spatial_categories[[category]]
  unique(spatial_long[tissue %in% tissues & is.finite(enrichment_score) & enrichment_score > 0,
                      .(test_family = "C3_spatial_proteome_enrichment",
                        resource = "spatial_proteome",
                        signature_id = paste0("category|", category),
                        signature_label = category, signature_group = "prespecified_tissue_category",
                        analysis_layer = "primary_spatial_category",
                        correction_family = "C3_primary_spatial_categories",
                        gene_symbol)])
}))
spatial_signatures <- unique(rbindlist(list(spatial_tissue_signatures, spatial_category_signatures)))

signature_rows <- rbindlist(list(age_signatures, marker_signatures, spatial_signatures), fill = TRUE)
signature_manifest <- unique(signature_rows[, .(
  test_family, resource, signature_id, signature_label, signature_group,
  analysis_layer, correction_family
)])
signature_manifest[, background_resource := fifelse(
  resource == "immune_age_reported_DEG", "immune_age",
  fifelse(resource == "immune_marker", "immune_marker", "spatial_proteome")
)]
signature_manifest[, interpretation_boundary := fifelse(
  resource == "immune_age_reported_DEG",
  "conditional concentration among reported age-DEG genes; not unbiased detectable-transcriptome enrichment",
  "cross-resource enrichment/support; not blinded independent replication"
)]
setorder(signature_manifest, test_family, resource, signature_id)
write_csv_utf8(signature_manifest, file.path(out_dir, "C3_signature_manifest.csv"))

make_balance_strata <- function(universe, target_genes) {
  x <- copy(universe)
  x[, is_target := gene_symbol %in% target_genes]
  x[, log_gene_length := log1p(as.numeric(gene_length_bp))]
  x[, log_magma_nsnps := log1p(as.numeric(magma_nsnps_median))]
  covars <- c("log_gene_length", "gene_density_plusminus_1mb",
              "log_magma_nsnps", "effective_parameter_ratio")
  for (v in covars) {
    z <- as.numeric(x[[v]])
    z[!is.finite(z)] <- median(z[is.finite(z)], na.rm = TRUE)
    med <- median(z, na.rm = TRUE)
    s <- 1.4826 * median(abs(z - med), na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- sd(z, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- 1
    x[[paste0(v, "_z")]] <- (z - med) / s
  }
  zcols <- paste0(covars, "_z")
  direction <- colMeans(x[is_target == TRUE, ..zcols]) - colMeans(x[is_target == FALSE, ..zcols])
  x[, balance_score := as.numeric(as.matrix(.SD) %*% direction), .SDcols = zcols]
  strata <- rbindlist(lapply(sort(unique(x$chr)), function(chromosome) {
    y <- x[chr == chromosome]
    if (!any(y$is_target)) return(NULL)
    y[, matching_fallback_to_whole_chromosome := FALSE]
    breaks <- unique(as.numeric(quantile(y$balance_score, probs = seq(0, 1, length.out = 6),
                                         na.rm = TRUE, type = 8)))
    if (length(breaks) < 2L) {
      y[, score_bin := 1L]
    } else {
      breaks[1] <- -Inf
      breaks[length(breaks)] <- Inf
      y[, score_bin := as.integer(cut(balance_score, breaks = breaks,
                                      include.lowest = TRUE, labels = FALSE))]
      counts <- y[, .(target_n = sum(is_target), control_n = sum(!is_target)), by = score_bin]
      if (any(counts[target_n > 0, control_n < target_n])) {
        y[, `:=`(score_bin = 1L, matching_fallback_to_whole_chromosome = TRUE)]
      }
    }
    y[, stratum := paste(chr, score_bin, sep = ":")]
    y
  }), fill = TRUE)
  attr(strata, "matching_covariates") <- zcols
  strata
}

build_null_counts <- function(universe, target_genes, hit_matrix) {
  matched <- make_balance_strata(universe, target_genes)
  balance_covariates <- attr(matched, "matching_covariates")
  target <- matched[is_target == TRUE]
  controls <- matched[is_target == FALSE]
  if (!nrow(target)) stop("No frozen-set genes intersect the resource background.")
  controls[, matrix_row := match(gene_symbol, rownames(hit_matrix))]
  if (anyNA(controls$matrix_row)) stop("Signature matrix is missing matched controls.")
  groups <- target[, .(target_n = .N), by = stratum]
  pools <- merge(groups, controls[, .(pool = list(matrix_row)), by = stratum],
                 by = "stratum", all.x = TRUE)
  if (any(lengths(pools$pool) < pools$target_n)) stop("Insufficient controls in a matching stratum.")
  null_counts <- matrix(0L, nrow = n_permutations, ncol = ncol(hit_matrix))
  colnames(null_counts) <- colnames(hit_matrix)
  zcols <- balance_covariates
  target_cov <- as.matrix(target[, ..zcols])
  control_cov <- as.matrix(controls[, ..zcols])
  target_mean <- colMeans(target_cov)
  target_var <- apply(target_cov, 2, var)
  pre_mean <- colMeans(control_cov)
  pre_var <- apply(control_cov, 2, var)
  pre_smd <- (target_mean - pre_mean) / sqrt((target_var + pre_var) / 2)
  post_smd <- matrix(NA_real_, nrow = n_permutations, ncol = length(zcols))
  colnames(post_smd) <- sub("_z$", "", zcols)
  lookup <- controls$matrix_row
  for (p in seq_len(n_permutations)) {
    selected <- unlist(Map(function(pool, k) sample(pool, k, replace = FALSE),
                           pools$pool, pools$target_n), use.names = FALSE)
    null_counts[p, ] <- colSums(hit_matrix[selected, , drop = FALSE])
    selected_cov <- control_cov[match(selected, lookup), , drop = FALSE]
    selected_var <- apply(selected_cov, 2, var)
    post_smd[p, ] <- (target_mean - colMeans(selected_cov)) /
      sqrt((target_var + selected_var) / 2)
  }
  pre_smd[!is.finite(pre_smd)] <- NA_real_
  post_smd[!is.finite(post_smd)] <- NA_real_
  list(
    null_counts = null_counts,
    target_genes = target$gene_symbol,
    balance = data.table(
      covariate = colnames(post_smd), pre_match_smd = as.numeric(pre_smd),
      pre_match_abs_smd = abs(as.numeric(pre_smd)),
      post_match_mean_smd = colMeans(post_smd, na.rm = TRUE),
      post_match_mean_abs_smd = colMeans(abs(post_smd), na.rm = TRUE),
      post_match_max_abs_smd = apply(abs(post_smd), 2, max, na.rm = TRUE),
      post_match_prop_abs_smd_lt_0_1 = colMeans(abs(post_smd) < 0.1, na.rm = TRUE),
      smd_estimable = is.finite(colMeans(abs(post_smd), na.rm = TRUE)),
      ideal_abs_smd_reference = 0.1
    ),
    diagnostics = data.table(
      target_genes_in_background = nrow(target), control_genes = nrow(controls),
      matching_strata = nrow(pools),
      whole_chromosome_fallback_count = uniqueN(matched[matching_fallback_to_whole_chromosome == TRUE, chr]),
      matching_failure_count = 0L, permutations = n_permutations
    )
  )
}

results <- list()
overlaps <- list()
balance_audit <- list()
matching_diagnostics <- list()
coverage_audit <- list()
failures <- list()
counter <- 0L

for (resource_name in unique(signature_manifest$resource)) {
  resource_key <- signature_manifest[resource == resource_name, unique(background_resource)]
  full_background <- backgrounds[[resource_key]]
  resource_signatures <- signature_manifest[resource == resource_name]
  resource_rows <- signature_rows[resource == resource_name]
  for (set_version in set_versions) {
    background_genes <- full_background
    if (set_version == "non_mhc") {
      mhc_genes <- covariates[is_extended_mhc == TRUE, gene_symbol]
      background_genes <- setdiff(background_genes, mhc_genes)
    }
    universe <- covariates[gene_symbol %in% background_genes]
    setorder(universe, gene_symbol)
    hit_matrix <- matrix(FALSE, nrow = nrow(universe), ncol = nrow(resource_signatures),
                         dimnames = list(universe$gene_symbol, resource_signatures$signature_id))
    for (j in seq_len(nrow(resource_signatures))) {
      genes <- resource_rows[signature_id == resource_signatures$signature_id[j], gene_symbol]
      hit_matrix[, j] <- rownames(hit_matrix) %in% genes
    }
    for (set_name in set_names) {
      counter <- counter + 1L
      target_genes <- intersect(read_frozen_set(set_name, set_version), universe$gene_symbol)
      coverage_n <- length(target_genes)
      coverage_tier <- if (coverage_n < 5L) {
        "descriptive_only_lt5"
      } else if (coverage_n < 10L) {
        "sparse_exploratory_5to9"
      } else {
        "inferential_ge10"
      }
      coverage_audit[[length(coverage_audit) + 1L]] <- data.table(
        resource = resource_name, set_name, set_version,
        background_n = nrow(universe), set_n_in_background = coverage_n,
        coverage_tier,
        fisher_estimable = coverage_n >= 5L,
        matched_permutation_estimable = coverage_n >= 10L
      )

      null <- NULL
      if (coverage_n >= 10L) {
        null <- tryCatch(build_null_counts(universe, target_genes, hit_matrix), error = function(e) e)
      }
      if (inherits(null, "error")) {
        failures[[length(failures) + 1L]] <- data.table(
          resource = resource_name, set_name, set_version, coverage_n,
          error_message = conditionMessage(null)
        )
        next
      }
      if (coverage_n >= 10L) {
        balance_audit[[length(balance_audit) + 1L]] <- cbind(
          data.table(resource = resource_name, set_name, set_version,
                     matching_status = "estimated", coverage_n), null$balance
        )
        matching_diagnostics[[length(matching_diagnostics) + 1L]] <- cbind(
          data.table(resource = resource_name, set_name, set_version,
                     matching_status = "estimated", coverage_n), null$diagnostics
        )
      } else {
        balance_audit[[length(balance_audit) + 1L]] <- data.table(
          resource = resource_name, set_name, set_version,
          matching_status = "not_estimable_due_to_target_coverage", coverage_n,
          covariate = c("log_gene_length", "gene_density_plusminus_1mb",
                        "log_magma_nsnps", "effective_parameter_ratio"),
          pre_match_smd = NA_real_, pre_match_abs_smd = NA_real_,
          post_match_mean_smd = NA_real_, post_match_mean_abs_smd = NA_real_,
          post_match_max_abs_smd = NA_real_, post_match_prop_abs_smd_lt_0_1 = NA_real_,
          smd_estimable = FALSE,
          ideal_abs_smd_reference = 0.1
        )
        matching_diagnostics[[length(matching_diagnostics) + 1L]] <- data.table(
          resource = resource_name, set_name, set_version,
          matching_status = "not_estimable_due_to_target_coverage", coverage_n,
          target_genes_in_background = coverage_n,
          control_genes = nrow(universe) - coverage_n,
          matching_strata = NA_integer_, whole_chromosome_fallback_count = NA_integer_,
          matching_failure_count = NA_integer_, permutations = 0L
        )
      }
      for (j in seq_len(ncol(hit_matrix))) {
        signature <- resource_signatures[signature_id == colnames(hit_matrix)[j]]
        target_hit <- sum(hit_matrix[target_genes, j])
        target_nonhit <- coverage_n - target_hit
        background_hit <- sum(hit_matrix[, j])
        non_target_hit <- background_hit - target_hit
        non_target_nonhit <- nrow(hit_matrix) - coverage_n - non_target_hit
        tab <- matrix(c(target_hit, target_nonhit, non_target_hit, non_target_nonhit),
                      nrow = 2, byrow = TRUE)
        if (coverage_n >= 5L) {
          fisher_greater <- fisher.test(tab, alternative = "greater")
          fisher_two_sided <- fisher.test(tab, alternative = "two.sided", conf.int = TRUE)
          fisher_or <- unname(fisher_two_sided$estimate)
          fisher_low <- fisher_two_sided$conf.int[1]
          fisher_high <- fisher_two_sided$conf.int[2]
          fisher_p <- fisher_greater$p.value
        } else {
          fisher_or <- fisher_low <- fisher_high <- fisher_p <- NA_real_
        }
        if (coverage_n >= 10L) {
          null_hits <- null$null_counts[, j]
          permutation_p <- (1 + sum(null_hits >= target_hit)) / (n_permutations + 1)
          null_mean <- mean(null_hits)
          null_sd <- sd(null_hits)
          permutations_used <- n_permutations
        } else {
          permutation_p <- null_mean <- null_sd <- NA_real_
          permutations_used <- 0L
        }
        overlap_genes <- target_genes[hit_matrix[target_genes, j]]
        results[[length(results) + 1L]] <- data.table(
          test_family = signature$test_family,
          resource = resource_name,
          signature_id = signature$signature_id,
          signature_label = signature$signature_label,
          signature_group = signature$signature_group,
          analysis_layer = signature$analysis_layer,
          correction_family = signature$correction_family,
          interpretation_boundary = signature$interpretation_boundary,
          set_name, set_version, coverage_tier,
          background_mhc_excluded = set_version == "non_mhc",
          background_n = nrow(hit_matrix), set_n_in_background = coverage_n,
          signature_n_in_background = background_hit, overlap_n = target_hit,
          fisher_odds_ratio = fisher_or,
          fisher_ci95_low = fisher_low,
          fisher_ci95_high = fisher_high,
          fisher_p_greater = fisher_p,
          matched_permutation_p = permutation_p,
          matched_null_mean_hits = null_mean, matched_null_sd_hits = null_sd,
          n_permutations = permutations_used,
          overlap_genes = paste(overlap_genes, collapse = ";")
        )
        if (length(overlap_genes)) {
          overlaps[[length(overlaps) + 1L]] <- data.table(
            test_family = signature$test_family, resource = resource_name,
            signature_id = signature$signature_id, signature_label = signature$signature_label,
            analysis_layer = signature$analysis_layer,
            set_name, set_version, coverage_tier, gene_symbol = overlap_genes
          )
        }
      }
      message(sprintf("[%d/18] %s | %s | %s", counter, resource_name, set_name, set_version))
    }
  }
}

result <- rbindlist(results, fill = TRUE)
failure_table <- rbindlist(failures, fill = TRUE)
if (!ncol(failure_table)) {
  failure_table <- data.table(resource = character(), set_name = character(),
                              set_version = character(), error_message = character())
}
overlap_table <- rbindlist(overlaps, fill = TRUE)
if (!ncol(overlap_table)) {
  overlap_table <- data.table(test_family = character(), resource = character(),
                              signature_id = character(), signature_label = character(),
                              set_name = character(), set_version = character(), gene_symbol = character())
}

result[, `:=`(fisher_bh_within_family_version = NA_real_,
              permutation_bh_within_family_version = NA_real_)]
result[coverage_tier != "descriptive_only_lt5",
       fisher_bh_within_family_version := bh_adjust(fisher_p_greater),
       by = .(correction_family, set_version)]
result[coverage_tier == "inferential_ge10",
       permutation_bh_within_family_version := bh_adjust(matched_permutation_p),
       by = .(correction_family, set_version)]
result[, `:=`(fisher_global_bh_full = NA_real_, permutation_global_bh_full = NA_real_)]
primary_rows <- result$set_version == "full" & result$coverage_tier == "inferential_ge10" &
  result$analysis_layer %in% c("primary_single_cell", "primary_spatial_category")
result[which(primary_rows), fisher_global_bh_full := bh_adjust(fisher_p_greater)]
result[which(primary_rows), permutation_global_bh_full := bh_adjust(matched_permutation_p)]
result[, both_methods_family_significant :=
         fisher_bh_within_family_version < 0.05 & permutation_bh_within_family_version < 0.05]
setorder(result, test_family, set_version, set_name, fisher_p_greater)

summary_table <- result[, .(
  tests = .N,
  fisher_estimable_n = sum(!is.na(fisher_p_greater)),
  permutation_estimable_n = sum(!is.na(matched_permutation_p)),
  fisher_family_significant_n = sum(fisher_bh_within_family_version < 0.05, na.rm = TRUE),
  permutation_family_significant_n = sum(permutation_bh_within_family_version < 0.05, na.rm = TRUE),
  both_methods_family_significant_n = sum(both_methods_family_significant, na.rm = TRUE),
  minimum_fisher_p = if (all(is.na(fisher_p_greater))) NA_real_ else min(fisher_p_greater, na.rm = TRUE),
  minimum_permutation_p = if (all(is.na(matched_permutation_p))) NA_real_ else min(matched_permutation_p, na.rm = TRUE)
), by = .(test_family, resource, analysis_layer, correction_family,
          set_name, set_version, coverage_tier)]

write_csv_utf8(result, file.path(out_dir, "C3_external_omics_enrichment_results.csv"))
write_csv_utf8(summary_table, file.path(out_dir, "C3_external_omics_summary.csv"))
write_csv_utf8(overlap_table, file.path(out_dir, "C3_gene_signature_overlaps.csv"))
coverage_table <- rbindlist(coverage_audit, fill = TRUE)
balance <- rbindlist(balance_audit, fill = TRUE)
matching_table <- rbindlist(matching_diagnostics, fill = TRUE)
write_csv_utf8(coverage_table, file.path(out_dir, "C3_coverage_and_estimability.csv"))
write_csv_utf8(balance, file.path(out_dir, "C3_matched_covariate_balance.csv"))
write_csv_utf8(rbindlist(matching_diagnostics, fill = TRUE), file.path(out_dir, "C3_matching_diagnostics.csv"))
write_csv_utf8(failure_table, file.path(out_dir, "C3_matching_failures.csv"))

combo_balance <- balance[matching_status == "estimated", .(
  max_post_match_mean_abs_smd = if (all(!is.finite(post_match_mean_abs_smd))) NA_real_ else
    max(post_match_mean_abs_smd[is.finite(post_match_mean_abs_smd)]),
  max_post_match_max_abs_smd = if (all(!is.finite(post_match_max_abs_smd))) NA_real_ else
    max(post_match_max_abs_smd[is.finite(post_match_max_abs_smd)]),
  any_smd_not_estimable = any(!smd_estimable)
), by = .(resource, set_name, set_version, coverage_n)]
combo_balance[, `:=`(
  imbalance_threshold = 0.1,
  requires_uniform_C3b_adjustment = any_smd_not_estimable |
    is.na(max_post_match_mean_abs_smd) | max_post_match_mean_abs_smd >= 0.1,
  C3b_scope = fifelse(any_smd_not_estimable | is.na(max_post_match_mean_abs_smd) |
                        max_post_match_mean_abs_smd >= 0.1,
                      "all signatures in this resource-set-version combination",
                      "not required")
)]
write_csv_utf8(combo_balance, file.path(out_dir, "C3_combination_balance_and_C3b_triggers.csv"))

input_paths <- c(covariates_path, background_paths, age_path, marker_path, spatial_m7_path,
                 file.path(frozen_dir, paste0(rep(set_names, each = 2), "_", rep(set_versions, 3), ".csv")))
input_audit <- data.table(
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv_utf8(input_audit, file.path(out_dir, "C3_input_provenance_sha256.csv"))

expected_tests <- nrow(signature_manifest) * length(set_names) * length(set_versions)
audit <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  status = "COMPLETE",
  frozen_sets = set_names,
  versions = set_versions,
  signature_count = nrow(signature_manifest),
  expected_tests = expected_tests,
  completed_tests = nrow(result),
  matching_failures = nrow(failure_table),
  permutations = n_permutations,
  seed = random_seed,
  full_background_counts = as.list(actual_background_n),
  coverage_rules = c("<5 descriptive overlap only", "5-9 sparse exploratory Fisher analysis",
                     ">=10 Fisher plus matched permutation"),
  non_mhc_rule = "MHC removed from both target sets and resource backgrounds",
  analysis_hierarchy = c("primary: sex-collapsed single-cell marker signatures and three prespecified spatial categories",
                         "secondary conditional: 61 reported age-DEG signatures",
                         "exploratory: 64 individual spatial tissues"),
  multiple_testing = "BH within prespecified correction family and set version; global BH across estimable full-set primary C3 tests",
  immune_age_limitation = "reported DEG universe only; conditional concentration, not unbiased detectable-transcriptome enrichment",
  external_evidence_label = "cross-resource validation / external orthogonal support; not blinded independent replication",
  max_post_match_mean_abs_smd = max(balance$post_match_mean_abs_smd, na.rm = TRUE),
  max_post_match_max_abs_smd = max(balance$post_match_max_abs_smd, na.rm = TRUE),
  combinations_requiring_uniform_C3b_adjustment = sum(combo_balance$requires_uniform_C3b_adjustment),
  frozen_gene_membership_changed = FALSE,
  d3_role = "D3 may summarize these C3 results after C3 QA; D3 must not re-score or re-rank frozen genes"
)
write_json(audit, file.path(out_dir, "C3_audit.json"), pretty = TRUE, auto_unbox = TRUE)

qa_bad_fisher <- result[coverage_tier != "descriptive_only_lt5" & is.na(fisher_p_greater), .N] > 0L
qa_bad_permutation <- result[coverage_tier == "inferential_ge10" & is.na(matched_permutation_p), .N] > 0L
qa_bad_descriptive <- result[coverage_tier == "descriptive_only_lt5" &
                               (!is.na(fisher_p_greater) | !is.na(matched_permutation_p)), .N] > 0L
if (nrow(failure_table) > 0L || nrow(result) != expected_tests ||
    qa_bad_fisher || qa_bad_permutation || qa_bad_descriptive) {
  stop("C3 failed QA. INCOMPLETE marker retained; inspect output diagnostics.")
}

write_json(
  list(stage = "C3_external_omics", status = "COMPLETE",
       completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
       signatures = nrow(signature_manifest), tests = nrow(result), failures = 0L),
  file.path(out_dir, "_STAGE_COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)
unlink(file.path(out_dir, "_STAGE_INCOMPLETE.json"))
if (any(combo_balance$requires_uniform_C3b_adjustment)) {
  write_json(
    list(status = "C3B_REQUIRED_BEFORE_D3",
         rule = "adjust every signature in each imbalanced resource-set-version combination; do not select only positive signatures",
         trigger_threshold = "maximum post-match mean absolute SMD >= 0.1",
         triggered_combinations = sum(combo_balance$requires_uniform_C3b_adjustment)),
    file.path(out_dir, "_C3B_REQUIRED.json"), pretty = TRUE, auto_unbox = TRUE
  )
}
message("Track C C3 completed: ", out_dir)
