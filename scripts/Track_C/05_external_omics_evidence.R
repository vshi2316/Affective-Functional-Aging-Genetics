## ============================================================================
## Track C D3: evidence audit after endpoint influence analysis and C3/C3b
## User-run only. This stage summarizes evidence; it never re-ranks frozen genes.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
trackc_root <- Sys.getenv("TRACKC_ROOT",
                      unset = file.path(base_dir, "analysis_ready_core", "trackC"))
c3_dir <- Sys.getenv("TRACKC_C3_DIR", unset = file.path(trackc_root, "06_C3_external_omics"))
c3b_dir <- Sys.getenv("TRACKC_C3B_DIR", unset = file.path(trackc_root, "07_C3b_uniform_adjustment"))
out_dir <- Sys.getenv("TRACKC_D3_OUT_DIR", unset = file.path(trackc_root, "08_D3_evidence_audit"))

trackc_root <- normalizePath(trackc_root, winslash = "/", mustWork = TRUE)
c3_dir <- normalizePath(c3_dir, winslash = "/", mustWork = TRUE)
c3b_dir <- normalizePath(c3b_dir, winslash = "/", mustWork = TRUE)
frozen_dir <- file.path(trackc_root, "00_frozen_gene_sets")
loco_dir <- file.path(trackc_root, "04_targeted_loco")
lolo_dir <- file.path(trackc_root, "05_targeted_lolo")

required_markers <- c(
  file.path(frozen_dir, "_STAGE_COMPLETE.json"),
  file.path(loco_dir, "_STAGE_COMPLETE.json"),
  file.path(lolo_dir, "_STAGE_COMPLETE.json"),
  file.path(c3_dir, "_STAGE_COMPLETE.json"),
  file.path(c3b_dir, "_STAGE_COMPLETE.json")
)
if (!all(file.exists(required_markers))) stop("Required upstream COMPLETE markers are missing.")
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Overwrite protection: output directory already contains files: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_json(list(stage = "D3_evidence_audit", status = "INCOMPLETE_OR_RUNNING",
                started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
           file.path(out_dir, "_STAGE_INCOMPLETE.json"), pretty = TRUE, auto_unbox = TRUE)

write_csv_utf8 <- function(x, path) fwrite(x, path, bom = TRUE, na = "")
sha256_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)

set_names <- c("cross_model_core", "affective_prioritized", "functional_prioritized")
read_set <- function(set_name, version) {
  unique(toupper(trimws(as.character(
    fread(file.path(frozen_dir, paste0(set_name, "_", version, ".csv")))$gene_symbol
  ))))
}

gene_universe <- data.table(gene_symbol = sort(unique(unlist(lapply(set_names, read_set, version = "full")))))
for (set_name in set_names) {
  gene_universe[, (paste0("in_", set_name)) := gene_symbol %in% read_set(set_name, "full")]
  gene_universe[, (paste0("in_", set_name, "_non_mhc")) := gene_symbol %in% read_set(set_name, "non_mhc")]
}

c3_path <- file.path(c3_dir, "C3_external_omics_enrichment_results.csv")
c3b_path <- file.path(c3b_dir, "C3b_adjusted_vs_C3_comparison.csv")
loco_path <- file.path(loco_dir, "phase1c_targeted_loco_summary.csv")
lolo_path <- file.path(lolo_dir, "phase1c_targeted_lolo_summary.csv")
c3 <- fread(c3_path)
c3b <- fread(c3b_path)
loco <- fread(loco_path, colClasses = list(character = c("phenocode", "most_influential_chromosome")))
lolo <- fread(lolo_path, colClasses = list(character = "phenocode"))

c3b[, all_three_family_significant :=
      adjusted_family_significant == TRUE & fisher_family_significant == TRUE &
      permutation_family_significant == TRUE]

primary <- c3b[analysis_layer %in% c("primary_single_cell", "primary_spatial_category")]
full_primary <- primary[set_version == "full"]
nonmhc_primary <- primary[set_version == "non_mhc",
  .(resource, signature_id, set_name,
    non_mhc_all_three_family_significant = all_three_family_significant,
    non_mhc_adjusted_or = adjusted_or,
    non_mhc_adjusted_bh = adjusted_bh_within_family_version)]
support_manifest <- merge(full_primary, nonmhc_primary,
                          by = c("resource", "signature_id", "set_name"),
                          all.x = TRUE, sort = FALSE)
support_manifest[, support_tier := fifelse(
  all_three_family_significant & adjusted_global_bh_primary_full < 0.05 &
    non_mhc_all_three_family_significant == TRUE,
  "tier1_global_primary_non_mhc_reinforced",
  fifelse(all_three_family_significant & non_mhc_all_three_family_significant == TRUE,
          "tier2_family_primary_non_mhc_reinforced",
          fifelse(all_three_family_significant,
                  "tier3_family_primary_full_only", "no_primary_support"))
)]
support_manifest[, manuscript_role := fifelse(
  support_tier == "tier1_global_primary_non_mhc_reinforced", "primary_external_spatial_support",
  fifelse(support_tier == "tier2_family_primary_non_mhc_reinforced", "secondary_external_spatial_support",
          fifelse(support_tier == "tier3_family_primary_full_only", "sensitivity_only", "not_supportive"))
)]
write_csv_utf8(support_manifest, file.path(out_dir, "D3_C3_primary_support_manifest.csv"))

## Map supported C3 signatures back to genes using the original frozen-set overlaps.
supported_keys <- support_manifest[support_tier != "no_primary_support",
  .(resource, signature_id, set_name, support_tier, manuscript_role)]
overlap_long <- list()
for (i in seq_len(nrow(supported_keys))) {
  key <- supported_keys[i]
  row <- c3[resource == key$resource & signature_id == key$signature_id &
              set_name == key$set_name & set_version == "full"]
  genes <- unlist(strsplit(row$overlap_genes, ";", fixed = TRUE))
  genes <- genes[nzchar(genes)]
  if (length(genes)) overlap_long[[length(overlap_long) + 1L]] <- data.table(
    gene_symbol = genes, set_name = key$set_name, resource = key$resource,
    signature_id = key$signature_id, signature_label = row$signature_label,
    support_tier = key$support_tier, manuscript_role = key$manuscript_role
  )
}
supported_gene_rows <- rbindlist(overlap_long, fill = TRUE)
if (!ncol(supported_gene_rows)) supported_gene_rows <- data.table(
  gene_symbol = character(), set_name = character(), resource = character(),
  signature_id = character(), signature_label = character(),
  support_tier = character(), manuscript_role = character())

gene_support <- supported_gene_rows[, .(
  supported_primary_signature_n = uniqueN(signature_id),
  supported_primary_signatures = paste(sort(unique(signature_label)), collapse = ";"),
  best_support_tier = if (any(support_tier == "tier1_global_primary_non_mhc_reinforced"))
    "tier1_global_primary_non_mhc_reinforced" else if
    (any(support_tier == "tier2_family_primary_non_mhc_reinforced"))
    "tier2_family_primary_non_mhc_reinforced" else "tier3_family_primary_full_only"
), by = gene_symbol]
gene_matrix <- merge(gene_universe, gene_support, by = "gene_symbol", all.x = TRUE)
gene_matrix[is.na(supported_primary_signature_n), `:=`(
  supported_primary_signature_n = 0L,
  supported_primary_signatures = "",
  best_support_tier = "no_primary_external_omics_support"
)]

## Descriptive overlaps are retained but cannot upgrade the support tier.
## c3 is signature-level; expand its overlap strings explicitly.
desc_rows <- c3[set_version == "full" & nzchar(overlap_genes),
  .(gene_symbol = unlist(strsplit(overlap_genes, ";", fixed = TRUE))),
  by = .(resource, signature_id)]
desc_summary <- desc_rows[nzchar(gene_symbol), .(
  descriptive_C3_signature_overlap_n = uniqueN(signature_id),
  descriptive_C3_resources = paste(sort(unique(resource)), collapse = ";")
), by = gene_symbol]
gene_matrix <- merge(gene_matrix, desc_summary, by = "gene_symbol", all.x = TRUE)
gene_matrix[is.na(descriptive_C3_signature_overlap_n), `:=`(
  descriptive_C3_signature_overlap_n = 0L, descriptive_C3_resources = ""
)]
setorder(gene_matrix, best_support_tier, gene_symbol)
write_csv_utf8(gene_matrix, file.path(out_dir, "D3_frozen_gene_external_omics_matrix.csv"))
write_csv_utf8(supported_gene_rows, file.path(out_dir, "D3_supported_gene_signature_long.csv"))

## Endpoint influence classification remains set/endpoint-level evidence.
triggered <- lolo[, .(target_id,
  lolo_class = fifelse(any_direction_reversal,
                       "dominant_locus_dependent_association",
                       "locus_concentrated_association"),
  most_influential_locus,
  maximum_lolo_log_or_attenuation = maximum_relative_log_or_attenuation)]
endpoint_audit <- merge(loco, triggered, by = "target_id", all.x = TRUE)
endpoint_audit[, endpoint_influence_class := fifelse(
  !is.na(lolo_class), lolo_class,
  fifelse(maximum_relative_log_or_attenuation >= 0.4,
          "moderately_chromosome_sensitive_association",
          "relatively_stable_no_major_single_chromosome_influence"))
]
write_csv_utf8(endpoint_audit, file.path(out_dir, "D3_endpoint_influence_classification.csv"))

set_summary <- rbindlist(lapply(set_names, function(sn) {
  genes <- read_set(sn, "full")
  nonmhc <- read_set(sn, "non_mhc")
  gm <- gene_matrix[gene_symbol %in% genes]
  data.table(
    set_name = sn, full_gene_n = length(genes), non_mhc_gene_n = length(nonmhc),
    tier1_supported_gene_n = sum(gm$best_support_tier == "tier1_global_primary_non_mhc_reinforced"),
    tier2_supported_gene_n = sum(gm$best_support_tier == "tier2_family_primary_non_mhc_reinforced"),
    tier3_supported_gene_n = sum(gm$best_support_tier == "tier3_family_primary_full_only"),
    any_primary_supported_gene_n = sum(gm$supported_primary_signature_n > 0),
    primary_single_cell_supported_signature_n = support_manifest[
      set_name == sn & analysis_layer == "primary_single_cell" &
        support_tier != "no_primary_support", .N],
    primary_spatial_supported_signature_n = support_manifest[
      set_name == sn & analysis_layer == "primary_spatial_category" &
        support_tier != "no_primary_support", .N]
  )
}))
write_csv_utf8(set_summary, file.path(out_dir, "D3_set_level_evidence_summary.csv"))

input_paths <- c(c3_path, c3b_path, loco_path, lolo_path,
                 file.path(frozen_dir, paste0(rep(set_names, each = 2), "_",
                                              rep(c("full", "non_mhc"), 3), ".csv")))
input_audit <- data.table(
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv_utf8(input_audit, file.path(out_dir, "D3_input_provenance_sha256.csv"))

audit <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  role = "evidence audit and classification only",
  frozen_unique_gene_n = nrow(gene_universe),
  primary_support_rule = "Fisher, matched permutation and uniform Firth adjustment all family-significant",
  tier1_rule = "primary full-set signal also passes adjusted global BH and is non-MHC reinforced",
  tier2_rule = "primary family-level signal is non-MHC reinforced but not globally significant",
  tier3_rule = "primary family-level full-set signal without non-MHC reinforcement",
  single_cell_conclusion = "no primary adjusted single-cell marker support",
  age_deg_conclusion = "no adjusted reported-DEG conditional support",
  spatial_conclusion = "limited category-level support; not broad neuroimmune validation",
  gene_reranking_performed = FALSE
)
write_json(audit, file.path(out_dir, "D3_audit.json"), pretty = TRUE, auto_unbox = TRUE)

required_outputs <- c(
  "D3_C3_primary_support_manifest.csv",
  "D3_frozen_gene_external_omics_matrix.csv",
  "D3_endpoint_influence_classification.csv",
  "D3_set_level_evidence_summary.csv"
)
if (!all(file.exists(file.path(out_dir, required_outputs))) || nrow(gene_matrix) != nrow(gene_universe)) {
  stop("D3 failed QA. INCOMPLETE marker retained.")
}
write_json(list(stage = "D3_evidence_audit", status = "COMPLETE",
                completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                frozen_unique_genes = nrow(gene_universe), gene_reranking = FALSE),
           file.path(out_dir, "_STAGE_COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE)
unlink(file.path(out_dir, "_STAGE_INCOMPLETE.json"))
message("Track C D3 evidence audit completed: ", out_dir)


