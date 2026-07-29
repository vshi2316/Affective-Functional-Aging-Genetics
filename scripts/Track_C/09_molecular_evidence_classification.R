options(stringsAsFactors = FALSE)

message("Track C final QTL closure and prespecified-rule reclassification")

required_packages <- c("data.table", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(data.table))

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
trackc_root <- Sys.getenv(
  "TRACKC_ROOT",
  unset = file.path(base_dir, "analysis_ready_core", "trackC")
)
qtl_root <- Sys.getenv(
  "TRACKC_QTL_ROOT",
  unset = file.path(trackc_root, "09_QTL_molecular_followup")
)
primary_dir <- Sys.getenv(
  "TRACKC_FINAL_QTL_PRIMARY_DIR",
  unset = file.path(qtl_root, "07_metabrain_three_locus_primary")
)
susie_dir <- Sys.getenv(
  "TRACKC_FINAL_QTL_SUSIE_DIR",
  unset = file.path(qtl_root, "08_frozen_three_locus_susie")
)
out_dir <- Sys.getenv(
  "TRACKC_FINAL_QTL_CLOSURE_OUT_DIR",
  unset = file.path(qtl_root, "09_final_QTL_closure_reclassification")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

require_file <- function(path, label = "file") {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

safe_ratio <- function(h4, h3) {
  h4 <- as.numeric(h4)
  h3 <- as.numeric(h3)
  fifelse(is.finite(h4) & is.finite(h3), h4 / pmax(h3, 1e-300), NA_real_)
}

require_file(
  file.path(primary_dir, "05_integrated_primary_evidence", "COMPLETE.txt"),
  "primary QTL completion marker"
)
require_file(file.path(susie_dir, "COMPLETE.txt"), "SuSiE completion marker")

frozen_spec_file <- require_file(
  file.path(primary_dir, "01_frozen_analysis_spec", "README.txt"),
  "frozen analysis specification"
)
frozen_spec_text <- readLines(frozen_spec_file, warn = FALSE, encoding = "UTF-8")
if (!any(grepl("PP.H4>=0.80 and H4/H3>=3", frozen_spec_text, fixed = TRUE))) {
  stop("Frozen PP.H4/H3 rule not found in specification")
}

susie <- fread(require_file(
  file.path(susie_dir, "all_frozen_susie_results.csv"),
  "all frozen SuSiE results"
))
trigger <- fread(require_file(
  file.path(susie_dir, "frozen_trigger_manifest_with_short_run_ids.csv"),
  "frozen trigger manifest"
))
candidate <- fread(require_file(
  file.path(
    primary_dir, "05_integrated_primary_evidence",
    "MetaBrain_candidate_summary.csv"
  ),
  "MetaBrain candidate summary"
))
failure <- fread(require_file(
  file.path(susie_dir, "frozen_susie_failures.csv"),
  "SuSiE failure audit"
))

if (nrow(trigger) != 91L) stop("Expected 91 frozen triggers; observed ", nrow(trigger))
if (nrow(susie) != 91L) stop("Expected 91 final SuSiE classifications; observed ", nrow(susie))
if (nrow(failure) != 0L) stop("Unresolved processing failures remain")
if (anyDuplicated(susie$run_id) || anyDuplicated(trigger$run_id)) stop("Duplicate run_id")
if (!setequal(susie$run_id, trigger$run_id)) stop("SuSiE results do not match frozen triggers")

# Match the two priors by the actual signal identifiers, never by row order.
pair_rows <- list()
for (i in seq_len(nrow(susie))) {
  meta <- susie[i]
  if (meta$susie_status != "completed") next
  run_dir <- file.path(susie_dir, "runs", meta$run_id)
  primary_file <- require_file(
    file.path(run_dir, "coloc_susie_primary_signal_pairs.csv"),
    paste0(meta$run_id, " primary signal pairs")
  )
  sensitivity_file <- require_file(
    file.path(run_dir, "coloc_susie_sensitivity_signal_pairs.csv"),
    paste0(meta$run_id, " sensitivity signal pairs")
  )
  p <- fread(primary_file)
  s <- fread(sensitivity_file)
  # A completed SuSiE fit can legitimately yield no credible signal pair
  # (for example, when one side has zero credible sets).  data.table writes
  # such an object as a zero-row, zero-column CSV.  Treat two empty files as
  # an analysed combination with no colocalising signal pair; only one-sided
  # emptiness is a primary/sensitivity inconsistency.
  if (!nrow(p) && !nrow(s)) next
  if (!nrow(p) || !nrow(s)) {
    stop("Primary/sensitivity signal-pair emptiness mismatch for ", meta$run_id)
  }
  required_pair <- c("hit1", "hit2", "PP.H3.abf", "PP.H4.abf")
  if (!all(required_pair %in% names(p)) || !all(required_pair %in% names(s))) {
    stop("Signal-pair columns missing for ", meta$run_id)
  }
  if (anyDuplicated(p, by = c("hit1", "hit2")) || anyDuplicated(s, by = c("hit1", "hit2"))) {
    stop("Duplicate hit1+hit2 signal pair for ", meta$run_id)
  }
  p <- p[, .(
    hit1, hit2,
    primary_h3 = as.numeric(PP.H3.abf),
    primary_h4 = as.numeric(PP.H4.abf)
  )]
  s <- s[, .(
    hit1, hit2,
    sensitivity_h3 = as.numeric(PP.H3.abf),
    sensitivity_h4 = as.numeric(PP.H4.abf)
  )]
  z <- merge(p, s, by = c("hit1", "hit2"), all = TRUE, sort = FALSE)
  if (nrow(z) != nrow(p) || nrow(z) != nrow(s) || anyNA(z)) {
    stop("Primary/sensitivity signal-pair mismatch for ", meta$run_id)
  }
  z[, `:=`(
    run_id = meta$run_id,
    combo_id = meta$combo_id,
    locus_label = meta$locus_label,
    factor_gwas = meta$factor_gwas,
    gene_symbol = meta$gene_symbol,
    ensembl_id = meta$ensembl_id,
    resource = meta$resource,
    tissue = meta$tissue
  )]
  z[, `:=`(
    primary_h4_h3_ratio = safe_ratio(primary_h4, primary_h3),
    sensitivity_h4_h3_ratio = safe_ratio(sensitivity_h4, sensitivity_h3)
  )]
  z[, `:=`(
    default_prior_positive = primary_h4 >= 0.80 & primary_h4_h3_ratio >= 3,
    strict_prior_stable = sensitivity_h4 >= 0.80 & sensitivity_h4_h3_ratio >= 3
  )]
  setcolorder(z, c(
    "run_id", "combo_id", "locus_label", "factor_gwas", "gene_symbol",
    "ensembl_id", "resource", "tissue", "hit1", "hit2",
    "primary_h3", "primary_h4", "primary_h4_h3_ratio",
    "sensitivity_h3", "sensitivity_h4", "sensitivity_h4_h3_ratio",
    "default_prior_positive", "strict_prior_stable"
  ))
  pair_rows[[length(pair_rows) + 1L]] <- z
}
pair_dt <- rbindlist(pair_rows, fill = TRUE)
fwrite(pair_dt, file.path(out_dir, "all_susie_signal_pairs_reclassified.csv"))

default_pairs <- pair_dt[default_prior_positive == TRUE]
strict_pairs <- pair_dt[strict_prior_stable == TRUE]
fwrite(default_pairs, file.path(out_dir, "default_prior_positive_signal_pairs.csv"))
fwrite(strict_pairs, file.path(out_dir, "strict_prior_stable_signal_pairs.csv"))

default_combo <- unique(default_pairs$run_id)
strict_combo <- unique(strict_pairs$run_id)
susie[, `:=`(
  original_abf_h4_h3_ratio = safe_ratio(original_abf_h4, original_abf_h3),
  default_positive_signal_pair_n = vapply(
    run_id, function(x) nrow(default_pairs[run_id == x]), integer(1)
  ),
  strict_stable_signal_pair_n = vapply(
    run_id, function(x) nrow(strict_pairs[run_id == x]), integer(1)
  )
)]
susie[, final_evidence_class := fcase(
  run_id %in% strict_combo & resource == "MetaBrain_2021_07_23" & tissue == "cortex_EUR",
    "prior_stable_brain_colocalization",
  run_id %in% default_combo,
    "primary_prior_supported_sensitivity_attenuated",
  susie_status != "completed" & original_abf_h4 >= 0.80 & original_abf_h4_h3_ratio >= 3,
    "single_signal_H4_support_multisignal_not_estimable",
  susie_status != "completed",
    "single_signal_predominantly_H3_multisignal_not_estimable",
  pmax(original_abf_h3, susie_max_h3, na.rm = TRUE) >= 0.80,
    "predominantly_distinct_signals",
  default = "molecularly_unresolved_after_susie"
)]
fwrite(susie, file.path(out_dir, "all_91_combinations_final_reclassification.csv"))

not_estimable <- susie[susie_status != "completed", .(
  run_id, combo_id, locus_label, factor_gwas, gene_symbol, ensembl_id,
  resource, tissue, original_abf_h3, original_abf_h4,
  original_abf_h4_h3_ratio, susie_status, final_evidence_class,
  permitted_interpretation = fifelse(
    final_evidence_class == "single_signal_H4_support_multisignal_not_estimable",
    "single-signal ABF H4 support; multi-signal analysis not estimable with the frozen external 1000 Genomes European LD reference",
    "single-signal result predominantly supports H3; multi-signal analysis not estimable with the frozen external 1000 Genomes European LD reference"
  )
)]
fwrite(not_estimable, file.path(out_dir, "multisignal_not_estimable_manifest.csv"))

gene_summary <- default_pairs[, .(
  default_positive_combination_n = uniqueN(run_id),
  default_positive_signal_pair_n = .N,
  strict_stable_combination_n = uniqueN(run_id[strict_prior_stable == TRUE]),
  strict_stable_signal_pair_n = sum(strict_prior_stable),
  resources = paste(sort(unique(resource)), collapse = ";"),
  tissues = paste(sort(unique(tissue)), collapse = ";"),
  strongest_primary_h4 = max(primary_h4),
  strongest_sensitivity_h4 = max(sensitivity_h4),
  final_gene_tier = if (any(strict_prior_stable & resource == "MetaBrain_2021_07_23" & tissue == "cortex_EUR")) {
    "prior_stable_brain_colocalization"
  } else "primary_prior_supported_sensitivity_attenuated"
), by = .(locus_label, factor_gwas, gene_symbol, ensembl_id)]
setorder(gene_summary, locus_label, final_gene_tier, gene_symbol)
fwrite(gene_summary, file.path(out_dir, "six_gene_three_locus_summary.csv"))

locus_summary <- default_pairs[, .(
  default_positive_combination_n = uniqueN(run_id),
  default_positive_signal_pair_n = .N,
  default_positive_gene_n = uniqueN(gene_symbol),
  strict_stable_combination_n = uniqueN(run_id[strict_prior_stable == TRUE]),
  strict_stable_signal_pair_n = sum(strict_prior_stable),
  strict_stable_gene_n = uniqueN(gene_symbol[strict_prior_stable == TRUE])
), by = .(locus_label, factor_gwas)]
fwrite(locus_summary, file.path(out_dir, "three_locus_final_summary.csv"))

candidate_final <- unique(candidate[, .(
  factor_gwas, gene_symbol, ensembl_id, resource, tissue,
  min_qtl_p, strong_cis_qtl, evidence_class,
  final_candidate_class = "no_candidate_level_MetaBrain_anchor"
)])
fwrite(candidate_final, file.path(out_dir, "candidate_MetaBrain_no_anchor_summary.csv"))

expected_strict <- c(
  "ARHGAP1|MetaBrain_2021_07_23|cortex_EUR",
  "NPIPB7|MetaBrain_2021_07_23|cortex_EUR"
)
observed_strict <- unique(paste(
  strict_pairs$gene_symbol, strict_pairs$resource, strict_pairs$tissue,
  sep = "|"
))
if (!setequal(expected_strict, observed_strict)) {
  stop("Strict-prior stable signal set differs from the audited expected set")
}
if (nrow(default_pairs) != 10L || uniqueN(default_pairs$run_id) != 9L ||
    uniqueN(default_pairs$gene_symbol) != 6L || uniqueN(default_pairs$locus_label) != 3L) {
  stop("Default-prior positive counts do not match audited results")
}
if (nrow(strict_pairs) != 2L || uniqueN(strict_pairs$run_id) != 2L) {
  stop("Strict-prior stable counts do not match audited results")
}
if (nrow(not_estimable) != 5L) stop("Expected five multi-signal not-estimable combinations")

audit <- list(
  stage = "final_QTL_closure_reclassification",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  statistical_models_recomputed = FALSE,
  source_susie_trigger_n = nrow(susie),
  susie_completed_n = sum(susie$susie_status == "completed"),
  multisignal_not_estimable_n = nrow(not_estimable),
  processing_failure_n = nrow(failure),
  default_prior_rule = "PP.H4>=0.80 and PP.H4/PP.H3>=3",
  default_prior_positive_signal_pair_n = nrow(default_pairs),
  default_prior_positive_combination_n = uniqueN(default_pairs$run_id),
  default_prior_positive_gene_n = uniqueN(default_pairs$gene_symbol),
  default_prior_positive_locus_n = uniqueN(default_pairs$locus_label),
  sensitivity_rule = "p12=1e-6: PP.H4>=0.80 and PP.H4/PP.H3>=3 on the same hit1+hit2 signal pair",
  strict_prior_stable_signal_pair_n = nrow(strict_pairs),
  strict_prior_stable_combination_n = uniqueN(strict_pairs$run_id),
  strict_prior_stable_genes = sort(unique(strict_pairs$gene_symbol)),
  regional_SMR_HEIDI_performed = FALSE,
  strict_regulatory_anchor_claim_authorised = FALSE,
  permitted_main_term = "prior-stable brain eQTL colocalization signal",
  qtl_stage_closed = TRUE,
  closure_reason = paste(
    "All three post-LOLO frozen loci were classified under the prespecified rule;",
    "no new gene, tissue, locus, resource, threshold or model was added."
  ),
  prohibited_claims = c(
    "expression mediation", "causal gene", "strict regulatory anchor",
    "independent factor-GWAS replication", "confirmed LD mismatch"
  )
)
jsonlite::write_json(
  audit, file.path(out_dir, "final_QTL_closure_audit.json"),
  pretty = TRUE, auto_unbox = TRUE
)

writeLines(
  c(
    "Final bounded Track C QTL closure.",
    "No ABF or SuSiE model was recomputed.",
    "Primary and p12=1e-6 results were matched by hit1+hit2.",
    "Default rule: PP.H4>=0.80 and PP.H4/PP.H3>=3.",
    "Prior-stable brain signals: ARHGAP1-MetaBrain cortex and NPIPB7-MetaBrain cortex.",
    "The other default-prior positive signals attenuated under p12=1e-6.",
    "Five combinations were not estimable with the frozen external 1000 Genomes European LD reference.",
    "Regional SMR/HEIDI was not performed; strict regulatory-anchor terminology is not authorised.",
    "Track C QTL analysis is closed. No additional gene, tissue, locus, resource or model is authorised."
  ),
  file.path(out_dir, "README.txt"),
  useBytes = TRUE
)
writeLines(
  c(
    "stage=final_QTL_closure_reclassification",
    paste0("completed_at=", audit$completed_at),
    "statistical_models_recomputed=false",
    paste0("default_positive_signal_pair_n=", nrow(default_pairs)),
    paste0("default_positive_combination_n=", uniqueN(default_pairs$run_id)),
    paste0("strict_prior_stable_signal_pair_n=", nrow(strict_pairs)),
    paste0("multisignal_not_estimable_n=", nrow(not_estimable)),
    "qtl_stage_closed=true"
  ),
  file.path(out_dir, "COMPLETE.txt"),
  useBytes = TRUE
)

message("Completed final QTL closure: ", out_dir)
message("Default-prior positive: 10 signal pairs / 9 combinations / 6 genes / 3 loci")
message("Strict-prior stable: ARHGAP1-MetaBrain cortex and NPIPB7-MetaBrain cortex")
message("Track C QTL stage closed: TRUE")


