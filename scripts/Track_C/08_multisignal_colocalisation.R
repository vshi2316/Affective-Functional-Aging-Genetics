options(stringsAsFactors = FALSE)

message("Track C frozen three-locus SuSiE diagnostic")

required_packages <- c("data.table", "bigsnpr", "coloc", "susieR", "jsonlite")
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
out_dir <- Sys.getenv(
  "TRACKC_FINAL_QTL_SUSIE_OUT_DIR",
  unset = file.path(qtl_root, "08_frozen_three_locus_susie")
)
bfile_prefix <- Sys.getenv(
  "TRACKC_LD_BFILE",
  unset = file.path(base_dir, "reference", "plink", "EUR")
)
ascii_work <- Sys.getenv(
  "TRACKC_FINAL_QTL_SUSIE_ASCII_WORK",
  unset = file.path(tempdir(), "trackC_multisignal_colocalisation")
)
maxit <- as.integer(Sys.getenv("TRACKC_SUSIE_MAXIT", unset = "1000"))
if (!is.finite(maxit) || maxit < 100L) stop("TRACKC_SUSIE_MAXIT must be >=100")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ascii_work, recursive = TRUE, showWarnings = FALSE)
run_root <- file.path(out_dir, "runs")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

require_file <- function(path, label = "file") {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

empty_failure <- function() {
  data.table(run_id = character(), combo_id = character(), error_message = character())
}

trigger_file <- require_file(
  file.path(
    primary_dir, "05_integrated_primary_evidence",
    "frozen_susie_trigger_manifest.csv"
  ),
  "frozen SuSiE trigger manifest"
)
primary_complete <- require_file(
  file.path(primary_dir, "05_integrated_primary_evidence", "COMPLETE.txt"),
  "primary integration completion marker"
)
trigger <- fread(trigger_file)
required_trigger <- c(
  "combo_id", "analysis_role", "locus_label", "factor_gwas", "gene_symbol",
  "ensembl_id", "resource", "tissue", "harmonized_file", "nsnps",
  "pp_h3", "pp_h4", "pp_h4_p12_1e6", "trigger_reason"
)
missing_trigger <- setdiff(required_trigger, names(trigger))
if (length(missing_trigger)) {
  stop("Trigger manifest missing columns: ", paste(missing_trigger, collapse = ", "))
}
if (nrow(trigger) != 91L) {
  stop("Frozen trigger count changed: expected 91, observed ", nrow(trigger))
}
if (anyDuplicated(trigger$combo_id)) stop("Duplicate combo_id in frozen trigger manifest")
if (any(trigger$analysis_role != "locus")) stop("Non-locus item found in frozen trigger manifest")
if (!all(trigger$locus_label %in% c(
  "chr3_functional", "chr11_affective", "chr16_functional"
))) stop("Unexpected locus in frozen trigger manifest")
if (any(trigger$nsnps < 100L)) stop("Frozen trigger with fewer than 100 SNPs")
if (any(!(trigger$pp_h3 >= 0.80 | trigger$pp_h4 >= 0.80))) {
  stop("Manifest contains an item that does not meet the frozen H3/H4 trigger")
}
trigger[, run_id := sprintf("S%03d", .I)]
fwrite(trigger, file.path(out_dir, "frozen_trigger_manifest_with_short_run_ids.csv"))

harmonized_required <- c(
  "SNP", "qtl_BP", "aligned_effect_allele", "aligned_other_allele",
  "eaf_qtl", "maf_qtl", "beta_gwas_aligned", "se_gwas", "p_gwas",
  "n_gwas", "beta_qtl", "se_qtl", "p_qtl", "n_qtl"
)

# Build the exact union of SNPs authorised by the frozen trigger files.
trigger_snps <- vector("list", nrow(trigger))
for (i in seq_len(nrow(trigger))) {
  f <- require_file(trigger$harmonized_file[i], paste0("harmonized input ", trigger$run_id[i]))
  h0 <- fread(f, select = "SNP")
  trigger_snps[[i]] <- unique(as.character(h0$SNP))
}
wanted_snps <- unique(unlist(trigger_snps, use.names = FALSE))

bim_file <- require_file(paste0(bfile_prefix, ".bim"), "1000G EUR BIM")
bed_file <- require_file(paste0(bfile_prefix, ".bed"), "1000G EUR BED")
require_file(paste0(bfile_prefix, ".fam"), "1000G EUR FAM")
bim <- fread(
  bim_file,
  col.names = c("CHR", "SNP", "CM", "BP_LD", "BIM_A1", "BIM_A2")
)
bim[, `:=`(
  SNP = as.character(SNP),
  BIM_A1 = toupper(as.character(BIM_A1)),
  BIM_A2 = toupper(as.character(BIM_A2))
)]
bim[, bim_row := .I]
bim <- bim[!duplicated(SNP)]
union_map <- bim[SNP %in% wanted_snps]
setorder(union_map, bim_row)
if (nrow(union_map) < 100L) stop("Fewer than 100 frozen SNPs occur in 1000G EUR")

# Read the union once.  This avoids 91 repeated BED scans while preserving each
# trigger's own LD matrix and allele-frequency QC downstream.
backing <- file.path(ascii_work, "frozen_trigger_union_1000G_EUR")
unlink(paste0(backing, c(".bk", ".rds")), force = TRUE)
bigsnpr::snp_readBed2(
  bedfile = bed_file,
  backingfile = backing,
  ind.col = union_map$bim_row,
  ncores = 1
)
obj <- bigsnpr::snp_attach(paste0(backing, ".rds"))
geno_union <- as.matrix(obj$genotypes[, ])
storage.mode(geno_union) <- "numeric"
geno_union[geno_union == 3] <- NA_real_
for (j in seq_len(ncol(geno_union))) {
  miss <- is.na(geno_union[, j])
  if (any(miss)) {
    mu <- mean(geno_union[, j], na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    geno_union[miss, j] <- mu
  }
}
union_map[, union_col := .I]
union_col_by_snp <- setNames(union_map$union_col, union_map$SNP)

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(!is.finite(x))) return(NA_real_)
  max(x[is.finite(x)])
}

write_ld_mismatch_result <- function(meta, error_message) {
  run_dir <- file.path(run_root, meta$run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  qa_file <- file.path(run_dir, "frequency_LD_QA.csv")
  qa <- if (file.exists(qa_file)) fread(qa_file) else data.table()
  ld_n <- if (nrow(qa) && "ld_reference_snp_n" %in% names(qa)) {
    as.integer(qa$ld_reference_snp_n[1L])
  } else NA_integer_
  retained_n <- if (nrow(qa) && "retained_for_susie_n" %in% names(qa)) {
    as.integer(qa$retained_for_susie_n[1L])
  } else NA_integer_
  ans <- data.table(
    run_id = meta$run_id, combo_id = meta$combo_id,
    locus_label = meta$locus_label, factor_gwas = meta$factor_gwas,
    gene_symbol = meta$gene_symbol, ensembl_id = meta$ensembl_id,
    resource = meta$resource, tissue = meta$tissue,
    original_abf_h3 = meta$pp_h3, original_abf_h4 = meta$pp_h4,
    harmonized_snp_n = as.integer(meta$nsnps),
    ld_reference_snp_n = ld_n, retained_snp_n = retained_n,
    robust_shared_signal_pair_n = 0L,
    susie_status = "not_estimable_external_LD_mismatch",
    decision = "molecularly_unresolved_due_external_LD_mismatch"
  )
  fwrite(ans, file.path(run_dir, "susie_summary.csv"))
  writeLines(
    c(
      "SuSiE was not estimable with the frozen external 1000G EUR LD panel.",
      "The single-signal ABF result remains reportable but cannot be upgraded to a multi-signal molecular anchor.",
      "The LD panel, thresholds, genes, tissues and resources were not changed.",
      paste0("Original error: ", gsub("[\r\n]+", " ", error_message))
    ),
    file.path(run_dir, "LD_MISMATCH_NOT_ESTIMABLE.txt")
  )
  writeLines(
    c(
      "stage=frozen_susie",
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      "status=not_estimable_external_LD_mismatch",
      "decision=molecularly_unresolved_due_external_LD_mismatch"
    ),
    file.path(run_dir, "COMPLETE.txt")
  )
  ans
}

# A resumed run may already contain the exact SuSiE diagnostic failures.  Reuse
# that audit evidence so that known LD-mismatch combinations are classified
# without repeating the futile automatic 100,000-iteration retry.
prior_failure_file <- file.path(out_dir, "frozen_susie_failures.csv")
known_ld_mismatch <- data.table()
if (file.exists(prior_failure_file)) {
  prior_failure <- fread(prior_failure_file)
  if (all(c("combo_id", "error_message") %in% names(prior_failure))) {
    known_ld_mismatch <- prior_failure[
      grepl("estimated prior variance is unreasonably large", error_message, fixed = TRUE)
    ]
  }
}

run_one <- function(meta) {
  run_dir <- file.path(run_root, meta$run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  summary_file <- file.path(run_dir, "susie_summary.csv")
  complete_file <- file.path(run_dir, "COMPLETE.txt")
  if (file.exists(summary_file) && file.exists(complete_file)) {
    old <- fread(summary_file)
    if (nrow(old) == 1L && identical(old$combo_id, meta$combo_id)) return(old)
  }

  h <- fread(require_file(meta$harmonized_file, paste0("harmonized input ", meta$run_id)))
  missing_h <- setdiff(harmonized_required, names(h))
  if (length(missing_h)) {
    stop("Harmonized input missing columns: ", paste(missing_h, collapse = ", "))
  }
  h <- h[!duplicated(SNP)]
  cols <- unname(union_col_by_snp[as.character(h$SNP)])
  in_ld <- is.finite(cols)
  h_ld <- h[in_ld]
  cols <- as.integer(cols[in_ld])
  map <- union_map[cols]
  geno <- geno_union[, cols, drop = FALSE]
  if (nrow(h_ld) < 100L) {
    ans <- data.table(
      run_id = meta$run_id, combo_id = meta$combo_id,
      locus_label = meta$locus_label, factor_gwas = meta$factor_gwas,
      gene_symbol = meta$gene_symbol, ensembl_id = meta$ensembl_id,
      resource = meta$resource, tissue = meta$tissue,
      original_abf_h3 = meta$pp_h3, original_abf_h4 = meta$pp_h4,
      harmonized_snp_n = nrow(h), ld_reference_snp_n = nrow(h_ld),
      retained_snp_n = 0L, susie_status = "insufficient_LD_overlap",
      robust_shared_signal_pair_n = 0L,
      decision = "not_estimable_after_frozen_LD_QC"
    )
    fwrite(ans, summary_file)
    writeLines(c("stage=frozen_susie", "status=not_estimable"), complete_file)
    return(ans)
  }

  dosage_freq <- colMeans(geno) / 2
  ea <- toupper(h_ld$aligned_effect_allele)
  oa <- toupper(h_ld$aligned_other_allele)
  same_a1 <- map$BIM_A1 == ea & map$BIM_A2 == oa
  same_a2 <- map$BIM_A2 == ea & map$BIM_A1 == oa
  allele_match <- same_a1 | same_a2
  expected_a1 <- ifelse(same_a1, h_ld$eaf_qtl, ifelse(same_a2, 1 - h_ld$eaf_qtl, NA_real_))
  expected_a2 <- ifelse(same_a2, h_ld$eaf_qtl, ifelse(same_a1, 1 - h_ld$eaf_qtl, NA_real_))
  diff_a1 <- median(abs(dosage_freq - expected_a1), na.rm = TRUE)
  diff_a2 <- median(abs(dosage_freq - expected_a2), na.rm = TRUE)
  counted <- if (is.finite(diff_a1) && (!is.finite(diff_a2) || diff_a1 <= diff_a2)) {
    "BIM_A1"
  } else "BIM_A2"
  expected_freq <- if (counted == "BIM_A1") expected_a1 else expected_a2
  freq_diff <- abs(dosage_freq - expected_freq)
  nonconstant <- apply(geno, 2L, stats::sd) > 0
  keep <- allele_match & is.finite(freq_diff) & freq_diff <= 0.20 & nonconstant

  qa <- data.table(
    run_id = meta$run_id, combo_id = meta$combo_id,
    harmonized_snp_n = nrow(h), ld_reference_snp_n = nrow(h_ld),
    exact_allele_match_n = sum(allele_match),
    frequency_difference_gt_0p20_n = sum(
      allele_match & is.finite(freq_diff) & freq_diff > 0.20
    ),
    retained_for_susie_n = sum(keep), dosage_counted_allele = counted,
    median_frequency_difference_if_BIM_A1_counted = diff_a1,
    median_frequency_difference_if_BIM_A2_counted = diff_a2
  )
  fwrite(qa, file.path(run_dir, "frequency_LD_QA.csv"))
  if (sum(keep) < 100L) {
    ans <- data.table(
      run_id = meta$run_id, combo_id = meta$combo_id,
      locus_label = meta$locus_label, factor_gwas = meta$factor_gwas,
      gene_symbol = meta$gene_symbol, ensembl_id = meta$ensembl_id,
      resource = meta$resource, tissue = meta$tissue,
      original_abf_h3 = meta$pp_h3, original_abf_h4 = meta$pp_h4,
      harmonized_snp_n = nrow(h), ld_reference_snp_n = nrow(h_ld),
      retained_snp_n = sum(keep), susie_status = "insufficient_after_frequency_LD_QC",
      robust_shared_signal_pair_n = 0L,
      decision = "not_estimable_after_frozen_LD_QC"
    )
    fwrite(ans, summary_file)
    writeLines(c("stage=frozen_susie", "status=not_estimable"), complete_file)
    return(ans)
  }

  h_use <- h_ld[keep]
  map_use <- map[keep]
  geno_use <- geno[, keep, drop = FALSE]
  counted_allele <- if (counted == "BIM_A1") map_use$BIM_A1 else map_use$BIM_A2
  flip <- counted_allele != toupper(h_use$aligned_effect_allele)
  h_use[, `:=`(
    beta_gwas_ld = ifelse(flip, -beta_gwas_aligned, beta_gwas_aligned),
    beta_qtl_ld = ifelse(flip, -beta_qtl, beta_qtl)
  )]
  ld <- stats::cor(geno_use)
  ld[!is.finite(ld)] <- 0
  diag(ld) <- 1
  rownames(ld) <- h_use$SNP
  colnames(ld) <- h_use$SNP

  qtl_sdY <- sqrt(stats::median(
    2 * h_use$n_qtl * h_use$maf_qtl * (1 - h_use$maf_qtl) * h_use$se_qtl^2,
    na.rm = TRUE
  ))
  if (!is.finite(qtl_sdY) || qtl_sdY <= 0) qtl_sdY <- 1
  d_gwas <- list(
    beta = h_use$beta_gwas_ld, varbeta = h_use$se_gwas^2,
    snp = h_use$SNP, position = as.integer(h_use$qtl_BP),
    MAF = h_use$maf_qtl, N = as.integer(round(median(h_use$n_gwas))),
    sdY = 1, type = "quant", LD = ld
  )
  d_qtl <- list(
    beta = h_use$beta_qtl_ld, varbeta = h_use$se_qtl^2,
    snp = h_use$SNP, position = as.integer(h_use$qtl_BP),
    MAF = h_use$maf_qtl, N = as.integer(round(median(h_use$n_qtl))),
    sdY = qtl_sdY, type = "quant", LD = ld
  )

  filtered_abf <- coloc::coloc.abf(
    within(d_gwas, rm(LD)), within(d_qtl, rm(LD)),
    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
  )
  susie_gwas <- coloc::runsusie(
    d_gwas, maxit = maxit, repeat_until_convergence = TRUE
  )
  susie_qtl <- coloc::runsusie(
    d_qtl, maxit = maxit, repeat_until_convergence = TRUE
  )
  coloc_primary <- coloc::coloc.susie(
    susie_gwas, susie_qtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
  )
  coloc_sensitivity <- coloc::coloc.susie(
    susie_gwas, susie_qtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-6
  )
  saveRDS(filtered_abf, file.path(run_dir, "frequency_QC_filtered_ABF.rds"))
  saveRDS(susie_gwas, file.path(run_dir, "susie_factor_GWAS.rds"))
  saveRDS(susie_qtl, file.path(run_dir, "susie_QTL.rds"))
  saveRDS(coloc_primary, file.path(run_dir, "coloc_susie_primary.rds"))
  saveRDS(coloc_sensitivity, file.path(run_dir, "coloc_susie_p12_1e6.rds"))

  ptab <- as.data.table(coloc_primary$summary)
  stab <- as.data.table(coloc_sensitivity$summary)
  ptab[, signal_pair_id := .I]
  stab[, signal_pair_id := .I]
  fwrite(ptab, file.path(run_dir, "coloc_susie_primary_signal_pairs.csv"))
  fwrite(stab, file.path(run_dir, "coloc_susie_sensitivity_signal_pairs.csv"))
  pair_n <- min(nrow(ptab), nrow(stab))
  robust_n <- 0L
  if (pair_n > 0L && "PP.H4.abf" %in% names(ptab) && "PP.H4.abf" %in% names(stab)) {
    robust_n <- sum(
      ptab$PP.H4.abf[seq_len(pair_n)] >= 0.80 &
        stab$PP.H4.abf[seq_len(pair_n)] >= 0.50,
      na.rm = TRUE
    )
  }
  fabf <- as.list(filtered_abf$summary)
  decision <- if (robust_n > 0L) {
    "robust_susie_shared_signal_requires_final_review"
  } else {
    "no_robust_susie_shared_signal"
  }
  ans <- data.table(
    run_id = meta$run_id, combo_id = meta$combo_id,
    locus_label = meta$locus_label, factor_gwas = meta$factor_gwas,
    gene_symbol = meta$gene_symbol, ensembl_id = meta$ensembl_id,
    resource = meta$resource, tissue = meta$tissue,
    original_abf_h3 = meta$pp_h3, original_abf_h4 = meta$pp_h4,
    harmonized_snp_n = nrow(h), ld_reference_snp_n = nrow(h_ld),
    retained_snp_n = sum(keep),
    filtered_abf_h3 = as.numeric(fabf[["PP.H3.abf"]]),
    filtered_abf_h4 = as.numeric(fabf[["PP.H4.abf"]]),
    susie_gwas_signal_n = length(susie_gwas$sets$cs),
    susie_qtl_signal_n = length(susie_qtl$sets$cs),
    susie_signal_pair_n = nrow(ptab),
    susie_max_h3 = if ("PP.H3.abf" %in% names(ptab)) safe_max(ptab$PP.H3.abf) else NA_real_,
    susie_max_h4 = if ("PP.H4.abf" %in% names(ptab)) safe_max(ptab$PP.H4.abf) else NA_real_,
    susie_max_h4_p12_1e6 = if ("PP.H4.abf" %in% names(stab)) safe_max(stab$PP.H4.abf) else NA_real_,
    robust_shared_signal_pair_n = robust_n,
    susie_status = "completed", decision = decision
  )
  fwrite(ans, summary_file)
  writeLines(
    c(
      "stage=frozen_susie",
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      paste0("decision=", decision)
    ),
    complete_file
  )
  ans
}

results <- list()
failures <- list()
for (i in seq_len(nrow(trigger))) {
  meta <- trigger[i]
  message(
    "SuSiE ", i, "/", nrow(trigger), ": ", meta$locus_label,
    " / ", meta$resource, " / ", meta$tissue, " / ", meta$gene_symbol
  )
  known <- known_ld_mismatch[combo_id == meta$combo_id]
  if (nrow(known)) {
    ans <- write_ld_mismatch_result(meta, known$error_message[1L])
  } else {
    ans <- tryCatch(run_one(meta), error = function(e) e)
  }
  if (inherits(ans, "error")) {
    error_text <- conditionMessage(ans)
    if (grepl("estimated prior variance is unreasonably large", error_text, fixed = TRUE)) {
      results[[length(results) + 1L]] <- write_ld_mismatch_result(meta, error_text)
    } else {
      failures[[length(failures) + 1L]] <- data.table(
        run_id = meta$run_id, combo_id = meta$combo_id,
        error_message = error_text
      )
    }
  } else {
    results[[length(results) + 1L]] <- ans
  }
}
result_dt <- if (length(results)) rbindlist(results, fill = TRUE) else data.table()
failure_dt <- if (length(failures)) rbindlist(failures, fill = TRUE) else empty_failure()
fwrite(result_dt, file.path(out_dir, "all_frozen_susie_results.csv"))
fwrite(failure_dt, file.path(out_dir, "frozen_susie_failures.csv"))

if (nrow(failure_dt)) {
  jsonlite::write_json(
    list(
      stage = "frozen_three_locus_susie", complete = FALSE,
      trigger_n = nrow(trigger), completed_n = nrow(result_dt),
      failure_n = nrow(failure_dt)
    ),
    file.path(out_dir, "frozen_susie_audit.json"),
    pretty = TRUE, auto_unbox = TRUE
  )
  stop("Frozen SuSiE stage has failures; inspect frozen_susie_failures.csv")
}

if (nrow(result_dt) != nrow(trigger)) stop("Not all frozen triggers produced a result")
robust <- result_dt[robust_shared_signal_pair_n > 0L]
not_estimable <- result_dt[susie_status != "completed"]
fwrite(robust, file.path(out_dir, "robust_susie_shared_signal_manifest.csv"))
summary_by_locus <- result_dt[, .(
  trigger_n = .N,
  completed_n = sum(susie_status == "completed"),
  not_estimable_n = sum(susie_status != "completed"),
  robust_shared_combination_n = sum(robust_shared_signal_pair_n > 0L),
  robust_shared_gene_n = uniqueN(gene_symbol[robust_shared_signal_pair_n > 0L])
), by = .(locus_label, factor_gwas)]
fwrite(summary_by_locus, file.path(out_dir, "frozen_susie_summary_by_locus.csv"))

stage_closed <- nrow(robust) == 0L
audit <- list(
  stage = "frozen_three_locus_susie",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  complete = TRUE,
  frozen_trigger_n = nrow(trigger),
  completed_n = sum(result_dt$susie_status == "completed"),
  not_estimable_n = nrow(not_estimable),
  robust_shared_combination_n = nrow(robust),
  qtl_stage_closed = stage_closed,
  next_step = if (stage_closed) {
    "Close bounded QTL stage and construct the final evidence matrix."
  } else {
    "Review only the robust shared-signal manifest; do not add genes, tissues, loci or resources."
  },
  interpretation_boundary = paste(
    "Post-LOLO targeted molecular resolution only;",
    "not independent factor-GWAS replication and not proof of mediation."
  ),
  prohibited_expansion = c(
    "new genes", "new tissues", "new loci", "new QTL resources",
    "threshold relaxation", "new plasma proteomics"
  )
)
jsonlite::write_json(
  audit, file.path(out_dir, "frozen_susie_audit.json"),
  pretty = TRUE, auto_unbox = TRUE
)
writeLines(
  c(
    "stage=frozen_three_locus_susie",
    paste0("completed_at=", audit$completed_at),
    paste0("frozen_trigger_n=", nrow(trigger)),
    paste0("robust_shared_combination_n=", nrow(robust)),
    paste0("qtl_stage_closed=", tolower(as.character(stage_closed)))
  ),
  file.path(out_dir, "COMPLETE.txt")
)
message("Completed all frozen SuSiE diagnostics: ", out_dir)
message("Robust shared-signal combinations: ", nrow(robust))
message("QTL stage closed: ", stage_closed)


