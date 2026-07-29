options(stringsAsFactors = FALSE)

message("Track C candidate QTL primary molecular analysis")

required_packages <- c("data.table", "coloc")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}
library(data.table)

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
trackc_root <- Sys.getenv(
  "TRACKC_ROOT",
  unset = file.path(base_dir, "analysis_ready_core", "trackC")
)
qtl_root <- Sys.getenv(
  "TRACKC_QTL_ROOT",
  unset = file.path(trackc_root, "09_QTL_molecular_followup")
)
out_dir <- Sys.getenv(
  "TRACKC_QTL_PRIMARY_OUT_DIR",
  unset = file.path(qtl_root, "03_candidate_primary_molecular_analysis")
)
stages <- unique(trimws(strsplit(
  Sys.getenv("TRACKC_QTL_PRIMARY_STAGES", unset = "prepare"),
  ",",
  fixed = TRUE
)[[1]]))
allowed_stages <- c("prepare", "abf", "smr", "integrate")
if (any(!stages %in% allowed_stages)) {
  stop("Unknown stage(s): ", paste(setdiff(stages, allowed_stages), collapse = ", "))
}

gtex_dir <- file.path(qtl_root, "01_GTEx_V8_remote_tabix_regions")
pqtl_dir <- file.path(qtl_root, "02_pQTL_raw")
spec_dir <- file.path(qtl_root, "00_frozen_region_spec")
prepare_dir <- file.path(out_dir, "00_harmonized_inputs")
abf_dir <- file.path(out_dir, "01_coloc_abf")
smr_dir <- file.path(out_dir, "02_smr_heidi")
integrated_dir <- file.path(out_dir, "03_integrated_evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

factor_files <- c(
  affective = file.path(
    base_dir, "analysis_ready_core", "trackB", "three_projection_locus_magma",
    "affective", "trackB_affective_fuma_input.tsv.gz"
  ),
  functional = file.path(
    base_dir, "analysis_ready_core", "trackB", "three_projection_locus_magma",
    "functional", "trackB_functional_fuma_input.tsv.gz"
  )
)

bfile_prefix <- Sys.getenv(
  "TRACKC_LD_BFILE",
  unset = file.path(base_dir, "reference", "plink", "EUR")
)
smr_exe <- Sys.getenv(
  "TRACKC_SMR_EXE",
  unset = Sys.getenv("SMR_EXECUTABLE", unset = "smr")
)
if (!file.exists(smr_exe)) {
  located_smr <- Sys.which(smr_exe)
  if (nzchar(located_smr)) smr_exe <- unname(located_smr)
}

require_file <- function(path, label = "file") {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_complete <- function(directory, stage) {
  writeLines(
    c(
      paste0("stage=", stage),
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    ),
    file.path(directory, "COMPLETE.txt"),
    useBytes = TRUE
  )
}

hash_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(path, algo = "sha256", file = TRUE)
}

pick_column <- function(x, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0L) {
    lower_names <- tolower(names(x))
    idx <- match(tolower(candidates), lower_names, nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) > 0L) hit <- names(x)[idx[1L]]
  }
  if (length(hit) == 0L) {
    if (required) stop("Missing column; expected one of: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  hit[1L]
}

is_palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")

genes <- fread(require_file(file.path(spec_dir, "candidate_genes_grch38.csv")))
factor_map <- fread(require_file(file.path(spec_dir, "candidate_factor_mapping.csv")))
tissue_spec <- fread(require_file(file.path(spec_dir, "frozen_GTEx_remote_tabix_datasets.csv")))

find_pqtl_text <- function(resource, gene) {
  directory <- file.path(pqtl_dir, resource, gene)
  if (!dir.exists(directory)) return(NA_character_)
  candidates <- list.files(directory, recursive = TRUE, full.names = TRUE)
  candidates <- candidates[
    grepl("\\.(txt|tsv|csv)(\\.gz)?$", candidates, ignore.case = TRUE) &
      !grepl("(manifest|audit|failure|log|matrix|flist)", basename(candidates), ignore.case = TRUE)
  ]
  if (length(candidates) == 0L) return(NA_character_)
  sizes <- file.info(candidates)$size
  candidates[order(sizes, decreasing = TRUE)][1L]
}

standardize_pqtl <- function(path, gene, resource, build) {
  x <- fread(path)
  snp <- pick_column(x, c("SNP", "snp", "rsid", "rsID"))
  chr <- pick_column(x, c("CHR", "chr", "chromosome"))
  bp <- pick_column(x, c("BP", "bp", "position", "pos"))
  ea <- pick_column(x, c("effect_allele", "A1", "a1", "ALT", "alt"))
  oa <- pick_column(x, c("other_allele", "A2", "a2", "REF", "ref"))
  beta <- pick_column(x, c("beta", "BETA", "b"))
  se <- pick_column(x, c("se", "SE", "stderr"))
  p <- pick_column(x, c("pval", "p", "P", "pvalue"))
  eaf <- pick_column(x, c("eaf", "EAF", "effect_allele_freq", "freq", "Freq"), FALSE)
  maf <- pick_column(x, c("MAF", "maf"), FALSE)
  n <- pick_column(x, c("N", "n", "samplesize", "sample_size"), FALSE)

  out <- data.table(
    gene_symbol = gene,
    source = resource,
    build = as.integer(build),
    SNP = trimws(as.character(x[[snp]])),
    CHR = as.integer(x[[chr]]),
    BP = as.integer(x[[bp]]),
    effect_allele = toupper(trimws(as.character(x[[ea]]))),
    other_allele = toupper(trimws(as.character(x[[oa]]))),
    beta = as.numeric(x[[beta]]),
    se = as.numeric(x[[se]]),
    pvalue = as.numeric(x[[p]]),
    eaf = if (!is.na(eaf)) as.numeric(x[[eaf]]) else NA_real_,
    MAF = if (!is.na(maf)) as.numeric(x[[maf]]) else NA_real_,
    N = if (!is.na(n)) as.numeric(x[[n]]) else NA_real_
  )
  out[!is.finite(MAF) & is.finite(eaf), MAF := pmin(eaf, 1 - eaf)]
  out <- out[
    grepl("^rs[0-9]+$", SNP) &
      nchar(effect_allele) == 1L & nchar(other_allele) == 1L &
      effect_allele %in% c("A", "C", "G", "T") &
      other_allele %in% c("A", "C", "G", "T") &
      !is_palindromic(effect_allele, other_allele) &
      is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(pvalue) & pvalue > 0 & pvalue <= 1 &
      is.finite(eaf) & eaf > 0 & eaf < 1
  ]
  setorder(out, SNP, pvalue)
  out[!duplicated(SNP)]
}

read_factor <- function(factor_name, wanted_snps) {
  path <- require_file(factor_files[[factor_name]], paste(factor_name, "factor GWAS"))
  x <- fread(path, select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"))
  x <- x[SNP %in% wanted_snps]
  x[, `:=`(
    SNP = trimws(as.character(SNP)),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    N = as.numeric(N)
  )]
  x[
    grepl("^rs[0-9]+$", SNP) &
      !is_palindromic(A1, A2) &
      is.finite(BETA) & is.finite(SE) & SE > 0 &
      is.finite(P) & P > 0 & P <= 1 & is.finite(N) & N > 0
  ][!duplicated(SNP)]
}

harmonize_factor_qtl <- function(gwas, qtl) {
  gwas <- copy(gwas)
  qtl <- copy(qtl)
  setnames(
    gwas,
    c("BETA", "SE", "P", "N"),
    c("gwas_beta_input", "gwas_se_input", "gwas_p_input", "gwas_n_input")
  )
  setnames(
    qtl,
    c("beta", "se", "pvalue", "N"),
    c("qtl_beta_input", "qtl_se_input", "qtl_p_input", "qtl_n_input")
  )
  merged <- merge(gwas, qtl, by = "SNP", all = FALSE, sort = FALSE)
  same <- merged$A1 == merged$effect_allele & merged$A2 == merged$other_allele
  reverse <- merged$A1 == merged$other_allele & merged$A2 == merged$effect_allele
  merged <- merged[same | reverse]
  reverse <- reverse[same | reverse]
  merged[, `:=`(
    beta_gwas_aligned = ifelse(reverse, -gwas_beta_input, gwas_beta_input),
    se_gwas = gwas_se_input,
    p_gwas = gwas_p_input,
    n_gwas = gwas_n_input,
    beta_qtl = qtl_beta_input,
    se_qtl = qtl_se_input,
    p_qtl = qtl_p_input,
    n_qtl = qtl_n_input,
    maf_qtl = MAF,
    eaf_qtl = eaf,
    aligned_effect_allele = effect_allele,
    aligned_other_allele = other_allele
  )]
  merged <- merged[
    is.finite(beta_gwas_aligned) & is.finite(se_gwas) & se_gwas > 0 &
      is.finite(beta_qtl) & is.finite(se_qtl) & se_qtl > 0 &
      is.finite(maf_qtl) & maf_qtl > 0 & maf_qtl < 0.5 &
      is.finite(eaf_qtl) & eaf_qtl > 0 & eaf_qtl < 1
  ]
  merged[!duplicated(SNP)]
}

build_combo_manifest <- function() {
  gtex_coverage <- fread(require_file(file.path(gtex_dir, "GTEx_remote_tabix_candidate_coverage.csv")))
  if (nrow(gtex_coverage) != 9L) stop("Expected 9 GTEx gene-tissue coverage rows")
  gtex_combos <- merge(
    factor_map,
    gtex_coverage[, .(gene_symbol, ensembl_id, tissue, dataset_id, sample_size, output_file)],
    by = c("gene_symbol", "ensembl_id"),
    allow.cartesian = TRUE
  )
  gtex_combos[, `:=`(
    qtl_class = "eQTL",
    resource = "GTEx_remote_tabix",
    source_build = 38L,
    qtl_file = output_file
  )]

  pqtl_base <- data.table(
    gene_symbol = c("GPX1", "TCF4", "TCF4", "GPX1", "TCF4", "TCF4"),
    ensembl_id = c(
      "ENSG00000233276", "ENSG00000196628", "ENSG00000196628",
      "ENSG00000233276", "ENSG00000196628", "ENSG00000196628"
    ),
    factor_gwas = c("functional", "affective", "functional", "functional", "affective", "functional"),
    resource = c("decode", "decode", "decode", "fenland", "fenland", "fenland"),
    source_build = c(38L, 38L, 38L, 37L, 37L, 37L)
  )
  pqtl_base[, qtl_file := vapply(
    seq_len(.N),
    function(i) find_pqtl_text(resource[i], gene_symbol[i]),
    character(1)
  )]
  pqtl_base[, `:=`(qtl_class = "pQTL", tissue = "circulating_protein", dataset_id = resource, sample_size = NA_real_)]
  if (any(is.na(pqtl_base$qtl_file))) {
    stop("Missing targeted pQTL text file(s); complete acquisition stage first")
  }

  combos <- rbindlist(list(
    gtex_combos[, .(gene_symbol, ensembl_id, factor_gwas, qtl_class, resource, tissue, dataset_id, source_build, sample_size, qtl_file)],
    pqtl_base[, .(gene_symbol, ensembl_id, factor_gwas, qtl_class, resource, tissue, dataset_id, source_build, sample_size, qtl_file)]
  ), fill = TRUE)
  combos[, combo_id := sprintf(
    "Q%02d_%s_%s_%s_%s",
    seq_len(.N), gene_symbol, factor_gwas, resource, tissue
  )]
  setcolorder(combos, c("combo_id", setdiff(names(combos), "combo_id")))
  combos
}

if ("prepare" %in% stages) {
  require_file(file.path(gtex_dir, "COMPLETE.txt"), "GTEx API completion marker")
  require_file(file.path(pqtl_dir, "COMPLETE.txt"), "pQTL completion marker")
  dir.create(prepare_dir, recursive = TRUE, showWarnings = FALSE)
  combos <- build_combo_manifest()

  qtl_cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(combos))) {
    row <- combos[i]
    key <- paste(row$resource, row$tissue, row$gene_symbol, sep = "|")
    if (!exists(key, qtl_cache, inherits = FALSE)) {
      qtl <- if (row$qtl_class == "eQTL") {
        fread(require_file(row$qtl_file))
      } else {
        standardize_pqtl(row$qtl_file, row$gene_symbol, row$resource, row$source_build)
      }
      assign(key, qtl, qtl_cache)
    }
  }
  wanted_snps <- unique(unlist(lapply(ls(qtl_cache), function(k) get(k, qtl_cache)$SNP)))
  factor_cache <- list(
    affective = read_factor("affective", wanted_snps),
    functional = read_factor("functional", wanted_snps)
  )

  qa <- list()
  failures <- list()
  for (i in seq_len(nrow(combos))) {
    row <- combos[i]
    key <- paste(row$resource, row$tissue, row$gene_symbol, sep = "|")
    qtl <- get(key, qtl_cache)
    gwas <- factor_cache[[row$factor_gwas]]
    harmonized <- harmonize_factor_qtl(gwas, qtl)
    combo_dir <- file.path(prepare_dir, row$combo_id)
    dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
    hfile <- file.path(combo_dir, "harmonized_factor_qtl.tsv.gz")
    qfile <- file.path(combo_dir, "standardized_qtl.tsv.gz")
    fwrite(harmonized, hfile, sep = "\t", compress = "gzip")
    fwrite(qtl, qfile, sep = "\t", compress = "gzip")
    eligible_abf <- nrow(harmonized) >= 50L
    eligible_smr <- eligible_abf && any(harmonized$p_qtl < 5e-8, na.rm = TRUE)
    eligible_coloc_interpretation <- eligible_smr &&
      any(harmonized$p_gwas < 1e-5, na.rm = TRUE)
    qa[[i]] <- cbind(
      row,
      data.table(
        qtl_rows = nrow(qtl),
        harmonized_rows = nrow(harmonized),
        qtl_min_p_raw = if (nrow(qtl) > 0L) min(qtl$pvalue) else NA_real_,
        qtl_min_p = if (nrow(harmonized) > 0L) min(harmonized$p_qtl) else NA_real_,
        qtl_p5e8_harmonized_n = if (nrow(harmonized) > 0L) sum(harmonized$p_qtl < 5e-8) else 0L,
        gwas_min_p = if (nrow(harmonized) > 0L) min(harmonized$p_gwas) else NA_real_,
        eligible_abf = eligible_abf,
        eligible_smr = eligible_smr,
        eligible_coloc_interpretation = eligible_coloc_interpretation,
        eligible_multisignal_followup = eligible_coloc_interpretation &&
          nrow(harmonized) >= 100L &&
          any(harmonized$p_qtl < 5e-8, na.rm = TRUE),
        harmonized_file = hfile,
        standardized_qtl_file = qfile,
        harmonized_sha256 = hash_file(hfile)
      )
    )
    if (!eligible_abf) {
      failures[[length(failures) + 1L]] <- data.table(
        combo_id = row$combo_id,
        stage = "prepare",
        reason = "fewer_than_50_harmonized_nonpalindromic_unique_rsid_snps"
      )
    }
  }
  qa_dt <- rbindlist(qa, fill = TRUE)
  failure_dt <- rbindlist(failures, fill = TRUE)
  if (ncol(failure_dt) == 0L) {
    failure_dt <- data.table(combo_id = character(), stage = character(), reason = character())
  }
  fwrite(combos, file.path(prepare_dir, "frozen_18_combination_manifest.csv"))
  fwrite(qa_dt, file.path(prepare_dir, "harmonization_QA.csv"))
  fwrite(failure_dt, file.path(prepare_dir, "preparation_noneligible_combinations.csv"))
  writeLines(
    c(
      "Candidate QTL harmonisation specification",
      "Factor GWAS coordinates are GRCh37; GTEx/deCODE QTL coordinates are GRCh38; Fenland is GRCh37.",
      "Cross-build matching uses stable unique rsIDs plus exact unordered allele pairs.",
      "Primary analysis is restricted to unique biallelic non-palindromic SNPs.",
      "GWAS effects are reoriented to the QTL effect allele before analysis.",
      "The GRCh37 1000G EUR LD panel is used by rsID; coordinates are not used to match cross-build variants.",
      "At least 50 harmonised SNPs are required for ABF colocalisation.",
      "A genome-wide significant cis-QTL is required for SMR/HEIDI eligibility.",
      "H3/H4 is inferentially interpreted only when cis-QTL P < 5e-8 and regional factor-GWAS P < 1e-5.",
      "ABF results outside that gate are retained as descriptive diagnostics only.",
      "No result may be interpreted as causal or as independent factor-GWAS replication."
    ),
    file.path(prepare_dir, "README_harmonization.txt"),
    useBytes = TRUE
  )
  if (nrow(qa_dt) != 18L) stop("Preparation did not produce all 18 prespecified combinations")
  write_complete(prepare_dir, "prepare")
  message("Completed preparation stage: ", prepare_dir)
}

run_abf <- function(row) {
  h <- fread(row$harmonized_file)
  if (!isTRUE(row$eligible_coloc_interpretation)) {
    return(data.table(
      combo_id = row$combo_id,
      status = "not_eligible_without_strong_cis_qtl_and_regional_gwas_signal"
    ))
  }
  if (nrow(h) < 50L) {
    return(data.table(combo_id = row$combo_id, status = "not_estimable_sparse_overlap"))
  }
  d_gwas <- list(
    beta = h$beta_gwas_aligned,
    varbeta = h$se_gwas^2,
    snp = h$SNP,
    type = "quant",
    sdY = 1
  )
  qtl_sdY <- sqrt(stats::median(
    2 * h$n_qtl * h$maf_qtl * (1 - h$maf_qtl) * h$se_qtl^2,
    na.rm = TRUE
  ))
  if (!is.finite(qtl_sdY) || qtl_sdY <= 0) qtl_sdY <- 1
  d_qtl <- list(
    beta = h$beta_qtl,
    varbeta = h$se_qtl^2,
    snp = h$SNP,
    type = "quant",
    sdY = qtl_sdY
  )
  primary <- coloc::coloc.abf(d_gwas, d_qtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  sensitivity <- coloc::coloc.abf(d_gwas, d_qtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-6)
  primary_summary <- as.list(primary$summary)
  sensitivity_summary <- as.list(sensitivity$summary)
  combo_dir <- file.path(abf_dir, row$combo_id)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(primary, file.path(combo_dir, "coloc_abf_primary.rds"))
  saveRDS(sensitivity, file.path(combo_dir, "coloc_abf_p12_1e6_sensitivity.rds"))
  fwrite(as.data.table(primary$results), file.path(combo_dir, "coloc_abf_snp_results.tsv.gz"), sep = "\t", compress = "gzip")
  data.table(
    combo_id = row$combo_id,
    status = "completed",
    nsnps = as.integer(primary_summary[["nsnps"]]),
    pp_h0 = as.numeric(primary_summary[["PP.H0.abf"]]),
    pp_h1 = as.numeric(primary_summary[["PP.H1.abf"]]),
    pp_h2 = as.numeric(primary_summary[["PP.H2.abf"]]),
    pp_h3 = as.numeric(primary_summary[["PP.H3.abf"]]),
    pp_h4 = as.numeric(primary_summary[["PP.H4.abf"]]),
    pp_h4_p12_1e6 = as.numeric(sensitivity_summary[["PP.H4.abf"]]),
    qtl_sdY_estimate = qtl_sdY
  )
}

if ("abf" %in% stages) {
  require_file(file.path(prepare_dir, "COMPLETE.txt"), "prepare completion marker")
  dir.create(abf_dir, recursive = TRUE, showWarnings = FALSE)
  qa <- fread(file.path(prepare_dir, "harmonization_QA.csv"))
  results <- list()
  failures <- list()
  for (i in seq_len(nrow(qa))) {
    row <- qa[i]
    result <- tryCatch(run_abf(row), error = function(e) e)
    if (inherits(result, "error")) {
      failures[[length(failures) + 1L]] <- data.table(
        combo_id = row$combo_id,
        error_message = conditionMessage(result)
      )
    } else {
      results[[length(results) + 1L]] <- result
    }
  }
  result_dt <- merge(qa, rbindlist(results, fill = TRUE), by = "combo_id", all.x = TRUE)
  failure_dt <- rbindlist(failures, fill = TRUE)
  if (ncol(failure_dt) == 0L) {
    failure_dt <- data.table(combo_id = character(), error_message = character())
  }
  fwrite(result_dt, file.path(abf_dir, "candidate_coloc_abf_summary.csv"))
  fwrite(failure_dt, file.path(abf_dir, "candidate_coloc_abf_failures.csv"))
  if (nrow(failure_dt) > 0L) stop("ABF stage has failures; inspect candidate_coloc_abf_failures.csv")
  write_complete(abf_dir, "abf")
  message("Completed ABF stage: ", abf_dir)
}

run_command <- function(executable, args, log_file) {
  output <- system2(executable, args = args, stdout = TRUE, stderr = TRUE)
  writeLines(output, log_file, useBytes = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("External command failed with status ", status, "; inspect ", log_file)
  }
  invisible(output)
}

make_besd <- function(qtl, row, combo_dir) {
  qtl <- qtl[
    grepl("^rs[0-9]+$", SNP) &
      is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(pvalue) & pvalue > 0 & pvalue <= 1 &
      is.finite(eaf) & eaf > 0 & eaf < 1
  ][!duplicated(SNP)]
  if (nrow(qtl) == 0L) stop("No valid QTL rows for BESD")
  matrix_file <- file.path(combo_dir, "qtl_matrix_eqtl_format.txt")
  prefix <- file.path(combo_dir, "qtl_bundle")
  matrix_dt <- qtl[, .(
    SNP,
    gene = row$gene_symbol,
    beta,
    `t-stat` = beta / se,
    `p-value` = pvalue,
    FDR = p.adjust(pvalue, method = "BH")
  )]
  fwrite(matrix_dt, matrix_file, sep = " ", quote = FALSE, na = "NA")
  run_command(
    smr_exe,
    c(
      "--eqtl-summary", normalizePath(matrix_file, winslash = "/", mustWork = TRUE),
      "--matrix-eqtl-format", "--make-besd",
      "--out", normalizePath(prefix, winslash = "/", mustWork = FALSE)
    ),
    file.path(combo_dir, "make_besd.log.txt")
  )
  epi_file <- paste0(prefix, ".epi")
  esi_file <- paste0(prefix, ".esi")
  besd_file <- paste0(prefix, ".besd")
  require_file(epi_file); require_file(esi_file); require_file(besd_file)

  gene_row <- genes[gene_symbol == row$gene_symbol]
  if (row$source_build == 38L) {
    gene_bp <- round((gene_row$start_grch38 + gene_row$end_grch38) / 2)
  } else {
    grch37 <- list(GPX1 = c(49394609L, 49396033L), TCF4 = c(52889562L, 53332018L))
    gene_bp <- round(mean(grch37[[row$gene_symbol]]))
  }
  epi <- fread(epi_file, header = FALSE)
  if (nrow(epi) != 1L || ncol(epi) < 6L) stop("Unexpected EPI structure")
  epi <- data.table(
    chr = as.integer(qtl$CHR[1]), probe = row$gene_symbol, cm = 0,
    bp = as.integer(gene_bp), gene = row$gene_symbol, strand = "NA"
  )
  fwrite(epi, epi_file, sep = "\t", col.names = FALSE, quote = FALSE)

  esi <- fread(esi_file, header = FALSE)
  if (ncol(esi) < 7L || nrow(esi) != nrow(qtl)) stop("Unexpected ESI row count/structure")
  setnames(esi, 1:7, c("chr", "SNP", "cm", "bp", "A1", "A2", "freq"))
  info <- qtl[, .(
    SNP,
    chr_new = as.integer(CHR), bp_new = as.integer(BP),
    A1_new = effect_allele, A2_new = other_allele, freq_new = eaf
  )]
  idx <- match(as.character(esi$SNP), info$SNP)
  if (anyNA(idx)) stop("ESI SNPs do not map one-to-one to standardized QTL")
  esi[, `:=`(
    chr = info$chr_new[idx], cm = 0, bp = info$bp_new[idx],
    A1 = info$A1_new[idx], A2 = info$A2_new[idx], freq = info$freq_new[idx]
  )]
  fwrite(esi[, .(chr, SNP, cm, bp, A1, A2, freq)], esi_file, sep = "\t", col.names = FALSE, quote = FALSE)
  prefix
}

parse_smr <- function(path, combo_id) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    return(data.table(combo_id = combo_id, smr_status = "no_result"))
  }
  x <- fread(path)
  if (nrow(x) == 0L) return(data.table(combo_id = combo_id, smr_status = "no_result"))
  x <- x[1L]
  p_smr <- pick_column(x, c("p_SMR", "p.SMR", "pSMR"), FALSE)
  p_heidi <- pick_column(x, c("p_HEIDI", "p.HEIDI", "pHEIDI"), FALSE)
  n_heidi <- pick_column(x, c("nsnp_HEIDI", "nsnp.HEIDI", "nsnpHEIDI"), FALSE)
  b_smr <- pick_column(x, c("b_SMR", "b.SMR", "bSMR"), FALSE)
  se_smr <- pick_column(x, c("se_SMR", "se.SMR", "seSMR"), FALSE)
  data.table(
    combo_id = combo_id,
    smr_status = "completed",
    smr_beta = if (!is.na(b_smr)) as.numeric(x[[b_smr]]) else NA_real_,
    smr_se = if (!is.na(se_smr)) as.numeric(x[[se_smr]]) else NA_real_,
    smr_p = if (!is.na(p_smr)) as.numeric(x[[p_smr]]) else NA_real_,
    heidi_p = if (!is.na(p_heidi)) as.numeric(x[[p_heidi]]) else NA_real_,
    heidi_n_snp = if (!is.na(n_heidi)) as.integer(x[[n_heidi]]) else NA_integer_
  )
}

run_smr_combo <- function(row) {
  if (!isTRUE(row$eligible_smr)) {
    return(data.table(combo_id = row$combo_id, smr_status = "not_eligible_no_strong_cis_qtl_or_sparse"))
  }
  combo_dir <- file.path(smr_dir, row$combo_id)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
  qtl <- fread(row$standardized_qtl_file)
  h <- fread(row$harmonized_file)
  prefix <- make_besd(qtl, row, combo_dir)
  gwas_ma <- h[, .(
    SNP,
    A1 = aligned_effect_allele,
    A2 = aligned_other_allele,
    freq = eaf_qtl,
    b = beta_gwas_aligned,
    se = se_gwas,
    p = p_gwas,
    N = as.integer(round(n_gwas))
  )]
  gwas_ma_file <- file.path(combo_dir, "factor_gwas_region.ma")
  fwrite(gwas_ma, gwas_ma_file, sep = "\t")
  out_prefix <- file.path(combo_dir, "smr_heidi")
  args <- c(
    "--bfile", normalizePath(bfile_prefix, winslash = "/", mustWork = FALSE),
    "--gwas-summary", normalizePath(gwas_ma_file, winslash = "/", mustWork = TRUE),
    "--beqtl-summary", normalizePath(prefix, winslash = "/", mustWork = FALSE),
    "--gene", row$gene_symbol,
    "--smr", "--heidi-mtd", "1",
    "--peqtl-smr", "5e-8", "--peqtl-heidi", "1.57e-3",
    "--ld-upper-limit", "0.9", "--ld-lower-limit", "0.05",
    "--heidi-min-m", "10", "--heidi-max-m", "20",
    "--diff-freq", "0.2", "--diff-freq-prop", "0.05",
    "--thread-num", "4",
    "--out", normalizePath(out_prefix, winslash = "/", mustWork = FALSE)
  )
  smr_log <- paste0(out_prefix, ".log.txt")
  command_result <- tryCatch(
    {
      run_command(smr_exe, args, smr_log)
      NULL
    },
    error = function(e) e
  )
  if (inherits(command_result, "error")) {
    log_lines <- if (file.exists(smr_log)) readLines(smr_log, warn = FALSE) else character()
    if (any(grepl("more than .* SNPs were removed by the allele frequency difference check", log_lines))) {
      return(data.table(
        combo_id = row$combo_id,
        smr_status = "not_estimable_prespecified_frequency_qc_failure",
        smr_qc_note = "More than 5% of SNPs differed in allele frequency by >0.20; thresholds were not relaxed."
      ))
    }
    stop(command_result)
  }
  parse_smr(paste0(out_prefix, ".smr"), row$combo_id)
}

if ("smr" %in% stages) {
  require_file(file.path(prepare_dir, "COMPLETE.txt"), "prepare completion marker")
  require_file(paste0(bfile_prefix, ".bed"), "LD BED")
  require_file(paste0(bfile_prefix, ".bim"), "LD BIM")
  require_file(paste0(bfile_prefix, ".fam"), "LD FAM")
  require_file(smr_exe, "SMR executable")
  dir.create(smr_dir, recursive = TRUE, showWarnings = FALSE)
  qa <- fread(file.path(prepare_dir, "harmonization_QA.csv"))
  results <- list()
  failures <- list()
  for (i in seq_len(nrow(qa))) {
    row <- qa[i]
    result <- tryCatch(run_smr_combo(row), error = function(e) e)
    if (inherits(result, "error")) {
      failures[[length(failures) + 1L]] <- data.table(
        combo_id = row$combo_id,
        error_message = conditionMessage(result)
      )
    } else {
      results[[length(results) + 1L]] <- result
    }
  }
  result_dt <- merge(qa, rbindlist(results, fill = TRUE), by = "combo_id", all.x = TRUE)
  for (column in c("smr_beta", "smr_se", "smr_p", "heidi_p")) {
    if (!column %in% names(result_dt)) result_dt[, (column) := NA_real_]
  }
  if (!"heidi_n_snp" %in% names(result_dt)) result_dt[, heidi_n_snp := NA_integer_]
  if (!"smr_qc_note" %in% names(result_dt)) result_dt[, smr_qc_note := NA_character_]
  failure_dt <- rbindlist(failures, fill = TRUE)
  if (ncol(failure_dt) == 0L) {
    failure_dt <- data.table(combo_id = character(), error_message = character())
  }
  fwrite(result_dt, file.path(smr_dir, "candidate_smr_heidi_summary.csv"))
  fwrite(failure_dt, file.path(smr_dir, "candidate_smr_heidi_failures.csv"))
  if (nrow(failure_dt) > 0L) stop("SMR stage has failures; inspect candidate_smr_heidi_failures.csv")
  write_complete(smr_dir, "smr")
  message("Completed SMR/HEIDI stage: ", smr_dir)
}

if ("integrate" %in% stages) {
  require_file(file.path(abf_dir, "COMPLETE.txt"), "ABF completion marker")
  require_file(file.path(smr_dir, "COMPLETE.txt"), "SMR completion marker")
  dir.create(integrated_dir, recursive = TRUE, showWarnings = FALSE)
  abf <- fread(file.path(abf_dir, "candidate_coloc_abf_summary.csv"))
  smr <- fread(file.path(smr_dir, "candidate_smr_heidi_summary.csv"))
  smr_keep <- smr[, .(combo_id, smr_status, smr_beta, smr_se, smr_p, heidi_p, heidi_n_snp)]
  out <- merge(abf, smr_keep, by = "combo_id", all.x = TRUE)
  out[, smr_fdr_within_qtl_class := p.adjust(smr_p, method = "BH"), by = qtl_class]
  out[, smr_fdr_global := p.adjust(smr_p, method = "BH")]
  out[, `:=`(
    heidi_pass_strict = is.finite(heidi_p) & heidi_p > 0.01 & heidi_n_snp >= 10,
    coloc_abf_strong = eligible_coloc_interpretation &
      is.finite(pp_h4) & pp_h4 >= 0.80 & pp_h4 > pp_h3,
    coloc_distinct_signal = eligible_coloc_interpretation &
      is.finite(pp_h3) & pp_h3 >= 0.80 & pp_h3 > pp_h4,
    smr_global_significant = is.finite(smr_fdr_global) & smr_fdr_global < 0.05
  )]
  out[, evidence_class_primary := fifelse(
    smr_global_significant & heidi_pass_strict & coloc_abf_strong,
    "provisional_strict_anchor_pending_multisignal_diagnostic",
    fifelse(
      coloc_distinct_signal,
      "distinct_signals",
      fifelse(
        (smr_global_significant & heidi_pass_strict) | coloc_abf_strong,
        "candidate_driven_supportive_evidence",
        fifelse(
          smr_global_significant & !heidi_pass_strict,
          "locus_shared_but_molecularly_unresolved",
          fifelse(
            !eligible_smr,
            "insufficient_qtl_coverage_or_no_strong_cis_qtl",
            "no_primary_molecular_support"
          )
        )
      )
    )
  )]
  out[, multisignal_followup_trigger :=
    eligible_multisignal_followup &
      (
        evidence_class_primary %in% c(
          "provisional_strict_anchor_pending_multisignal_diagnostic",
          "candidate_driven_supportive_evidence",
          "distinct_signals",
          "locus_shared_but_molecularly_unresolved"
        ) |
          (is.finite(pp_h3) & pp_h3 >= 0.50) |
          (is.finite(pp_h4) & pp_h4 >= 0.50)
      )
  ]
  out[, permitted_interpretation := fifelse(
    evidence_class_primary == "provisional_strict_anchor_pending_multisignal_diagnostic",
    "strong convergent molecular evidence; not causal; requires the prespecified multi-signal diagnostic",
    "candidate-driven molecular follow-up; no causal or therapeutic-target claim"
  )]
  fwrite(out, file.path(integrated_dir, "candidate_molecular_evidence_summary.csv"))
  fwrite(
    out[multisignal_followup_trigger == TRUE],
    file.path(integrated_dir, "frozen_multisignal_followup_trigger_manifest.csv")
  )
  writeLines(
    c(
      "Primary candidate molecular evidence integration",
      "18 combinations were prespecified: 12 GTEx eQTL and 6 circulating pQTL combinations.",
      "SMR FDR is reported within eQTL/pQTL families and globally.",
      "Strict HEIDI support requires P_HEIDI > 0.01 and at least 10 HEIDI SNPs.",
      "ABF strong colocalisation requires PP.H4 >= 0.80 and PP.H4 > PP.H3.",
      "No result is called causal.",
      "Any provisional signal remains provisional until the prespecified multi-signal stage is complete.",
      "deCODE and Fenland are different cohorts but share SomaScan assays and are not cross-platform replication.",
      "No additional cohorts were analysed."
    ),
    file.path(integrated_dir, "README_interpretation.txt"),
    useBytes = TRUE
  )
  write_complete(integrated_dir, "integrate")
  message("Completed integrated evidence stage: ", integrated_dir)
}

message("Track C candidate QTL primary analysis finished")



