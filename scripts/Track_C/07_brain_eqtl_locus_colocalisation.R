options(stringsAsFactors = FALSE)

message("Track C final bounded molecular follow-up: MetaBrain candidates and three frozen loci")

required_packages <- c(
  "data.table", "Rsamtools", "GenomicRanges", "IRanges",
  "coloc", "digest", "jsonlite"
)
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
input_root <- Sys.getenv(
  "TRACKC_FINAL_QTL_INPUT_ROOT",
  unset = file.path(qtl_root, "06_inputs")
)
metabrain_dir <- Sys.getenv(
  "TRACKC_METABRAIN_DIR",
  unset = file.path(input_root, "MetaBrain")
)
gtex_dir <- Sys.getenv(
  "TRACKC_EQTL_CATALOGUE_DIR",
  unset = file.path(input_root, "eQTL Catalogue")
)
out_dir <- Sys.getenv(
  "TRACKC_FINAL_QTL_OUT_DIR",
  unset = file.path(qtl_root, "07_metabrain_three_locus_primary")
)
stages <- unique(trimws(strsplit(
  Sys.getenv(
    "TRACKC_FINAL_QTL_STAGES",
    unset = "index,spec,extract,prepare,abf,integrate"
  ),
  ",", fixed = TRUE
)[[1L]]))
allowed_stages <- c("index", "spec", "extract", "prepare", "abf", "integrate")
if (any(!stages %in% allowed_stages)) {
  stop("Unknown stage(s): ", paste(setdiff(stages, allowed_stages), collapse = ", "))
}

index_dir <- file.path(out_dir, "00_input_QA_and_indices")
spec_dir <- file.path(out_dir, "01_frozen_analysis_spec")
extract_dir <- file.path(out_dir, "02_standardized_QTL")
prepare_dir <- file.path(out_dir, "03_harmonized_factor_QTL")
abf_dir <- file.path(out_dir, "04_coloc_ABF")
integrate_dir <- file.path(out_dir, "05_integrated_primary_evidence")
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

require_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop(label, " is missing: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_complete <- function(directory, stage, extra = character()) {
  writeLines(
    c(
      paste0("stage=", stage),
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      extra
    ),
    file.path(directory, "COMPLETE.txt"),
    useBytes = TRUE
  )
}

safe_id <- function(...) {
  x <- paste(..., sep = "__")
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

is_palindromic <- function(a1, a2) {
  paste0(toupper(a1), toupper(a2)) %in% c("AT", "TA", "CG", "GC")
}

empty_failure <- function() {
  data.table(stage = character(), item_id = character(), error_message = character())
}

meta_columns <- c(
  "Gene", "GeneChr", "GenePos", "GeneStrand", "GeneSymbol", "SNP",
  "SNPChr", "SNPPos", "SNPAlleles", "SNPEffectAllele",
  "SNPEffectAlleleFreq", "MetaP", "MetaPN", "MetaPZ", "MetaBeta", "MetaSE",
  "NrDatasets", "DatasetCorrelationCoefficients", "DatasetZScores", "DatasetSampleSizes"
)
gtex_columns <- c(
  "variant", "r2", "pvalue", "molecular_trait_object_id", "molecular_trait_id",
  "maf", "gene_id", "median_tpm", "beta", "se", "an", "ac", "chromosome",
  "position", "ref", "alt", "type", "rsid"
)

candidate_spec <- data.table(
  gene_symbol = c("PHF2", "GPX1", "TCF4", "TCF4"),
  ensembl_id = c(
    "ENSG00000197724", "ENSG00000233276",
    "ENSG00000196628", "ENSG00000196628"
  ),
  factor_gwas = c("affective", "functional", "affective", "functional"),
  chr_grch38 = c(9L, 3L, 18L, 18L),
  gene_start_grch38 = c(93576366L, 49357174L, 55222185L, 55222185L),
  gene_end_grch38 = c(93679599L, 49358605L, 55664787L, 55664787L),
  analysis_role = "prespecified_candidate_MetaBrain_cortex_followup"
)

# The GRCh38 bounds were mapped from the already frozen GRCh37 bounds with
# Ensembl REST assembly mapping on 2026-07-16. The broad interval is retained
# from the minimum to maximum mapped coordinate; no result informed the map.
locus_spec <- data.table(
  locus_label = c("chr3_functional", "chr11_affective", "chr16_functional"),
  factor_gwas = c("functional", "affective", "functional"),
  chr = c(3L, 11L, 16L),
  grch37_start = c(48052921L, 46376411L, 27477974L),
  grch37_end = c(51233949L, 48870107L, 30002104L),
  grch38_start = c(48011431L, 46354861L, 27466653L),
  grch38_end = c(51196518L, 48848555L, 29990783L),
  ensembl_mapping_segments = c(41L, 3L, 1L),
  evidence_origin = "post_LOLO_prospectively_specified_targeted_followup",
  replication_status = "not_independent_replication"
)

tissue_spec <- data.table(
  resource = c("GTEx_V8", "GTEx_V8", "GTEx_V8"),
  tissue = c("brain_cortex", "spleen", "skeletal_muscle"),
  filename = c("Brain_Cortex.tsv.gz", "Spleen.tsv.gz", "Muscle_Skeletal.tsv.gz"),
  sample_size = c(205L, 227L, 706L),
  role = "prespecified_tissue_sensitivity_for_MetaBrain_defined_locus_genes"
)

meta_file_for_chr <- function(chr) {
  file.path(metabrain_dir, paste0("2021-07-23-cortex-EUR-80PCs-chr", chr, ".txt.gz"))
}

parse_md5_file <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  hit <- regexec("^([0-9a-fA-F]{32})[[:space:]]+(.+)$", lines)
  pieces <- regmatches(lines, hit)
  pieces <- pieces[lengths(pieces) == 3L]
  if (!length(pieces)) return(data.table(filename = character(), expected_md5 = character()))
  rbindlist(lapply(pieces, function(z) data.table(
    filename = basename(trimws(z[3L])), expected_md5 = tolower(z[2L])
  )))
}

scan_local_tabix <- function(path, index, chr, start, end) {
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(chr),
    ranges = IRanges::IRanges(start = as.integer(start), end = as.integer(end))
  )
  tbx <- Rsamtools::TabixFile(path, index = index)
  lines <- Rsamtools::scanTabix(tbx, param = gr)[[1L]]
  lines
}

lines_to_dt <- function(lines, columns) {
  if (!length(lines)) return(as.data.table(setNames(replicate(
    length(columns), character(), simplify = FALSE
  ), columns)))
  x <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE)
  if (ncol(x) != length(columns)) {
    stop("Unexpected column count: ", ncol(x), " rather than ", length(columns))
  }
  setnames(x, columns)
  x
}

standardize_metabrain <- function(x, gene_symbol = NULL, ensembl_id = NULL) {
  if (!nrow(x)) return(data.table())
  x[, `:=`(
    Gene = sub("\\..*$", "", as.character(Gene)),
    GeneChr = as.integer(GeneChr),
    GenePos = as.integer(GenePos),
    GeneSymbol = toupper(trimws(as.character(GeneSymbol))),
    SNPChr = as.integer(SNPChr),
    SNPPos = as.integer(SNPPos),
    SNPEffectAllele = toupper(trimws(as.character(SNPEffectAllele))),
    SNPEffectAlleleFreq = as.numeric(SNPEffectAlleleFreq),
    MetaP = as.numeric(MetaP),
    MetaPN = as.numeric(MetaPN),
    MetaBeta = as.numeric(MetaBeta),
    MetaSE = as.numeric(MetaSE)
  )]
  if (!is.null(gene_symbol)) x <- x[GeneSymbol == toupper(gene_symbol)]
  if (!is.null(ensembl_id)) x <- x[Gene == sub("\\..*$", "", ensembl_id)]
  if (!nrow(x)) return(data.table())
  x[, rsid := sub("^[^:]+:[^:]+:(rs[0-9]+):.*$", "\\1", as.character(SNP))]
  x[!grepl("^rs[0-9]+$", rsid), rsid := NA_character_]
  allele_parts <- tstrsplit(as.character(x$SNPAlleles), "/", fixed = TRUE)
  x[, `:=`(allele_1 = toupper(allele_parts[[1L]]), allele_2 = toupper(allele_parts[[2L]]))]
  x[, other_allele := fifelse(
    SNPEffectAllele == allele_1, allele_2,
    fifelse(SNPEffectAllele == allele_2, allele_1, NA_character_)
  )]
  x[, variant_key := paste(SNPChr, SNPPos, allele_1, allele_2, sep = ":")]
  x[, rsid_variant_n := uniqueN(variant_key), by = rsid]
  out <- x[
    !is.na(rsid) & rsid_variant_n == 1L &
      nchar(SNPEffectAllele) == 1L & nchar(other_allele) == 1L &
      SNPEffectAllele %in% c("A", "C", "G", "T") &
      other_allele %in% c("A", "C", "G", "T") &
      !is_palindromic(SNPEffectAllele, other_allele) &
      is.finite(MetaBeta) & is.finite(MetaSE) & MetaSE > 0 &
      is.finite(MetaP) & MetaP > 0 & MetaP <= 1 &
      is.finite(SNPEffectAlleleFreq) & SNPEffectAlleleFreq > 0 & SNPEffectAlleleFreq < 1,
    .(
      gene_symbol = GeneSymbol,
      ensembl_id = Gene,
      gene_chr = GeneChr,
      gene_tss_grch38 = GenePos,
      resource = "MetaBrain_2021_07_23",
      tissue = "cortex_EUR",
      build = 38L,
      SNP = rsid,
      CHR = SNPChr,
      BP = SNPPos,
      effect_allele = SNPEffectAllele,
      other_allele = other_allele,
      eaf = SNPEffectAlleleFreq,
      MAF = pmin(SNPEffectAlleleFreq, 1 - SNPEffectAlleleFreq),
      beta = MetaBeta,
      se = MetaSE,
      pvalue = MetaP,
      N = MetaPN
    )
  ]
  setorder(out, SNP, pvalue)
  out[!duplicated(SNP)]
}

standardize_gtex <- function(x, genes, tissue, sample_size, locus_label) {
  if (!nrow(x)) return(data.table())
  x[, `:=`(
    gene_id = sub("\\..*$", "", as.character(gene_id)),
    chromosome = as.integer(chromosome),
    position = as.integer(position),
    ref = toupper(as.character(ref)),
    alt = toupper(as.character(alt)),
    rsid = trimws(as.character(rsid)),
    beta = as.numeric(beta),
    se = as.numeric(se),
    pvalue = as.numeric(pvalue),
    maf = as.numeric(maf),
    ac = as.numeric(ac),
    an = as.numeric(an)
  )]
  x <- x[gene_id %in% genes$ensembl_id]
  if (!nrow(x)) return(data.table())
  x[, eaf := ac / an]
  x[, variant_key := paste(chromosome, position, ref, alt, sep = ":")]
  x[, rsid_variant_n := uniqueN(variant_key), by = rsid]
  x <- merge(
    x, genes[, .(ensembl_id, gene_symbol, gene_tss_grch38)],
    by.x = "gene_id", by.y = "ensembl_id", all.x = TRUE, sort = FALSE
  )
  out <- x[
    grepl("^rs[0-9]+$", rsid) & rsid_variant_n == 1L &
      nchar(ref) == 1L & nchar(alt) == 1L &
      ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T") &
      !is_palindromic(alt, ref) &
      is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(pvalue) & pvalue > 0 & pvalue <= 1 &
      is.finite(eaf) & eaf > 0 & eaf < 1,
    .(
      gene_symbol, ensembl_id = gene_id,
      gene_chr = chromosome, gene_tss_grch38,
      resource = "eQTL_Catalogue_imported_GTEx_V8",
      tissue = tissue, build = 38L,
      SNP = rsid, CHR = chromosome, BP = position,
      effect_allele = alt, other_allele = ref,
      eaf, MAF = pmin(eaf, 1 - eaf),
      beta, se, pvalue, N = as.numeric(sample_size),
      locus_label
    )
  ]
  setorder(out, ensembl_id, SNP, pvalue)
  out[!duplicated(out, by = c("ensembl_id", "SNP"))]
}

read_factor <- function(factor_name, wanted_snps) {
  path <- require_file(factor_files[[factor_name]], paste(factor_name, "factor GWAS"))
  x <- fread(path, select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"))
  x <- x[SNP %in% wanted_snps]
  x[, `:=`(
    SNP = as.character(SNP), CHR = as.integer(CHR), BP = as.integer(BP),
    A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)),
    BETA = as.numeric(BETA), SE = as.numeric(SE),
    P = as.numeric(P), N = as.numeric(N)
  )]
  x[
    grepl("^rs[0-9]+$", SNP) &
      nchar(A1) == 1L & nchar(A2) == 1L &
      !is_palindromic(A1, A2) &
      is.finite(BETA) & is.finite(SE) & SE > 0 &
      is.finite(P) & P > 0 & P <= 1 & is.finite(N) & N > 0
  ][!duplicated(SNP)]
}

harmonize_factor_qtl <- function(gwas, qtl) {
  if (!nrow(gwas) || !nrow(qtl)) return(data.table())
  g <- copy(gwas)
  q <- copy(qtl)
  setnames(g, c("BETA", "SE", "P", "N"), c(
    "gwas_beta_input", "gwas_se_input", "gwas_p_input", "gwas_n_input"
  ))
  setnames(q, c("beta", "se", "pvalue", "N", "BP", "CHR"), c(
    "qtl_beta_input", "qtl_se_input", "qtl_p_input", "qtl_n_input",
    "qtl_BP", "qtl_CHR"
  ))
  z <- merge(g, q, by = "SNP", all = FALSE, sort = FALSE)
  same <- z$A1 == z$effect_allele & z$A2 == z$other_allele
  reverse <- z$A1 == z$other_allele & z$A2 == z$effect_allele
  z <- z[same | reverse]
  reverse <- reverse[same | reverse]
  z[, `:=`(
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
  z <- z[
    is.finite(beta_gwas_aligned) & is.finite(se_gwas) & se_gwas > 0 &
      is.finite(beta_qtl) & is.finite(se_qtl) & se_qtl > 0 &
      is.finite(maf_qtl) & maf_qtl > 0 & maf_qtl < 0.5 &
      is.finite(eaf_qtl) & eaf_qtl > 0 & eaf_qtl < 1
  ]
  z[!duplicated(SNP)]
}

if ("index" %in% stages) {
  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  md5_file <- require_file(file.path(metabrain_dir, "MD5 sums.txt"), "MetaBrain MD5 list")
  md5 <- parse_md5_file(md5_file)
  expected_chr <- c(3L, 9L, 11L, 16L, 18L)
  qa_rows <- list()
  for (chr in expected_chr) {
    f <- require_file(meta_file_for_chr(chr), paste0("MetaBrain chr", chr))
    expected <- md5[filename == basename(f), expected_md5][1L]
    actual <- tolower(digest::digest(f, algo = "md5", file = TRUE))
    if (!length(expected) || is.na(expected)) stop("No official MD5 found for ", basename(f))
    if (!identical(expected, actual)) stop("MD5 mismatch for ", basename(f))
    idx <- paste0(f, ".tbi")
    generated <- FALSE
    if (!file.exists(idx)) {
      message("Creating local MetaBrain gene-position index: ", basename(f))
      Rsamtools::indexTabix(
        f, seq = 2L, start = 3L, end = 3L,
        skip = 1L, comment = "#", zeroBased = FALSE
      )
      generated <- TRUE
    }
    require_file(idx, paste0("MetaBrain chr", chr, " index"))
    qa_rows[[length(qa_rows) + 1L]] <- data.table(
      resource = "MetaBrain", chromosome = chr, file = f,
      bytes = as.numeric(file.info(f)$size), expected_md5 = expected,
      observed_md5 = actual, md5_match = TRUE,
      index_file = normalizePath(idx, winslash = "/", mustWork = TRUE),
      index_generated_locally = generated
    )
  }
  for (i in seq_len(nrow(tissue_spec))) {
    f <- require_file(file.path(gtex_dir, tissue_spec$filename[i]), tissue_spec$tissue[i])
    idx <- require_file(paste0(f, ".tbi"), paste0(tissue_spec$tissue[i], " index"))
    qa_rows[[length(qa_rows) + 1L]] <- data.table(
      resource = "GTEx_V8", chromosome = NA_integer_, file = f,
      bytes = as.numeric(file.info(f)$size), expected_md5 = NA_character_,
      observed_md5 = NA_character_, md5_match = NA,
      index_file = idx, index_generated_locally = FALSE
    )
  }
  qa <- rbindlist(qa_rows, fill = TRUE)
  fwrite(qa, file.path(index_dir, "input_file_and_index_QA.csv"))
  write_complete(index_dir, "index")
  message("Completed index/input-QA stage: ", index_dir)
}

if ("spec" %in% stages) {
  dir.create(spec_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(candidate_spec, file.path(spec_dir, "frozen_MetaBrain_candidate_scope.csv"))
  fwrite(locus_spec, file.path(spec_dir, "frozen_three_locus_coordinates_both_builds.csv"))
  fwrite(tissue_spec, file.path(spec_dir, "frozen_GTEx_tissue_sensitivity_scope.csv"))
  writeLines(
    c(
      "Final bounded Track C molecular follow-up.",
      "MetaBrain 2021-07-23 Cortex-EUR is the primary brain eQTL resource.",
      "Candidate scope is restricted to PHF2-affective, GPX1-functional and TCF4-affective/functional.",
      "Locus scope is restricted to the three post-LOLO frozen regions.",
      "All MetaBrain genes with indexed gene position inside each GRCh38 region are retained in the coverage audit.",
      "GTEx cortex, spleen and skeletal muscle are sensitivity resources for the same MetaBrain-defined locus-gene universe.",
      "Formal ABF interpretation requires cis-eQTL P<5e-8, regional factor-GWAS P<1e-5 and >=50 harmonized SNPs.",
      "Coloc posterior probabilities are not BH-adjusted.",
      "Primary shared-signal rule: PP.H4>=0.80 and H4/H3>=3; p12=1e-6 is reported as prior sensitivity.",
      "SuSiE is not run mechanically. A frozen trigger manifest is generated for H3>=0.80 or H4>=0.80 results with >=100 harmonized SNPs.",
      "No new gene, tissue, locus, disease endpoint or QTL resource may be added.",
      "No additional cohorts or proteomic resources were analysed."
    ),
    file.path(spec_dir, "README.txt"), useBytes = TRUE
  )
  write_complete(spec_dir, "spec")
  message("Completed specification stage: ", spec_dir)
}

if ("extract" %in% stages) {
  require_file(file.path(index_dir, "COMPLETE.txt"), "index completion marker")
  require_file(file.path(spec_dir, "COMPLETE.txt"), "spec completion marker")
  meta_out <- file.path(extract_dir, "MetaBrain")
  gtex_out <- file.path(extract_dir, "GTEx_sensitivity")
  dir.create(meta_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(gtex_out, recursive = TRUE, showWarnings = FALSE)
  failures <- list()
  meta_manifest <- list()

  # Candidate coverage is recorded even when the gene is absent or lacks a strong cis-eQTL.
  candidate_unique <- unique(candidate_spec[, .(
    gene_symbol, ensembl_id, chr_grch38, gene_start_grch38, gene_end_grch38
  )])
  for (i in seq_len(nrow(candidate_unique))) {
    row <- candidate_unique[i]
    result <- tryCatch({
      f <- meta_file_for_chr(row$chr_grch38)
      lines <- scan_local_tabix(
        f, paste0(f, ".tbi"), row$chr_grch38,
        row$gene_start_grch38, row$gene_end_grch38
      )
      raw <- lines_to_dt(lines, meta_columns)
      q <- standardize_metabrain(raw, row$gene_symbol, row$ensembl_id)
      outfile <- file.path(meta_out, paste0("candidate__", row$gene_symbol, ".tsv.gz"))
      if (nrow(q)) fwrite(q, outfile, sep = "\t", compress = "gzip")
      data.table(
        analysis_role = "candidate", locus_label = NA_character_,
        factor_gwas = NA_character_, gene_symbol = row$gene_symbol,
        ensembl_id = row$ensembl_id, gene_tss_grch38 = if (nrow(q)) q$gene_tss_grch38[1L] else NA_integer_,
        resource = "MetaBrain_2021_07_23", tissue = "cortex_EUR",
        qtl_file = if (nrow(q)) normalizePath(outfile, winslash = "/", mustWork = TRUE) else NA_character_,
        qtl_row_n = nrow(q), unique_snp_n = if (nrow(q)) uniqueN(q$SNP) else 0L,
        min_qtl_p = if (nrow(q)) min(q$pvalue) else NA_real_,
        max_qtl_n = if (nrow(q)) max(q$N) else NA_real_
      )
    }, error = function(e) e)
    if (inherits(result, "error")) {
      failures[[length(failures) + 1L]] <- data.table(
        stage = "extract_MetaBrain_candidate", item_id = row$gene_symbol,
        error_message = conditionMessage(result)
      )
    } else meta_manifest[[length(meta_manifest) + 1L]] <- result
  }

  # MetaBrain defines the complete primary gene universe by gene position.
  for (i in seq_len(nrow(locus_spec))) {
    loc <- locus_spec[i]
    result <- tryCatch({
      f <- meta_file_for_chr(loc$chr)
      lines <- scan_local_tabix(
        f, paste0(f, ".tbi"), loc$chr, loc$grch38_start, loc$grch38_end
      )
      raw <- lines_to_dt(lines, meta_columns)
      q_all <- standardize_metabrain(raw)
      if (nrow(q_all)) {
        q_all <- q_all[
          gene_tss_grch38 >= loc$grch38_start & gene_tss_grch38 <= loc$grch38_end
        ]
      }
      if (!nrow(q_all)) stop("No standardized MetaBrain genes in frozen locus")
      for (ens in unique(q_all$ensembl_id)) {
        q <- q_all[ensembl_id == ens]
        sym <- q$gene_symbol[1L]
        outfile <- file.path(meta_out, paste0(
          "locus__", loc$locus_label, "__", safe_id(sym, ens), ".tsv.gz"
        ))
        fwrite(q, outfile, sep = "\t", compress = "gzip")
        meta_manifest[[length(meta_manifest) + 1L]] <- data.table(
          analysis_role = "locus", locus_label = loc$locus_label,
          factor_gwas = loc$factor_gwas, gene_symbol = sym,
          ensembl_id = ens, gene_tss_grch38 = q$gene_tss_grch38[1L],
          resource = "MetaBrain_2021_07_23", tissue = "cortex_EUR",
          qtl_file = normalizePath(outfile, winslash = "/", mustWork = TRUE),
          qtl_row_n = nrow(q), unique_snp_n = uniqueN(q$SNP),
          min_qtl_p = min(q$pvalue), max_qtl_n = max(q$N)
        )
      }
      TRUE
    }, error = function(e) e)
    if (inherits(result, "error")) failures[[length(failures) + 1L]] <- data.table(
      stage = "extract_MetaBrain_locus", item_id = loc$locus_label,
      error_message = conditionMessage(result)
    )
  }

  meta_dt <- rbindlist(meta_manifest, fill = TRUE)
  locus_genes <- unique(meta_dt[analysis_role == "locus", .(
    locus_label, factor_gwas, gene_symbol, ensembl_id, gene_tss_grch38
  )])
  fwrite(meta_dt, file.path(extract_dir, "MetaBrain_candidate_and_locus_gene_manifest.csv"))
  fwrite(locus_genes, file.path(extract_dir, "MetaBrain_defined_locus_gene_universe.csv"))

  # GTEx region queries include +/-1 Mb so that all cis variants for TSS-inside genes are available.
  gtex_manifest <- list()
  for (ti in seq_len(nrow(tissue_spec))) {
    tissue <- tissue_spec[ti]
    f <- file.path(gtex_dir, tissue$filename)
    idx <- paste0(f, ".tbi")
    for (li in seq_len(nrow(locus_spec))) {
      loc <- locus_spec[li]
      genes <- locus_genes[locus_label == loc$locus_label]
      item <- paste(tissue$tissue, loc$locus_label, sep = "__")
      result <- tryCatch({
        lines <- scan_local_tabix(
          f, idx, loc$chr,
          max(1L, loc$grch38_start - 1000000L), loc$grch38_end + 1000000L
        )
        raw <- lines_to_dt(lines, gtex_columns)
        q_all <- standardize_gtex(
          raw, genes, tissue$tissue, tissue$sample_size, loc$locus_label
        )
        for (ens in genes$ensembl_id) {
          q <- q_all[ensembl_id == ens]
          gene_row <- genes[ensembl_id == ens][1L]
          outfile <- file.path(gtex_out, paste0(
            "locus__", loc$locus_label, "__", tissue$tissue, "__",
            safe_id(gene_row$gene_symbol, ens), ".tsv.gz"
          ))
          if (nrow(q)) fwrite(q, outfile, sep = "\t", compress = "gzip")
          gtex_manifest[[length(gtex_manifest) + 1L]] <- data.table(
            analysis_role = "locus", locus_label = loc$locus_label,
            factor_gwas = loc$factor_gwas, gene_symbol = gene_row$gene_symbol,
            ensembl_id = ens, gene_tss_grch38 = gene_row$gene_tss_grch38,
            resource = "eQTL_Catalogue_imported_GTEx_V8", tissue = tissue$tissue,
            qtl_file = if (nrow(q)) normalizePath(outfile, winslash = "/", mustWork = TRUE) else NA_character_,
            qtl_row_n = nrow(q), unique_snp_n = if (nrow(q)) uniqueN(q$SNP) else 0L,
            min_qtl_p = if (nrow(q)) min(q$pvalue) else NA_real_,
            max_qtl_n = if (nrow(q)) max(q$N) else as.numeric(tissue$sample_size)
          )
        }
        TRUE
      }, error = function(e) e)
      if (inherits(result, "error")) failures[[length(failures) + 1L]] <- data.table(
        stage = "extract_GTEx_locus", item_id = item,
        error_message = conditionMessage(result)
      )
    }
  }
  gtex_dt <- rbindlist(gtex_manifest, fill = TRUE)
  fwrite(gtex_dt, file.path(extract_dir, "GTEx_locus_gene_coverage_manifest.csv"))
  failure_dt <- if (length(failures)) rbindlist(failures, fill = TRUE) else empty_failure()
  fwrite(failure_dt, file.path(extract_dir, "QTL_extraction_failures.csv"))
  if (nrow(failure_dt)) stop("QTL extraction has failures; inspect QTL_extraction_failures.csv")
  write_complete(
    extract_dir, "extract",
    c(
      paste0("metabrain_locus_gene_n=", uniqueN(locus_genes$ensembl_id)),
      paste0("metabrain_strong_locus_gene_n=", uniqueN(meta_dt[analysis_role == "locus" & min_qtl_p < 5e-8, ensembl_id]))
    )
  )
  message("Completed QTL extraction stage: ", extract_dir)
}

if ("prepare" %in% stages) {
  require_file(file.path(extract_dir, "COMPLETE.txt"), "extract completion marker")
  dir.create(prepare_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- fread(file.path(extract_dir, "MetaBrain_candidate_and_locus_gene_manifest.csv"))
  gtex <- fread(file.path(extract_dir, "GTEx_locus_gene_coverage_manifest.csv"))

  meta_candidate <- merge(
    candidate_spec,
    meta[analysis_role == "candidate", .(
      gene_symbol, ensembl_id, resource, tissue, qtl_file, qtl_row_n,
      unique_snp_n, min_qtl_p, max_qtl_n
    )],
    by = c("gene_symbol", "ensembl_id"), all.x = TRUE, allow.cartesian = TRUE
  )
  meta_candidate[, `:=`(
    locus_label = NA_character_, analysis_role = "candidate",
    gene_tss_grch38 = NA_integer_
  )]
  meta_locus <- meta[analysis_role == "locus"]
  combos <- rbindlist(list(
    meta_candidate[, .(
      analysis_role, locus_label, factor_gwas, gene_symbol, ensembl_id,
      gene_tss_grch38, resource, tissue, qtl_file, qtl_row_n,
      unique_snp_n, min_qtl_p, max_qtl_n
    )],
    meta_locus,
    gtex
  ), fill = TRUE)
  combos[, combo_id := safe_id(
    analysis_role,
    fifelse(is.na(locus_label), "no_locus", locus_label),
    factor_gwas, resource, tissue, gene_symbol, ensembl_id
  )]
  if (anyDuplicated(combos$combo_id)) stop("Duplicate combination IDs")
  # Keep the descriptive combination identifier in all audit tables, but never
  # use it as a Windows file or directory name.  Some GTEx resource/tissue/gene
  # combinations exceed the traditional Windows MAX_PATH limit.  The stable,
  # sequential file_id is deliberately assigned only after the complete frozen
  # combination table has been constructed.
  combos[, file_id := sprintf("C%04d", .I)]
  combos[, strong_cis_qtl := is.finite(min_qtl_p) & min_qtl_p < 5e-8]
  fwrite(combos, file.path(prepare_dir, "complete_candidate_and_locus_combination_manifest.csv"))
  fwrite(
    combos[, .(
      file_id, combo_id, analysis_role, locus_label, factor_gwas,
      resource, tissue, gene_symbol, ensembl_id
    )],
    file.path(prepare_dir, "short_file_id_to_full_combination_map.csv")
  )

  strong_files <- unique(combos[strong_cis_qtl == TRUE & !is.na(qtl_file), qtl_file])
  wanted <- unique(unlist(lapply(strong_files, function(f) fread(f, select = "SNP")$SNP)))
  factor_cache <- list(
    affective = read_factor("affective", wanted),
    functional = read_factor("functional", wanted)
  )
  qa_rows <- list()
  failures <- list()
  for (i in seq_len(nrow(combos))) {
    row <- combos[i]
    result <- tryCatch({
      if (is.na(row$qtl_file) || !file.exists(row$qtl_file) || !isTRUE(row$strong_cis_qtl)) {
        data.table(
          combo_id = row$combo_id,
          harmonized_file = NA_character_, harmonized_snp_n = 0L,
          min_harmonized_qtl_p = row$min_qtl_p,
          min_harmonized_gwas_p = NA_real_, regional_gwas_p1e5 = FALSE,
          eligible_coloc_interpretation = FALSE,
          prepare_status = if (is.na(row$qtl_file)) "no_qtl_coverage" else "no_strong_cis_qtl"
        )
      } else {
        q <- fread(row$qtl_file)
        g <- factor_cache[[row$factor_gwas]]
        h <- harmonize_factor_qtl(g, q)
        # Use a new short-name subdirectory so that files left by an interrupted
        # pre-fix run cannot be mistaken for current formal inputs.
        hfile <- file.path(prepare_dir, "harmonized_short", paste0(row$file_id, ".tsv.gz"))
        dir.create(dirname(hfile), recursive = TRUE, showWarnings = FALSE)
        if (nrow(h)) fwrite(h, hfile, sep = "\t", compress = "gzip")
        min_q <- if (nrow(h)) min(h$p_qtl) else NA_real_
        min_g <- if (nrow(h)) min(h$p_gwas) else NA_real_
        eligible <- nrow(h) >= 50L && is.finite(min_q) && min_q < 5e-8 &&
          is.finite(min_g) && min_g < 1e-5
        data.table(
          combo_id = row$combo_id,
          harmonized_file = if (nrow(h)) normalizePath(hfile, winslash = "/", mustWork = TRUE) else NA_character_,
          harmonized_snp_n = nrow(h),
          min_harmonized_qtl_p = min_q,
          min_harmonized_gwas_p = min_g,
          regional_gwas_p1e5 = is.finite(min_g) && min_g < 1e-5,
          eligible_coloc_interpretation = eligible,
          prepare_status = if (!nrow(h)) "no_harmonized_snps" else if (nrow(h) < 50L) {
            "sparse_harmonized_overlap"
          } else if (!is.finite(min_g) || min_g >= 1e-5) {
            "no_regional_factor_gwas_signal"
          } else "eligible"
        )
      }
    }, error = function(e) e)
    if (inherits(result, "error")) failures[[length(failures) + 1L]] <- data.table(
      stage = "prepare", item_id = row$combo_id, error_message = conditionMessage(result)
    ) else qa_rows[[length(qa_rows) + 1L]] <- result
  }
  qa <- merge(combos, rbindlist(qa_rows, fill = TRUE), by = "combo_id", all.x = TRUE)
  failure_dt <- if (length(failures)) rbindlist(failures, fill = TRUE) else empty_failure()
  fwrite(qa, file.path(prepare_dir, "harmonization_and_eligibility_QA.csv"))
  fwrite(failure_dt, file.path(prepare_dir, "harmonization_failures.csv"))
  if (nrow(failure_dt)) stop("Preparation has failures; inspect harmonization_failures.csv")
  write_complete(
    prepare_dir, "prepare",
    c(
      paste0("combination_n=", nrow(qa)),
      paste0("eligible_coloc_n=", sum(qa$eligible_coloc_interpretation, na.rm = TRUE))
    )
  )
  message("Completed preparation/harmonization stage: ", prepare_dir)
}

run_abf <- function(row) {
  if (!isTRUE(row$eligible_coloc_interpretation)) {
    return(data.table(combo_id = row$combo_id, abf_status = "not_eligible"))
  }
  h <- fread(row$harmonized_file)
  if (nrow(h) < 50L) return(data.table(combo_id = row$combo_id, abf_status = "sparse"))
  d_gwas <- list(
    beta = h$beta_gwas_aligned, varbeta = h$se_gwas^2,
    snp = h$SNP, type = "quant", sdY = 1
  )
  qtl_sdY <- sqrt(stats::median(
    2 * h$n_qtl * h$maf_qtl * (1 - h$maf_qtl) * h$se_qtl^2,
    na.rm = TRUE
  ))
  if (!is.finite(qtl_sdY) || qtl_sdY <= 0) qtl_sdY <- 1
  d_qtl <- list(
    beta = h$beta_qtl, varbeta = h$se_qtl^2,
    snp = h$SNP, type = "quant", sdY = qtl_sdY
  )
  primary <- coloc::coloc.abf(d_gwas, d_qtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  sensitivity <- coloc::coloc.abf(d_gwas, d_qtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-6)
  combo_dir <- file.path(abf_dir, row$file_id)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(primary, file.path(combo_dir, "coloc_abf_primary.rds"))
  saveRDS(sensitivity, file.path(combo_dir, "coloc_abf_p12_1e6.rds"))
  fwrite(
    as.data.table(primary$results),
    file.path(combo_dir, "coloc_abf_snp_results.tsv.gz"),
    sep = "\t", compress = "gzip"
  )
  p <- as.list(primary$summary)
  s <- as.list(sensitivity$summary)
  data.table(
    combo_id = row$combo_id, abf_status = "completed",
    nsnps = as.integer(p[["nsnps"]]),
    pp_h0 = as.numeric(p[["PP.H0.abf"]]),
    pp_h1 = as.numeric(p[["PP.H1.abf"]]),
    pp_h2 = as.numeric(p[["PP.H2.abf"]]),
    pp_h3 = as.numeric(p[["PP.H3.abf"]]),
    pp_h4 = as.numeric(p[["PP.H4.abf"]]),
    pp_h3_p12_1e6 = as.numeric(s[["PP.H3.abf"]]),
    pp_h4_p12_1e6 = as.numeric(s[["PP.H4.abf"]]),
    qtl_sdY_estimate = qtl_sdY
  )
}

if ("abf" %in% stages) {
  require_file(file.path(prepare_dir, "COMPLETE.txt"), "prepare completion marker")
  dir.create(abf_dir, recursive = TRUE, showWarnings = FALSE)
  qa <- fread(file.path(prepare_dir, "harmonization_and_eligibility_QA.csv"))
  results <- list()
  failures <- list()
  for (i in seq_len(nrow(qa))) {
    row <- qa[i]
    result <- tryCatch(run_abf(row), error = function(e) e)
    if (inherits(result, "error")) failures[[length(failures) + 1L]] <- data.table(
      stage = "abf", item_id = row$combo_id, error_message = conditionMessage(result)
    ) else results[[length(results) + 1L]] <- result
  }
  result_dt <- rbindlist(results, fill = TRUE)
  failure_dt <- if (length(failures)) rbindlist(failures, fill = TRUE) else empty_failure()
  fwrite(result_dt, file.path(abf_dir, "all_coloc_abf_results.csv"))
  fwrite(failure_dt, file.path(abf_dir, "coloc_abf_failures.csv"))
  if (nrow(failure_dt)) stop("ABF stage has failures; inspect coloc_abf_failures.csv")
  write_complete(
    abf_dir, "abf",
    paste0("completed_abf_n=", sum(result_dt$abf_status == "completed", na.rm = TRUE))
  )
  message("Completed ABF stage: ", abf_dir)
}

if ("integrate" %in% stages) {
  require_file(file.path(abf_dir, "COMPLETE.txt"), "ABF completion marker")
  dir.create(integrate_dir, recursive = TRUE, showWarnings = FALSE)
  qa <- fread(file.path(prepare_dir, "harmonization_and_eligibility_QA.csv"))
  abf <- fread(file.path(abf_dir, "all_coloc_abf_results.csv"))
  final <- merge(qa, abf, by = "combo_id", all.x = TRUE)
  final[, h4_h3_ratio := fifelse(
    is.finite(pp_h4) & is.finite(pp_h3), pp_h4 / pmax(pp_h3, 1e-300), NA_real_
  )]
  final[, `:=`(
    primary_strong_shared = abf_status == "completed" &
      pp_h4 >= 0.80 & h4_h3_ratio >= 3,
    prior_sensitivity_support = abf_status == "completed" & pp_h4_p12_1e6 >= 0.50,
    distinct_signal_support = abf_status == "completed" & pp_h3 >= 0.80
  )]
  final[, evidence_class := fcase(
    is.na(qtl_file), "no_qtl_coverage",
    strong_cis_qtl == FALSE, "no_strong_cis_qtl",
    prepare_status == "no_regional_factor_gwas_signal", "strong_qtl_without_regional_factor_gwas_signal",
    prepare_status == "sparse_harmonized_overlap", "insufficient_harmonized_overlap",
    abf_status == "completed" & primary_strong_shared & prior_sensitivity_support,
      "strong_shared_regulatory_signal_pending_multisignal_review",
    abf_status == "completed" & distinct_signal_support, "distinct_local_signals",
    abf_status == "completed", "molecularly_unresolved",
    default = "not_estimable"
  )]
  final[, susie_trigger := abf_status == "completed" & nsnps >= 100L &
    (pp_h3 >= 0.80 | pp_h4 >= 0.80)]

  triggers <- final[susie_trigger == TRUE, .(
    combo_id, analysis_role, locus_label, factor_gwas, gene_symbol, ensembl_id,
    resource, tissue, harmonized_file, nsnps, pp_h3, pp_h4, pp_h4_p12_1e6,
    trigger_reason = fifelse(pp_h4 >= 0.80, "single_signal_H4_ge_0.80", "single_signal_H3_ge_0.80")
  )]
  fwrite(final, file.path(integrate_dir, "final_primary_candidate_and_locus_evidence.csv"))
  fwrite(final[analysis_role == "candidate"], file.path(integrate_dir, "MetaBrain_candidate_summary.csv"))
  fwrite(final[analysis_role == "locus"], file.path(integrate_dir, "three_locus_all_gene_results.csv"))
  fwrite(triggers, file.path(integrate_dir, "frozen_susie_trigger_manifest.csv"))
  summary <- final[, .(
    combination_n = .N,
    qtl_covered_n = sum(!is.na(qtl_file)),
    strong_cis_qtl_n = sum(strong_cis_qtl, na.rm = TRUE),
    eligible_coloc_n = sum(eligible_coloc_interpretation, na.rm = TRUE),
    completed_abf_n = sum(abf_status == "completed", na.rm = TRUE),
    strong_shared_n = sum(evidence_class == "strong_shared_regulatory_signal_pending_multisignal_review"),
    distinct_signal_n = sum(evidence_class == "distinct_local_signals"),
    susie_trigger_n = sum(susie_trigger, na.rm = TRUE)
  ), by = .(analysis_role, locus_label, resource, tissue, factor_gwas)]
  fwrite(summary, file.path(integrate_dir, "primary_evidence_summary_by_family.csv"))
  audit <- list(
    stage = "MetaBrain_and_three_locus_primary",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    primary_analysis_complete = TRUE,
    susie_trigger_n = nrow(triggers),
    qtl_stage_closed = nrow(triggers) == 0L,
    interpretation = if (nrow(triggers)) {
      "Run only the frozen SuSiE triggers before final locus classification."
    } else {
      "No SuSiE trigger; final bounded QTL stage can be closed."
    },
    prohibited_expansion = c(
      "new genes", "new tissues", "new loci", "new QTL resources",
      "threshold relaxation", "new cellular-aging plasma proteomics"
    )
  )
  jsonlite::write_json(
    audit, file.path(integrate_dir, "primary_molecular_followup_audit.json"),
    pretty = TRUE, auto_unbox = TRUE
  )
  write_complete(
    integrate_dir, "integrate",
    c(
      "primary_analysis_complete=true",
      paste0("susie_trigger_n=", nrow(triggers)),
      paste0("qtl_stage_closed=", tolower(as.character(nrow(triggers) == 0L)))
    )
  )
  message("Completed primary molecular integration: ", integrate_dir)
  message("Frozen SuSiE trigger count: ", nrow(triggers))
}

message("Track C final bounded primary molecular follow-up finished.")


