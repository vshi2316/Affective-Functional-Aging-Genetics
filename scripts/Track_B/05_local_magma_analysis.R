options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir) || !dir.exists(base_dir)) {
  stop("Set DCV_BASE_DIR to the project directory containing analysis_ready_core.")
}
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)

trackb_dir <- file.path(base_dir, "analysis_ready_core", "trackB")
one_file <- file.path(trackb_dir, "factor_gwas_qsnp", "trackB_afffunc_factor_gwas_qsnp.csv")
two_file <- file.path(trackb_dir, "twofactor_gwas_qsnp", "trackB_twofactor_gwas_qsnp.csv")
out_dir <- file.path(trackb_dir, "three_projection_locus_magma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

resolve_executable <- function(path, label) {
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  located <- Sys.which(path)
  if (nzchar(located)) return(unname(located))
  stop(label, " was not found. Set its environment variable in config/config.R.")
}

plink_exe <- resolve_executable(Sys.getenv("PLINK_EXECUTABLE", unset = "plink"), "PLINK")
magma_exe <- resolve_executable(Sys.getenv("MAGMA_EXECUTABLE", unset = "magma"), "MAGMA")
reference_root <- Sys.getenv("MAGMA_REFERENCE_DIR", unset = file.path(base_dir, "reference", "magma"))
bfile <- file.path(reference_root, "g1000_eur", "g1000_eur")
gene_loc <- file.path(reference_root, "ENSGv110.coding.genes.txt")
set_annot <- file.path(reference_root, "MSigDB_20231Hs_MAGMA.txt")
gtex_covar <- file.path(reference_root, "gtex_ts_avg_log2TPM.txt")

required <- c(
  one_file, two_file, plink_exe, magma_exe,
  paste0(bfile, ".bed"), paste0(bfile, ".bim"), paste0(bfile, ".fam"),
  gene_loc, set_annot, gtex_covar
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) stop("Missing required file(s): ", paste(missing, collapse = "; "))

one <- fread(one_file, select = c("CHR", "BP", "SNP", "A1", "A2", "MAF", "beta_factor", "se_factor", "p_factor"))
two <- fread(two_file, select = c(
  "CHR", "BP", "SNP", "A1", "A2", "MAF",
  "beta_affective", "se_affective", "p_affective",
  "beta_functional", "se_functional", "p_functional"
))

analyses <- list(
  common = one[, .(CHR, BP, SNP, A1, A2, MAF, BETA = beta_factor, SE = se_factor, P = p_factor)],
  affective = two[, .(CHR, BP, SNP, A1, A2, MAF, BETA = beta_affective, SE = se_affective, P = p_affective)],
  functional = two[, .(CHR, BP, SNP, A1, A2, MAF, BETA = beta_functional, SE = se_functional, P = p_functional)]
)

run_command <- function(executable, arguments, log_file) {
  output <- system2(executable, args = arguments, stdout = TRUE, stderr = TRUE)
  writeLines(output, log_file)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) stop("External command failed. See: ", log_file)
  invisible(output)
}

read_clumped <- function(path, analysis_name) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.table())
  result <- tryCatch(fread(path, fill = TRUE), error = function(e) data.table())
  if (nrow(result) == 0) return(result)
  result[, analysis := analysis_name]
  result[, in_extended_mhc := CHR == 6 & BP >= 25000000 & BP <= 34000000]
  result
}

summary_rows <- list()
clump_rows <- list()

for (analysis_name in names(analyses)) {
  dat <- copy(analyses[[analysis_name]])
  dat <- dat[
    is.finite(CHR) & is.finite(BP) & !is.na(SNP) & nzchar(SNP) &
      is.finite(P) & P > 0 & P <= 1 & is.finite(BETA) & is.finite(SE) & SE > 0 &
      is.finite(MAF) & MAF > 0 & MAF < 1
  ]
  dat[, N_EFF := 1 / (2 * MAF * (1 - MAF) * SE^2)]
  plausible_n <- dat[N_EFF >= 1000 & N_EFF <= 5000000, N_EFF]
  median_n <- if (length(plausible_n) > 0) round(median(plausible_n, na.rm = TRUE)) else NA_real_
  dat[, N := median_n]
  dat[, in_extended_mhc := CHR == 6 & BP >= 25000000 & BP <= 34000000]

  analysis_dir <- file.path(out_dir, analysis_name)
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

  fuma_file <- file.path(analysis_dir, paste0("trackB_", analysis_name, "_fuma_input.tsv.gz"))
  magma_file <- file.path(analysis_dir, paste0("trackB_", analysis_name, "_magma_input.txt"))
  clump_input <- file.path(analysis_dir, paste0("trackB_", analysis_name, "_plink_clump_input.txt"))
  fwrite(dat[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, N)], fuma_file, sep = "\t")
  fwrite(dat[, .(SNP, P, N)], magma_file, sep = "\t")
  fwrite(dat[, .(SNP, P)], clump_input, sep = "\t")

  clump_prefix <- file.path(analysis_dir, paste0("trackB_", analysis_name, "_plink"))
  run_command(
    plink_exe,
    c(
      "--bfile", shQuote(bfile),
      "--clump", shQuote(clump_input),
      "--clump-snp-field", "SNP",
      "--clump-field", "P",
      "--clump-p1", "5e-8",
      "--clump-p2", "1e-5",
      "--clump-r2", "0.1",
      "--clump-kb", "500",
      "--out", shQuote(clump_prefix)
    ),
    file.path(analysis_dir, paste0("trackB_", analysis_name, "_plink.log.txt"))
  )
  clumped <- read_clumped(paste0(clump_prefix, ".clumped"), analysis_name)
  if (nrow(clumped) > 0) {
    fwrite(clumped, file.path(analysis_dir, paste0("trackB_", analysis_name, "_independent_lead_snps.csv")))
    clump_rows[[analysis_name]] <- clumped
  }

  magma_prefix <- file.path(analysis_dir, paste0("trackB_", analysis_name))
  run_command(
    magma_exe,
    c(
      "--annotate",
      "--snp-loc", shQuote(paste0(bfile, ".bim")),
      "--gene-loc", shQuote(gene_loc),
      "--out", shQuote(magma_prefix)
    ),
    file.path(analysis_dir, paste0("trackB_", analysis_name, "_magma_annotate.log.txt"))
  )
  run_command(
    magma_exe,
    c(
      "--bfile", shQuote(bfile),
      "--pval", shQuote(magma_file), "use=1,2", "ncol=N",
      "--gene-annot", shQuote(paste0(magma_prefix, ".genes.annot")),
      "--out", shQuote(magma_prefix)
    ),
    file.path(analysis_dir, paste0("trackB_", analysis_name, "_magma_gene.log.txt"))
  )
  run_command(
    magma_exe,
    c(
      "--gene-results", shQuote(paste0(magma_prefix, ".genes.raw")),
      "--set-annot", shQuote(set_annot),
      "--out", shQuote(paste0(magma_prefix, "_geneset"))
    ),
    file.path(analysis_dir, paste0("trackB_", analysis_name, "_magma_geneset.log.txt"))
  )
  run_command(
    magma_exe,
    c(
      "--gene-results", shQuote(paste0(magma_prefix, ".genes.raw")),
      "--gene-covar", shQuote(gtex_covar),
      "--model", "direction=positive", "condition-hide=Average",
      "--out", shQuote(paste0(magma_prefix, "_tissue"))
    ),
    file.path(analysis_dir, paste0("trackB_", analysis_name, "_magma_tissue.log.txt"))
  )

  summary_rows[[analysis_name]] <- data.table(
    analysis = analysis_name,
    n_variants = nrow(dat),
    n_gws_variants = dat[P < 5e-8, .N],
    n_gws_mhc_variants = dat[P < 5e-8 & in_extended_mhc, .N],
    n_gws_non_mhc_variants = dat[P < 5e-8 & !in_extended_mhc, .N],
    median_effective_n = median_n,
    n_independent_lead_snps = nrow(clumped),
    n_independent_lead_snps_mhc = if (nrow(clumped) > 0) sum(clumped$in_extended_mhc) else 0,
    n_independent_lead_snps_non_mhc = if (nrow(clumped) > 0) sum(!clumped$in_extended_mhc) else 0
  )
}

summary_table <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
fwrite(summary_table, file.path(out_dir, "trackB_three_projection_locus_summary.csv"))
if (length(clump_rows) > 0) {
  fwrite(
    rbindlist(clump_rows, use.names = TRUE, fill = TRUE),
    file.path(out_dir, "trackB_three_projection_independent_lead_snps.csv")
  )
}

writeLines(
  c(
    "Track B three-projection locus and MAGMA analysis",
    "PLINK clumping: European 1000 Genomes reference, P1 5e-8, P2 1e-5, r2 0.1, 500 kb.",
    "Independent lead SNPs are PLINK clumps and must not be labelled FUMA genomic risk loci.",
    "Extended MHC sensitivity: chromosome 6, 25 to 34 Mb.",
    "FUMA inputs are provided separately for common, affective and functional projections.",
    "FUMA SNP2GENE must be run on the FUMA platform with GRCh37 and the same European reference population.",
    "MAGMA gene, gene-set and GTEx tissue analyses are run locally using the same reference files for all projections."
  ),
  file.path(out_dir, "trackB_three_projection_locus_magma_note.txt")
)

print(summary_table)
message("Completed. Outputs written to: ", out_dir)


