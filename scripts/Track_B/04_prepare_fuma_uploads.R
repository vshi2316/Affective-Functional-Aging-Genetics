options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir) || !dir.exists(base_dir)) stop("Set DCV_BASE_DIR.")

input_file <- file.path(
  base_dir, "analysis_ready_core", "trackB", "twofactor_gwas_qsnp",
  "trackB_twofactor_gwas_qsnp.csv"
)
fuma_root <- file.path(base_dir, "FUMA", "FUMA_twofactor_domain_GWAS_uploads")
dir.create(fuma_root, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing two-factor GWAS output: ", input_file)

dat <- fread(input_file, select = c(
  "SNP", "CHR", "BP", "A1", "A2", "MAF",
  "beta_affective", "se_affective", "p_affective",
  "beta_functional", "se_functional", "p_functional"
))

jobs <- list(
  FUMA_JOB06_Affective_Factor_GWAS = list(
    analysis = "affective_factor_GWAS",
    beta = "beta_affective", se = "se_affective", p = "p_affective",
    N = as.integer(Sys.getenv("TRACKB_AFFECTIVE_FUMA_N", unset = "500199")),
    filename = "FUMA_JOB06_affective_factor_GWAS_full_SNP2GENE_GRCh37.tsv.gz"
  ),
  FUMA_JOB07_Functional_Factor_GWAS = list(
    analysis = "functional_factor_GWAS",
    beta = "beta_functional", se = "se_functional", p = "p_functional",
    N = as.integer(Sys.getenv("TRACKB_FUNCTIONAL_FUMA_N", unset = "460844")),
    filename = "FUMA_JOB07_functional_factor_GWAS_full_SNP2GENE_GRCh37.tsv.gz"
  )
)

manifest <- list()
for (job_name in names(jobs)) {
  specification <- jobs[[job_name]]
  job_dir <- file.path(fuma_root, job_name)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  upload <- dat[, .(
    SNP, CHR, BP, A1, A2,
    P = get(specification$p),
    BETA = get(specification$beta),
    SE = get(specification$se),
    N = specification$N,
    MAF
  )]
  upload <- upload[
    !is.na(SNP) & nzchar(SNP) & is.finite(CHR) & is.finite(BP) &
      is.finite(P) & P > 0 & P <= 1 & is.finite(BETA) & is.finite(SE) & SE > 0 &
      is.finite(MAF) & MAF > 0 & MAF < 1
  ]
  output_file <- file.path(job_dir, specification$filename)
  fwrite(upload, output_file, sep = "\t", quote = FALSE, na = "")

  preview_file <- file.path(job_dir, sub("full_SNP2GENE_GRCh37.tsv.gz", "preview_first100.tsv", specification$filename))
  fwrite(upload[1:min(.N, 100)], preview_file, sep = "\t", quote = FALSE, na = "")

  writeLines(c(
    paste0("FUMA upload package: ", job_name),
    paste0("Mandatory upload file: ", specification$filename),
    "Genome build: GRCh37/hg19",
    "Reference population: 1000G Phase 3 EUR",
    "Column mapping: SNP=SNP; chromosome=CHR; position=BP; effect allele=A1; non-effect allele=A2; P=P; beta=BETA; standard error=SE; sample size=N; minor allele frequency=MAF.",
    "Maximum P-value of lead SNPs: 5e-8. FUMA accepts values no larger than 1e-5; do not use 0.05.",
    "Maximum P-value of candidate SNPs: 0.05.",
    "Independent significant SNP r2: 0.6.",
    "Lead SNP r2: 0.1.",
    "Merge distance: 250 kb.",
    "Use standard FUMA lead-SNP discovery. Leave predefined lead SNP file empty.",
    "Run positional mapping and MAGMA. Keep eQTL and chromatin interaction mapping disabled for this stage to avoid circular candidate selection."
  ), file.path(job_dir, "README_FUMA_UPLOAD_SETTINGS.txt"))

  manifest[[job_name]] <- data.table(
    job = job_name,
    analysis = specification$analysis,
    upload_file = output_file,
    rows = nrow(upload),
    gws_variants = upload[P < 5e-8, .N],
    min_p = min(upload$P),
    sample_size_N = specification$N,
    build = "GRCh37/hg19"
  )
}

fwrite(rbindlist(manifest), file.path(fuma_root, "FUMA_twofactor_upload_manifest.csv"))
message("FUMA upload packages written to: ", fuma_root)
