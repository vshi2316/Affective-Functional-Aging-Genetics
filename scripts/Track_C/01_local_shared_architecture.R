# ==============================================================================
# Track C local shared-architecture analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# Component: run_trackC1_lava_conjfdr.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(LAVA)
  library(cfdr.pleio)
})

project_dir <- normalizePath(Sys.getenv("DCV_PROJECT_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
data_root <- normalizePath(Sys.getenv("DCV_BASE_DIR", unset = project_dir), winslash = "/", mustWork = TRUE)
trait_file <- Sys.getenv("TRACKC1_TRAIT_MANIFEST", unset = file.path(project_dir, "config", "trackC1_trait_manifest.csv"))
pair_file <- Sys.getenv("TRACKC1_PAIR_MANIFEST", unset = file.path(project_dir, "config", "trackC1_pair_manifest.csv"))
ref_prefix <- Sys.getenv("LAVA_REFERENCE_PREFIX", unset = file.path(data_root, "reference", "lava", "lava-ukb"))
loci_file <- Sys.getenv("LAVA_LOCI_FILE", unset = file.path(data_root, "reference", "lava", "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
conj_ref_dir <- Sys.getenv("CONJFDR_REFERENCE_DIR", unset = file.path(data_root, "reference", "conjfdr"))
out_dir <- file.path(data_root, "analysis_ready_core", "trackC", "local_shared_architecture")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- fread(trait_file)
pairs <- fread(pair_file)
required_traits <- c("trait_id", "trait_name", "lava_file", "conjfdr_file", "cases", "controls")
required_pairs <- c("trait_1", "trait_2", "overlap_rho")
if (length(setdiff(required_traits, names(traits)))) stop("Track C1 trait manifest has missing columns.")
if (length(setdiff(required_pairs, names(pairs)))) stop("Track C1 pair manifest has missing columns.")
if (!file.exists(loci_file)) stop("Missing LAVA loci file: ", basename(loci_file))
if (!dir.exists(conj_ref_dir)) stop("Missing conjunctional FDR reference directory.")

resolve_input <- function(x) {
  if (file.exists(x)) return(normalizePath(x, winslash = "/"))
  normalizePath(file.path(data_root, "GWAS", "formatted", x), winslash = "/", mustWork = FALSE)
}
traits[, lava_file := vapply(lava_file, resolve_input, character(1))]
traits[, conjfdr_file := vapply(conjfdr_file, resolve_input, character(1))]
if (any(!file.exists(traits$lava_file)) || any(!file.exists(traits$conjfdr_file))) stop("One or more Track C1 GWAS inputs are missing.")

run_lava_pair <- function(pair, rows, pair_dir) {
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  info <- rows[, .(phenotype = trait_name, cases, controls, filename = lava_file)]
  info_file <- file.path(pair_dir, "lava_input_info.txt")
  fwrite(info, info_file, sep = "\t", na = "NA")

  overlap_file <- NULL
  if (is.finite(pair$overlap_rho) && pair$overlap_rho != 0) {
    overlap <- matrix(c(1, pair$overlap_rho, pair$overlap_rho, 1), 2, 2,
                      dimnames = list(rows$trait_name, rows$trait_name))
    overlap_file <- file.path(pair_dir, "lava_sample_overlap.txt")
    write.table(overlap, overlap_file, quote = FALSE, sep = "\t", col.names = NA)
  }

  input <- LAVA::process.input(
    input.info.file = info_file,
    sample.overlap.file = overlap_file,
    ref.prefix = ref_prefix,
    phenos = rows$trait_name
  )
  loci <- LAVA::read.loci(loci_file)
  univ <- vector("list", nrow(loci))
  bivar <- vector("list", nrow(loci))
  failures <- vector("list", nrow(loci))

  set.seed(154226)
  for (i in seq_len(nrow(loci))) {
    locus <- tryCatch(LAVA::process.locus(loci[i, ], input), error = function(e) e)
    if (inherits(locus, "error") || is.null(locus)) {
      failures[[i]] <- data.table(locus_index = i, message = if (inherits(locus, "error")) conditionMessage(locus) else "null_locus")
      next
    }
    result <- tryCatch(LAVA::run.univ.bivar(locus, univ.thresh = 0.05 / nrow(loci)), error = function(e) e)
    if (inherits(result, "error")) {
      failures[[i]] <- data.table(locus_index = i, message = conditionMessage(result))
      next
    }
    loc <- data.table(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n_snps = locus$n.snps, n_pcs = locus$K)
    if (!is.null(result$univ)) univ[[i]] <- cbind(loc, as.data.table(result$univ))
    if (!is.null(result$bivar)) bivar[[i]] <- cbind(loc, as.data.table(result$bivar))
  }
  fwrite(rbindlist(univ, fill = TRUE), file.path(pair_dir, "lava_univariate.csv"))
  fwrite(rbindlist(bivar, fill = TRUE), file.path(pair_dir, "lava_bivariate.csv"))
  fwrite(rbindlist(failures, fill = TRUE), file.path(pair_dir, "lava_skipped_loci.csv"))
}

read_cfdr <- function(path) {
  x <- fread(path)
  snp <- intersect(c("SNP", "rsid", "RSID"), names(x))[1]
  beta <- intersect(c("BETA", "beta", "b"), names(x))[1]
  pval <- intersect(c("PVAL", "P", "p", "pval"), names(x))[1]
  if (any(is.na(c(snp, beta, pval)))) stop("conjFDR input requires SNP, BETA and PVAL columns: ", basename(path))
  setnames(x, c(snp, beta, pval), c("SNP", "BETA", "PVAL"))
  unique(x[!is.na(SNP) & is.finite(BETA) & is.finite(PVAL) & PVAL > 0 & PVAL <= 1, .(SNP, BETA, PVAL)], by = "SNP")
}

run_conjfdr_pair <- function(pair, rows, pair_dir) {
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  trait1 <- read_cfdr(rows$conjfdr_file[1])
  trait2 <- read_cfdr(rows$conjfdr_file[2])
  refdat <- cfdr.pleio::refdata_location(conj_ref_dir)
  local_ref <- file.path(pair_dir, "local_reference")
  obj <- cfdr.pleio::cfdr_pleio$new()
  obj$init_data(
    trait1 = trait1,
    trait2 = trait2,
    trait_names = rows$trait_name,
    refdat = refdat,
    local_refdat_path = local_ref,
    correct_GC = TRUE,
    correct_SO = FALSE,
    filter_maf_min = 0.005,
    verbose = TRUE
  )
  obj$initialize_pruning_index(n_iter = 50, seed = 154226, verbose = TRUE)
  obj$calculate_cond_fdr(fdr_trait = 1, verbose = TRUE)
  obj$calculate_cond_fdr(fdr_trait = 2, verbose = TRUE)
  fwrite(obj$get_trait_results(), file.path(pair_dir, "conjfdr_results.csv"))
  saveRDS(obj, file.path(pair_dir, "conjfdr_analysis_object.rds"))
}

run_log <- list()
for (i in seq_len(nrow(pairs))) {
  pair <- pairs[i]
  rows <- traits[match(c(pair$trait_1, pair$trait_2), trait_id)]
  if (any(is.na(rows$trait_id))) stop("Unknown trait ID in Track C1 pair manifest.")
  pair_id <- paste(pair$trait_1, pair$trait_2, sprintf("rho%.2f", pair$overlap_rho), sep = "__")
  pair_dir <- file.path(out_dir, pair_id)
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  lava_status <- tryCatch({ run_lava_pair(pair, rows, file.path(pair_dir, "LAVA")); "completed" }, error = function(e) paste("failed", conditionMessage(e), sep = ": "))
  if (!identical(lava_status, "completed")) stop("LAVA failed for ", pair_id, ": ", lava_status)
  conj_status <- tryCatch({ run_conjfdr_pair(pair, rows, file.path(pair_dir, "conjFDR")); "completed" }, error = function(e) paste("failed", conditionMessage(e), sep = ": "))
  if (!identical(conj_status, "completed")) stop("conjFDR failed for ", pair_id, ": ", conj_status)
  run_log[[i]] <- data.table(pair_id, lava_status, conj_status)
}
fwrite(rbindlist(run_log), file.path(out_dir, "trackC1_run_log.csv"))
message("Track C1 completed: ", out_dir)

})


