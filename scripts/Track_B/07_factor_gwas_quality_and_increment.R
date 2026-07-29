options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

# Track B final, scope-locked quality and incremental-value audit.
#
# This script does not refit GenomicSEM, does not redefine the frozen
# 17/41/259 gene sets, and does not authorize rerunning Track C or Track D.
# The original factor GWAS remains primary. Q-clean results are sensitivity
# analyses for residual SNP heterogeneity only.

resolve_existing_dir <- function(path, label) {
  if (!nzchar(path) || !dir.exists(path)) stop(label, " directory is missing: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

resolve_existing_file <- function(path, label) {
  if (!nzchar(path) || !file.exists(path)) stop(label, " file is missing: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

resolve_executable <- function(path, label) {
  if (nzchar(path) && file.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  located <- Sys.which(path)
  if (nzchar(located)) return(unname(located))
  stop(label, " was not found. Set the corresponding environment variable.")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (default) "true" else "false")))
  value %in% c("1", "true", "yes", "y")
}

split_stages <- function(x) {
  unique(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_min <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

write_stage_marker <- function(stage_dir, stage, details = character()) {
  writeLines(
    c(
      paste0("stage=", stage),
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      details
    ),
    file.path(stage_dir, "COMPLETE.txt")
  )
}

require_stage <- function(stage_dir, stage) {
  marker <- file.path(stage_dir, "COMPLETE.txt")
  if (!file.exists(marker)) stop("Required stage is incomplete: ", stage, ". Missing: ", marker)
}

run_external <- function(executable, arguments, log_file) {
  output <- system2(executable, args = arguments, stdout = TRUE, stderr = TRUE)
  writeLines(output, log_file)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop("External command failed with status ", status, ". Inspect: ", log_file)
  }
  invisible(output)
}

first_existing_file <- function(paths, label = "output") {
  paths <- unique(paths[nzchar(paths)])
  hit <- paths[file.exists(paths) & file.info(paths)$size > 0]
  if (!length(hit)) {
    stop(label, " was not generated. Checked: ", paste(paths, collapse = "; "))
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

find_existing_file <- function(paths) {
  paths <- unique(paths[nzchar(paths)])
  hit <- paths[file.exists(paths) & file.info(paths)$size > 0]
  if (!length(hit)) return(NA_character_)
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

read_magma_gene_output <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(trimws(lines), "#")]
  if (!length(lines)) stop("MAGMA gene output is empty: ", path)
  fread(text = paste(lines, collapse = "\n"), header = TRUE)
}

read_fuma_loci <- function(path, factor_name) {
  x <- fread(path)
  required <- c("GenomicLocus", "rsID", "chr", "pos", "p", "start", "end")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("FUMA loci missing columns: ", paste(missing, collapse = ", "))
  x[, `:=`(
    factor = factor_name,
    GenomicLocus = as.integer(GenomicLocus),
    chr = as.integer(chr),
    pos = as.integer(pos),
    start = as.integer(start),
    end = as.integer(end),
    p = safe_num(p)
  )]
  x
}

read_fuma_leads <- function(path, factor_name) {
  x <- fread(path)
  required <- c("GenomicLocus", "rsID", "chr", "pos", "p")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("FUMA lead file missing columns: ", paste(missing, collapse = ", "))
  x[, `:=`(
    factor = factor_name,
    GenomicLocus = as.integer(GenomicLocus),
    chr = as.integer(chr),
    pos = as.integer(pos),
    p = safe_num(p)
  )]
  x
}

make_qq_plot <- function(p_values, title, output_file, max_points = 200000L) {
  p_values <- safe_num(p_values)
  p_values <- sort(p_values[is.finite(p_values) & p_values > 0 & p_values <= 1])
  n <- length(p_values)
  if (n < 100L) return(invisible(FALSE))
  index <- unique(as.integer(round(seq(1, n, length.out = min(n, max_points)))))
  plot_data <- data.table(
    expected = -log10(ppoints(n)[index]),
    observed = -log10(p_values[index])
  )
  plot <- ggplot(plot_data, aes(expected, observed)) +
    geom_point(size = 0.45, alpha = 0.55, colour = "#2C6E9E") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#B23A48") +
    labs(title = title, x = "Expected -log10(P)", y = "Observed -log10(P)") +
    theme_classic(base_size = 11)
  ggsave(output_file, plot, width = 5.8, height = 5.4, dpi = 600)
  invisible(TRUE)
}

two_sided_p <- function(beta, se) {
  beta <- safe_num(beta)
  se <- safe_num(se)
  z <- beta / se
  2 * pnorm(abs(z), lower.tail = FALSE)
}

base_dir <- resolve_existing_dir(
  Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd())),
  "DCV_BASE_DIR"
)
trackb_dir <- file.path(base_dir, "analysis_ready_core", "trackB")
if (!dir.exists(trackb_dir)) stop("Track B directory is missing: ", trackb_dir)

out_dir <- Sys.getenv(
  "TRACKB_FINAL_AUDIT_OUT_DIR",
  unset = file.path(trackb_dir, "final_qc_increment_audit")
)
out_dir <- ensure_dir(out_dir)

requested_stages <- split_stages(Sys.getenv(
  "TRACKB_FINAL_AUDIT_STAGES",
  unset = "prepare,ldsc,qsnp,magma,enabled,integrate"
))
valid_stages <- c("prepare", "ldsc", "qsnp", "magma", "enabled", "integrate")
unknown_stages <- setdiff(requested_stages, valid_stages)
if (length(unknown_stages)) stop("Unknown stage(s): ", paste(unknown_stages, collapse = ", "))

stage_dirs <- list(
  prepare = ensure_dir(file.path(out_dir, "00_prepare_and_specification")),
  ldsc = ensure_dir(file.path(out_dir, "01_factor_gwas_ldsc_qc")),
  qsnp = ensure_dir(file.path(out_dir, "02_qsnp_locus_sensitivity")),
  magma = ensure_dir(file.path(out_dir, "03_qclean_magma_sensitivity")),
  enabled = ensure_dir(file.path(out_dir, "04_factor_enabled_loci")),
  integrate = ensure_dir(file.path(out_dir, "05_integrated_closure"))
)

inputs <- list(
  twofactor = file.path(trackb_dir, "twofactor_gwas_qsnp", "trackB_twofactor_gwas_qsnp.csv"),
  five_trait_snp = file.path(trackb_dir, "factor_gwas_qsnp", "trackB_afffunc_sumstats.rds"),
  job06_dir = file.path(base_dir, "FUMA", "JOB06"),
  job07_dir = file.path(base_dir, "FUMA", "JOB07"),
  frozen_dir = file.path(trackb_dir, "trackC_twofactor_frozen_gene_sets", "tables"),
  fuma_upload_manifest = file.path(base_dir, "FUMA", "FUMA_twofactor_domain_GWAS_uploads", "FUMA_twofactor_upload_manifest.csv")
)

inputs$job06_loci <- file.path(inputs$job06_dir, "GenomicRiskLoci.txt")
inputs$job07_loci <- file.path(inputs$job07_dir, "GenomicRiskLoci.txt")
inputs$job06_leads <- file.path(inputs$job06_dir, "leadSNPs.txt")
inputs$job07_leads <- file.path(inputs$job07_dir, "leadSNPs.txt")
inputs$job06_magma <- file.path(inputs$job06_dir, "magma.genes.out")
inputs$job07_magma <- file.path(inputs$job07_dir, "magma.genes.out")
inputs$job06_params <- file.path(inputs$job06_dir, "params.config")
inputs$job07_params <- file.path(inputs$job07_dir, "params.config")
inputs$frozen_membership <- file.path(inputs$frozen_dir, "trackC_master_gene_membership_matrix.csv")
inputs$frozen_manifest <- file.path(inputs$frozen_dir, "trackC_frozen_gene_set_manifest.csv")

magma_exe <- resolve_executable(
  Sys.getenv("TRACKB_MAGMA_EXE", unset = Sys.getenv("MAGMA_EXECUTABLE", unset = "magma")),
  "MAGMA executable"
)

magma_bfile <- Sys.getenv(
  "TRACKB_MAGMA_BFILE",
  unset = file.path(base_dir, "reference", "magma", "g1000_eur", "g1000_eur")
)
for (extension in c(".bed", ".bim", ".fam")) {
  resolve_existing_file(paste0(magma_bfile, extension), paste0("MAGMA reference ", extension))
}
magma_bfile <- normalizePath(magma_bfile, winslash = "/", mustWork = FALSE)

ldsc_ld_dir <- resolve_existing_dir(Sys.getenv(
  "TRACKB_LDSC_LD_DIR",
  unset = file.path(base_dir, "reference", "ldsc", "eur_w_ld_chr")
), "LDSC LD-score")
ldsc_wld_dir <- resolve_existing_dir(Sys.getenv(
  "TRACKB_LDSC_WLD_DIR",
  unset = file.path(base_dir, "reference", "ldsc", "weights_hm3_no_hla")
), "LDSC regression-weight")
hm3_file <- resolve_existing_file(Sys.getenv(
  "TRACKB_LDSC_HM3_FILE",
  unset = file.path(base_dir, "reference", "ldsc", "w_hm3.snplist")
), "HapMap3 SNP list")

q_threshold <- safe_num(Sys.getenv("TRACKB_QSNP_THRESHOLD", unset = "5e-8"))
factor_threshold <- safe_num(Sys.getenv("TRACKB_FACTOR_GWS_THRESHOLD", unset = "5e-8"))
mhc_chr <- 6L
mhc_start <- 25000000L
mhc_end <- 34000000L

core_rds <- file.path(stage_dirs$prepare, "trackB_twofactor_compact_core.rds")
gene_loc_file <- file.path(stage_dirs$prepare, "fuma_18516_gene_locations.txt")
spec_file <- file.path(stage_dirs$prepare, "final_trackB_audit_specification.json")

message("Track B final quality and incremental-value audit")
message("Output root: ", out_dir)
message("Requested stages: ", paste(requested_stages, collapse = ", "))

if ("prepare" %in% requested_stages) {
  required_inputs <- unlist(inputs, use.names = FALSE)
  missing_inputs <- required_inputs[!file.exists(required_inputs) & !dir.exists(required_inputs)]
  if (length(missing_inputs)) stop("Missing required inputs: ", paste(missing_inputs, collapse = "; "))

  params06 <- readLines(inputs$job06_params, warn = FALSE)
  params07 <- readLines(inputs$job07_params, warn = FALSE)
  required_param_patterns <- c(
    "FUMA = v1.8.2", "snp2genegrch38 = 0",
    "leadP = 5e-8", "r2 = 0.6", "r2_2 = 0.1", "pop = EUR",
    "mergeDist = 250", "ensembl = v102", "genetype = protein_coding",
    "magma_window = 0"
  )
  for (pattern in required_param_patterns) {
    if (!any(grepl(pattern, params06, fixed = TRUE)) || !any(grepl(pattern, params07, fixed = TRUE))) {
      stop("FUMA parameter mismatch for frozen pattern: ", pattern)
    }
  }

  gene06 <- fread(inputs$job06_magma, select = c("GENE", "CHR", "START", "STOP", "SYMBOL"))
  gene07 <- fread(inputs$job07_magma, select = c("GENE", "CHR", "START", "STOP", "SYMBOL"))
  if (nrow(gene06) != 18516L || nrow(gene07) != 18516L) {
    stop("Expected 18,516 tested FUMA MAGMA genes in both JOB06 and JOB07.")
  }
  setorder(gene06, GENE)
  setorder(gene07, GENE)
  if (!identical(gene06[, .(GENE, CHR, START, STOP)], gene07[, .(GENE, CHR, START, STOP)])) {
    stop("JOB06 and JOB07 do not share identical FUMA MAGMA gene coordinates.")
  }
  fwrite(gene06[, .(GENE, CHR, START, STOP)], gene_loc_file, sep = "\t", col.names = FALSE)
  fwrite(gene06, file.path(stage_dirs$prepare, "fuma_18516_gene_universe.csv"))

  selected_columns <- c(
    "CHR", "BP", "SNP", "A1", "A2", "MAF",
    "beta_affective", "se_affective", "p_affective",
    "beta_functional", "se_functional", "p_functional",
    "q_snp", "q_snp_df", "q_snp_p"
  )
  core <- fread(inputs$twofactor, select = selected_columns, showProgress = TRUE)
  core[, CHR := as.integer(CHR)]
  core[, BP := as.integer(BP)]
  numeric_columns <- setdiff(selected_columns, c("SNP", "A1", "A2", "CHR", "BP"))
  core[, (numeric_columns) := lapply(.SD, safe_num), .SDcols = numeric_columns]
  if (nrow(core) != 6735089L) warning("Two-factor row count differs from the prior audit: ", nrow(core))
  saveRDS(core, core_rds, compress = FALSE)

  manifest <- rbindlist(lapply(names(inputs), function(name) {
    path <- inputs[[name]]
    info <- file.info(path)
    data.table(
      role = name,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      is_directory = isTRUE(info$isdir),
      bytes = if (isTRUE(info$isdir)) NA_real_ else info$size,
      file_mtime = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
      sha256 = if (isTRUE(info$isdir) || info$size > 100000000) NA_character_ else sha256_file(path)
    )
  }), fill = TRUE)
  fwrite(manifest, file.path(stage_dirs$prepare, "input_manifest.csv"))

  specification <- list(
    audit_name = "Track B final quality and incremental-value audit",
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    primary_factor_gwas_unchanged = TRUE,
    genomicsem_refit_authorized = FALSE,
    q_snp_threshold = q_threshold,
    factor_gws_threshold = factor_threshold,
    qclean_role = "sensitivity analysis only",
    frozen_gene_sets_changed = FALSE,
    downstream_trackC_trackD_rerun_authorized = FALSE,
    factor_enabled_definition = "FUMA factor interval in which none of the contributing univariate GWAS reaches P < 5e-8 within the same interval",
    mhc_rule = "GRCh37 chr6:25,000,000-34,000,000",
    magma_software_rule = "Use the configured MAGMA executable with the frozen gene universe and analysis settings",
    qclean_gene_universe = 18516L,
    qclean_gene_coordinates = "reconstructed exactly from JOB06/JOB07 FUMA final outputs",
    effective_sample_size_claim_authorized = FALSE,
    ldsc_scaling_note = "LDSC h2 is reported under explicit computational N scalings and is not a traditional cohort sample size"
  )
  write_json(specification, spec_file, pretty = TRUE, auto_unbox = TRUE, na = "null")
  write_stage_marker(stage_dirs$prepare, "prepare", c(
    paste0("twofactor_variants=", nrow(core)),
    "fuma_gene_universe=18516",
    "frozen_sets_changed=false"
  ))
  rm(core)
  gc()
  message("Completed prepare stage: ", stage_dirs$prepare)
}

if ("ldsc" %in% requested_stages) {
  require_stage(stage_dirs$prepare, "prepare")
  if (!requireNamespace("ldscr", quietly = TRUE)) stop("The ldscr R package is required for the ldsc stage.")
  core <- readRDS(core_rds)
  hm3 <- fread(hm3_file, select = "SNP")
  hm3_snps <- unique(hm3$SNP)

  fuma_upload <- fread(inputs$fuma_upload_manifest)
  scaling_reference <- data.table(
    factor = c("affective", "functional"),
    fuma_upload_N = c(
      fuma_upload[analysis == "affective_factor_GWAS", sample_size_N][1],
      fuma_upload[analysis == "functional_factor_GWAS", sample_size_N][1]
    )
  )

  loci <- rbindlist(list(
    read_fuma_loci(inputs$job06_loci, "affective"),
    read_fuma_loci(inputs$job07_loci, "functional")
  ), fill = TRUE)
  top_loci <- loci[order(p), .SD[1], by = factor]

  n_audit <- list()
  ldsc_rows <- list()
  row_id <- 0L

  for (factor_name in c("affective", "functional")) {
    beta_col <- paste0("beta_", factor_name)
    se_col <- paste0("se_", factor_name)
    p_col <- paste0("p_", factor_name)
    z <- core[[beta_col]] / core[[se_col]]
    snpwise_n <- 1 / (2 * core$MAF * (1 - core$MAF) * core[[se_col]]^2)
    n_subset <- snpwise_n[
      is.finite(snpwise_n) & is.finite(core$MAF) & core$MAF >= 0.1 & core$MAF <= 0.4
    ]
    formula_mean_n <- mean(n_subset)
    formula_median_n <- median(n_subset)
    fuma_n <- scaling_reference[factor == factor_name, fuma_upload_N][1]
    n_audit[[factor_name]] <- data.table(
      factor = factor_name,
      n_snps_maf_0p1_0p4 = length(n_subset),
      formula_mean_N = formula_mean_n,
      formula_median_N = formula_median_n,
      formula_q05_N = as.numeric(quantile(n_subset, 0.05)),
      formula_q95_N = as.numeric(quantile(n_subset, 0.95)),
      fuma_upload_scaling_N = fuma_n,
      traditional_cohort_N_authorized = FALSE
    )

    top <- top_loci[factor == factor_name][1]
    scenarios <- list(
      full = rep(TRUE, nrow(core)),
      non_mhc = !(core$CHR == mhc_chr & core$BP >= mhc_start & core$BP <= mhc_end),
      top_fuma_locus_excluded = !(core$CHR == top$chr & core$BP >= top$start & core$BP <= top$end)
    )
    scaling_methods <- list(
      formula_mean_maf_0p1_0p4 = formula_mean_n,
      fuma_upload_scaling = fuma_n
    )

    for (scenario_name in names(scenarios)) {
      for (scaling_name in names(scaling_methods)) {
        if (scenario_name != "full" && scaling_name != "formula_mean_maf_0p1_0p4") next
        keep <- scenarios[[scenario_name]] & is.finite(z) & !is.na(core$SNP) & core$SNP %chin% hm3_snps
        n_scale <- scaling_methods[[scaling_name]]
        munged <- data.table(
          SNP = core$SNP[keep],
          N = rep(n_scale, sum(keep)),
          Z = z[keep],
          A1 = core$A1[keep],
          A2 = core$A2[keep]
        )
        chisq_max <- max(0.001 * n_scale, 80)
        n_after_chisq <- sum(munged$Z^2 <= chisq_max, na.rm = TRUE)
        row_id <- row_id + 1L
        result <- tryCatch(
          ldscr::ldsc_h2(
            munged_sumstats = munged,
            ld = ldsc_ld_dir,
            wld = ldsc_wld_dir,
            n_blocks = 200,
            chisq_max = NA
          ),
          error = function(e) e
        )
        if (inherits(result, "error")) {
          ldsc_rows[[row_id]] <- data.table(
            factor = factor_name, scenario = scenario_name, N_scaling_method = scaling_name,
            N_scaling_value = n_scale, n_hm3_input = nrow(munged),
            chisq_max = chisq_max, n_after_chisq_filter_expected = n_after_chisq,
            status = "failed", error_message = conditionMessage(result)
          )
        } else {
          result <- as.data.table(result)
          result[, `:=`(
            factor = factor_name,
            scenario = scenario_name,
            N_scaling_method = scaling_name,
            N_scaling_value = n_scale,
            n_hm3_input = nrow(munged),
            chisq_max = chisq_max,
            n_after_chisq_filter_expected = n_after_chisq,
            ratio_interpretable = is.finite(mean_chisq) & mean_chisq > 1.02,
            status = "completed",
            error_message = NA_character_
          )]
          ldsc_rows[[row_id]] <- result
        }
        rm(munged, result)
        gc()
      }
    }

    make_qq_plot(
      core[[p_col]],
      paste0("", tools::toTitleCase(factor_name), " factor GWAS"),
      file.path(stage_dirs$ldsc, paste0("qq_", factor_name, "_factor_gwas.png"))
    )
  }

  n_audit_dt <- rbindlist(n_audit, fill = TRUE)
  ldsc_dt <- rbindlist(ldsc_rows, fill = TRUE)
  fwrite(n_audit_dt, file.path(stage_dirs$ldsc, "factor_gwas_N_scaling_audit.csv"))
  fwrite(ldsc_dt, file.path(stage_dirs$ldsc, "factor_gwas_ldsc_results.csv"))
  fwrite(top_loci, file.path(stage_dirs$ldsc, "top_fuma_loci_used_for_sensitivity.csv"))

  if (any(ldsc_dt$status != "completed")) stop("One or more LDSC models failed; inspect factor_gwas_ldsc_results.csv")
  write_stage_marker(stage_dirs$ldsc, "ldsc", c(
    paste0("models_completed=", nrow(ldsc_dt)),
    "traditional_effective_sample_size_claim_authorized=false",
    "primary_N_scaling=formula_mean_maf_0p1_0p4",
    "sensitivity_N_scaling=fuma_upload_scaling"
  ))
  rm(core)
  gc()
  message("Completed LDSC stage: ", stage_dirs$ldsc)
}

if ("qsnp" %in% requested_stages) {
  require_stage(stage_dirs$prepare, "prepare")
  core <- readRDS(core_rds)
  core[, q_gws := is.finite(q_snp_p) & q_snp_p < q_threshold]
  core[, affective_gws := is.finite(p_affective) & p_affective < factor_threshold]
  core[, functional_gws := is.finite(p_functional) & p_functional < factor_threshold]

  direct_summary <- rbindlist(lapply(c("affective", "functional"), function(factor_name) {
    factor_gws <- core[[paste0(factor_name, "_gws")]]
    data.table(
      factor = factor_name,
      n_all_variants = nrow(core),
      n_factor_gws_snps = sum(factor_gws, na.rm = TRUE),
      n_qsnp_gws_snps = sum(core$q_gws, na.rm = TRUE),
      n_direct_overlap_snps = sum(factor_gws & core$q_gws, na.rm = TRUE),
      proportion_factor_gws_snps_qsnp = sum(factor_gws & core$q_gws, na.rm = TRUE) / sum(factor_gws, na.rm = TRUE),
      interpretation_unit = "correlated SNP count; not an independent locus proportion"
    )
  }))
  fwrite(direct_summary, file.path(stage_dirs$qsnp, "factor_qsnp_direct_overlap_summary.csv"))

  leads <- rbindlist(list(
    read_fuma_leads(inputs$job06_leads, "affective"),
    read_fuma_leads(inputs$job07_leads, "functional")
  ), fill = TRUE)
  lead_core <- core[, .(SNP, q_snp_p, q_gws, p_affective, p_functional)]
  lead_audit <- merge(leads, lead_core, by.x = "rsID", by.y = "SNP", all.x = TRUE)
  lead_audit[, original_lead_exactly_qclean := !is.na(q_gws) & !q_gws]
  fwrite(lead_audit, file.path(stage_dirs$qsnp, "fuma_lead_snp_qsnp_audit.csv"))
  lead_summary <- lead_audit[, .(
    n_original_fuma_leads = .N,
    n_original_fuma_leads_qclean = sum(original_lead_exactly_qclean, na.rm = TRUE),
    any_original_fuma_lead_qclean = any(original_lead_exactly_qclean, na.rm = TRUE),
    all_original_fuma_leads_qclean = all(original_lead_exactly_qclean)
  ), by = .(factor, GenomicLocus)]

  loci <- rbindlist(list(
    read_fuma_loci(inputs$job06_loci, "affective"),
    read_fuma_loci(inputs$job07_loci, "functional")
  ), fill = TRUE)
  setkey(core, CHR, BP)
  locus_rows <- vector("list", nrow(loci))
  for (i in seq_len(nrow(loci))) {
    locus <- loci[i]
    region <- core[CHR == locus$chr & BP >= locus$start & BP <= locus$end]
    p_col <- paste0("p_", locus$factor)
    factor_sig <- is.finite(region[[p_col]]) & region[[p_col]] < factor_threshold
    qclean_factor_sig <- factor_sig & !region$q_gws
    clean_region <- region[qclean_factor_sig]
    top_clean <- if (nrow(clean_region)) clean_region[which.min(get(p_col))] else NULL
    original_top <- region[SNP == locus$rsID]
    locus_rows[[i]] <- data.table(
      factor = locus$factor,
      GenomicLocus = locus$GenomicLocus,
      chr = locus$chr,
      start = locus$start,
      end = locus$end,
      original_top_snp = locus$rsID,
      original_top_p = locus$p,
      original_top_qsnp_p = if (nrow(original_top)) safe_num(original_top$q_snp_p[1]) else NA_real_,
      original_top_exactly_qclean = if (nrow(original_top)) !isTRUE(original_top$q_gws[1]) else NA,
      n_variants_in_interval = nrow(region),
      n_qsnp_gws_in_interval = sum(region$q_gws, na.rm = TRUE),
      contains_qsnp_gws = any(region$q_gws, na.rm = TRUE),
      n_factor_gws_in_interval = sum(factor_sig, na.rm = TRUE),
      n_factor_gws_qsnp_overlap = sum(factor_sig & region$q_gws, na.rm = TRUE),
      n_qclean_factor_gws_in_interval = sum(qclean_factor_sig, na.rm = TRUE),
      qclean_interval_retained = any(qclean_factor_sig, na.rm = TRUE),
      qclean_top_snp = if (nrow(clean_region)) as.character(top_clean$SNP) else NA_character_,
      qclean_top_p = if (nrow(clean_region)) safe_num(top_clean[[p_col]]) else NA_real_
    )
  }
  locus_audit <- rbindlist(locus_rows, fill = TRUE)
  locus_audit <- merge(
    locus_audit,
    lead_summary,
    by = c("factor", "GenomicLocus"), all.x = TRUE
  )
  locus_audit[, locus_retention_class := fifelse(
    !qclean_interval_retained, "not_retained_after_qclean",
    fifelse(original_top_exactly_qclean, "retained_with_original_top_snp", "retained_with_alternative_qclean_signal")
  )]
  fwrite(locus_audit, file.path(stage_dirs$qsnp, "fuma_locus_qsnp_qclean_audit.csv"))

  fuma_upload <- fread(inputs$fuma_upload_manifest)
  qclean_manifest <- list()
  for (factor_name in c("affective", "functional")) {
    p_col <- paste0("p_", factor_name)
    fuma_n <- if (factor_name == "affective") {
      fuma_upload[analysis == "affective_factor_GWAS", sample_size_N][1]
    } else {
      fuma_upload[analysis == "functional_factor_GWAS", sample_size_N][1]
    }
    qclean <- core[!q_gws & is.finite(get(p_col)) & get(p_col) > 0 & get(p_col) <= 1,
      .(SNP, P = get(p_col), N = as.integer(fuma_n))
    ]
    qclean_file <- file.path(stage_dirs$qsnp, paste0(factor_name, "_qclean_magma_input.txt"))
    fwrite(qclean, qclean_file, sep = "\t")
    qclean_manifest[[factor_name]] <- data.table(
      factor = factor_name,
      file = normalizePath(qclean_file, winslash = "/", mustWork = TRUE),
      n_variants = nrow(qclean),
      n_removed_qsnp_gws = sum(core$q_gws, na.rm = TRUE),
      n_qclean_factor_gws = qclean[P < factor_threshold, .N],
      magma_scaling_N_inherited_from_fuma_upload = fuma_n
    )
  }
  fwrite(rbindlist(qclean_manifest), file.path(stage_dirs$qsnp, "qclean_magma_input_manifest.csv"))

  write_stage_marker(stage_dirs$qsnp, "qsnp", c(
    paste0("fuma_loci_audited=", nrow(locus_audit)),
    "original_factor_gwas_primary=true",
    "qclean_role=sensitivity_only"
  ))
  rm(core)
  gc()
  message("Completed Q_SNP stage: ", stage_dirs$qsnp)
}

if ("magma" %in% requested_stages) {
  require_stage(stage_dirs$prepare, "prepare")
  require_stage(stage_dirs$qsnp, "qsnp")

  qclean_manifest <- fread(file.path(stage_dirs$qsnp, "qclean_magma_input_manifest.csv"))
  annotation_prefix <- file.path(stage_dirs$magma, "fuma_exact_coordinates")
  annotation_candidates <- c(
    paste0(annotation_prefix, ".genes.annot"),
    paste0(annotation_prefix, ".genes.annot.txt")
  )
  annotation_file <- find_existing_file(annotation_candidates)
  if (is.na(annotation_file)) {
    run_external(
      magma_exe,
      c(
        "--annotate",
        "--snp-loc", shQuote(paste0(magma_bfile, ".bim")),
        "--gene-loc", shQuote(gene_loc_file),
        "--out", shQuote(annotation_prefix)
      ),
      file.path(stage_dirs$magma, "magma_annotation.log")
    )
    annotation_file <- first_existing_file(
      annotation_candidates,
      "MAGMA gene annotation output"
    )
  } else {
    message("Reusing completed MAGMA annotation: ", annotation_file)
  }

  gene_results <- list()
  for (factor_name in c("affective", "functional")) {
    qclean_file <- qclean_manifest[factor == factor_name, file][1]
    output_prefix <- file.path(stage_dirs$magma, paste0(factor_name, "_qclean"))
    result_candidates <- c(
      paste0(output_prefix, ".genes.out"),
      paste0(output_prefix, ".genes.out.txt")
    )
    result_file <- find_existing_file(result_candidates)
    if (is.na(result_file)) {
      run_external(
        magma_exe,
        c(
          "--bfile", shQuote(magma_bfile),
          "--pval", shQuote(qclean_file), "use=1,2", "ncol=N",
          "--gene-annot", shQuote(annotation_file),
          "--out", shQuote(output_prefix)
        ),
        file.path(stage_dirs$magma, paste0("magma_", factor_name, "_qclean.log"))
      )
      result_file <- first_existing_file(
        result_candidates,
        paste0("MAGMA ", factor_name, " Q-clean gene output")
      )
    } else {
      message("Reusing completed ", factor_name, " Q-clean MAGMA result: ", result_file)
    }
    result <- read_magma_gene_output(result_file)
    result[, factor := factor_name]
    gene_results[[factor_name]] <- result
    fwrite(result, file.path(stage_dirs$magma, paste0(factor_name, "_qclean_magma_genes.csv")))
  }

  original_gene <- rbindlist(list(
    fread(inputs$job06_magma)[, factor := "affective"],
    fread(inputs$job07_magma)[, factor := "functional"]
  ), fill = TRUE)
  original_gene[, P := safe_num(P)]
  original_gene[, SYMBOL := toupper(as.character(SYMBOL))]
  threshold_fixed <- 0.05 / 18516

  comparison_rows <- list()
  qclean_symbol_rows <- list()
  for (factor_name in c("affective", "functional")) {
    q <- gene_results[[factor_name]]
    if (!"GENE" %in% names(q) || !"P" %in% names(q)) stop("Unexpected MAGMA output columns for ", factor_name)
    q[, P := safe_num(P)]
    original <- original_gene[factor == factor_name, .(GENE, SYMBOL, original_P = P)]
    comparison <- merge(original, q[, .(GENE, qclean_P = P)], by = "GENE", all.x = TRUE)
    comparison[, `:=`(
      factor = factor_name,
      original_significant = is.finite(original_P) & original_P < threshold_fixed,
      qclean_significant = is.finite(qclean_P) & qclean_P < threshold_fixed,
      fixed_bonferroni_threshold = threshold_fixed
    )]
    comparison_rows[[factor_name]] <- comparison
    qclean_symbol_rows[[factor_name]] <- comparison[, .(
      qclean_P = safe_min(qclean_P),
      qclean_significant = any(qclean_significant, na.rm = TRUE)
    ), by = .(factor, SYMBOL)]
  }
  comparison_all <- rbindlist(comparison_rows, fill = TRUE)
  qclean_symbol <- rbindlist(qclean_symbol_rows, fill = TRUE)
  fwrite(comparison_all, file.path(stage_dirs$magma, "original_vs_qclean_magma_gene_comparison.csv"))

  frozen <- fread(inputs$frozen_membership)
  frozen[, symbol_upper := toupper(symbol)]
  aff_q <- qclean_symbol[factor == "affective", .(
    symbol_upper = SYMBOL,
    qclean_p_affective = qclean_P,
    qclean_sig_affective = qclean_significant
  )]
  fun_q <- qclean_symbol[factor == "functional", .(
    symbol_upper = SYMBOL,
    qclean_p_functional = qclean_P,
    qclean_sig_functional = qclean_significant
  )]
  frozen_q <- merge(frozen, aff_q, by = "symbol_upper", all.x = TRUE)
  frozen_q <- merge(frozen_q, fun_q, by = "symbol_upper", all.x = TRUE)
  frozen_q[, qclean_original_driver_retained := fifelse(
    primary_set == "cross_model_core", qclean_sig_affective & qclean_sig_functional,
    fifelse(primary_set == "affective_prioritized", qclean_sig_affective,
      fifelse(primary_set == "functional_prioritized", qclean_sig_functional, NA)
    )
  )]
  frozen_q[, qclean_original_pattern_retained := fifelse(
    primary_set == "cross_model_core", qclean_sig_affective & qclean_sig_functional,
    fifelse(primary_set == "affective_prioritized", qclean_sig_affective & !qclean_sig_functional,
      fifelse(primary_set == "functional_prioritized", qclean_sig_functional & !qclean_sig_affective, NA)
    )
  )]
  frozen_q[, interpretation := "Sensitivity status only; frozen membership is unchanged"]
  fwrite(frozen_q, file.path(stage_dirs$magma, "frozen_17_41_259_qclean_retention.csv"))
  frozen_summary <- frozen_q[primary_set %chin% c("cross_model_core", "affective_prioritized", "functional_prioritized"), .(
    frozen_n = .N,
    qclean_driver_retained_n = sum(qclean_original_driver_retained, na.rm = TRUE),
    qclean_pattern_retained_n = sum(qclean_original_pattern_retained, na.rm = TRUE),
    qclean_driver_retained_proportion = mean(qclean_original_driver_retained, na.rm = TRUE),
    qclean_pattern_retained_proportion = mean(qclean_original_pattern_retained, na.rm = TRUE)
  ), by = primary_set]
  fwrite(frozen_summary, file.path(stage_dirs$magma, "frozen_set_qclean_retention_summary.csv"))

  qclean_tested_counts <- vapply(gene_results, nrow, integer(1))
  parameter_audit <- list(
    annotation_gene_universe_n = 18516L,
    qclean_tested_gene_counts = as.list(qclean_tested_counts),
    same_frozen_annotation_universe_used = TRUE,
    exact_gene_coordinates_source = "FUMA JOB06/JOB07 Ensembl final outputs",
    gene_window = 0,
    reference_panel = magma_bfile,
    interpretation = "Q-clean local MAGMA sensitivity; not a replacement for primary FUMA JOB06/JOB07",
    frozen_sets_changed = FALSE
  )
  write_json(parameter_audit, file.path(stage_dirs$magma, "qclean_magma_parameter_audit.json"), pretty = TRUE, auto_unbox = TRUE)

  if (any(qclean_tested_counts < 0.95 * 18516L)) {
    stop("Q-clean MAGMA tested fewer than 95% of the frozen 18,516-gene annotation universe.")
  }
  write_stage_marker(stage_dirs$magma, "magma", c(
    "same_frozen_annotation_universe_used=true",
    "frozen_sets_changed=false"
  ))
  message("Completed Q-clean MAGMA stage: ", stage_dirs$magma)
}

if ("enabled" %in% requested_stages) {
  require_stage(stage_dirs$prepare, "prepare")
  require_stage(stage_dirs$qsnp, "qsnp")
  raw <- as.data.table(readRDS(inputs$five_trait_snp))
  raw[, CHR := as.integer(CHR)]
  raw[, BP := as.integer(BP)]
  setkey(raw, CHR, BP)

  loci <- rbindlist(list(
    read_fuma_loci(inputs$job06_loci, "affective"),
    read_fuma_loci(inputs$job07_loci, "functional")
  ), fill = TRUE)
  locus_qclean <- fread(file.path(stage_dirs$qsnp, "fuma_locus_qsnp_qclean_audit.csv"))
  enabled_rows <- vector("list", nrow(loci))

  factor_traits <- list(
    affective = c("mdd", "lonely"),
    functional = c("frailty", "srh", "walkpace")
  )
  for (i in seq_len(nrow(loci))) {
    locus <- loci[i]
    region <- raw[CHR == locus$chr & BP >= locus$start & BP <= locus$end]
    traits <- factor_traits[[locus$factor]]
    trait_min <- setNames(rep(NA_real_, length(traits)), traits)
    for (trait in traits) {
      p <- two_sided_p(region[[paste0("beta.", trait)]], region[[paste0("se.", trait)]])
      trait_min[[trait]] <- safe_min(p)
    }
    enabled_rows[[i]] <- data.table(
      factor = locus$factor,
      GenomicLocus = locus$GenomicLocus,
      chr = locus$chr,
      start = locus$start,
      end = locus$end,
      factor_top_snp = locus$rsID,
      factor_top_p = locus$p,
      n_harmonized_input_snps_in_interval = nrow(region),
      min_p_mdd = if ("mdd" %in% traits) trait_min[["mdd"]] else NA_real_,
      min_p_lonely = if ("lonely" %in% traits) trait_min[["lonely"]] else NA_real_,
      min_p_frailty = if ("frailty" %in% traits) trait_min[["frailty"]] else NA_real_,
      min_p_srh = if ("srh" %in% traits) trait_min[["srh"]] else NA_real_,
      min_p_walkpace = if ("walkpace" %in% traits) trait_min[["walkpace"]] else NA_real_,
      minimum_contributing_trait_p = safe_min(trait_min),
      any_contributing_trait_gws = any(is.finite(trait_min) & trait_min < factor_threshold),
      factor_enabled_locus = !any(is.finite(trait_min) & trait_min < factor_threshold),
      terminology = "factor-enabled locus; not a novel locus claim"
    )
  }
  enabled <- rbindlist(enabled_rows, fill = TRUE)
  enabled <- merge(
    enabled,
    locus_qclean[, .(factor, GenomicLocus, qclean_interval_retained, locus_retention_class)],
    by = c("factor", "GenomicLocus"), all.x = TRUE
  )
  fwrite(enabled, file.path(stage_dirs$enabled, "factor_enabled_fuma_loci.csv"))
  enabled_summary <- enabled[, .(
    n_fuma_loci = .N,
    n_factor_enabled_loci = sum(factor_enabled_locus, na.rm = TRUE),
    n_factor_enabled_qclean_retained = sum(factor_enabled_locus & qclean_interval_retained, na.rm = TRUE),
    proportion_factor_enabled = mean(factor_enabled_locus, na.rm = TRUE)
  ), by = factor]
  fwrite(enabled_summary, file.path(stage_dirs$enabled, "factor_enabled_loci_summary.csv"))

  write_stage_marker(stage_dirs$enabled, "enabled", c(
    paste0("fuma_loci_evaluated=", nrow(enabled)),
    "novel_locus_claim_authorized=false"
  ))
  rm(raw)
  gc()
  message("Completed factor-enabled loci stage: ", stage_dirs$enabled)
}

if ("integrate" %in% requested_stages) {
  for (stage in c("prepare", "ldsc", "qsnp", "magma", "enabled")) {
    require_stage(stage_dirs[[stage]], stage)
  }
  ldsc_results <- fread(file.path(stage_dirs$ldsc, "factor_gwas_ldsc_results.csv"))
  direct <- fread(file.path(stage_dirs$qsnp, "factor_qsnp_direct_overlap_summary.csv"))
  locus <- fread(file.path(stage_dirs$qsnp, "fuma_locus_qsnp_qclean_audit.csv"))
  frozen_summary <- fread(file.path(stage_dirs$magma, "frozen_set_qclean_retention_summary.csv"))
  enabled_summary <- fread(file.path(stage_dirs$enabled, "factor_enabled_loci_summary.csv"))
  magma_audit <- read_json(file.path(stage_dirs$magma, "qclean_magma_parameter_audit.json"), simplifyVector = TRUE)

  integrated <- list(
    audit_name = "Track B final quality and incremental-value audit",
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    statistical_scope_closed = TRUE,
    genomicsem_refitted = FALSE,
    primary_factor_gwas_changed = FALSE,
    qclean_is_sensitivity_only = TRUE,
    frozen_17_41_259_changed = FALSE,
    trackC_trackD_rerun_authorized = FALSE,
    factor_effective_sample_size_claim_authorized = FALSE,
    ldsc_models_planned = nrow(ldsc_results),
    ldsc_models_completed = sum(ldsc_results$status == "completed"),
    direct_snp_overlap = direct,
    fuma_locus_counts = locus[, .(
      n_loci = .N,
      n_containing_qsnp_gws = sum(contains_qsnp_gws, na.rm = TRUE),
      n_qclean_retained = sum(qclean_interval_retained, na.rm = TRUE),
      n_not_retained = sum(!qclean_interval_retained, na.rm = TRUE)
    ), by = factor],
    frozen_set_qclean_retention = frozen_summary,
    factor_enabled_loci = enabled_summary,
    qclean_magma_same_frozen_annotation_universe = magma_audit$same_frozen_annotation_universe_used,
    authorized_terms = c(
      "Q_SNP sensitivity analysis",
      "factor-enabled locus",
      "locus not genome-wide significant in contributing univariate GWAS"
    ),
    prohibited_terms = c(
      "novel locus without GWAS Catalog audit",
      "independent factor-GWAS replication",
      "Q_SNP proves the factor association is invalid",
      "Q-clean re-frozen gene set",
      "traditional effective sample size for the latent factor"
    ),
    stopping_rule = "No further statistical expansion after this audit, irrespective of positive or negative sensitivity results"
  )
  write_json(integrated, file.path(stage_dirs$integrate, "final_trackB_audit_closure.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

  overview <- rbindlist(list(
    data.table(section = "ldsc", metric = "models_completed", value = as.character(sum(ldsc_results$status == "completed"))),
    direct[, .(section = "qsnp_snp", metric = paste0(factor, "_direct_overlap"), value = as.character(n_direct_overlap_snps))],
    locus[, .(section = "qsnp_locus", metric = paste0(factor[1], "_qclean_retained_loci"), value = as.character(sum(qclean_interval_retained))), by = factor][, factor := NULL],
    frozen_summary[, .(section = "magma", metric = paste0(primary_set, "_driver_retained"), value = as.character(qclean_driver_retained_n))],
    enabled_summary[, .(section = "factor_enabled", metric = paste0(factor, "_qclean_retained"), value = as.character(n_factor_enabled_qclean_retained))]
  ), fill = TRUE)
  fwrite(overview, file.path(stage_dirs$integrate, "final_trackB_audit_overview.csv"))
  write_stage_marker(stage_dirs$integrate, "integrate", c(
    "statistical_scope_closed=true",
    "primary_factor_gwas_changed=false",
    "frozen_17_41_259_changed=false",
    "trackC_trackD_rerun_authorized=false"
  ))
  writeLines(
    c(
      "Track B final audit completed.",
      "The original factor GWAS remains primary.",
      "Q-clean analyses are sensitivity analyses only.",
      "The frozen 17/41/259 gene sets were not changed.",
      "No Track C or Track D rerun is authorized by this audit.",
      "All statistical expansion stops here irrespective of the results."
    ),
    file.path(out_dir, "FINAL_AUDIT_COMPLETE.txt")
  )
  message("Completed integrated closure stage: ", stage_dirs$integrate)
}

message("Requested Track B final-audit stages finished.")



