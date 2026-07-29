options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(GenomicSEM)
})

resolve_base_dir <- function() {
  explicit_base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (nzchar(explicit_base_dir) && dir.exists(explicit_base_dir)) {
    return(normalizePath(explicit_base_dir, winslash = "/", mustWork = TRUE))
  }

  project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
  if (dir.exists(project_dir)) {
    return(normalizePath(project_dir, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate the data root. Set DCV_BASE_DIR before source().")
}

prepare_trait <- function(path, trait_name, trait_type = c("binary", "continuous")) {
  trait_type <- match.arg(trait_type)
  dat <- read_tsv(path, show_col_types = FALSE, progress = FALSE)

  if (!"MAF" %in% names(dat)) {
    if ("eaf" %in% names(dat)) {
      dat$MAF <- pmin(suppressWarnings(as.numeric(dat$eaf)), 1 - suppressWarnings(as.numeric(dat$eaf)))
    } else if ("EAF" %in% names(dat)) {
      dat$MAF <- pmin(suppressWarnings(as.numeric(dat$EAF)), 1 - suppressWarnings(as.numeric(dat$EAF)))
    } else {
      stop("Missing both MAF and eaf/EAF columns in input file: ", path)
    }
  }

  dat <- dat %>%
    transmute(
      CHR = suppressWarnings(as.numeric(CHR)),
      BP = suppressWarnings(as.numeric(BP)),
      SNP = as.character(SNP),
      A1 = as.character(effect_allele),
      A2 = as.character(other_allele),
      MAF = suppressWarnings(as.numeric(MAF)),
      N = suppressWarnings(as.numeric(N)),
      beta_raw = suppressWarnings(as.numeric(beta)),
      se_raw = suppressWarnings(as.numeric(se)),
      p_raw = suppressWarnings(as.numeric(pval))
    ) %>%
    filter(!is.na(SNP), nzchar(SNP), !is.na(A1), !is.na(A2), !is.na(MAF), !is.na(N), !is.na(beta_raw))

  z_val <- ifelse(!is.na(dat$se_raw) & dat$se_raw > 0, dat$beta_raw / dat$se_raw, NA_real_)
  z_from_p <- sign(dat$beta_raw) * qnorm(dat$p_raw / 2, lower.tail = FALSE)
  dat$Z <- ifelse(is.finite(z_val), z_val, z_from_p)
  varSNP <- 2 * dat$MAF * (1 - dat$MAF)

  if (trait_type == "binary") {
    denom <- sqrt((dat$beta_raw^2) * varSNP + (pi^2) / 3)
    beta_std <- dat$beta_raw / denom
    se_std <- dat$se_raw / denom
  } else {
    beta_std <- dat$Z / sqrt(dat$N * varSNP)
    se_std <- 1 / sqrt(dat$N * varSNP)
  }

  dat %>%
    transmute(
      CHR, BP, SNP, A1, A2, MAF,
      !!paste0("beta.", trait_name) := beta_std,
      !!paste0("se.", trait_name) := se_std
    ) %>%
    filter(is.finite(.data[[paste0("beta.", trait_name)]]), is.finite(.data[[paste0("se.", trait_name)]]))
}

merge_trait <- function(base_df, add_df) {
  merged <- inner_join(base_df, add_df, by = "SNP", suffix = c("", ".new"))
  same <- merged$A1 == merged$A1.new & merged$A2 == merged$A2.new
  flip <- merged$A1 == merged$A2.new & merged$A2 == merged$A1.new
  keep <- same | flip
  merged <- merged[keep, , drop = FALSE]

  beta_new <- grep("^beta\\..*\\.new$", names(merged), value = TRUE)
  se_new <- grep("^se\\..*\\.new$", names(merged), value = TRUE)
  if (length(beta_new) == 1) {
    merged[[beta_new]] <- ifelse(flip[keep], -merged[[beta_new]], merged[[beta_new]])
  }

  merged %>%
    transmute(
      CHR = CHR,
      BP = BP,
      SNP = SNP,
      A1 = A1,
      A2 = A2,
      MAF = MAF,
      across(matches("^beta\\.|^se\\."), identity),
      across(all_of(c(beta_new, se_new)), identity, .names = "{.col}")
    ) %>%
    rename_with(~ gsub("\\.new$", "", .x))
}

base_dir <- resolve_base_dir()
trackb_dir <- file.path(base_dir, "analysis_ready_core", "trackB")
model_dir <- file.path(trackb_dir, "genomicsem_model_comparison")
out_dir <- file.path(trackb_dir, "factor_gwas_qsnp")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

raw_files <- c(
  mdd = file.path(base_dir, "GWAS", "formatted", "ieu-b-102_mdd_METAL.txt"),
  lonely = file.path(base_dir, "GWAS", "formatted", "ukb-b-8476_loneliness_METAL.txt"),
  frailty = file.path(base_dir, "GWAS", "formatted", "ebi-a-GCST90020053_frailty_MTAG.txt"),
  srh = file.path(base_dir, "GWAS", "formatted", "ukb-b-6306_self_rated_health_MTAG.txt"),
  walkpace = file.path(base_dir, "GWAS", "formatted", "ukb-b-4711_walking_pace_MTAG.txt")
)

missing_required <- raw_files[!file.exists(raw_files)]
if (length(missing_required) > 0) {
  stop("Missing required input file(s): ", paste(missing_required, collapse = "; "))
}

# Build a 5-trait SNP-level summary object compatible with userGWASa.
mdd_dat <- prepare_trait(raw_files[["mdd"]], "mdd", "binary")
lonely_dat <- prepare_trait(raw_files[["lonely"]], "lonely", "binary")
frailty_dat <- prepare_trait(raw_files[["frailty"]], "frailty", "continuous")
srh_dat <- prepare_trait(raw_files[["srh"]], "srh", "continuous")
walk_dat <- prepare_trait(raw_files[["walkpace"]], "walkpace", "continuous")

sumstats_af <- mdd_dat %>%
  merge_trait(lonely_dat) %>%
  merge_trait(frailty_dat) %>%
  merge_trait(srh_dat) %>%
  merge_trait(walk_dat)

saveRDS(sumstats_af, file.path(out_dir, "trackB_afffunc_sumstats.rds"))

# Recompute a matching 5-trait LDSC structure rather than reusing the 7-trait object.
hm3_file <- normalizePath(Sys.getenv("TRACKB_LDSC_HM3_FILE", unset = file.path(base_dir, "reference", "ldsc", "w_hm3.snplist")), winslash = "/", mustWork = TRUE)
ld_dir <- normalizePath(Sys.getenv("TRACKB_LDSC_LD_DIR", unset = file.path(base_dir, "reference", "ldsc", "eur_w_ld_chr")), winslash = "/", mustWork = TRUE)
wld_dir <- normalizePath(Sys.getenv("TRACKB_LDSC_WLD_DIR", unset = file.path(base_dir, "reference", "ldsc", "weights_hm3_no_hla")), winslash = "/", mustWork = TRUE)
munged_dir <- file.path(model_dir, "munged_sumstats")
dir.create(munged_dir, recursive = TRUE, showWarnings = FALSE)
trait_names <- c("mdd", "lonely", "frailty", "srh", "walkpace")
munged_files <- file.path(munged_dir, paste0(trait_names, ".sumstats.gz"))
if (!all(file.exists(munged_files))) {
  previous_wd <- getwd()
  setwd(munged_dir)
  tryCatch(
    GenomicSEM::munge(
      files = unname(raw_files),
      hm3 = hm3_file,
      trait.names = trait_names,
      info.filter = 0.9,
      maf.filter = 0.01,
      parallel = FALSE,
      overwrite = TRUE
    ),
    finally = setwd(previous_wd)
  )
}
if (!all(file.exists(munged_files))) {
  stop("GenomicSEM::munge did not create all five expected .sumstats.gz files.")
}

ldsc_af <- GenomicSEM::ldsc(
  traits = munged_files,
  sample.prev = rep(NA_real_, 5),
  population.prev = rep(NA_real_, 5),
  ld = ld_dir,
  wld = wld_dir,
  trait.names = trait_names,
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = file.path(out_dir, "trackB_afffunc_ldsc"),
  stand = TRUE,
  select = FALSE,
  chisq.max = NA
)
saveRDS(ldsc_af, file.path(out_dir, "trackB_afffunc_ldsc_covstruc.rds"))

model_afffunc_snp <- "
aff_func =~ mdd + lonely + frailty + srh + walkpace
aff_func ~ SNP
"

gwas_afffunc <- GenomicSEM::userGWAS(
  covstruc = ldsc_af,
  SNPs = sumstats_af,
  estimation = "DWLS",
  model = model_afffunc_snp,
  printwarn = TRUE,
  sub = FALSE,
  parallel = FALSE,
  std.lv = TRUE,
  fix_measurement = TRUE,
  Q_SNP = TRUE,
  analytic = TRUE,
  batch_size = 50000
)

factor_results <- as_tibble(gwas_afffunc) %>%
  rename_with(~ gsub("^beta_aff_func$", "beta_factor", .x)) %>%
  rename_with(~ gsub("^SE_aff_func$", "se_factor", .x)) %>%
  rename_with(~ gsub("^Z_beta_aff_func$", "z_factor", .x)) %>%
  rename_with(~ gsub("^p_val_aff_func$", "p_factor", .x)) %>%
  rename_with(~ gsub("^Q_omnibus$", "q_snp", .x)) %>%
  rename_with(~ gsub("^Q_omnibus_df$", "q_snp_df", .x)) %>%
  rename_with(~ gsub("^Q_omnibus_pval$", "q_snp_p", .x))

write_csv(factor_results, file.path(out_dir, "trackB_afffunc_factor_gwas_qsnp.csv"), na = "")

summary_tbl <- tibble(
  n_snps = nrow(factor_results),
  n_factor_gws = sum(!is.na(factor_results$p_factor) & factor_results$p_factor < 5e-8),
  n_qsnp_gws = sum(!is.na(factor_results$q_snp_p) & factor_results$q_snp_p < 5e-8),
  min_p_factor = suppressWarnings(min(factor_results$p_factor, na.rm = TRUE)),
  min_p_qsnp = suppressWarnings(min(factor_results$q_snp_p, na.rm = TRUE))
)
write_csv(summary_tbl, file.path(out_dir, "trackB_afffunc_factor_gwas_summary.csv"), na = "")

manhattan_prep <- factor_results %>%
  mutate(
    CHR = suppressWarnings(as.numeric(CHR)),
    BP = suppressWarnings(as.numeric(BP)),
    neglog10_p_factor = -log10(p_factor),
    neglog10_p_qsnp = -log10(q_snp_p)
  ) %>%
  filter(!is.na(CHR), !is.na(BP))

write_csv(manhattan_prep, file.path(out_dir, "trackB_afffunc_manhattan_input.csv"), na = "")

plot_manhattan <- function(df, y_col, title_text, out_file) {
  if (!y_col %in% names(df)) return(invisible(NULL))
  p <- ggplot(df, aes(x = BP, y = .data[[y_col]], color = factor(CHR %% 2))) +
    geom_point(size = 0.35, alpha = 0.7) +
    facet_wrap(~ CHR, scales = "free_x", nrow = 1) +
    geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "firebrick") +
    labs(title = title_text, x = "Genomic position within chromosome", y = "-log10(P)", color = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", strip.text = element_text(size = 7), panel.grid.minor = element_blank())
  ggsave(out_file, p, width = 14, height = 4.8, dpi = 300)
}

plot_manhattan(
  manhattan_prep,
  "neglog10_p_factor",
  "Track B Affective-Functional Factor GWAS",
  file.path(fig_dir, "figure_trackB_afffunc_factor_gwas_manhattan.png")
)

plot_manhattan(
  manhattan_prep,
  "neglog10_p_qsnp",
  "Track B Affective-Functional Q_SNP",
  file.path(fig_dir, "figure_trackB_afffunc_qsnp_manhattan.png")
)

qq_from_p <- function(pvals, title_text, out_file) {
  pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 100) return(invisible(NULL))
  obs <- -log10(sort(pvals))
  exp <- -log10(ppoints(length(pvals)))
  qq_df <- tibble(expected = exp, observed = obs)
  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 0.5, alpha = 0.6, color = "#2c7fb8") +
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = 2) +
    labs(title = title_text, x = "Expected -log10(P)", y = "Observed -log10(P)") +
    theme_minimal(base_size = 11)
  ggsave(out_file, p, width = 5.5, height = 5.5, dpi = 300)
}

qq_from_p(
  factor_results$p_factor,
  "QQ Plot: Affective-Functional Factor GWAS",
  file.path(fig_dir, "figure_trackB_afffunc_factor_gwas_qq.png")
)

qq_from_p(
  factor_results$q_snp_p,
  "QQ Plot: Affective-Functional Q_SNP",
  file.path(fig_dir, "figure_trackB_afffunc_qsnp_qq.png")
)

lead_hits <- factor_results %>%
  mutate(CHR = suppressWarnings(as.numeric(CHR)), BP = suppressWarnings(as.numeric(BP))) %>%
  filter((!is.na(p_factor) & p_factor < 5e-8) | (!is.na(q_snp_p) & q_snp_p < 5e-8)) %>%
  arrange(p_factor, q_snp_p)

write_csv(lead_hits, file.path(out_dir, "trackB_afffunc_significant_hits.csv"), na = "")

notes <- c(
  "Track B Factor GWAS + Q_SNP v1",
  "Primary downstream latent target: affective-functional factor.",
  "This script manually constructs a 5-trait SNP-level summary object because GenomicSEM::sumstats() should not be rerun on already munged .sumstats.gz files.",
  "A fresh 5-trait LDSC covariance structure is recomputed to match the SNP-level summary object."
)
writeLines(notes, file.path(out_dir, "trackB_factor_gwas_qsnp_note.txt"))

message("Done.")
message("Outputs written to: ", out_dir)
