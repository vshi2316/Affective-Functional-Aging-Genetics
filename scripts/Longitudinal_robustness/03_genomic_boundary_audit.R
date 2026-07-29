## Genomic structural boundary audit.
## Existing five-GWAS input only. Writes to a dedicated frozen output directory.

options(stringsAsFactors = FALSE)
base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir)) stop("Set DCV_BASE_DIR before running this script.")
if (!dir.exists(base_dir)) stop("DCV_BASE_DIR does not exist: ", base_dir)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
trackb <- file.path(base_dir, "analysis_ready_core", "trackB_v1")
cov_file <- file.path(trackb, "factor_gwas_qsnp_v1", "trackB_afffunc_ldsc_covstruc_v1.rds")
out_dir <- file.path(base_dir, "robustness_analysis", "02_genomicsem_boundary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("GenomicSEM", quietly = TRUE)) stop("Install GenomicSEM before running this script.")
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble); library(GenomicSEM)})
if (!file.exists(cov_file)) stop("Missing covstruc: ", cov_file)
covstruc <- readRDS(cov_file)
free_syntax <- "affective_domain =~ mdd + lonely\nfunctional_domain =~ frailty + srh + walkpace\naffective_domain ~~ functional_domain"
boundary_syntax <- "affective_domain =~ mdd + lonely\nfunctional_domain =~ frailty + srh + walkpace\naffective_domain ~~ 1*functional_domain"
safe <- function(s) tryCatch(GenomicSEM::usermodel(covstruc = covstruc, estimation = "DWLS", model = s, CFIcalc = TRUE, std.lv = TRUE), error = function(e) e)
free <- safe(free_syntax); boundary <- safe(boundary_syntax)
extract_fit <- function(x, label) {
  if (inherits(x, "error") || is.null(x$modelfit)) return(tibble(model = label, status = "failed", error = if (inherits(x, "error")) conditionMessage(x) else "modelfit unavailable"))
  m <- as.data.frame(x$modelfit)
  get1 <- function(nm) if (nm %in% names(m)) suppressWarnings(as.numeric(m[[nm]][1])) else NA_real_
  tibble(model = label, status = "estimated", error = NA_character_, chisq = get1("chisq"),
    df = get1("df"), p_chisq = get1("p_chisq"), CFI = get1("CFI"), SRMR = get1("SRMR"))
}
fit <- bind_rows(extract_fit(free, "correlated_two_factor_free"), extract_fit(boundary, "rho_equal_1_boundary"))
write_csv(fit, file.path(out_dir, "genomicsem_free_vs_rho1_fit.csv"), na = "")
free_par <- if (!inherits(free, "error") && !is.null(free$results)) as_tibble(free$results) else tibble()
if (nrow(free_par)) {
  se_col <- intersect(c("STD_All_SE", "STD_Genotype_SE", "Unstand_SE"), names(free_par))[1]
  est_col <- intersect(c("STD_All", "STD_Genotype", "Unstand_Est"), names(free_par))[1]
  rho <- free_par %>% filter(lhs %in% c("affective_domain", "functional_domain"), rhs %in% c("affective_domain", "functional_domain"), op == "~~", lhs != rhs)
  rho$estimate <- if (!is.na(est_col)) suppressWarnings(as.numeric(rho[[est_col]])) else NA_real_
  rho$se <- if (!is.na(se_col)) suppressWarnings(as.numeric(rho[[se_col]])) else NA_real_
  p_col <- if ("p_value" %in% names(rho)) "p_value" else if ("pval" %in% names(rho)) "pval" else NA_character_
  rho <- rho %>% transmute(lhs, op, rhs, estimate, se, p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(.data[[p_col]])) else NA_real_)
} else rho <- tibble()
write_csv(rho, file.path(out_dir, "free_factor_correlation_parameter.csv"), na = "")
writeLines(c(
  "Interpretation rule:",
  "A converged rho=1 boundary model and a valid GenomicSEM-compatible constraint comparison may support statistical separation within this GWAS covariance model.",
  "A failed or inadmissible boundary model is a boundary limitation, not proof of one factor or two factors.",
  "Never call rho<1 biological independence or external replication."
), file.path(out_dir, "genomicsem_readme.txt"))
message("GenomicSEM boundary audit complete: ", out_dir)
