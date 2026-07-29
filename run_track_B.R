project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
stages <- unique(trimws(strsplit(Sys.getenv("TRACKB_RUN_STAGES", unset = "pre_fuma"), ",", fixed = TRUE)[[1]]))
valid <- c("pre_fuma", "post_fuma", "audit")
if (any(!stages %in% valid)) stop("Unknown Track B stage: ", paste(setdiff(stages, valid), collapse = ", "))

source_track_b <- function(number, name) {
  source(file.path(project_dir, "scripts", "Track_B", sprintf("%02d_%s.R", number, name)), encoding = "UTF-8", chdir = FALSE)
}
if ("pre_fuma" %in% stages) {
  source_track_b(1, "prepare_genetic_covariance")
  source_track_b(2, "compare_genomicsem_models")
  source_track_b(3, "factor_gwas_and_qsnp")
  source_track_b(4, "prepare_fuma_uploads")
  source_track_b(5, "local_magma_analysis")
}
if ("post_fuma" %in% stages) source_track_b(6, "import_fuma_results")
if ("audit" %in% stages) source_track_b(7, "factor_gwas_quality_and_increment")
