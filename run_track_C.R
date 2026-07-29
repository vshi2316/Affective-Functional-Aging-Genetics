project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
stages <- unique(trimws(strsplit(Sys.getenv("TRACKC_RUN_STAGES", unset = "architecture,gene_sets"), ",", fixed = TRUE)[[1]]))
valid <- c("architecture", "gene_sets", "omics", "molecular")
if (any(!stages %in% valid)) stop("Unknown Track C stage: ", paste(setdiff(stages, valid), collapse = ", "))
src <- function(number, name) source(file.path(project_dir, "scripts", "Track_C", sprintf("%02d_%s.R", number, name)), encoding = "UTF-8", chdir = FALSE)
if ("architecture" %in% stages) src(1, "local_shared_architecture")
if ("gene_sets" %in% stages) src(2, "freeze_gene_sets")
if ("omics" %in% stages) {
  src(3, "external_omics")
  src(4, "external_omics_adjustment")
  src(5, "external_omics_evidence")
}
if ("molecular" %in% stages) {
  src(6, "candidate_molecular_followup")
  src(7, "brain_eqtl_locus_colocalisation")
  src(8, "multisignal_colocalisation")
  src(9, "molecular_evidence_classification")
}
