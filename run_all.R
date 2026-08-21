project_dir <- Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
pipeline <- unique(trimws(strsplit(Sys.getenv("DCV_RUN_PIPELINE", unset = ""), ",", fixed = TRUE)[[1]]))
pipeline <- pipeline[nzchar(pipeline)]
valid <- c(
  "track_a", "longitudinal_robustness", "mortality_bridge", "dementia_bridge",
  "track_b_pre_fuma", "track_b_post_fuma", "track_d",
  "track_c_architecture", "track_c_gene_sets", "track_c_omics",
  "track_c_molecular", "track_b_audit"
)
if (!length(pipeline)) stop("Set DCV_RUN_PIPELINE to one or more documented stages before source('run_all.R').")
if (any(!pipeline %in% valid)) stop("Unknown pipeline stage: ", paste(setdiff(pipeline, valid), collapse = ", "))

for (stage in pipeline) {
  message("Running ", stage)
  if (stage == "track_a") source(file.path(project_dir, "run_track_A.R"), encoding = "UTF-8")
  if (stage == "longitudinal_robustness") {
    source(file.path(project_dir, "run_longitudinal_robustness.R"), encoding = "UTF-8")
  }
  if (stage == "mortality_bridge") {
    source(file.path(project_dir, "run_mortality_bridge.R"), encoding = "UTF-8")
  }
  if (stage == "dementia_bridge") {
    source(file.path(project_dir, "run_dementia_bridge.R"), encoding = "UTF-8")
  }
  if (stage == "track_b_pre_fuma") {
    Sys.setenv(TRACKB_RUN_STAGES = "pre_fuma"); source(file.path(project_dir, "run_track_B.R"), encoding = "UTF-8")
  }
  if (stage == "track_b_post_fuma") {
    Sys.setenv(TRACKB_RUN_STAGES = "post_fuma"); source(file.path(project_dir, "run_track_B.R"), encoding = "UTF-8")
  }
  if (stage == "track_d") source(file.path(project_dir, "run_track_D.R"), encoding = "UTF-8")
  if (stage == "track_c_architecture") {
    Sys.setenv(TRACKC_RUN_STAGES = "architecture"); source(file.path(project_dir, "run_track_C.R"), encoding = "UTF-8")
  }
  if (stage == "track_c_gene_sets") {
    Sys.setenv(TRACKC_RUN_STAGES = "gene_sets"); source(file.path(project_dir, "run_track_C.R"), encoding = "UTF-8")
  }
  if (stage == "track_c_omics") {
    Sys.setenv(TRACKC_RUN_STAGES = "omics"); source(file.path(project_dir, "run_track_C.R"), encoding = "UTF-8")
  }
  if (stage == "track_c_molecular") {
    Sys.setenv(TRACKC_RUN_STAGES = "molecular"); source(file.path(project_dir, "run_track_C.R"), encoding = "UTF-8")
  }
  if (stage == "track_b_audit") {
    Sys.setenv(TRACKB_RUN_STAGES = "audit"); source(file.path(project_dir, "run_track_B.R"), encoding = "UTF-8")
  }
}
