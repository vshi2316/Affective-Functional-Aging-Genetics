## Copy this file to config/config.R and edit only local paths.

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_root <- "/path/to/private/data-root"

Sys.setenv(
  DCV_PROJECT_DIR = project_dir,
  DCV_BASE_DIR = data_root,
  TRACKC_ROOT = file.path(data_root, "analysis_ready_core", "trackC"),
  MAGMA_EXECUTABLE = "magma",
  SMR_EXECUTABLE = "smr",
  PLINK_EXECUTABLE = "plink",
  TRACKB_LDSC_HM3_FILE = file.path(data_root, "reference", "ldsc", "w_hm3.snplist"),
  TRACKB_LDSC_LD_DIR = file.path(data_root, "reference", "ldsc", "eur_w_ld_chr"),
  TRACKB_LDSC_WLD_DIR = file.path(data_root, "reference", "ldsc", "weights_hm3_no_hla"),
  TRACKB_MAGMA_BFILE = file.path(data_root, "reference", "magma", "g1000_eur", "g1000_eur"),
  TRACKB_AFFECTIVE_FUMA_N = "500199",
  TRACKB_FUNCTIONAL_FUMA_N = "460844",
  TRACKC_LD_BFILE = file.path(data_root, "reference", "plink", "EUR"),
  MORTALITY_LONGITUDINAL_INPUT =
    file.path(data_root, "mortality_bridge", "longitudinal_input.rds"),
  MORTALITY_EXACT_VISIT_DATE_INPUT =
    file.path(data_root, "mortality_bridge", "exact_visit_dates.csv"),
  MORTALITY_WAVE_STATUS_INPUT =
    file.path(data_root, "mortality_bridge", "wave_interview_status.csv"),
  MORTALITY_ENDPOINT_MAIN_INPUT =
    file.path(data_root, "mortality_bridge", "mortality_endpoints_main.rds"),
  MORTALITY_ENDPOINT_SUPPORT_INPUT =
    file.path(data_root, "mortality_bridge", "mortality_endpoints_support.rds"),
  MORTALITY_OUTPUT_ROOT =
    file.path(project_dir, "results", "mortality_bridge"),
  DEMENTIA_DYNAMIC_CANDIDATE_INPUT =
    file.path(project_dir, "results", "mortality_bridge", "candidate", "dynamic_candidate_exact.rds"),
  DEMENTIA_ID_CROSSWALK_INPUT =
    file.path(data_root, "dementia_bridge", "dementia_id_crosswalk.rds"),
  DEMENTIA_ELSA_OUTCOME_INPUT =
    file.path(data_root, "dementia_bridge", "elsa_harmonized_panel.dta"),
  DEMENTIA_SHARE_OUTCOME_INPUT =
    file.path(data_root, "dementia_bridge", "share_harmonized_panel.dta"),
  DEMENTIA_ELSA_COGNITION_INPUT =
    file.path(data_root, "dementia_bridge", "elsa_analysis_ready_long.rds"),
  DEMENTIA_SHARE_COGNITION_INPUT =
    file.path(data_root, "dementia_bridge", "share_analysis_ready_long.rds"),
  DEMENTIA_OUTPUT_ROOT =
    file.path(project_dir, "results", "dementia_bridge")
)
