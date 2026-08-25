
#!/usr/bin/env Rscript
# Build the locked ADNI factor-PGS CDR-SB analysis dataset.
# PGS weights and all sample rules are fixed before examining CDR-SB results.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

required <- c("dplyr", "readr", "tibble")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))

project_dir <- normalizePath(Sys.getenv("DCV_PROJECT_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
clinical_file <- Sys.getenv("ADNI_CLINICAL_FILE", unset = "")
pgs_file <- Sys.getenv("ADNI_FACTOR_PGS_FILE", unset = "")
output_root <- Sys.getenv("ADNI_OUTPUT_ROOT", unset = file.path(project_dir, "results", "adni_validation"))
out_dir <- file.path(output_root, "dataset")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(clinical_file) || !file.exists(clinical_file)) stop("Set ADNI_CLINICAL_FILE to the governed ADNIMERGE CSV.")
if (!nzchar(pgs_file) || !file.exists(pgs_file)) stop("Set ADNI_FACTOR_PGS_FILE to the frozen factor-PGS TSV.")

clean_id <- function(x) trimws(gsub('^"|"$', "", as.character(x)))
pgs <- utils::read.delim(pgs_file, check.names = FALSE, stringsAsFactors = FALSE)
clinical <- utils::read.csv(clinical_file, check.names = FALSE, stringsAsFactors = FALSE,
                            na.strings = c("", "NA", "N/A", "."))
if ("IID" %in% names(pgs) && !"PTID" %in% names(pgs)) names(pgs)[names(pgs) == "IID"] <- "PTID"
pgs$PTID <- clean_id(pgs$PTID)
clinical$PTID <- clean_id(clinical$PTID)

pgs_columns <- c("PTID", "PGS_affective_z", "PGS_functional_z", "primary_ancestry_proxy", paste0("PC", 1:10))
clinical_columns <- c("PTID", "EXAMDATE", "CDRSB", "AGE", "PTGENDER", "PTEDUCAT", "DX_bl", "APOE4")
if (!all(pgs_columns %in% names(pgs))) stop("Frozen PGS file is missing: ", paste(setdiff(pgs_columns, names(pgs)), collapse = ", "))
if (!all(clinical_columns %in% names(clinical))) stop("ADNIMERGE file is missing: ", paste(setdiff(clinical_columns, names(clinical)), collapse = ", "))

dat <- dplyr::transmute(
  clinical,
  PTID = .data$PTID,
  visit_date = as.Date(.data$EXAMDATE),
  CDRSB = suppressWarnings(as.numeric(.data$CDRSB)),
  age = suppressWarnings(as.numeric(.data$AGE)),
  sex = as.character(.data$PTGENDER),
  education = suppressWarnings(as.numeric(.data$PTEDUCAT)),
  diagnosis_baseline = as.character(.data$DX_bl),
  APOE4 = suppressWarnings(as.numeric(.data$APOE4))
) |>
  dplyr::inner_join(dplyr::select(pgs, dplyr::all_of(pgs_columns)), by = "PTID") |>
  dplyr::filter(.data$primary_ancestry_proxy %in% TRUE) |>
  dplyr::filter(!is.na(.data$visit_date), is.finite(.data$CDRSB), is.finite(.data$age), !is.na(.data$sex)) |>
  dplyr::group_by(.data$PTID) |>
  dplyr::arrange(.data$visit_date, .by_group = TRUE) |>
  dplyr::mutate(
    baseline_date = dplyr::first(.data$visit_date),
    time_years = as.numeric(.data$visit_date - .data$baseline_date) / 365.25,
    baseline_age = dplyr::first(.data$age),
    baseline_diagnosis = dplyr::first(.data$diagnosis_baseline),
    baseline_education = dplyr::first(.data$education),
    baseline_APOE4 = dplyr::first(.data$APOE4),
    baseline_sex = dplyr::first(.data$sex)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    baseline_age_c = as.numeric(scale(.data$baseline_age)),
    education_c = as.numeric(scale(.data$baseline_education)),
    PGS_affective_z = as.numeric(scale(.data$PGS_affective_z)),
    PGS_functional_z = as.numeric(scale(.data$PGS_functional_z)),
    sex = factor(.data$baseline_sex),
    baseline_diagnosis = factor(.data$baseline_diagnosis)
  ) |>
  dplyr::filter(is.finite(.data$time_years), .data$time_years >= 0, is.finite(.data$baseline_age_c),
                is.finite(.data$education_c), !is.na(.data$baseline_diagnosis), !is.na(.data$sex))

if (dplyr::n_distinct(dat$PTID) < 200L || nrow(dat) < 500L) stop("Insufficient linked ADNI sample after prespecified exclusions.")

protocol <- tibble::tribble(
  ~field, ~value,
  "analysis_status", "prespecified independent ADNI clinical validation",
  "primary_outcome", "longitudinal CDR-SB change",
  "primary_exposures", "frozen affective-factor and functional-factor PGS",
  "primary_sample", "primary ancestry-proxy sample",
  "prohibited", "ADNI-based SNP selection, PGS retuning, outcome switching, or time-window selection"
)
flow <- tibble::tibble(
  stage = c("clinical_rows", "linked_primary_ancestry_rows", "linked_primary_ancestry_participants", "participants_with_2plus_CDRSB_visits"),
  n = c(nrow(clinical), nrow(dat), dplyr::n_distinct(dat$PTID), sum(table(dat$PTID) >= 2L))
)
readr::write_csv(protocol, file.path(out_dir, "01_frozen_protocol.csv"), na = "")
readr::write_csv(flow, file.path(out_dir, "02_sample_flow.csv"), na = "")
readr::write_csv(dat, file.path(out_dir, "03_analysis_dataset_locked.csv"), na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
message("Locked ADNI analysis dataset written to: ", out_dir)
