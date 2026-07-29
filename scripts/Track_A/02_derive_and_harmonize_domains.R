# ==============================================================================
# Track A compact analysis script 2
# Final longitudinal analysis components.
# ==============================================================================

# ------------------------------------------------------------------------------
# Component: derive_domain_variables.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_core4_dir <- file.path(base_dir, "analysis_ready_core", "ready_core4")
input_hrs_dir <- file.path(base_dir, "analysis_ready_core", "ready_hrs")
out_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  file.path(input_core4_dir, "share_analysis_ready_long_core4.rds"),
  file.path(input_core4_dir, "elsa_analysis_ready_long_core4.rds"),
  file.path(input_core4_dir, "mhas_analysis_ready_long_core4.rds"),
  file.path(input_core4_dir, "charls_analysis_ready_long_core4.rds"),
  file.path(input_hrs_dir, "hrs_analysis_ready_long.rds")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

z_score <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

z_by_cohort <- function(x, cohort) {
  ave(x, cohort, FUN = z_score)
}

mean_na <- function(...) {
  mat <- cbind(...)
  out <- rowMeans(mat, na.rm = TRUE)
  all_missing <- apply(is.na(mat), 1, all)
  out[all_missing] <- NA_real_
  out
}

count_nonmissing <- function(...) {
  mat <- cbind(...)
  rowSums(!is.na(mat))
}

read_inputs <- function() {
  files <- list(
    SHARE = file.path(input_core4_dir, "share_analysis_ready_long_core4.rds"),
    ELSA = file.path(input_core4_dir, "elsa_analysis_ready_long_core4.rds"),
    MHAS = file.path(input_core4_dir, "mhas_analysis_ready_long_core4.rds"),
    CHARLS = file.path(input_core4_dir, "charls_analysis_ready_long_core4.rds"),
    HRS = file.path(input_hrs_dir, "hrs_analysis_ready_long.rds")
  )
  out <- purrr::imap(files, ~ readRDS(.x) %>% mutate(cohort = .y))
  dplyr::bind_rows(out)
}

derive_domains <- function(df) {
  df %>%
    mutate(
      depression_z = z_by_cohort(depression_score, cohort),
      loneliness_z = z_by_cohort(loneliness_score, cohort),
      cognition_z = z_by_cohort(cognition_score, cohort),
      adl_z = z_by_cohort(adl_count, cohort),
      iadl_z = z_by_cohort(iadl_count, cohort),
      srh_z = z_by_cohort(self_rated_health, cohort),
      social_contact_z = z_by_cohort(social_contact, cohort),
      social_participation_z = z_by_cohort(social_participation, cohort),
      internet_z = z_by_cohort(internet_use, cohort),
      email_z = z_by_cohort(email_use, cohort),
      digital_index_z = z_by_cohort(digital_index, cohort),

      affective_score = mean_na(depression_z, loneliness_z),
      affective_component_n = count_nonmissing(depression_score, loneliness_score),

      functional_score = mean_na(adl_z, iadl_z, srh_z),
      functional_component_n = count_nonmissing(adl_count, iadl_count, self_rated_health),
      functional_any_limitation = ifelse(
        is.na(adl_any) & is.na(iadl_any),
        NA_real_,
        pmax(adl_any, iadl_any, na.rm = TRUE)
      ),

      basic_cognition_score = cognition_z,
      basic_cognition_component_n = ifelse(is.na(cognition_score), 0, 1),

      social_exposure_patch = mean_na(social_contact_z, social_participation_z),
      social_patch_component_n = count_nonmissing(social_contact, social_participation),

      digital_exposure_patch = mean_na(internet_z, email_z, digital_index_z),
      digital_patch_component_n = count_nonmissing(internet_use, email_use, digital_index),

      affective_ready = ifelse(!is.na(depression_score), 1, 0),
      functional_ready = ifelse(functional_component_n >= 2, 1, 0),
      cognition_ready = ifelse(!is.na(cognition_score), 1, 0),
      social_patch_ready = ifelse(social_patch_component_n >= 1, 1, 0),
      digital_patch_ready = ifelse(digital_patch_component_n >= 1, 1, 0)
    )
}

build_dictionary <- function(df) {
  tibble(
    variable = names(df),
    class = purrr::map_chr(df, ~ paste(class(.x), collapse = ";")),
    non_missing_n = purrr::map_int(df, ~ sum(!is.na(.x))),
    example = purrr::map_chr(df, ~ {
      vals <- unique(.x[!is.na(.x)])
      if (length(vals) == 0) return("")
      paste(head(as.character(vals), 3), collapse = " | ")
    })
  )
}

save_split_outputs <- function(df) {
  split(df, df$cohort) |>
    purrr::iwalk(function(x, nm) {
      saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_domain_ready.rds")))
    })
}

analysis_ready <- read_inputs()
domain_ready <- derive_domains(analysis_ready)

save_split_outputs(domain_ready)
saveRDS(domain_ready, file.path(out_dir, "pooled_domain_ready.rds"))
write_csv(build_dictionary(domain_ready), file.path(out_dir, "domain_variable_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: build_domain_coverage_report.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
out_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")

pooled_file <- file.path(input_dir, "pooled_domain_ready.rds")
if (!file.exists(pooled_file)) {
  stop("Missing required input file: ", pooled_file)
}

domain_ready <- readRDS(pooled_file)

score_status <- function(domain_type, pct_non_missing) {
  if (domain_type == "mature") {
    if (is.na(pct_non_missing)) return("unknown")
    if (pct_non_missing >= 0.70) return("ready_for_harmonization")
    if (pct_non_missing >= 0.40) return("usable_with_caution")
    return("needs_targeted_fix")
  }

  if (is.na(pct_non_missing)) return("unknown")
  if (pct_non_missing >= 0.50) return("strong_patch_layer")
  if (pct_non_missing >= 0.20) return("partial_patch_layer")
  "major_gap"
}

domain_specs <- tribble(
  ~domain_group, ~domain_type, ~ready_var, ~score_var, ~component_var,
  "affective", "mature", "affective_ready", "affective_score", "affective_component_n",
  "functional", "mature", "functional_ready", "functional_score", "functional_component_n",
  "basic_cognition", "mature", "cognition_ready", "basic_cognition_score", "basic_cognition_component_n",
  "social_exposure_patch", "patch", "social_patch_ready", "social_exposure_patch", "social_patch_component_n",
  "digital_exposure_patch", "patch", "digital_patch_ready", "digital_exposure_patch", "digital_patch_component_n"
)

domain_coverage_report <- purrr::pmap_dfr(
  domain_specs,
  function(domain_group, domain_type, ready_var, score_var, component_var) {
    domain_group_value <- domain_group
    domain_type_value <- domain_type

    domain_ready %>%
      group_by(cohort) %>%
      summarise(
        n_rows = n(),
        ready_n = sum(.data[[ready_var]] == 1, na.rm = TRUE),
        ready_pct = ready_n / n_rows,
        score_non_missing_n = sum(!is.na(.data[[score_var]])),
        score_non_missing_pct = score_non_missing_n / n_rows,
        mean_component_n = mean(.data[[component_var]], na.rm = TRUE),
        median_component_n = stats::median(.data[[component_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        domain_group = domain_group_value,
        domain_type = domain_type_value,
        status = vapply(
          score_non_missing_pct,
          function(x) score_status(domain_type_value, x),
          character(1)
        )
      ) %>%
      relocate(domain_group, domain_type, .before = cohort)
  }
)

cohort_variable_coverage <- domain_ready %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    across(
      c(
        depression_score, loneliness_score, cognition_score,
        adl_count, iadl_count, self_rated_health,
        social_contact, social_participation,
        internet_use, email_use, digital_index
      ),
      ~ sum(!is.na(.x)),
      .names = "{.col}_non_missing_n"
    ),
    .groups = "drop"
  )

cohort_readiness_summary <- domain_coverage_report %>%
  filter(domain_type == "mature") %>%
  group_by(cohort) %>%
  summarise(
    n_mature_domains_ready = sum(status == "ready_for_harmonization"),
    n_mature_domains_usable = sum(status %in% c("ready_for_harmonization", "usable_with_caution")),
    .groups = "drop"
  ) %>%
  mutate(
    next_step = case_when(
      n_mature_domains_ready == 3 ~ "enter_harmonization_now",
      n_mature_domains_usable == 3 ~ "harmonization_with_documented_caution",
      TRUE ~ "targeted_domain_fix_first"
    )
  )

patch_priority_report <- domain_coverage_report %>%
  filter(domain_type == "patch") %>%
  arrange(cohort, score_non_missing_pct)

write_csv(domain_coverage_report, file.path(out_dir, "domain_coverage_report.csv"), na = "")
write_csv(cohort_variable_coverage, file.path(out_dir, "cohort_variable_coverage.csv"), na = "")
write_csv(cohort_readiness_summary, file.path(out_dir, "cohort_readiness_summary.csv"), na = "")
write_csv(patch_priority_report, file.path(out_dir, "patch_priority_report.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: build_harmonization_input_affective_functional.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
out_dir <- file.path(base_dir, "analysis_ready_core", "harmonization_input")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pooled_file <- file.path(input_dir, "pooled_domain_ready.rds")
if (!file.exists(pooled_file)) {
  stop("Missing required input file: ", pooled_file)
}

domain_ready <- readRDS(pooled_file)

count_nonmissing <- function(...) {
  mat <- cbind(...)
  rowSums(!is.na(mat))
}

harmonization_input <- domain_ready %>%
  transmute(
    cohort,
    source_file,
    source_type,
    id,
    wave,
    year,
    country,
    age,
    sex,
    education_years,
    education_level,
    marital_status,
    partnered,
    race_ethnicity,
    weight_person,
    weight_household,
    weight_longitudinal,

    depression_score,
    loneliness_score,
    affective_score,
    affective_component_n,
    affective_ready,

    adl_count,
    iadl_count,
    adl_any,
    iadl_any,
    self_rated_health,
    functional_any_limitation,
    functional_score,
    functional_component_n,
    functional_ready,

    affective_functional_joint_ready =
      ifelse(affective_ready == 1 & functional_ready == 1, 1, 0),
    affective_functional_complete_case =
      ifelse(
        count_nonmissing(
          affective_score, functional_score,
          age, sex, education_level, marital_status
        ) == 6,
        1, 0
      ),
    affective_source_profile = case_when(
      !is.na(depression_score) & !is.na(loneliness_score) ~ "depression_plus_loneliness",
      !is.na(depression_score) ~ "depression_only",
      !is.na(loneliness_score) ~ "loneliness_only",
      TRUE ~ NA_character_
    ),
    functional_source_profile = case_when(
      !is.na(adl_count) & !is.na(iadl_count) & !is.na(self_rated_health) ~ "adl_plus_iadl_plus_srh",
      !is.na(adl_count) & !is.na(iadl_count) ~ "adl_plus_iadl",
      !is.na(adl_count) & !is.na(self_rated_health) ~ "adl_plus_srh",
      !is.na(iadl_count) & !is.na(self_rated_health) ~ "iadl_plus_srh",
      !is.na(adl_count) ~ "adl_only",
      !is.na(iadl_count) ~ "iadl_only",
      !is.na(self_rated_health) ~ "srh_only",
      TRUE ~ NA_character_
    )
  )

split(harmonization_input, harmonization_input$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_harmonization_input_affective_functional.rds")))
  })

harmonization_summary <- harmonization_input %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    affective_ready_n = sum(affective_ready == 1, na.rm = TRUE),
    functional_ready_n = sum(functional_ready == 1, na.rm = TRUE),
    joint_ready_n = sum(affective_functional_joint_ready == 1, na.rm = TRUE),
    complete_case_n = sum(affective_functional_complete_case == 1, na.rm = TRUE),
    affective_score_non_missing_n = sum(!is.na(affective_score)),
    functional_score_non_missing_n = sum(!is.na(functional_score)),
    age_non_missing_n = sum(!is.na(age)),
    sex_non_missing_n = sum(!is.na(sex)),
    education_level_non_missing_n = sum(!is.na(education_level)),
    marital_status_non_missing_n = sum(!is.na(marital_status)),
    .groups = "drop"
  )

harmonization_dictionary <- tibble(
  variable = names(harmonization_input),
  class = purrr::map_chr(harmonization_input, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(harmonization_input, ~ sum(!is.na(.x))),
  example = purrr::map_chr(harmonization_input, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(harmonization_input, file.path(out_dir, "pooled_harmonization_input_affective_functional.rds"))
write_csv(harmonization_summary, file.path(out_dir, "harmonization_input_affective_functional_summary.csv"), na = "")
write_csv(harmonization_dictionary, file.path(out_dir, "harmonization_input_affective_functional_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: patch_charls_cognition.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
domain_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
patch_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready_patches")
dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)

charls_domain_file <- file.path(domain_dir, "charls_domain_ready.rds")
pooled_domain_file <- file.path(domain_dir, "pooled_domain_ready.rds")
charls_raw_file <- file.path(base_dir, "CHARLS", "02_harmonized", "Regular Waves", "H_CHARLS_D_Data", "H_CHARLS_D_Data.dta")

for (f in c(charls_domain_file, pooled_domain_file, charls_raw_file)) {
  if (!file.exists(f)) stop("Missing required input file: ", f)
}

safe_pull <- function(df, nm, type = c("numeric", "character")) {
  type <- match.arg(type)
  if (!nm %in% names(df)) {
    if (type == "character") return(rep(NA_character_, nrow(df)))
    return(rep(NA_real_, nrow(df)))
  }
  x <- df[[nm]]
  if (inherits(x, "labelled")) x <- haven::zap_labels(x)
  if (type == "character") return(as.character(x))
  suppressWarnings(as.numeric(x))
}

na_codes_to_na <- function(x, extra = c(-1, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

row_sum_or_na <- function(...) {
  mat <- cbind(...)
  out <- rowSums(mat, na.rm = TRUE)
  out[apply(is.na(mat), 1, all)] <- NA_real_
  out
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

charls_domain <- readRDS(charls_domain_file)
pooled_domain <- readRDS(pooled_domain_file)
charls_raw <- read_dta(charls_raw_file)

wave_to_year <- c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018)

charls_cognition_patch <- bind_rows(lapply(1:4, function(w) {
  imrc <- na_codes_to_na(safe_pull(charls_raw, paste0("r", w, "imrc")))
  orient <- na_codes_to_na(safe_pull(charls_raw, paste0("r", w, "orient")))
  tr20 <- na_codes_to_na(safe_pull(charls_raw, paste0("r", w, "tr20")))

  patch_raw <- row_sum_or_na(imrc, orient, tr20)

  tibble(
    id = safe_pull(charls_raw, "ID", "character"),
    wave = w,
    year = unname(wave_to_year[as.character(w)]),
    charls_imrc_patch = imrc,
    charls_orient_patch = orient,
    charls_tr20_patch = tr20,
    charls_cognition_patch_raw = patch_raw,
    charls_cognition_patch_component_n =
      rowSums(!is.na(cbind(imrc, orient, tr20))),
    charls_cognition_patch_source = ifelse(
      !is.na(patch_raw),
      "charls_harmonized_imrc_orient_tr20",
      NA_character_
    )
  )
}))

charls_patched <- charls_domain %>%
  left_join(charls_cognition_patch, by = c("id", "wave", "year")) %>%
  mutate(
    cognition_score_before_patch = cognition_score,
    cognition_score_patched = dplyr::coalesce(
      cognition_score,
      charls_cognition_patch_raw
    ),
    basic_cognition_score_before_patch = basic_cognition_score,
    basic_cognition_score_patched = z_score(cognition_score_patched),
    cognition_ready_before_patch = cognition_ready,
    cognition_ready_patched = ifelse(!is.na(cognition_score_patched), 1, 0),
    basic_cognition_component_n_patched = ifelse(
      !is.na(cognition_score_patched),
      1,
      0
    ),
    charls_cognition_patch_applied = ifelse(
      is.na(cognition_score_before_patch) & !is.na(charls_cognition_patch_raw),
      1,
      0
    )
  )

pooled_patched <- pooled_domain %>%
  filter(cohort != "CHARLS") %>%
  bind_rows(charls_patched)

patch_summary <- charls_patched %>%
  group_by(wave, year) %>%
  summarise(
    n_rows = n(),
    cognition_before_non_missing_n = sum(!is.na(cognition_score_before_patch)),
    cognition_patch_non_missing_n = sum(!is.na(charls_cognition_patch_raw)),
    cognition_after_non_missing_n = sum(!is.na(cognition_score_patched)),
    patch_applied_n = sum(charls_cognition_patch_applied == 1, na.rm = TRUE),
    .groups = "drop"
  )

patch_dictionary <- tibble(
  variable = names(charls_patched),
  class = purrr::map_chr(charls_patched, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(charls_patched, ~ sum(!is.na(.x))),
  example = purrr::map_chr(charls_patched, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(charls_patched, file.path(patch_dir, "charls_domain_ready_patched_cognition.rds"))
saveRDS(pooled_patched, file.path(patch_dir, "pooled_domain_ready_patched_charls_cognition.rds"))
write_csv(patch_summary, file.path(patch_dir, "charls_cognition_patch_summary.csv"), na = "")
write_csv(patch_dictionary, file.path(patch_dir, "charls_cognition_patch_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", patch_dir)

})

# ------------------------------------------------------------------------------
# Component: patch_social_digital_exposure.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
domain_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
patch_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready_patches")
dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)

pooled_domain_file <- file.path(domain_dir, "pooled_domain_ready.rds")
if (!file.exists(pooled_domain_file)) {
  stop("Missing required input file: ", pooled_domain_file)
}

domain_ready <- readRDS(pooled_domain_file)

safe_pull <- function(df, nm, type = c("numeric", "character")) {
  type <- match.arg(type)
  if (!nm %in% names(df)) {
    if (type == "character") return(rep(NA_character_, nrow(df)))
    return(rep(NA_real_, nrow(df)))
  }
  x <- df[[nm]]
  if (inherits(x, "labelled")) x <- haven::zap_labels(x)
  if (type == "character") return(as.character(x))
  suppressWarnings(as.numeric(x))
}

na_codes_to_na <- function(x, extra = c(-1, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

mean_na <- function(...) {
  mat <- cbind(...)
  out <- rowMeans(mat, na.rm = TRUE)
  out[apply(is.na(mat), 1, all)] <- NA_real_
  out
}

label_text <- function(df) {
  vapply(df, function(col) {
    lbl <- attr(col, "label")
    if (is.null(lbl)) "" else as.character(lbl)
  }, character(1))
}

matches_any <- function(text, patterns) {
  if (length(text) == 0 || is.na(text) || !nzchar(text)) return(FALSE)
  any(vapply(patterns, function(p) grepl(p, text, ignore.case = TRUE, perl = TRUE), logical(1)))
}

find_candidate_vars <- function(df, patterns, prefix = NULL) {
  nms <- names(df)
  lbls <- label_text(df)
  keep <- vapply(seq_along(nms), function(i) {
    name_ok <- if (is.null(prefix)) TRUE else startsWith(nms[[i]], prefix)
    if (!name_ok) return(FALSE)
    matches_any(nms[[i]], patterns) || matches_any(lbls[[i]], patterns)
  }, logical(1))
  nms[keep]
}

binary_from_candidates <- function(df, vars) {
  if (length(vars) == 0) return(rep(NA_real_, nrow(df)))
  mat <- sapply(vars, function(v) na_codes_to_na(safe_pull(df, v)))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
  pos <- mat > 0
  pos[is.na(pos)] <- FALSE
  nonmiss <- !is.na(mat)
  out <- ifelse(rowSums(nonmiss) == 0, NA_real_, ifelse(rowSums(pos) > 0, 1, 0))
  out
}

extract_patch_wide <- function(df, id_var, waves, year_map, social_patterns, digital_patterns) {
  bind_rows(lapply(waves, function(w) {
    prefix <- paste0("r", w)
    social_vars <- find_candidate_vars(df, social_patterns, prefix = prefix)
    digital_vars <- find_candidate_vars(df, digital_patterns, prefix = prefix)

    tibble(
      id = safe_pull(df, id_var, "character"),
      wave = w,
      year = unname(year_map[as.character(w)]),
      social_participation_patch_raw = binary_from_candidates(df, social_vars),
      digital_exposure_patch_raw = binary_from_candidates(df, digital_vars),
      social_patch_source = ifelse(length(social_vars) > 0, paste(social_vars, collapse = "|"), NA_character_),
      digital_patch_source = ifelse(length(digital_vars) > 0, paste(digital_vars, collapse = "|"), NA_character_)
    )
  }))
}

extract_patch_long <- function(df, id_var, wave_var, year_var = NULL, social_patterns, digital_patterns) {
  social_vars <- find_candidate_vars(df, social_patterns)
  digital_vars <- find_candidate_vars(df, digital_patterns)

  tibble(
    id = safe_pull(df, id_var, "character"),
    wave = safe_pull(df, wave_var),
    year = if (is.null(year_var)) rep(NA_real_, nrow(df)) else safe_pull(df, year_var),
    social_participation_patch_raw = binary_from_candidates(df, social_vars),
    digital_exposure_patch_raw = binary_from_candidates(df, digital_vars),
    social_patch_source = ifelse(length(social_vars) > 0, paste(social_vars, collapse = "|"), NA_character_),
    digital_patch_source = ifelse(length(digital_vars) > 0, paste(digital_vars, collapse = "|"), NA_character_)
  )
}

social_patterns <- c(
  "club", "volunt", "activity", "social", "group", "organization",
  "course", "community", "friend", "caclub", "asiste_club"
)

digital_patterns <- c(
  "internet", "email", "online", "web", "computer", "digital",
  "mobile", "smart", "tablet", "comu_telef_comp"
)

share_raw_file <- file.path(base_dir, "SHARE", "easySHARE_rel9-0-0_stata", "easySHARE_rel9-0-0.dta")
elsa_raw_file <- file.path(base_dir, "ELSA", "5050stata_A06F12276F837D82374CA07FCEA89BDC3945439B39E962E7B03B39CEBABE46BC_V1", "UKDA-5050-stata", "stata", "stata13_se", "wave_8_elsa_data_eul.dta")
mhas_raw_file <- file.path(base_dir, "MHAS", "03_constructed", "simplemhas", "simpleMHAS.dta")
charls_raw_file <- file.path(base_dir, "CHARLS", "02_harmonized", "Regular Waves", "H_CHARLS_D_Data", "H_CHARLS_D_Data.dta")

share_patch <- if (file.exists(share_raw_file)) {
  extract_patch_long(read_dta(share_raw_file), "mergeid", "wave", "int_year", social_patterns, digital_patterns)
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

elsa_patch <- if (file.exists(elsa_raw_file)) {
  x <- read_dta(elsa_raw_file)
  tibble(
    id = safe_pull(x, "idauniq", "character"),
    wave = 8,
    year = 2016,
    social_participation_patch_raw = binary_from_candidates(x, find_candidate_vars(x, social_patterns)),
    digital_exposure_patch_raw = binary_from_candidates(x, find_candidate_vars(x, digital_patterns)),
    social_patch_source = paste(find_candidate_vars(x, social_patterns), collapse = "|"),
    digital_patch_source = paste(find_candidate_vars(x, digital_patterns), collapse = "|")
  )
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

mhas_patch <- if (file.exists(mhas_raw_file)) {
  x <- read_dta(mhas_raw_file)
  tibble(
    id = safe_pull(x, "np", "character"),
    wave = safe_pull(x, "ronda"),
    year = safe_pull(x, "a_o_ent"),
    social_participation_patch_raw = binary_from_candidates(x, c("asiste_club", "voluntario")),
    digital_exposure_patch_raw = binary_from_candidates(x, c("comu_telef_comp")),
    social_patch_source = "asiste_club|voluntario",
    digital_patch_source = "comu_telef_comp"
  )
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

charls_patch <- if (file.exists(charls_raw_file)) {
  extract_patch_wide(
    read_dta(charls_raw_file),
    id_var = "ID",
    waves = 1:4,
    year_map = c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018),
    social_patterns = social_patterns,
    digital_patterns = digital_patterns
  )
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

patch_map <- bind_rows(
  share_patch %>% mutate(cohort = "SHARE"),
  elsa_patch %>% mutate(cohort = "ELSA"),
  mhas_patch %>% mutate(cohort = "MHAS"),
  charls_patch %>% mutate(cohort = "CHARLS")
) %>%
  distinct(cohort, id, wave, year, .keep_all = TRUE)

patched <- domain_ready %>%
  left_join(patch_map, by = c("cohort", "id", "wave", "year")) %>%
  mutate(
    social_participation_before_patch = social_participation,
    digital_before_patch = digital_exposure_patch,

    social_participation_after_patch = dplyr::coalesce(
      social_participation,
      social_participation_patch_raw
    ),
    digital_exposure_after_patch_raw = dplyr::coalesce(
      digital_index,
      internet_use,
      email_use,
      digital_exposure_patch_raw
    ),

    social_participation_after_patch_z =
      ave(social_participation_after_patch, cohort, FUN = z_score),
    digital_exposure_after_patch_z =
      ave(digital_exposure_after_patch_raw, cohort, FUN = z_score),

    social_exposure_patch = mean_na(social_contact_z, social_participation_after_patch_z),
    digital_exposure_patch = digital_exposure_after_patch_z,

    social_patch_component_n =
      rowSums(!is.na(cbind(social_contact, social_participation_after_patch))),
    digital_patch_component_n =
      rowSums(!is.na(cbind(digital_exposure_after_patch_raw))),

    social_patch_ready = ifelse(social_patch_component_n >= 1, 1, 0),
    digital_patch_ready = ifelse(digital_patch_component_n >= 1, 1, 0),

    social_patch_applied = ifelse(
      is.na(social_participation_before_patch) & !is.na(social_participation_patch_raw),
      1,
      0
    ),
    digital_patch_applied = ifelse(
      is.na(digital_before_patch) & !is.na(digital_exposure_patch_raw),
      1,
      0
    )
  )

split(patched, patched$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(patch_dir, paste0(tolower(nm), "_domain_ready_patched_social_digital.rds")))
  })

patch_summary <- patched %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    social_before_non_missing_n = sum(!is.na(social_exposure_patch)),
    social_after_non_missing_n = sum(!is.na(social_exposure_patch)),
    social_participation_patch_non_missing_n = sum(!is.na(social_participation_patch_raw)),
    social_patch_applied_n = sum(social_patch_applied == 1, na.rm = TRUE),
    digital_before_non_missing_n = sum(!is.na(digital_exposure_patch)),
    digital_after_non_missing_n = sum(!is.na(digital_exposure_patch)),
    digital_patch_non_missing_n = sum(!is.na(digital_exposure_patch_raw)),
    digital_patch_applied_n = sum(digital_patch_applied == 1, na.rm = TRUE),
    social_patch_source_example = dplyr::first(na.omit(social_patch_source)),
    digital_patch_source_example = dplyr::first(na.omit(digital_patch_source)),
    .groups = "drop"
  )

patch_dictionary <- tibble(
  variable = names(patched),
  class = purrr::map_chr(patched, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(patched, ~ sum(!is.na(.x))),
  example = purrr::map_chr(patched, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(patched, file.path(patch_dir, "pooled_domain_ready_patched_social_digital.rds"))
write_csv(patch_summary, file.path(patch_dir, "social_digital_patch_summary.csv"), na = "")
write_csv(patch_dictionary, file.path(patch_dir, "social_digital_patch_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", patch_dir)

})

# ------------------------------------------------------------------------------
# Component: harmonize_affective_functional.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "harmonization_input")
out_dir <- file.path(base_dir, "analysis_ready_core", "harmonized_affective_functional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "pooled_harmonization_input_affective_functional.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

harm_input <- readRDS(input_file)

na_codes_to_na <- function(x, extra = c(-1, -3, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

z_by_group <- function(x, grp) {
  ave(x, grp, FUN = z_score)
}

affective_functional_harmonized <- harm_input %>%
  mutate(
    across(
      c(age, sex, education_years, education_level, marital_status, partnered,
        depression_score, loneliness_score, adl_count, iadl_count, self_rated_health,
        affective_score, functional_score),
      na_codes_to_na
    ),
    sex_female = case_when(
      cohort == "SHARE" & sex %in% c(0, 1) ~ sex,
      sex == 2 ~ 1,
      sex == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    education_level_clean = education_level,
    marital_status_clean = marital_status,
    affective_harmonized = z_by_group(affective_score, cohort),
    functional_harmonized = z_by_group(functional_score, cohort),
    affective_harmonized_wave = z_by_group(
      affective_score,
      interaction(cohort, wave, drop = TRUE)
    ),
    functional_harmonized_wave = z_by_group(
      functional_score,
      interaction(cohort, wave, drop = TRUE)
    ),
    affective_function_joint_harmonized = rowMeans(
      cbind(affective_harmonized, functional_harmonized),
      na.rm = TRUE
    ),
    affective_function_joint_harmonized = ifelse(
      is.nan(affective_function_joint_harmonized),
      NA_real_,
      affective_function_joint_harmonized
    ),
    harmonization_eligibility = ifelse(
      affective_functional_joint_ready == 1 &
        !is.na(age) &
        !is.na(sex_female) &
        !is.na(education_level_clean),
      1,
      0
    )
  )

split(affective_functional_harmonized, affective_functional_harmonized$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_affective_functional_harmonized.rds")))
  })

harmonization_summary <- affective_functional_harmonized %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    affective_non_missing_n = sum(!is.na(affective_harmonized)),
    functional_non_missing_n = sum(!is.na(functional_harmonized)),
    joint_non_missing_n = sum(!is.na(affective_function_joint_harmonized)),
    eligible_n = sum(harmonization_eligibility == 1, na.rm = TRUE),
    sex_female_non_missing_n = sum(!is.na(sex_female)),
    education_non_missing_n = sum(!is.na(education_level_clean)),
    marital_non_missing_n = sum(!is.na(marital_status_clean)),
    .groups = "drop"
  )

harmonization_dictionary <- tibble(
  variable = names(affective_functional_harmonized),
  class = purrr::map_chr(affective_functional_harmonized, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(affective_functional_harmonized, ~ sum(!is.na(.x))),
  example = purrr::map_chr(affective_functional_harmonized, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(affective_functional_harmonized, file.path(out_dir, "pooled_affective_functional_harmonized.rds"))
write_csv(harmonization_summary, file.path(out_dir, "affective_functional_harmonization_summary.csv"), na = "")
write_csv(harmonization_dictionary, file.path(out_dir, "affective_functional_harmonization_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: integrate_charls_cognition_patch.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
domain_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
patch_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready_patches")
out_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready_integrated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pooled_domain_file <- file.path(domain_dir, "pooled_domain_ready.rds")
charls_patch_file <- file.path(patch_dir, "charls_domain_ready_patched_cognition.rds")

for (f in c(pooled_domain_file, charls_patch_file)) {
  if (!file.exists(f)) stop("Missing required input file: ", f)
}

pooled_domain <- readRDS(pooled_domain_file)
charls_patch <- readRDS(charls_patch_file)

integrated <- pooled_domain %>%
  filter(cohort != "CHARLS") %>%
  bind_rows(
    charls_patch %>%
      mutate(
        cognition_score = cognition_score_patched,
        basic_cognition_score = basic_cognition_score_patched,
        basic_cognition_component_n = basic_cognition_component_n_patched,
        cognition_ready = cognition_ready_patched
      )
  ) %>%
  arrange(cohort, id, wave, year)

integration_summary <- integrated %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    cognition_score_non_missing_n = sum(!is.na(cognition_score)),
    basic_cognition_non_missing_n = sum(!is.na(basic_cognition_score)),
    cognition_ready_n = sum(cognition_ready == 1, na.rm = TRUE),
    .groups = "drop"
  )

charls_delta_summary <- charls_patch %>%
  summarise(
    n_rows = n(),
    cognition_before_non_missing_n = sum(!is.na(cognition_score_before_patch)),
    cognition_after_non_missing_n = sum(!is.na(cognition_score_patched)),
    basic_cognition_before_non_missing_n = sum(!is.na(basic_cognition_score_before_patch)),
    basic_cognition_after_non_missing_n = sum(!is.na(basic_cognition_score_patched)),
    patch_applied_n = sum(charls_cognition_patch_applied == 1, na.rm = TRUE)
  )

integration_dictionary <- tibble(
  variable = names(integrated),
  class = purrr::map_chr(integrated, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(integrated, ~ sum(!is.na(.x))),
  example = purrr::map_chr(integrated, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(integrated, file.path(out_dir, "pooled_domain_ready_integrated_charls_cognition.rds"))
write_csv(integration_summary, file.path(out_dir, "integrated_charls_cognition_summary.csv"), na = "")
write_csv(charls_delta_summary, file.path(out_dir, "charls_cognition_integration_delta.csv"), na = "")
write_csv(integration_dictionary, file.path(out_dir, "integrated_charls_cognition_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: build_basic_cognition_harmonization_input.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready_integrated")
out_dir <- file.path(base_dir, "analysis_ready_core", "harmonization_input")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "pooled_domain_ready_integrated_charls_cognition.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

cog_input_base <- readRDS(input_file)

na_codes_to_na <- function(x, extra = c(-1, -3, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

z_by_group <- function(x, grp) {
  ave(x, grp, FUN = z_score)
}

count_nonmissing <- function(...) {
  mat <- cbind(...)
  rowSums(!is.na(mat))
}

basic_cognition_harmonization_input <- cog_input_base %>%
  mutate(
    across(
      c(age, sex, education_years, education_level, marital_status,
        cognition_score, basic_cognition_score),
      na_codes_to_na
    ),
    sex_female = case_when(
      cohort == "SHARE" & sex %in% c(0, 1) ~ sex,
      sex == 2 ~ 1,
      sex == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    education_level_clean = education_level,
    marital_status_clean = marital_status,
    cognition_harmonized = z_by_group(cognition_score, cohort),
    cognition_harmonized_wave = z_by_group(
      cognition_score,
      interaction(cohort, wave, drop = TRUE)
    ),
    cognition_harmonization_eligibility = ifelse(
      cognition_ready == 1 &
        !is.na(age) &
        !is.na(sex_female) &
        !is.na(education_level_clean),
      1,
      0
    ),
    cognition_complete_case = ifelse(
      count_nonmissing(
        cognition_score, age, sex_female,
        education_level_clean, marital_status_clean
      ) == 5,
      1,
      0
    ),
    cognition_source_profile = case_when(
      cohort == "CHARLS" & !is.na(charls_cognition_patch_source) ~ charls_cognition_patch_source,
      TRUE ~ "cohort_native_cognition_field"
    )
  ) %>%
  transmute(
    cohort,
    source_file,
    source_type,
    id,
    wave,
    year,
    country,
    age,
    sex,
    sex_female,
    education_years,
    education_level,
    education_level_clean,
    marital_status,
    marital_status_clean,
    partnered,
    race_ethnicity,
    weight_person,
    weight_household,
    weight_longitudinal,
    cognition_score,
    basic_cognition_score,
    cognition_harmonized,
    cognition_harmonized_wave,
    basic_cognition_component_n,
    cognition_ready,
    cognition_harmonization_eligibility,
    cognition_complete_case,
    cognition_source_profile
  )

split(basic_cognition_harmonization_input, basic_cognition_harmonization_input$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_basic_cognition_harmonization_input.rds")))
  })

cognition_summary <- basic_cognition_harmonization_input %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    cognition_non_missing_n = sum(!is.na(cognition_harmonized)),
    eligibility_n = sum(cognition_harmonization_eligibility == 1, na.rm = TRUE),
    complete_case_n = sum(cognition_complete_case == 1, na.rm = TRUE),
    age_non_missing_n = sum(!is.na(age)),
    sex_non_missing_n = sum(!is.na(sex_female)),
    education_non_missing_n = sum(!is.na(education_level_clean)),
    marital_non_missing_n = sum(!is.na(marital_status_clean)),
    .groups = "drop"
  )

cognition_dictionary <- tibble(
  variable = names(basic_cognition_harmonization_input),
  class = purrr::map_chr(basic_cognition_harmonization_input, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(basic_cognition_harmonization_input, ~ sum(!is.na(.x))),
  example = purrr::map_chr(basic_cognition_harmonization_input, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(basic_cognition_harmonization_input, file.path(out_dir, "pooled_basic_cognition_harmonization_input.rds"))
write_csv(cognition_summary, file.path(out_dir, "basic_cognition_harmonization_summary.csv"), na = "")
write_csv(cognition_dictionary, file.path(out_dir, "basic_cognition_harmonization_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: patch_social_digital_exposure.R
# ------------------------------------------------------------------------------
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  root <- Sys.getenv("DCV_BASE_DIR", unset = Sys.getenv("DCV_PROJECT_DIR", unset = getwd()))
  if (!dir.exists(root)) stop("DCV_BASE_DIR does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
domain_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready")
patch_dir <- file.path(base_dir, "analysis_ready_core", "domain_ready_patches")
dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)

pooled_domain_file <- file.path(domain_dir, "pooled_domain_ready.rds")
if (!file.exists(pooled_domain_file)) {
  stop("Missing required input file: ", pooled_domain_file)
}

domain_ready <- readRDS(pooled_domain_file)

safe_pull <- function(df, nm, type = c("numeric", "character")) {
  type <- match.arg(type)
  if (!nm %in% names(df)) {
    if (type == "character") return(rep(NA_character_, nrow(df)))
    return(rep(NA_real_, nrow(df)))
  }
  x <- df[[nm]]
  if (inherits(x, "labelled")) x <- haven::zap_labels(x)
  if (type == "character") return(as.character(x))
  suppressWarnings(as.numeric(x))
}

na_codes_to_na <- function(x, extra = c(-1, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

mean_na <- function(...) {
  mat <- cbind(...)
  out <- rowMeans(mat, na.rm = TRUE)
  out[apply(is.na(mat), 1, all)] <- NA_real_
  out
}

first_nonmissing_string <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

label_text <- function(df) {
  vapply(df, function(col) {
    lbl <- attr(col, "label")
    if (is.null(lbl)) "" else as.character(lbl)
  }, character(1))
}

matches_any <- function(text, patterns) {
  if (length(text) == 0 || is.na(text) || !nzchar(text)) return(FALSE)
  any(vapply(patterns, function(p) grepl(p, text, ignore.case = TRUE, perl = TRUE), logical(1)))
}

find_candidate_vars <- function(df, patterns, prefix = NULL) {
  nms <- names(df)
  lbls <- label_text(df)
  keep <- vapply(seq_along(nms), function(i) {
    name_ok <- if (is.null(prefix)) TRUE else startsWith(nms[[i]], prefix)
    if (!name_ok) return(FALSE)
    matches_any(nms[[i]], patterns) || matches_any(lbls[[i]], patterns)
  }, logical(1))
  nms[keep]
}

binary_from_candidates <- function(df, vars) {
  if (length(vars) == 0) return(rep(NA_real_, nrow(df)))
  mat <- sapply(vars, function(v) na_codes_to_na(safe_pull(df, v)))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
  pos <- mat > 0
  pos[is.na(pos)] <- FALSE
  nonmiss <- !is.na(mat)
  ifelse(rowSums(nonmiss) == 0, NA_real_, ifelse(rowSums(pos) > 0, 1, 0))
}

extract_patch_wide <- function(df, id_var, waves, year_map, social_patterns, digital_patterns) {
  bind_rows(lapply(waves, function(w) {
    prefix <- paste0("r", w)
    social_vars <- find_candidate_vars(df, social_patterns, prefix = prefix)
    digital_vars <- find_candidate_vars(df, digital_patterns, prefix = prefix)

    tibble(
      id = safe_pull(df, id_var, "character"),
      wave = w,
      year = unname(year_map[as.character(w)]),
      social_participation_patch_raw = binary_from_candidates(df, social_vars),
      digital_exposure_patch_raw = binary_from_candidates(df, digital_vars),
      social_patch_source = ifelse(length(social_vars) > 0, paste(social_vars, collapse = "|"), NA_character_),
      digital_patch_source = ifelse(length(digital_vars) > 0, paste(digital_vars, collapse = "|"), NA_character_)
    )
  }))
}

extract_patch_long <- function(df, id_var, wave_var, year_var = NULL, social_patterns, digital_patterns) {
  social_vars <- find_candidate_vars(df, social_patterns)
  digital_vars <- find_candidate_vars(df, digital_patterns)

  tibble(
    id = safe_pull(df, id_var, "character"),
    wave = safe_pull(df, wave_var),
    year = if (is.null(year_var)) rep(NA_real_, nrow(df)) else safe_pull(df, year_var),
    social_participation_patch_raw = binary_from_candidates(df, social_vars),
    digital_exposure_patch_raw = binary_from_candidates(df, digital_vars),
    social_patch_source = ifelse(length(social_vars) > 0, paste(social_vars, collapse = "|"), NA_character_),
    digital_patch_source = ifelse(length(digital_vars) > 0, paste(digital_vars, collapse = "|"), NA_character_)
  )
}

social_patterns <- c(
  "club", "volunt", "activity", "social", "group", "organization",
  "course", "community", "friend", "caclub", "asiste_club", "alone", "mmalone"
)

digital_patterns <- c(
  "internet", "email", "online", "web", "computer", "digital",
  "mobile", "smart", "tablet", "comu_telef_comp"
)

share_raw_file <- file.path(base_dir, "SHARE", "easySHARE_rel9-0-0_stata", "easySHARE_rel9-0-0.dta")
elsa_raw_file <- file.path(base_dir, "ELSA", "5050stata_A06F12276F837D82374CA07FCEA89BDC3945439B39E962E7B03B39CEBABE46BC_V1", "UKDA-5050-stata", "stata", "stata13_se", "wave_8_elsa_data_eul.dta")
mhas_raw_file <- file.path(base_dir, "MHAS", "03_constructed", "simplemhas", "simpleMHAS.dta")
charls_raw_file <- file.path(base_dir, "CHARLS", "02_harmonized", "Regular Waves", "H_CHARLS_D_Data", "H_CHARLS_D_Data.dta")

share_patch <- if (file.exists(share_raw_file)) {
  extract_patch_long(read_dta(share_raw_file), "mergeid", "wave", "int_year", social_patterns, digital_patterns)
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

elsa_patch <- if (file.exists(elsa_raw_file)) {
  x <- read_dta(elsa_raw_file)
  tibble(
    id = safe_pull(x, "idauniq", "character"),
    wave = 8,
    year = 2016,
    social_participation_patch_raw = binary_from_candidates(x, find_candidate_vars(x, social_patterns)),
    digital_exposure_patch_raw = binary_from_candidates(x, find_candidate_vars(x, digital_patterns)),
    social_patch_source = first_nonmissing_string(find_candidate_vars(x, social_patterns)),
    digital_patch_source = first_nonmissing_string(find_candidate_vars(x, digital_patterns))
  )
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

mhas_patch <- if (file.exists(mhas_raw_file)) {
  x <- read_dta(mhas_raw_file)
  tibble(
    id = safe_pull(x, "np", "character"),
    wave = safe_pull(x, "ronda"),
    year = safe_pull(x, "a_o_ent"),
    social_participation_patch_raw = binary_from_candidates(x, c("asiste_club", "voluntario")),
    digital_exposure_patch_raw = binary_from_candidates(x, c("comu_telef_comp")),
    social_patch_source = "asiste_club|voluntario",
    digital_patch_source = "comu_telef_comp"
  )
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

charls_patch <- if (file.exists(charls_raw_file)) {
  extract_patch_wide(
    read_dta(charls_raw_file),
    id_var = "ID",
    waves = 1:4,
    year_map = c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018),
    social_patterns = social_patterns,
    digital_patterns = digital_patterns
  )
} else {
  tibble(id = character(), wave = numeric(), year = numeric(),
         social_participation_patch_raw = numeric(), digital_exposure_patch_raw = numeric(),
         social_patch_source = character(), digital_patch_source = character())
}

patch_map <- bind_rows(
  share_patch %>% mutate(cohort = "SHARE"),
  elsa_patch %>% mutate(cohort = "ELSA"),
  mhas_patch %>% mutate(cohort = "MHAS"),
  charls_patch %>% mutate(cohort = "CHARLS")
) %>%
  distinct(cohort, id, wave, year, .keep_all = TRUE)

patched <- domain_ready %>%
  left_join(patch_map, by = c("cohort", "id", "wave", "year")) %>%
  mutate(
    social_participation_before_patch = social_participation,
    digital_before_patch = digital_exposure_patch,

    social_participation_after_patch = dplyr::coalesce(
      social_participation,
      social_participation_patch_raw
    ),
    digital_exposure_after_patch_raw = dplyr::coalesce(
      digital_index,
      internet_use,
      email_use,
      digital_exposure_patch_raw
    ),

    social_participation_after_patch_z =
      ave(social_participation_after_patch, cohort, FUN = z_score),
    digital_exposure_after_patch_z =
      ave(digital_exposure_after_patch_raw, cohort, FUN = z_score),

    social_exposure_patch = mean_na(social_contact_z, social_participation_after_patch_z),
    digital_exposure_patch = digital_exposure_after_patch_z,

    social_patch_component_n =
      rowSums(!is.na(cbind(social_contact, social_participation_after_patch))),
    digital_patch_component_n =
      rowSums(!is.na(cbind(digital_exposure_after_patch_raw))),

    social_patch_ready = ifelse(social_patch_component_n >= 1, 1, 0),
    digital_patch_ready = ifelse(digital_patch_component_n >= 1, 1, 0),

    social_patch_applied = ifelse(
      is.na(social_participation_before_patch) & !is.na(social_participation_patch_raw),
      1,
      0
    ),
    digital_patch_applied = ifelse(
      is.na(digital_before_patch) & !is.na(digital_exposure_patch_raw),
      1,
      0
    )
  )

split(patched, patched$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(patch_dir, paste0(tolower(nm), "_domain_ready_patched_social_digital.rds")))
  })

patch_summary <- patched %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    social_before_non_missing_n = sum(!is.na(social_exposure_patch)),
    social_after_non_missing_n = sum(!is.na(social_exposure_patch)),
    social_participation_patch_non_missing_n = sum(!is.na(social_participation_patch_raw)),
    social_patch_applied_n = sum(social_patch_applied == 1, na.rm = TRUE),
    digital_before_non_missing_n = sum(!is.na(digital_exposure_patch)),
    digital_after_non_missing_n = sum(!is.na(digital_exposure_patch)),
    digital_patch_non_missing_n = sum(!is.na(digital_exposure_patch_raw)),
    digital_patch_applied_n = sum(digital_patch_applied == 1, na.rm = TRUE),
    social_patch_source_example = first_nonmissing_string(social_patch_source),
    digital_patch_source_example = first_nonmissing_string(digital_patch_source),
    .groups = "drop"
  )

patch_dictionary <- tibble(
  variable = names(patched),
  class = purrr::map_chr(patched, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(patched, ~ sum(!is.na(.x))),
  example = purrr::map_chr(patched, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(patched, file.path(patch_dir, "pooled_domain_ready_patched_social_digital.rds"))
write_csv(patch_summary, file.path(patch_dir, "social_digital_patch_summary.csv"), na = "")
write_csv(patch_dictionary, file.path(patch_dir, "social_digital_patch_dictionary.csv"), na = "")

message("Done.")
message("Outputs written to: ", patch_dir)

})


