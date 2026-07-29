# ==============================================================================
# Track A compact analysis script 1
# Final longitudinal analysis components.
# ==============================================================================

# ------------------------------------------------------------------------------
# Component: build_minimal_inputs_long_core4.R
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
out_dir <- file.path(base_dir, "analysis_ready_core", "ready_core4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

schema_names <- c(
  "cohort", "source_file", "source_type", "id", "wave", "year",
  "country", "age", "sex", "education_years", "education_level",
  "marital_status", "partnered", "race_ethnicity",
  "depression_score", "loneliness_score",
  "cognition_score", "adl_count", "iadl_count",
  "adl_any", "iadl_any", "self_rated_health",
  "chronic_count", "social_participation", "social_contact",
  "internet_use", "email_use", "digital_index",
  "weight_person", "weight_household", "weight_longitudinal",
  "death_year", "death_month", "death_flag", "dropout_flag",
  "inw_flag"
)

make_na <- function(n, mode = c("numeric", "character")) {
  mode <- match.arg(mode)
  if (mode == "character") rep(NA_character_, n) else rep(NA_real_, n)
}

safe_pull <- function(df, nm, type = c("numeric", "character")) {
  type <- match.arg(type)
  if (!nm %in% names(df)) return(make_na(nrow(df), type))
  x <- df[[nm]]
  if (inherits(x, "labelled")) x <- haven::zap_labels(x)
  if (type == "character") return(as.character(x))
  suppressWarnings(as.numeric(x))
}

safe_pull_any <- function(df, nms, type = c("numeric", "character")) {
  type <- match.arg(type)
  present <- intersect(nms, names(df))
  if (length(present) == 0) return(make_na(nrow(df), type))
  out <- safe_pull(df, present[[1]], type)
  if (length(present) > 1) {
    for (nm in present[-1]) out <- dplyr::coalesce(out, safe_pull(df, nm, type))
  }
  out
}

na_codes_to_na <- function(x, extra = c(-1, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

derive_binary_any <- function(x) {
  x <- na_codes_to_na(x)
  ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0))
}

derive_digital_index <- function(internet_use, email_use) {
  iu <- na_codes_to_na(internet_use)
  eu <- na_codes_to_na(email_use)
  out <- rep(NA_real_, length(iu))
  out[!is.na(iu)] <- iu[!is.na(iu)]
  out[is.na(out) & !is.na(eu)] <- eu[is.na(out) & !is.na(eu)]
  out
}

finalize_long <- function(df) {
  df %>%
    mutate(
      across(where(is.numeric), na_codes_to_na),
      adl_any = derive_binary_any(adl_count),
      iadl_any = derive_binary_any(iadl_count),
      digital_index = derive_digital_index(internet_use, email_use)
    ) %>%
    select(all_of(schema_names)) %>%
    distinct()
}

build_dictionary <- function(df, cohort_name) {
  tibble(
    cohort = cohort_name,
    analysis_name = names(df),
    non_missing_n = map_int(df, ~ sum(!is.na(.x))),
    class = map_chr(df, ~ paste(class(.x), collapse = ";")),
    example = map_chr(df, ~ {
      vals <- unique(.x[!is.na(.x)])
      if (length(vals) == 0) return("")
      paste(head(as.character(vals), 3), collapse = " | ")
    })
  )
}

build_share_long <- function() {
  path <- file.path(base_dir, "SHARE", "easySHARE_rel9-0-0_stata", "easySHARE_rel9-0-0.dta")
  x <- read_dta(path)
  tibble(
    cohort = "SHARE",
    source_file = path,
    source_type = "easyshare_long",
    id = safe_pull(x, "mergeid", "character"),
    wave = safe_pull(x, "wave"),
    year = safe_pull_any(x, c("int_year")),
    country = safe_pull_any(x, c("country"), "character"),
    age = safe_pull_any(x, c("age")),
    sex = safe_pull_any(x, c("female")),
    education_years = safe_pull_any(x, c("eduyears_mod")),
    education_level = safe_pull_any(x, c("isced1997_r")),
    marital_status = safe_pull_any(x, c("mar_stat")),
    partnered = safe_pull_any(x, c("partnerinhh")),
    race_ethnicity = make_na(nrow(x), "character"),
    depression_score = safe_pull_any(x, c("eurod")),
    loneliness_score = make_na(nrow(x)),
    cognition_score = rowSums(cbind(
      safe_pull_any(x, c("recall_1")),
      safe_pull_any(x, c("recall_2")),
      safe_pull_any(x, c("orienti")),
      safe_pull_any(x, c("numeracy_1")),
      safe_pull_any(x, c("numeracy_2"))
    ), na.rm = TRUE),
    adl_count = safe_pull_any(x, c("adlwa", "adla")),
    iadl_count = safe_pull_any(x, c("iadlza", "iadla")),
    self_rated_health = safe_pull_any(x, c("sphus")),
    chronic_count = safe_pull_any(x, c("chronic_mod")),
    social_participation = make_na(nrow(x)),
    social_contact = safe_pull_any(x, c("sp008_", "sp009_1_mod")),
    internet_use = make_na(nrow(x)),
    email_use = make_na(nrow(x)),
    weight_person = make_na(nrow(x)),
    weight_household = make_na(nrow(x)),
    weight_longitudinal = make_na(nrow(x)),
    death_year = make_na(nrow(x)),
    death_month = make_na(nrow(x)),
    death_flag = make_na(nrow(x)),
    dropout_flag = make_na(nrow(x)),
    inw_flag = rep(1, nrow(x))
  ) %>%
    mutate(cognition_score = ifelse(cognition_score == 0, NA_real_, cognition_score)) %>%
    finalize_long()
}

build_mhas_long <- function() {
  path <- file.path(base_dir, "MHAS", "04_gateway_harmonized", "H_MHAS_d.dta")
  x <- read_dta(path)
  waves <- 1:6
  years <- c(`1` = 2001, `2` = 2003, `3` = 2012, `4` = 2015, `5` = 2018, `6` = 2021)

  bind_rows(lapply(waves, function(w) {
    tibble(
      cohort = "MHAS",
      source_file = path,
      source_type = "gateway_harmonized_wide",
      id = safe_pull_any(x, c("rahhidnp", "unhhidnp"), "character"),
      wave = w,
      year = safe_pull_any(x, c(paste0("r", w, "iwy"))),
      country = rep("Mexico", nrow(x)),
      age = safe_pull_any(x, c(paste0("r", w, "agey"))),
      sex = safe_pull_any(x, c("ragender")),
      education_years = safe_pull_any(x, c("raedyrs")),
      education_level = safe_pull_any(x, c("raeducl", "raedisced")),
      marital_status = safe_pull_any(x, c(paste0("r", w, "mstat"))),
      partnered = safe_pull_any(x, c(paste0("r", w, "mpart"))),
      race_ethnicity = make_na(nrow(x), "character"),
      depression_score = safe_pull_any(x, c(paste0("r", w, "depres"))),
      loneliness_score = make_na(nrow(x)),
      cognition_score = rowSums(cbind(
        safe_pull_any(x, c(paste0("r", w, "slfmem"))),
        safe_pull_any(x, c(paste0("r", w, "pstmem"))),
        safe_pull_any(x, c(paste0("r", w, "prmem")))
      ), na.rm = TRUE),
      adl_count = safe_pull_any(x, c(paste0("r", w, "adltot6"), paste0("r", w, "adlfive"), paste0("r", w, "adla"))),
      iadl_count = safe_pull_any(x, c(paste0("r", w, "iadlfour"))),
      self_rated_health = safe_pull_any(x, c(paste0("r", w, "shlt"))),
      chronic_count = make_na(nrow(x)),
      social_participation = make_na(nrow(x)),
      social_contact = safe_pull_any(x, c(paste0("r", w, "pcnt"), paste0("r", w, "kcnt"))),
      internet_use = make_na(nrow(x)),
      email_use = make_na(nrow(x)),
      weight_person = safe_pull_any(x, c(paste0("r", w, "wtresp"))),
      weight_household = safe_pull_any(x, c(paste0("r", w, "wthh"))),
      weight_longitudinal = make_na(nrow(x)),
      death_year = safe_pull_any(x, c("radyear")),
      death_month = safe_pull_any(x, c("radmonth")),
      death_flag = safe_pull_any(x, c("radyear")),
      dropout_flag = safe_pull_any(x, c(paste0("r", w, "iwstat"))),
      inw_flag = safe_pull_any(x, c(paste0("inw", w)))
    )
  })) %>%
    mutate(
      year = ifelse(is.na(year), unname(years[as.character(wave)]), year),
      cognition_score = ifelse(cognition_score == 0, NA_real_, cognition_score)
    ) %>%
    finalize_long()
}

build_charls_long <- function() {
  path <- file.path(base_dir, "CHARLS", "02_harmonized", "Regular Waves", "H_CHARLS_D_Data", "H_CHARLS_D_Data.dta")
  x <- read_dta(path)
  waves <- 1:4
  years <- c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018)

  bind_rows(lapply(waves, function(w) {
    tibble(
      cohort = "CHARLS",
      source_file = path,
      source_type = "gateway_harmonized_wide",
      id = safe_pull_any(x, c("ID"), "character"),
      wave = w,
      year = safe_pull_any(x, c(paste0("r", w, "iwy"))),
      country = rep("China", nrow(x)),
      age = safe_pull_any(x, c(paste0("r", w, "agey"))),
      sex = safe_pull_any(x, c("ragender")),
      education_years = safe_pull_any(x, c("raedyrs")),
      education_level = safe_pull_any(x, c("raeducl")),
      marital_status = safe_pull_any(x, c(paste0("r", w, "mstat"))),
      partnered = make_na(nrow(x)),
      race_ethnicity = make_na(nrow(x), "character"),
      depression_score = safe_pull_any(x, c(paste0("r", w, "cesd10"), paste0("r", w, "depresl"))),
      loneliness_score = safe_pull_any(x, c(paste0("r", w, "flonel"))),
      cognition_score = rowSums(cbind(
        safe_pull_any(x, c(paste0("r", w, "memrye"))),
        safe_pull_any(x, c(paste0("r", w, "memryf"))),
        safe_pull_any(x, c(paste0("r", w, "rxmemry")))
      ), na.rm = TRUE),
      adl_count = safe_pull_any(x, c(paste0("r", w, "adlwa"), paste0("r", w, "adlfive"))),
      iadl_count = safe_pull_any(x, c(paste0("r", w, "iadlza"))),
      self_rated_health = safe_pull_any(x, c(paste0("r", w, "shlt"))),
      chronic_count = make_na(nrow(x)),
      social_participation = make_na(nrow(x)),
      social_contact = make_na(nrow(x)),
      internet_use = make_na(nrow(x)),
      email_use = make_na(nrow(x)),
      weight_person = safe_pull_any(x, c(paste0("r", w, "wtresp"))),
      weight_household = safe_pull_any(x, c(paste0("r", w, "wthh"))),
      weight_longitudinal = safe_pull_any(x, c(paste0("r", w, "wthhl"), paste0("r", w, "wtrespl"))),
      death_year = safe_pull_any(x, c("radyear")),
      death_month = safe_pull_any(x, c("radmonth")),
      death_flag = safe_pull_any(x, c("radyear")),
      dropout_flag = safe_pull_any(x, c(paste0("r", w, "iwstat"))),
      inw_flag = safe_pull_any(x, c(paste0("inw", w)))
    )
  })) %>%
    mutate(
      year = ifelse(is.na(year), unname(years[as.character(wave)]), year),
      cognition_score = ifelse(cognition_score == 0, NA_real_, cognition_score)
    ) %>%
    finalize_long()
}

build_elsa_long <- function() {
  path <- file.path(base_dir, "ELSA", "5050stata_A06F12276F837D82374CA07FCEA89BDC3945439B39E962E7B03B39CEBABE46BC_V1", "UKDA-5050-stata", "stata", "stata13_se", "gh_elsa_h.dta")
  x <- read_dta(path)
  waves <- 1:10
  years <- c(`1` = 2002, `2` = 2004, `3` = 2006, `4` = 2008, `5` = 2010, `6` = 2012, `7` = 2014, `8` = 2016, `9` = 2018, `10` = 2021)

  bind_rows(lapply(waves, function(w) {
    tibble(
      cohort = "ELSA",
      source_file = path,
      source_type = "gateway_harmonized_wide",
      id = safe_pull_any(x, c("idauniq"), "character"),
      wave = w,
      year = safe_pull_any(x, c(paste0("r", w, "iwy"))),
      country = rep("England", nrow(x)),
      age = safe_pull_any(x, c(paste0("r", w, "agey"))),
      sex = safe_pull_any(x, c("ragender")),
      education_years = safe_pull_any(x, c("raedyrs")),
      education_level = safe_pull_any(x, c("raeducl")),
      marital_status = safe_pull_any(x, c(paste0("r", w, "mstat"))),
      partnered = safe_pull_any(x, c(paste0("r", w, "mpart"))),
      race_ethnicity = make_na(nrow(x), "character"),
      depression_score = safe_pull_any(x, c(paste0("r", w, "cesd"), paste0("r", w, "depres"))),
      loneliness_score = safe_pull_any(x, c(paste0("r", w, "lonely"), paste0("r", w, "lonela"))),
      cognition_score = rowSums(cbind(
        safe_pull_any(x, c(paste0("r", w, "imrc"))),
        safe_pull_any(x, c(paste0("r", w, "dlrc"))),
        safe_pull_any(x, c(paste0("r", w, "tr20"))),
        safe_pull_any(x, c(paste0("r", w, "orient"))),
        safe_pull_any(x, c(paste0("r", w, "numer_e")))
      ), na.rm = TRUE),
      adl_count = safe_pull_any(x, c(paste0("r", w, "adlwa"), paste0("r", w, "adlfive"))),
      iadl_count = safe_pull_any(x, c(paste0("r", w, "iadlfour"), paste0("r", w, "iadla"))),
      self_rated_health = safe_pull_any(x, c(paste0("r", w, "shlt"))),
      chronic_count = make_na(nrow(x)),
      social_participation = make_na(nrow(x)),
      social_contact = safe_pull_any(x, c(paste0("r", w, "kcnt"))),
      internet_use = make_na(nrow(x)),
      email_use = make_na(nrow(x)),
      weight_person = safe_pull_any(x, c(paste0("r", w, "wtresp"))),
      weight_household = make_na(nrow(x)),
      weight_longitudinal = safe_pull_any(x, c(paste0("r", w, "lwtresp"))),
      death_year = make_na(nrow(x)),
      death_month = make_na(nrow(x)),
      death_flag = make_na(nrow(x)),
      dropout_flag = safe_pull_any(x, c(paste0("r", w, "iwstat"))),
      inw_flag = safe_pull_any(x, c(paste0("r", w, "iwstat")))
    )
  })) %>%
    mutate(
      year = ifelse(is.na(year), unname(years[as.character(wave)]), year),
      cognition_score = ifelse(cognition_score == 0, NA_real_, cognition_score)
    ) %>%
    finalize_long()
}

builders <- list(
  SHARE = build_share_long,
  MHAS = build_mhas_long,
  ELSA = build_elsa_long,
  CHARLS = build_charls_long
)

analysis_ready <- list()
dictionary_parts <- list()

for (nm in names(builders)) {
  message("Building ", nm, " ...")
  dat <- builders[[nm]]()
  analysis_ready[[nm]] <- dat
  dictionary_parts[[nm]] <- build_dictionary(dat, nm)
  saveRDS(dat, file.path(out_dir, paste0(tolower(nm), "_analysis_ready_long_core4.rds")))
}

pooled_core4 <- bind_rows(analysis_ready) %>%
  mutate(
    cohort = factor(cohort, levels = c("ELSA", "SHARE", "CHARLS", "MHAS")),
    pooled_id = paste(cohort, id, wave, sep = "::")
  ) %>%
  relocate(pooled_id, .before = id)

dictionary <- bind_rows(dictionary_parts) %>%
  arrange(cohort, analysis_name)

saveRDS(pooled_core4, file.path(out_dir, "pooled_core4_long.rds"))
write_csv(dictionary, file.path(out_dir, "dcv_variable_dictionary_long_core4.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})

# ------------------------------------------------------------------------------
# Component: build_hrs_long.R
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
out_dir <- file.path(base_dir, "analysis_ready_core", "ready_hrs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

schema_names <- c(
  "cohort", "source_file", "source_type", "id", "wave", "year",
  "country", "age", "sex", "education_years", "education_level",
  "marital_status", "partnered", "race_ethnicity",
  "depression_score", "loneliness_score",
  "cognition_score", "adl_count", "iadl_count",
  "adl_any", "iadl_any", "self_rated_health",
  "chronic_count", "social_participation", "social_contact",
  "internet_use", "email_use", "digital_index",
  "weight_person", "weight_household", "weight_longitudinal",
  "death_year", "death_month", "death_flag", "dropout_flag",
  "inw_flag"
)

make_na <- function(n, mode = c("numeric", "character")) {
  mode <- match.arg(mode)
  if (mode == "character") rep(NA_character_, n) else rep(NA_real_, n)
}

safe_pull <- function(df, nm, type = c("numeric", "character")) {
  type <- match.arg(type)
  if (!nm %in% names(df)) return(make_na(nrow(df), type))
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

derive_binary_any <- function(x) {
  x <- na_codes_to_na(x)
  ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0))
}

derive_digital_index <- function(internet_use, email_use) {
  iu <- na_codes_to_na(internet_use)
  eu <- na_codes_to_na(email_use)
  out <- rep(NA_real_, length(iu))
  out[!is.na(iu)] <- iu[!is.na(iu)]
  out[is.na(out) & !is.na(eu)] <- eu[is.na(out) & !is.na(eu)]
  out
}

finalize_long <- function(df) {
  df %>%
    mutate(
      across(where(is.numeric), na_codes_to_na),
      adl_any = derive_binary_any(adl_count),
      iadl_any = derive_binary_any(iadl_count),
      digital_index = derive_digital_index(internet_use, email_use)
    ) %>%
    select(all_of(schema_names)) %>%
    distinct()
}

build_dictionary <- function(df, cohort_name) {
  tibble(
    cohort = cohort_name,
    analysis_name = names(df),
    non_missing_n = map_int(df, ~ sum(!is.na(.x))),
    class = map_chr(df, ~ paste(class(.x), collapse = ";")),
    example = map_chr(df, ~ {
      vals <- unique(.x[!is.na(.x)])
      if (length(vals) == 0) return("")
      paste(head(as.character(vals), 3), collapse = " | ")
    })
  )
}

hrs_h_path <- file.path(base_dir, "HRS", "_processed", "rds", "Gateway Harmonized HRS", "H_HRS_d_stata", "H_HRS_d.rds")
rand_path <- file.path(base_dir, "HRS", "_processed", "rds", "RAND HRS Longitudinal File 2022", "randhrs1992_2022v1_STATA", "randhrs1992_2022v1.rds")
h <- readRDS(hrs_h_path)
r <- readRDS(rand_path)

dedupe_on_hhidpn <- function(df) {
  df %>%
    mutate(hhidpn = as.character(hhidpn)) %>%
    distinct(hhidpn, .keep_all = TRUE)
}

r <- dedupe_on_hhidpn(r)
h <- dedupe_on_hhidpn(h)

r_join <- r %>% rename_with(~ paste0(.x, "__r"), -hhidpn)
h_join <- h %>% rename_with(~ paste0(.x, "__h"), -hhidpn)
base <- full_join(r_join, h_join, by = "hhidpn")

waves <- 1:15
years <- c(`1` = 1992, `2` = 1994, `3` = 1996, `4` = 1998, `5` = 2000, `6` = 2002, `7` = 2004,
           `8` = 2006, `9` = 2008, `10` = 2010, `11` = 2012, `12` = 2014, `13` = 2016, `14` = 2018, `15` = 2020)

hrs_long <- bind_rows(lapply(waves, function(w) {
  tibble(
    cohort = "HRS",
    source_file = paste(rand_path, hrs_h_path, sep = " + "),
    source_type = "rand_wide + gateway_harmonized_wide",
    id = base$hhidpn,
    wave = w,
    year = unname(years[as.character(w)]),
    country = make_na(nrow(base), "character"),
    age = safe_pull(base, paste0("r", w, "agey_e__r")),
    sex = safe_pull(base, "ragender__r"),
    education_years = safe_pull(base, "raeduc__r"),
    education_level = safe_pull(base, "raeducl__h"),
    marital_status = safe_pull(base, paste0("r", w, "mstat__r")),
    partnered = make_na(nrow(base)),
    race_ethnicity = safe_pull(base, "raracem__r", "character"),
    depression_score = safe_pull(base, paste0("r", w, "cesd__r")),
    loneliness_score = safe_pull(base, paste0("r", w, "ydlonely__h")),
    cognition_score = rowSums(cbind(
      safe_pull(base, paste0("r", w, "orient__h")),
      safe_pull(base, paste0("r", w, "numer__h")),
      safe_pull(base, paste0("r", w, "rxmemry__h"))
    ), na.rm = TRUE),
    adl_count = safe_pull(base, paste0("r", w, "adltot6__h")),
    iadl_count = safe_pull(base, paste0("r", w, "iadltot_h__h")),
    self_rated_health = safe_pull(base, paste0("r", w, "shlt__r")),
    chronic_count = make_na(nrow(base)),
    social_participation = make_na(nrow(base)),
    social_contact = safe_pull(base, paste0("r", w, "kcnt__h")),
    internet_use = safe_pull(base, paste0("r", w, "email__h")),
    email_use = safe_pull(base, paste0("r", w, "email__h")),
    weight_person = safe_pull(base, paste0("r", w, "nwtresp__h")),
    weight_household = make_na(nrow(base)),
    weight_longitudinal = make_na(nrow(base)),
    death_year = make_na(nrow(base)),
    death_month = make_na(nrow(base)),
    death_flag = make_na(nrow(base)),
    dropout_flag = make_na(nrow(base)),
    inw_flag = make_na(nrow(base))
  )
})) %>%
  mutate(cognition_score = ifelse(cognition_score == 0, NA_real_, cognition_score)) %>%
  finalize_long()

saveRDS(hrs_long, file.path(out_dir, "hrs_analysis_ready_long.rds"))
write_csv(build_dictionary(hrs_long, "HRS"), file.path(out_dir, "dcv_variable_dictionary_long_hrs.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)

})


