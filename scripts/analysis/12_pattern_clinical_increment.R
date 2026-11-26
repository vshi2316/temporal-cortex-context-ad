# Lean/full adjustment and complete-case pattern comparisons
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

project <- study_path("project")
main_results <- study_path("project", "results")
aligned_dir <- study_path("project", "results/13H1_aligned_regional_slope_map_release")
clinical_file <- file.path(main_results, "13F_R1_clinical_closure_confirmatory", "10_two_MRI_landmark_frame.csv")
key_file <- study_path("project", "input/analysis_definitions/scan_keys.csv")
map_file <- file.path(REPO_ROOT, "definitions/bilateral_region_mapping.csv")
aligned_file <- file.path(aligned_dir, "01_aligned_34_region_contrasts.csv")
fs_file <- study_path("cohorts", "ADNI/data/MRI/MRI_FreeSurfer.csv")
required <- c(clinical_file, key_file, map_file, aligned_file, fs_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing frozen inputs: ", paste(missing, collapse = "; "))

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
zscore <- function(x) {
  x <- num(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

temporal <- c("entorhinal", "fusiform", "inferiortemporal", "middletemporal",
              "temporalpole", "parahippocampal")
aligned <- fread(aligned_file)
pairs <- fread(map_file)
if (nrow(aligned) != 34L || nrow(pairs) != 34L ||
    !setequal(aligned$region, pairs$region) ||
    !all(aligned$status == "SUCCESS"))
  stop("The frozen 34-region genetic contrast map is incomplete.")
regions <- pairs$region
template <- pairs[, .(region)]
template[, genetic_slope_difference := aligned$estimate[match(region, aligned$region)]]
template[, weight := -genetic_slope_difference]
template[, weight := weight - mean(weight)]
template[, weight := weight / sqrt(sum(weight^2))]

keys <- fread(key_file)[in_main_model == TRUE,
                        .(RID = as.character(RID), mri_date = as.IDate(mri_date),
                          main_imageuid = as.character(main_imageuid))]
if (nrow(keys) != 4659L || uniqueN(keys$RID) != 695L)
  stop("The frozen 13K scan frame changed.")
fields <- unique(c(pairs$left_field, pairs$right_field))
fs <- fread(fs_file, select = c("RID", "EXAMDATE", "OVERALLQC", "IMAGEUID", fields),
            showProgress = FALSE)
fs[, `:=`(RID = as.character(RID), mri_date = as.IDate(as.character(EXAMDATE)))]
fs <- fs[!is.na(RID) & !is.na(mri_date)]
fs <- fs[!grepl("fail|exclude|poor|unusable", tolower(trimws(as.character(OVERALLQC))))]
setorderv(fs, c("RID", "mri_date", "IMAGEUID"), order = c(1L, 1L, -1L), na.last = TRUE)
fs <- fs[, .SD[1L], by = .(RID, mri_date)]
fs <- merge(fs, keys, by = c("RID", "mri_date"), sort = FALSE)
if (nrow(fs) != nrow(keys) || any(as.character(fs$IMAGEUID) != fs$main_imageuid))
  stop("FreeSurfer scan keys or Image IDs differ from 13K.")

for (i in seq_len(nrow(pairs))) {
  left <- pairs$left_field[i]
  right <- pairs$right_field[i]
  fs[, (regions[i]) := rowMeans(.SD, na.rm = TRUE), .SDcols = c(left, right)]
  strict_field <- paste0("strict__", regions[i])
  fs[, (strict_field) := rowMeans(.SD, na.rm = FALSE), .SDcols = c(left, right)]
}
setorder(fs, RID, mri_date)
first_scan <- fs[, .SD[1L], by = RID]
template[, baseline_sd_mm := vapply(region, function(v) sd(first_scan[[v]], na.rm = TRUE), numeric(1))]
if (any(!is.finite(template$baseline_sd_mm) | template$baseline_sd_mm <= 0))
  stop("At least one regional baseline thickness SD is invalid.")
fs[, temporal_thickness_mm := rowMeans(.SD, na.rm = FALSE), .SDcols = temporal]
strict_regions <- paste0("strict__", regions)
scan <- fs[, c("RID", "mri_date", "temporal_thickness_mm", "IMAGEUID", regions,
               strict_regions), with = FALSE]
setnames(scan, "IMAGEUID", "imageuid")
index <- copy(scan)
setnames(index, c("mri_date", "temporal_thickness_mm", "imageuid", regions, strict_regions),
         c("index_mri_date", "index_thickness_mm", "index_imageuid", paste0("index__", regions),
           paste0("index_strict__", regions)))
landmark <- copy(scan)
setnames(landmark, c("mri_date", "temporal_thickness_mm", "imageuid", regions, strict_regions),
         c("landmark_mri_date", "landmark_thickness_mm", "landmark_imageuid", paste0("landmark__", regions),
           paste0("landmark_strict__", regions)))

clinical <- fread(clinical_file)
if (!"scanner" %in% names(clinical)) {
  if ("site" %in% names(clinical)) {
    clinical[, scanner := site]
  } else {
    stop("Clinical frame needs scanner or site.")
  }
}
clinical[, `:=`(RID = as.character(RID), index_mri_date = as.IDate(index_mri_date),
                landmark_mri_date = as.IDate(landmark_mri_date), baseline_dx_date = as.IDate(baseline_dx_date),
                event_date = as.IDate(event_date), censor_date = as.IDate(censor_date))]
if (nrow(clinical) != 616L || uniqueN(clinical$RID) != 616L)
  stop("The frozen two-MRI clinical landmark frame changed.")
x <- merge(clinical, index, by = c("RID", "index_mri_date"), all.x = TRUE, sort = FALSE)
x <- merge(x, landmark, by = c("RID", "landmark_mri_date"), all.x = TRUE, sort = FALSE)
setDT(x)
if (nrow(x) != nrow(clinical)) stop("Clinical-to-MRI linkage is not one-to-one.")
x[, regional_thinning_mm_per_year := (index_thickness_mm - landmark_thickness_mm) /
    (as.numeric(landmark_mri_date - index_mri_date) / 365.25)]
equivalence_error <- max(abs(x$regional_thinning_mm_per_year - num(x$atrophy_rate)), na.rm = TRUE)
if (!is.finite(equivalence_error) || equivalence_error > 1e-8)
  stop("The six-region mean no longer reproduces the original atrophy rate.")
interval <- as.numeric(x$landmark_mri_date - x$index_mri_date) / 365.25
regional_rates <- sapply(seq_along(regions), function(i) {
  (num(x[[paste0("index__", regions[i])]]) - num(x[[paste0("landmark__", regions[i])]])) /
    interval / template$baseline_sd_mm[i]
})
x[, spatial_pattern_score := as.vector(regional_rates %*% template$weight)]
strict_regional_rates <- sapply(seq_along(regions), function(i) {
  (num(x[[paste0("index_strict__", regions[i])]]) -
     num(x[[paste0("landmark_strict__", regions[i])]])) /
    interval / template$baseline_sd_mm[i]
})
x[, spatial_pattern_score_strict_bilateral := as.vector(strict_regional_rates %*% template$weight)]
x[, strict_bilateral_complete := apply(is.finite(strict_regional_rates), 1L, all)]
x[, risk_start_strict := landmark_mri_date]
x[!is.na(baseline_dx_date) & baseline_dx_date > risk_start_strict, risk_start_strict := baseline_dx_date]
x[, `:=`(event_before_risk_start = !is.na(event_date) & event_date <= risk_start_strict,
         censor_before_risk_start = !is.na(censor_date) & censor_date <= risk_start_strict,
         event_strict = as.integer(!is.na(event_date)),
         followup_strict = (fifelse(!is.na(event_date), as.numeric(event_date), as.numeric(censor_date)) -
                            as.numeric(risk_start_strict)) / 365.25)]
flow <- data.table(stage = "frozen_two_MRI_landmark", people = nrow(x), events = sum(x$event_strict))
add_flow <- function(label, dat) data.table(stage = label, people = nrow(dat), events = sum(dat$event_strict))
linked <- x[!is.na(index_imageuid) & !is.na(landmark_imageuid)]
setDT(linked)
if (!nrow(linked)) stop("No clinical rows linked to both FreeSurfer scans.")
flow <- rbind(flow, add_flow("both_FreeSurfer_scans_matched", linked))
linked <- linked[is.finite(spatial_pattern_score)]
if (!nrow(linked)) stop("No linked rows have a finite 34-region spatial score.")
flow <- rbind(flow, add_flow("all_34_pattern_regions_complete", linked))
strict_linked <- x[strict_bilateral_complete == TRUE]
flow <- rbind(flow, add_flow("all_34_pattern_regions_strict_bilateral", strict_linked))
linked <- linked[!is.na(baseline_dx_date) & event_before_risk_start == FALSE &
                   censor_before_risk_start == FALSE & is.finite(followup_strict) & followup_strict > 0]
flow <- rbind(flow, add_flow("strict_post_landmark_risk_start", linked))

pc <- paste0("pc", 1:10, "_z")
raw_model_vars <- c("gAD", "gMAV", "baseline_structure_vulnerability_z", "atrophy_rate_z",
                    "baseline_mPACC", "baseline_CDRSB", "baseline_FAQ", "age_index_z",
                    "education_z", "apoe4_z", "sex", "scanner", "baseline_dx", pc)
missing_raw <- setdiff(raw_model_vars, names(linked))
if (length(missing_raw))
  stop("Clinical frame is missing required variables: ", paste(missing_raw, collapse = ", "))
linked[, spatial_pattern_z := zscore(spatial_pattern_score)]
linked[, event := as.integer(event_strict)]
linked[, followup_years := num(followup_strict)]
linked[, gAD := num(gAD)]
linked[, gMAV := num(gMAV)]
linked[, baseline_structure_vulnerability_z := num(baseline_structure_vulnerability_z)]
linked[, atrophy_rate_z := num(atrophy_rate_z)]
linked[, baseline_mPACC_z := zscore(baseline_mPACC)]
linked[, baseline_CDRSB_z := zscore(baseline_CDRSB)]
linked[, baseline_FAQ_z := zscore(baseline_FAQ)]
linked[, age_index_z := num(age_index_z)]
linked[, education_z := num(education_z)]
linked[, apoe4_z := num(apoe4_z)]
linked[, sex := factor(sex)]
linked[, scanner := factor(scanner)]
linked[, baseline_dx := factor(baseline_dx)]
for (v in pc) linked[, (v) := num(get(v))]
minimal <- c("age_index_z", "sex", "education_z", "apoe4_z", "scanner", "baseline_dx", pc,
             "gAD", "gMAV", "baseline_mPACC_z", "baseline_CDRSB_z", "baseline_FAQ_z")
baseline <- c(minimal, "baseline_structure_vulnerability_z")
temporal_model <- c(baseline, "atrophy_rate_z")
spatial_model <- c(temporal_model, "spatial_pattern_z")
all_vars <- c("followup_years", "event", spatial_model)
complete_rows <- function(dat, vars) {
  tab <- as.data.frame(dat)
  missing_vars <- setdiff(vars, colnames(tab))
  if (length(missing_vars))
    stop("Missing model variables: ", paste(missing_vars, collapse = ", "),
         "; input class=", paste(class(dat), collapse = "/"),
         "; input dimensions=", paste(dim(tab), collapse = "x"))
  keep <- tab[, vars, drop = FALSE]
  stats::complete.cases(keep)
}
keep_common <- complete_rows(linked, all_vars)
dat_common <- droplevels(linked[keep_common, ])
flow <- rbind(flow, add_flow("common_complete_model_frame", dat_common))
lean_terms <- c("age_index_z", "sex", "education_z", "apoe4_z", "scanner",
                "baseline_dx", "gAD", "gMAV")
keep_lean <- complete_rows(linked, c("followup_years", "event", lean_terms,
                                    "spatial_pattern_z"))
lean_dat <- droplevels(linked[keep_lean, ])
flow <- rbind(flow, add_flow("lean_complete_model_frame", lean_dat))

out <- file.path(project, "results", paste0("13H2_clinical_increment_audit_", "release"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(flow, file.path(out, "01_sample_flow.csv"), bom = TRUE)

bilateral <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  l <- pairs$left_field[i]
  r <- pairs$right_field[i]
  data.table(region = regions[i], left_nonmissing = sum(is.finite(num(fs[[l]])), na.rm = TRUE),
             right_nonmissing = sum(is.finite(num(fs[[r]])), na.rm = TRUE),
             both_nonmissing = sum(is.finite(num(fs[[l]])) & is.finite(num(fs[[r]])), na.rm = TRUE),
             left_only = sum(is.finite(num(fs[[l]])) & !is.finite(num(fs[[r]])), na.rm = TRUE),
             right_only = sum(!is.finite(num(fs[[l]])) & is.finite(num(fs[[r]])), na.rm = TRUE))
}))
fwrite(bilateral, file.path(out, "02_bilateral_region_completeness.csv"), bom = TRUE)
missing_table <- data.table(variable = spatial_model,
                            missing_people = vapply(spatial_model, function(v) sum(is.na(linked[[v]])), integer(1)),
                            missing_events = vapply(spatial_model, function(v) sum(is.na(linked[[v]]) &
                                                                               linked$event == 1L), integer(1)))
fwrite(missing_table, file.path(out, "02b_model_covariate_missingness.csv"), bom = TRUE)

fit_one <- function(dat, terms) {
  form <- as.formula(paste("Surv(followup_years, event) ~", paste(terms, collapse = " + ")))
  fit <- coxph(form, data = dat, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  bad <- names(coef(fit))[!is.finite(coef(fit))]
  if (length(bad)) stop("Non-finite coefficients in model: ", paste(bad, collapse = ", "))
  fit
}

# Directly isolate the two variables responsible for the 584-person reduction.
missing_sensitive <- c("baseline_CDRSB_z", "baseline_FAQ_z")
reduced_temporal <- setdiff(temporal_model, missing_sensitive)
reduced_spatial <- c(reduced_temporal, "spatial_pattern_z")
reduced_vars <- c("followup_years", "event", reduced_spatial)
keep_reduced <- complete_rows(linked, reduced_vars)
dat_reduced_594 <- droplevels(linked[keep_reduced, ])
keep_reduced_584 <- keep_reduced & complete_rows(linked, missing_sensitive)
dat_reduced_584 <- droplevels(linked[keep_reduced_584, ])
fit_reduced_594_base <- fit_one(dat_reduced_594, reduced_temporal)
fit_reduced_594_spatial <- fit_one(dat_reduced_594, reduced_spatial)
fit_reduced_584_base <- fit_one(dat_reduced_584, reduced_temporal)
fit_reduced_584_spatial <- fit_one(dat_reduced_584, reduced_spatial)
reduced_beta <- c(coef(fit_reduced_594_spatial)[["spatial_pattern_z"]],
                  coef(fit_reduced_584_spatial)[["spatial_pattern_z"]])
reduced_se <- c(sqrt(vcov(fit_reduced_594_spatial)["spatial_pattern_z", "spatial_pattern_z"]),
                sqrt(vcov(fit_reduced_584_spatial)["spatial_pattern_z", "spatial_pattern_z"]))
missingness_audit <- data.table(
  analysis = c("reduced_model_594", "reduced_model_584"),
  n = c(nrow(dat_reduced_594), nrow(dat_reduced_584)),
  events = c(sum(dat_reduced_594$event), sum(dat_reduced_584$event)),
  spatial_hr = exp(reduced_beta),
  spatial_lower_95 = exp(reduced_beta - 1.96 * reduced_se),
  spatial_upper_95 = exp(reduced_beta + 1.96 * reduced_se),
  spatial_wald_p = 2 * pnorm(-abs(reduced_beta / reduced_se)),
  spatial_lr_p = pchisq(c(
    2 * (as.numeric(logLik(fit_reduced_594_spatial)) - as.numeric(logLik(fit_reduced_594_base))),
    2 * (as.numeric(logLik(fit_reduced_584_spatial)) - as.numeric(logLik(fit_reduced_584_base)))
  ), df = 1, lower.tail = FALSE)
)
fwrite(missingness_audit, file.path(out, "03c_missing_covariate_case_sensitivity.csv"), bom = TRUE)
models <- list(minimal_clinical = minimal, plus_baseline_structure = baseline,
               plus_temporal_mean = temporal_model, plus_spatial_pattern = spatial_model)
fits <- lapply(models, function(v) fit_one(dat_common, v))
names(fits) <- names(models)
model_table <- rbindlist(lapply(seq_along(fits), function(i) {
  fit <- fits[[i]]
  term_name <- if (i == 1L) NA_character_ else names(models)[i]
  data.table(model = names(fits)[i], n = nrow(dat_common), events = sum(dat_common$event),
             terms = length(models[[i]]), logLik = as.numeric(logLik(fit)),
             added_term = term_name)
}))
model_table[, `:=`(lr_chisq = c(NA_real_, 2 * diff(logLik)),
                   lr_p_value = c(NA_real_, pchisq(2 * diff(logLik), df = 1, lower.tail = FALSE)))]
fwrite(model_table, file.path(out, "03_model_comparison_common_frame.csv"), bom = TRUE)
full_fit <- fits[["plus_spatial_pattern"]]
full_beta <- coef(full_fit)[["spatial_pattern_z"]]
full_se <- sqrt(vcov(full_fit)["spatial_pattern_z", "spatial_pattern_z"])
lean_base <- fit_one(lean_dat, lean_terms)
lean_aug <- fit_one(lean_dat, c(lean_terms, "spatial_pattern_z"))
lean_beta <- coef(lean_aug)[["spatial_pattern_z"]]
lean_se <- sqrt(vcov(lean_aug)["spatial_pattern_z", "spatial_pattern_z"])
increment <- data.table(
  analysis = c("strong_common_frame", "lean_available_frame"),
  n = c(nrow(dat_common), nrow(lean_dat)),
  events = c(sum(dat_common$event), sum(lean_dat$event)),
  hr_per_sd = exp(c(full_beta, lean_beta)),
  lower_95 = exp(c(full_beta - 1.96 * full_se, lean_beta - 1.96 * lean_se)),
  upper_95 = exp(c(full_beta + 1.96 * full_se, lean_beta + 1.96 * lean_se)),
  wald_p = 2 * pnorm(-abs(c(full_beta / full_se, lean_beta / lean_se))),
  lr_p_value = c(model_table[model == "plus_spatial_pattern", lr_p_value],
                 pchisq(2 * (as.numeric(logLik(lean_aug)) - as.numeric(logLik(lean_base))),
                        df = 1, lower.tail = FALSE))
)
fwrite(increment, file.path(out, "03a_complete_case_sensitivity.csv"), bom = TRUE)

strict_risk <- x[is.finite(spatial_pattern_score_strict_bilateral) &
                   !is.na(baseline_dx_date) & event_before_risk_start == FALSE &
                   censor_before_risk_start == FALSE & is.finite(followup_strict) & followup_strict > 0]
setDT(strict_risk)
strict_risk[, spatial_pattern_z := zscore(spatial_pattern_score)]
strict_risk[, spatial_pattern_strict_z := zscore(spatial_pattern_score_strict_bilateral)]
strict_risk[, event := as.integer(event_strict)]
strict_risk[, followup_years := num(followup_strict)]
strict_risk[, gAD := num(gAD)]
strict_risk[, gMAV := num(gMAV)]
strict_risk[, baseline_structure_vulnerability_z := num(baseline_structure_vulnerability_z)]
strict_risk[, atrophy_rate_z := num(atrophy_rate_z)]
strict_risk[, baseline_mPACC_z := zscore(baseline_mPACC)]
strict_risk[, baseline_CDRSB_z := zscore(baseline_CDRSB)]
strict_risk[, baseline_FAQ_z := zscore(baseline_FAQ)]
strict_risk[, age_index_z := num(age_index_z)]
strict_risk[, education_z := num(education_z)]
strict_risk[, apoe4_z := num(apoe4_z)]
strict_risk[, sex := factor(sex)]
strict_risk[, scanner := factor(scanner)]
strict_risk[, baseline_dx := factor(baseline_dx)]
for (v in pc) strict_risk[, (v) := num(get(v))]
keep_strict <- complete_rows(strict_risk,
                             c("followup_years", "event", spatial_model,
                               "spatial_pattern_strict_z"))
strict_dat <- droplevels(strict_risk[keep_strict, ])
strict_base <- fit_one(strict_dat, temporal_model)
strict_current <- fit_one(strict_dat, spatial_model)
strict_aug <- fit_one(strict_dat, c(temporal_model, "spatial_pattern_strict_z"))
strict_beta <- coef(strict_aug)[["spatial_pattern_strict_z"]]
strict_se <- sqrt(vcov(strict_aug)["spatial_pattern_strict_z", "spatial_pattern_strict_z"])
strict_lr <- 2 * (as.numeric(logLik(strict_aug)) - as.numeric(logLik(strict_base)))
current_beta <- coef(strict_current)[["spatial_pattern_z"]]
current_se <- sqrt(vcov(strict_current)["spatial_pattern_z", "spatial_pattern_z"])
current_lr <- 2 * (as.numeric(logLik(strict_current)) - as.numeric(logLik(strict_base)))
strict_result <- data.table(
  analysis = c("current_score_on_strict_subset", "strict_bilateral_34_region_score"),
  n = nrow(strict_dat), events = sum(strict_dat$event),
  hr_per_sd = exp(c(current_beta, strict_beta)),
  lower_95 = exp(c(current_beta - 1.96 * current_se, strict_beta - 1.96 * strict_se)),
  upper_95 = exp(c(current_beta + 1.96 * current_se, strict_beta + 1.96 * strict_se)),
  wald_p = 2 * pnorm(-abs(c(current_beta / current_se, strict_beta / strict_se))),
  lr_chisq = c(current_lr, strict_lr),
  lr_p_value = pchisq(c(current_lr, strict_lr), df = 1, lower.tail = FALSE),
  score_correlation_with_current = cor(strict_dat$spatial_pattern_z, strict_dat$spatial_pattern_strict_z)
)
fwrite(strict_result, file.path(out, "03b_strict_bilateral_sensitivity.csv"), bom = TRUE)

cor_table <- data.table(
  comparison = c("spatial_vs_atrophy_rate", "spatial_vs_baseline_structure", "spatial_vs_gAD", "spatial_vs_gMAV"),
  estimate = c(cor(dat_common$spatial_pattern_z, dat_common$atrophy_rate_z),
               cor(dat_common$spatial_pattern_z, dat_common$baseline_structure_vulnerability_z),
               cor(dat_common$spatial_pattern_z, dat_common$gAD),
               cor(dat_common$spatial_pattern_z, dat_common$gMAV)),
  n = nrow(dat_common))
fwrite(cor_table, file.path(out, "04_score_correlation_audit.csv"), bom = TRUE)
fwrite(dat_common[, .(RID, followup_years, event, spatial_pattern_z, atrophy_rate_z,
                      baseline_structure_vulnerability_z, gAD, gMAV)],
       file.path(out, "05_common_model_frame.csv"), bom = TRUE)

status <- c(
  "STATUS=CLINICAL_INCREMENT_AUDIT_COMPLETE",
  paste0("LINKED_PEOPLE=", flow[stage == "both_FreeSurfer_scans_matched", people]),
  paste0("RISK_PEOPLE=", flow[stage == "strict_post_landmark_risk_start", people]),
  paste0("COMMON_COMPLETE_PEOPLE=", nrow(dat_common)),
  paste0("COMMON_COMPLETE_EVENTS=", sum(dat_common$event)),
  paste0("LEAN_COMPLETE_PEOPLE=", nrow(lean_dat)),
  paste0("LEAN_COMPLETE_EVENTS=", sum(lean_dat$event)),
  paste0("SPATIAL_INCREMENT_LR_P=", signif(model_table[model == "plus_spatial_pattern", lr_p_value], 5)),
  paste0("TEMPORAL_INCREMENT_LR_P=", signif(model_table[model == "plus_temporal_mean", lr_p_value], 5)),
  paste0("STRICT_BILATERAL_PEOPLE=", nrow(strict_dat)),
  paste0("STRICT_BILATERAL_EVENTS=", sum(strict_dat$event)),
  paste0("STRICT_BILATERAL_LR_P=", signif(strict_result[analysis == "strict_bilateral_34_region_score", lr_p_value], 5)),
  paste0("REDUCED_594_SPATIAL_LR_P=", signif(missingness_audit[analysis == "reduced_model_594", spatial_lr_p], 5)),
  paste0("REDUCED_584_SPATIAL_LR_P=", signif(missingness_audit[analysis == "reduced_model_584", spatial_lr_p], 5)),
  paste0("SIX_REGION_MEAN_DUPLICATE_MAX_ERROR=", signif(equivalence_error, 5)),
  "NOTE=All four Cox models use the identical common complete-case frame.",
  "NOTE=The comparison separates clinical, baseline-structure, temporal-mean, and 34-region spatial increments.",

  paste0("OUTPUT_DIR=", out)
)
writeLines(status, file.path(out, "00_status.txt"), useBytes = TRUE)
message(paste(status, collapse = "\n"))
