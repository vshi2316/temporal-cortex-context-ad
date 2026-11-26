# Covariate-adjusted structural models and diagnostics
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

# Longitudinal genetic associations with the fixed temporal thickness index.
pick_root <- function(ascii, native) if (file.exists(ascii)) ascii else native
project_dir <- study_path("project")
project_dir_native <- study_path("project")
adni_root <- study_path("cohorts")
result_root <- file.path(project_dir, "results")
prep_dir <- file.path(result_root, "13F_R_only_structural_analysis")
out_dir <- file.path(project_dir_native, "results", "13F_R1_final_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  path
}
read_csv <- function(path) fread(file = need(path), showProgress = FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
finite <- function(x) is.finite(num(x))
zscore <- function(x) {
  x <- num(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
has_variation <- function(x) {
  if (is.factor(x)) return(nlevels(droplevels(x)) >= 2L)
  length(unique(x[!is.na(x) & is.finite(num(x))])) >= 2L
}

d <- as.data.table(readRDS(need(file.path(prep_dir, "04_ADNI_R_only_longitudinal_structural_frame.rds"))))
d[, RID := as.character(RID)]
d[, mri_date := as.IDate(as.character(mri_date))]
d[, `:=`(
  gMAV = num(clinical_gMAV_noAPOE_z_eur),
  gAD = num(clinical_gAD_noAPOE_z_eur),
  age_mri = num(age_mri),
  qc_explicit_fail = as.logical(qc_explicit_fail)
)]

# Obtain baseline covariates from the outcome-independent PGS/clinical file.
pgs_file <- file.path(project_dir, "results", "10A_ADNI_gMAV_PGS_assembly", "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
if (file.exists(pgs_file)) {
  pgs <- fread(file = pgs_file, sep = "\t", showProgress = FALSE)
  pgs[, RID := as.character(RID)]
  pgs[, exam_date__ := as.IDate(as.character(EXAMDATE))]
  setorder(pgs, RID, exam_date__)
  first_nonmissing <- function(x, numeric = FALSE) {
    if (numeric) {
      z <- num(x); return(if (any(is.finite(z))) z[which(is.finite(z))[1L]] else NA_real_)
    }
    z <- as.character(x); z <- z[!is.na(z) & nzchar(z)]
    if (length(z)) z[1L] else NA_character_
  }
  base_cov <- pgs[, .(
    baseline_age = first_nonmissing(age_bl, TRUE),
    baseline_sex = first_nonmissing(sex),
    baseline_education = first_nonmissing(education, TRUE),
    baseline_diagnosis = first_nonmissing(dx_bl)
  ), by = RID]
  d <- merge(d, base_cov, by = "RID", all.x = TRUE, sort = FALSE)
} else {
  d[, `:=`(baseline_age = age_mri, baseline_sex = NA_character_,
           baseline_education = NA_real_, baseline_diagnosis = NA_character_)]
}
d[, `:=`(
  baseline_age = num(baseline_age),
  sex = factor(baseline_sex),
  education = num(baseline_education),
  diagnosis = factor(baseline_diagnosis),
  scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_"))
)]
pc_cols <- paste0("clinical_PC", 1:10)
for (v in pc_cols) if (v %in% names(d)) d[, (v) := num(get(v))]

# Conservative participant-level raw-image anomaly flag for sensitivity only.
flagged_file <- file.path(result_root, "13F_automated_raw_T1_QC_and_processing_manifest", "04_automatic_QC_flagged_images.csv")
flagged_rids <- character()
raw_qc_available <- file.exists(flagged_file)
if (raw_qc_available) {
  flagged <- fread(file = flagged_file, showProgress = FALSE)
  flagged_rids <- unique(as.character(flagged[cohort == "ADNI", RID]))
}
if (!raw_qc_available) message("Raw-image QC file absent; raw-QC exclusion sensitivities are skipped.")
d[, raw_auto_flagged_participant := RID %in% flagged_rids]

# One record per participant and MRI date. The primary endpoint alone defines
# the tie-break so missing secondary endpoints cannot change the primary sample.
d[, primary_complete := finite(ad_vulnerable_index)]
setorderv(d, c("RID", "mri_date", "primary_complete", "IMAGEUID"), order = c(1, 1, -1, -1), na.last = TRUE)
d <- d[, .SD[1L], by = .(RID, mri_date)]
setorder(d, RID, mri_date)
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
d[, age_mri := baseline_age + time_years]

# Standardize outcomes on the deduplicated PGS-complete frame before applying
# the spacing rule; this preserves a fixed scale across primary/sensitivity runs.
d[, `:=`(
  ad_vulnerable_z = zscore(ad_vulnerable_index),
  global_cortical_z = zscore(global_cortical_thickness),
  hippocampus_z = zscore(hip_hippocampus_icv_ratio),
  wmh_z = zscore(wmh_wmh_brain_ratio),
  mPACCdigit_z = zscore(clinical_mPACCdigit),
  mPACCtrailsB_z = zscore(clinical_mPACCtrailsB)
)]
d <- d[finite(ad_vulnerable_z) & finite(gMAV) & finite(gAD) & !is.na(mri_date)]

# Preserve all PGS-complete records after deterministic deduplication, then
# require at least two distinct MRI dates for every longitudinal model.
d_all <- copy(d)
d_no270 <- d_all[RID %in% d_all[, .(n_dates = uniqueN(mri_date)), by = RID][n_dates >= 2L, RID]]

# Apply the frozen >=270-day rule to create the primary longitudinal cohort.
elig <- d_no270[, .(
  n_dates = uniqueN(mri_date),
  second_date_270d = {
    z <- sort(unique(mri_date)); z <- z[as.integer(z - z[1L]) >= 270L]
    if (length(z)) z[1L] else as.IDate(NA)
  }
), by = RID]
elig[, primary_longitudinal := n_dates >= 2L & !is.na(second_date_270d)]
d <- d_no270[RID %in% elig[primary_longitudinal == TRUE, RID]]

covariate_names <- c("baseline_age", "sex", "education", "diagnosis", "scanner", pc_cols)
usable_covariates <- function(x, omit = character()) {
  out <- character()
  for (v in setdiff(covariate_names, omit)) {
    if (!v %in% names(x)) next
    if (is.factor(x[[v]])) {
      if (nlevels(droplevels(x[[v]])) >= 2L) out <- c(out, v)
    } else if (has_variation(x[[v]])) out <- c(out, v)
  }
  unique(out)
}

fit_mixed <- function(x, outcome, rhs, model_id) {
  vars <- unique(c(outcome, "RID", "time_years", all.vars(as.formula(paste("~", rhs)))))
  vars <- vars[vars %in% names(x)]
  a <- x[complete.cases(x[, vars, with = FALSE])]
  np <- uniqueN(a$RID); nr <- nrow(a); nrep <- a[, .N, by = RID][N >= 2L, .N]
  if (np < 100L || nr < 200L) {
    return(list(
      fit = NULL,
      coefficients = data.table(model_id = model_id, outcome = outcome, term = NA_character_, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p_value = NA_real_, status = "NOT_RUN_INSUFFICIENT_DATA", n_people = np, n_rows = nr, n_repeated = nrep, random_structure = NA_character_, singular = NA),
      diagnostics = data.table(model_id = model_id, outcome = outcome, status = "NOT_RUN_INSUFFICIENT_DATA", n_people = np, n_rows = nr, n_repeated = nrep, random_structure = NA_character_, singular = NA, AIC = NA_real_, BIC = NA_real_, logLik = NA_real_, convergence_message = "")
    ))
  }
  a <- droplevels(a)
  f_slope <- as.formula(paste0(outcome, " ~ ", rhs, " + (1 + time_years | RID)"))
  f_int <- as.formula(paste0(outcome, " ~ ", rhs, " + (1 | RID)"))
  fit <- tryCatch(lmer(f_slope, data = a, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))), error = function(e) e)
  structure <- "random_intercept_and_time_slope"
  if (inherits(fit, "error") || isSingular(fit, tol = 1e-4)) {
    fit <- tryCatch(lmer(f_int, data = a, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))), error = function(e) e)
    structure <- "random_intercept_fallback"
  }
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      coefficients = data.table(model_id = model_id, outcome = outcome, term = NA_character_, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p_value = NA_real_, status = "FAILED", n_people = np, n_rows = nr, n_repeated = nrep, random_structure = structure, singular = NA, message = conditionMessage(fit)),
      diagnostics = data.table(model_id = model_id, outcome = outcome, status = "FAILED", n_people = np, n_rows = nr, n_repeated = nrep, random_structure = structure, singular = NA, AIC = NA_real_, BIC = NA_real_, logLik = NA_real_, convergence_message = conditionMessage(fit))
    ))
  }
  sm <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(sm, c("Estimate", "Std. Error", "t value"), c("estimate", "std_error", "statistic"), skip_absent = TRUE)
  sm[, `:=`(model_id = model_id, outcome = outcome, status = "SUCCESS", p_value = 2 * pnorm(-abs(statistic)), n_people = np, n_rows = nr, n_repeated = nrep, random_structure = structure)]
  conv <- fit@optinfo$conv$lme4$messages
  list(
    fit = fit,
    coefficients = sm[, .(model_id, outcome, term, estimate, std_error, statistic, p_value, status, n_people, n_rows, n_repeated, random_structure)],
    diagnostics = data.table(model_id = model_id, outcome = outcome, status = "SUCCESS", n_people = np, n_rows = nr, n_repeated = nrep, random_structure = structure, singular = isSingular(fit, tol = 1e-4), AIC = AIC(fit), BIC = BIC(fit), logLik = as.numeric(logLik(fit)), convergence_message = if (length(conv)) paste(conv, collapse = " | ") else "")
  )
}

# -----------------------------------------------------------------------------
# ADNI primary and sensitivity models
# -----------------------------------------------------------------------------
model_list <- list(); diagnostic_list <- list(); fits <- list(); k <- 0L
run_structural <- function(x, outcome, variant, spec = "time_interaction", omit = character()) {
  covars <- usable_covariates(x, omit = omit)
  rhs <- if (spec == "main") paste(c("time_years", "gMAV", "gAD", covars), collapse = " + ") else paste(c("time_years * gMAV", "time_years * gAD", covars), collapse = " + ")
  id <- paste(variant, outcome, spec, sep = "__")
  ans <- fit_mixed(x, outcome, rhs, id)
  model_list[[length(model_list) + 1L]] <<- ans$coefficients
  diagnostic_list[[length(diagnostic_list) + 1L]] <<- ans$diagnostics
  if (variant == "PRIMARY" && outcome == "ad_vulnerable_z" && spec == "time_interaction") fits[["primary_ad_vulnerable"]] <<- ans$fit
  invisible(ans)
}

# Primary endpoint: full model plus two prespecified robustness variants.
for (variant in c("PRIMARY", "EXCLUDE_RAW_AUTO_FLAGGED", "WITHOUT_BASELINE_DIAGNOSIS")) {
  if (variant == "EXCLUDE_RAW_AUTO_FLAGGED" && !raw_qc_available) next
  x <- d
  omit <- if (variant == "WITHOUT_BASELINE_DIAGNOSIS") "diagnosis" else character()
  if (variant == "EXCLUDE_RAW_AUTO_FLAGGED") x <- x[raw_auto_flagged_participant == FALSE]
  run_structural(x, "ad_vulnerable_z", variant, "main", omit)
  run_structural(x, "ad_vulnerable_z", variant, "time_interaction", omit)
}

# Secondary structural endpoints use the primary covariate specification.
for (outcome in c("global_cortical_z", "hippocampus_z", "wmh_z")) {
  if (!outcome %in% names(d)) next
  run_structural(d, outcome, "PRIMARY", "time_interaction")
}

# No-270-day sensitivity (at least two distinct dates) for the primary endpoint.
for (variant in c("NO270", "NO270_EXCLUDE_RAW_AUTO_FLAGGED")) {
  if (variant == "NO270_EXCLUDE_RAW_AUTO_FLAGGED" && !raw_qc_available) next
  x <- d_no270
  if (variant == "NO270_EXCLUDE_RAW_AUTO_FLAGGED") x <- x[raw_auto_flagged_participant == FALSE]
  run_structural(x, "ad_vulnerable_z", variant, "time_interaction")
}

# Baseline structure -> subsequent cognitive slope bridge.
baseline_structure <- d[finite(ad_vulnerable_z), .SD[1L], by = RID][, .(RID, baseline_ad_vulnerable_z = ad_vulnerable_z)]
d <- merge(d, baseline_structure, by = "RID", all.x = TRUE, sort = FALSE)
for (outcome in c("mPACCdigit_z", "mPACCtrailsB_z")) {
  if (!outcome %in% names(d)) next
  for (variant in c("PRIMARY", "EXCLUDE_RAW_AUTO_FLAGGED")) {
    if (variant == "EXCLUDE_RAW_AUTO_FLAGGED" && !raw_qc_available) next
    x <- d[finite(get(outcome)) & finite(baseline_ad_vulnerable_z)]
    if (variant == "EXCLUDE_RAW_AUTO_FLAGGED") x <- x[raw_auto_flagged_participant == FALSE]
    covars <- usable_covariates(x)
    rhs <- paste(c("time_years * baseline_ad_vulnerable_z", covars), collapse = " + ")
    id <- paste(variant, "COGNITIVE_BRIDGE", outcome, sep = "__")
    ans <- fit_mixed(x, outcome, rhs, id)
    model_list[[length(model_list) + 1L]] <- ans$coefficients
    diagnostic_list[[length(diagnostic_list) + 1L]] <- ans$diagnostics
  }
}

coefficients <- rbindlist(model_list, fill = TRUE)
diagnostics <- rbindlist(diagnostic_list, fill = TRUE)
coefficients[, family := fifelse(grepl("^PRIMARY__ad_vulnerable_z__time_interaction$", model_id), "primary", fifelse(grepl("COGNITIVE_BRIDGE", model_id), "cognitive_bridge", "secondary_or_sensitivity"))]
coefficients[, fdr := NA_real_]
primary_ix <- coefficients$family == "primary" & coefficients$term %in% c("time_years:gMAV", "time_years:gAD")
if (sum(primary_ix) > 0L) coefficients[primary_ix, fdr := p.adjust(p_value, method = "holm")]
secondary_ix <- coefficients$family == "secondary_or_sensitivity" & coefficients$term %in% c("time_years:gMAV", "time_years:gAD") & grepl("^PRIMARY__", coefficients$model_id)
if (sum(secondary_ix) > 0L) coefficients[secondary_ix, fdr := p.adjust(p_value, method = "BH")]
cognitive_ix <- coefficients$family == "cognitive_bridge" & coefficients$term == "time_years:baseline_ad_vulnerable_z"
if (sum(cognitive_ix) > 0L) coefficients[cognitive_ix, fdr := p.adjust(p_value, method = "BH")]

person_primary <- d[, .(n_dates = uniqueN(mri_date)), by = RID]
sample_flow <- rbindlist(list(
  data.table(set = "PGS_complete_after_dedup", rows = nrow(d_all), people = uniqueN(d_all$RID)),
  data.table(set = "at_least_2_dates_no_spacing_threshold", rows = nrow(d_no270), people = uniqueN(d_no270$RID)),
  data.table(set = "primary_>=270d_spacing", rows = nrow(d), people = uniqueN(d$RID)),
  data.table(set = "primary_ad_vulnerable_complete", rows = sum(finite(d$ad_vulnerable_z)), people = uniqueN(d[finite(ad_vulnerable_z), RID])),
  data.table(set = "global_cortical_available", rows = sum(finite(d$global_cortical_z)), people = uniqueN(d[finite(global_cortical_z), RID])),
  data.table(set = "hippocampus_available", rows = sum(finite(d$hippocampus_z)), people = uniqueN(d[finite(hippocampus_z), RID])),
  data.table(set = "WMH_available", rows = sum(finite(d$wmh_z)), people = uniqueN(d[finite(wmh_z), RID])),
  data.table(set = "mPACCdigit_available", rows = sum(finite(d$mPACCdigit_z)), people = uniqueN(d[finite(mPACCdigit_z), RID])),
  data.table(set = "mPACCtrailsB_available", rows = sum(finite(d$mPACCtrailsB_z)), people = uniqueN(d[finite(mPACCtrailsB_z), RID]))
), fill = TRUE)

# Diagnostics for the primary fit.
primary_fit <- fits[["primary_ad_vulnerable"]]
if (!is.null(primary_fit)) {
  diag_dt <- data.table(fitted = fitted(primary_fit), residual = residuals(primary_fit), standardized_residual = as.numeric(scale(residuals(primary_fit))))
  fwrite(diag_dt, file.path(out_dir, "07_primary_residuals.csv"), bom = TRUE)
}

fwrite(coefficients, file.path(out_dir, "01_model_coefficients.csv"), bom = TRUE)
fwrite(diagnostics, file.path(out_dir, "02_model_diagnostics.csv"), bom = TRUE)
fwrite(sample_flow, file.path(out_dir, "03_sample_flow.csv"), bom = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)

message("13F-R1 final analysis completed: ", out_dir)
