# Common and diagnosis-specific time trends with technical specifications
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

project <- study_path("project")
main_results <- study_path("project", "results")
structural_file <- file.path(
  main_results, "13F_R_only_structural_analysis",
  "04_ADNI_R_only_longitudinal_structural_frame.rds"
)
pgs_file <- file.path(
  main_results, "10A_ADNI_gMAV_PGS_assembly",
  "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"
)
fs_file <- study_path("cohorts", "ADNI/data/MRI/MRI_FreeSurfer.csv")

out_dir <- file.path(project, "results", paste0(
  "13H44_ADNI_diagnosis_time_specification_audit_",
  "release"
))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  path
}
invisible(lapply(c(structural_file, pgs_file), need))

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
first_num <- function(x) {
  z <- num(x)
  if (any(is.finite(z))) z[which(is.finite(z))[1L]] else NA_real_
}
first_chr <- function(x) {
  z <- trimws(as.character(x))
  z <- z[!is.na(z) & nzchar(z)]
  if (length(z)) z[1L] else NA_character_
}
has_variation <- function(x) {
  if (is.factor(x)) return(nlevels(droplevels(x)) >= 2L)
  if (is.character(x)) {
    z <- trimws(x)
    z <- z[!is.na(z) & nzchar(z)]
    return(length(unique(z)) >= 2L)
  }
  length(unique(x[is.finite(num(x))])) >= 2L
}
write_csv <- function(x, name) {
  fwrite(as.data.table(x), file.path(out_dir, name), bom = TRUE, na = "NA")
}

# Rebuild the frozen 13K analysis frame without changing its participant/date rule.
d <- as.data.table(readRDS(structural_file))
pgs <- fread(pgs_file, sep = "\t", showProgress = FALSE)
d[, `:=`(
  RID = as.character(RID),
  mri_date = as.IDate(as.character(mri_date)),
  gMAV = num(clinical_gMAV_noAPOE_z_eur),
  gAD = num(clinical_gAD_noAPOE_z_eur)
)]
pgs[, `:=`(
  RID = as.character(RID),
  exam_date__ = as.IDate(as.character(EXAMDATE))
)]
setorder(pgs, RID, exam_date__)

base_cov <- pgs[, .(
  baseline_age = first_num(age_bl),
  baseline_sex = first_chr(sex),
  baseline_education = first_num(education),
  baseline_diagnosis = first_chr(dx_bl)
), by = RID]
d <- merge(d, base_cov, by = "RID", all.x = TRUE, sort = FALSE)
d[, `:=`(
  baseline_age = num(baseline_age),
  sex = factor(baseline_sex),
  education = num(baseline_education),
  diagnosis = factor(baseline_diagnosis),
  scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_")),
  primary_complete = is.finite(num(ad_vulnerable_index))
)]

pc_cols <- paste0("clinical_PC", 1:10)
for (v in pc_cols) d[, (v) := num(get(v))]
setorderv(d, c("RID", "mri_date", "primary_complete", "IMAGEUID"),
         c(1L, 1L, -1L, -1L), na.last = TRUE)
d <- d[, .SD[1L], by = .(RID, mri_date)]
setorder(d, RID, mri_date)
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
d[, ad_vulnerable_z := {
  z <- num(ad_vulnerable_index)
  (z - mean(z, na.rm = TRUE)) / sd(z, na.rm = TRUE)
}]

d <- d[
  is.finite(ad_vulnerable_z) & is.finite(gMAV) & is.finite(gAD) &
    !is.na(mri_date)
]
elig <- d[, {
  dates <- sort(unique(mri_date))
  later <- dates[as.integer(dates - dates[1L]) >= 270L]
  list(n_dates = length(dates), has_270d_second = length(later) > 0L)
}, by = RID]
d <- d[RID %in% elig[n_dates >= 2L & has_270d_second == TRUE, RID]]

candidate_covars <- c(
  "baseline_age", "sex", "education", "diagnosis", "scanner", pc_cols
)
covars <- candidate_covars[vapply(candidate_covars, function(v) {
  v %in% names(d) && has_variation(d[[v]])
}, logical(1))]
model_vars <- unique(c(
  "ad_vulnerable_z", "RID", "time_years", "gMAV", "gAD", covars
))
d <- droplevels(d[complete.cases(d[, ..model_vars])])
if (uniqueN(d$RID) != 695L || nrow(d) != 4659L) {
  stop("Frozen 13K frame is not 695 people / 4659 MRI rows.", call. = FALSE)
}

# MRI phase is added for the technical sensitivity models.
d[, IMAGEUID := as.character(IMAGEUID)]
phase_available <- FALSE
if (file.exists(fs_file) && "IMAGEUID" %in% names(d)) {
  fs <- fread(fs_file, select = c(
    "RID", "IMAGEUID", "EXAMDATE", "PHASE", "FIELD_STRENGTH", "FSVER"
  ), showProgress = FALSE)
  fs[, `:=`(
    RID = as.character(RID), IMAGEUID = as.character(IMAGEUID),
    mri_date = as.IDate(as.character(EXAMDATE)),
    phase = as.character(PHASE)
  )]
  fs <- unique(fs[, .(RID, IMAGEUID, mri_date, phase)])
  if (anyDuplicated(fs, by = c("RID", "IMAGEUID", "mri_date"))) {
    stop("MRI_FreeSurfer contains conflicting phase rows.", call. = FALSE)
  }
  phase_join <- merge(
    d[, .(RID, IMAGEUID, mri_date)], fs,
    by = c("RID", "IMAGEUID", "mri_date"), all.x = TRUE, sort = FALSE
  )
  if (nrow(phase_join) != nrow(d)) {
    stop("Phase join changed the frozen MRI row count.", call. = FALSE)
  }
  d <- merge(d, phase_join[, .(RID, IMAGEUID, mri_date, phase)],
             by = c("RID", "IMAGEUID", "mri_date"), all.x = TRUE,
             sort = FALSE)
  phase_available <- any(!is.na(d$phase)) && has_variation(d$phase)
}
if (!"phase" %in% names(d)) d[, phase := NA_character_]
setorder(d, RID, mri_date)

if (uniqueN(d$RID) != 695L || nrow(d) != 4659L) {
  stop("Frame changed after phase audit.", call. = FALSE)
}

write_csv(d[, .(RID, mri_date, time_years, diagnosis, gAD, gMAV,
                FIELD_STRENGTH, FSVER, scanner, phase)],
          "01_analysis_frame_audit.csv")
write_csv(d[, .(people = uniqueN(RID), scans = .N), by = diagnosis],
          "02_diagnosis_distribution.csv")
write_csv(d[, .(people = uniqueN(RID), scans = .N),
            by = .(diagnosis, FIELD_STRENGTH, FSVER, scanner)],
          "03_technical_distribution.csv")

control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

fit_one <- function(dat, model_id, added_terms, technical = FALSE) {
  dat <- droplevels(dat)
  if (technical) {
    covars_model <- c("baseline_age", "sex", "education", "diagnosis",
                     "FIELD_STRENGTH", "FSVER", pc_cols)
    if (phase_available && "phase" %in% names(dat) && has_variation(dat$phase)) {
      covars_model <- c(covars_model, "phase")
    }
  } else {
    covars_model <- c("baseline_age", "sex", "education", "diagnosis",
                     "scanner", pc_cols)
  }
  covars_model <- covars_model[vapply(covars_model,
    function(v) has_variation(dat[[v]]), logical(1))]
  rhs <- c("time_years * gMAV", "time_years * gAD", covars_model)
  if ("diagnosis_time" %in% added_terms) {
    rhs <- c(rhs, "time_years * diagnosis")
  }
  formula_for <- function(random) as.formula(paste0(
    "ad_vulnerable_z ~ ", paste(rhs, collapse = " + "),
    if (random == "slope") " + (1 + time_years | RID)" else " + (1 | RID)"
  ))
  fit <- tryCatch(
    lmer(formula_for("slope"), data = dat, REML = FALSE, control = control),
    error = function(e) e
  )
  random_structure <- "random_intercept_and_time_slope"
  if (inherits(fit, "error") || isSingular(fit, tol = 1e-4)) {
    fit <- tryCatch(
      lmer(formula_for("intercept"), data = dat, REML = FALSE, control = control),
      error = function(e) e
    )
    random_structure <- "random_intercept_fallback"
  }
  if (inherits(fit, "error")) {
    stop(model_id, " failed: ", conditionMessage(fit), call. = FALSE)
  }
  b <- fixef(fit)
  V <- as.matrix(vcov(fit))
  required <- c("time_years:gAD", "time_years:gMAV")
  if (!all(required %in% names(b))) {
    stop(model_id, " lacks the two genetic slope terms.", call. = FALSE)
  }
  L <- setNames(rep(0, length(b)), names(b))
  L["time_years:gAD"] <- 1
  L["time_years:gMAV"] <- -1
  est <- unname(sum(L * b))
  se <- sqrt(as.numeric(t(L) %*% V %*% L))
  contrast <- data.table(
    model_id = model_id,
    added_terms = paste(added_terms, collapse = ";"),
    technical_covariates = paste(covars_model, collapse = ";"),
    estimate = est,
    standard_error = se,
    conf_low = est - 1.96 * se,
    conf_high = est + 1.96 * se,
    z_value = est / se,
    p_value = 2 * pnorm(-abs(est / se)),
    n_people = uniqueN(dat$RID), n_scans = nrow(dat),
    random_structure = random_structure
  )
  co <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(co, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
           c("estimate", "standard_error", "statistic", "p_value"),
           skip_absent = TRUE)
  co[, `:=`(model_id = model_id, random_structure = random_structure,
            n_people = uniqueN(dat$RID), n_scans = nrow(dat))]
  list(contrast = contrast, coefficients = co, fit = fit)
}

fits <- list(original_13K = fit_one(
  d, "original_13K", "none", technical = FALSE
))
frozen <- fread(file.path(
  main_results, "13K_core_claim_completion",
  "01_gAD_vs_gMAV_direct_slope_contrast.csv"
), showProgress = FALSE)
frozen_est <- frozen$estimate_difference_SD_per_year[1L]
frozen_se <- frozen$standard_error[1L]
if (abs(fits$original_13K$contrast$estimate - frozen_est) > 1e-6 ||
    abs(fits$original_13K$contrast$standard_error - frozen_se) > 1e-6) {
  stop("Original model does not reproduce 13K: current=",
       fits$original_13K$contrast$estimate, ", frozen=", frozen_est,
       call. = FALSE)
}
write_csv(data.table(
  quantity = c("direct_contrast", "standard_error"),
  frozen_13K = c(frozen_est, frozen_se),
  current = c(fits$original_13K$contrast$estimate,
              fits$original_13K$contrast$standard_error)
), "09_original_13K_reproduction.csv")
fits$diagnosis_time <- fit_one(
  d, "diagnosis_time", "diagnosis_time", technical = FALSE
)
fits$technical_separate <- fit_one(
  d, "technical_separate", "technical", technical = TRUE
)
if (phase_available) {
  fits$diagnosis_time_technical <- fit_one(
    d, "diagnosis_time_technical", c("diagnosis_time", "technical"),
    technical = TRUE
  )
}

contrasts <- rbindlist(lapply(fits, `[[`, "contrast"), fill = TRUE)
coefficients <- rbindlist(lapply(fits, `[[`, "coefficients"), fill = TRUE)
write_csv(contrasts, "04_model_comparison.csv")
write_csv(coefficients, "05_model_coefficients.csv")

original_est <- contrasts[model_id == "original_13K", estimate]
contrasts[, change_from_original := estimate - original_est]
contrasts[, percent_change_from_original :=
  100 * (estimate - original_est) / abs(original_est)]
write_csv(contrasts, "06_model_comparison_with_change.csv")

diag_terms <- coefficients[grepl("time_years:diagnosis", term, fixed = TRUE)]
write_csv(diag_terms, "07_diagnosis_time_terms.csv")

flow <- data.table(
  item = c("people", "MRI_rows", "diagnosis_levels", "phase_available",
           "models_fit", "new_data_downloaded"),
  value = c(uniqueN(d$RID), nrow(d), uniqueN(d$diagnosis), phase_available,
            length(fits), FALSE)
)
write_csv(flow, "08_audit_flow.csv")

status <- c(
  "STATUS=DIAGNOSIS_TIME_SPECIFICATION_AUDIT_COMPLETE",
  "SCRIPT_REVISION=20261027_13K_REPRODUCTION_FIX",
  "ORIGINAL_13K_REPRODUCTION=PASS",
  "PRIMARY_FRAME=695_people_4659_MRI_rows",
  "PRIMARY_MODEL=original_13K_reproduced",
  "AUDIT_MODEL=original_plus_time_years_by_baseline_diagnosis",
  paste0("PHASE_AVAILABLE=", phase_available),
  paste0("MODELS_FIT=", length(fits)),
  "NEW_DATA_DOWNLOAD=FALSE",
  paste0("OUTPUT_DIR=", out_dir)
)
writeLines(status, file.path(out_dir, "00_status.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"),
           useBytes = TRUE)
message(paste(status, collapse = "\n"))
