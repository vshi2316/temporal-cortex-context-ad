# Longitudinal genetic models and spacing/QC sensitivities
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
suppressPackageStartupMessages({ library(data.table); library(lme4) })

pick_root <- function(ascii, native) if (file.exists(ascii)) ascii else native
project_dir <- pick_root(study_path("project"), study_path("project"))
result_root <- file.path(project_dir, "results")
prep_dir <- file.path(result_root, "13F_R_only_structural_analysis")
out_dir <- file.path(result_root, "13F_R_only_longitudinal_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
need <- function(path) { if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE); path }
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
finite <- function(x) is.finite(num(x))
zscore <- function(x) { x <- num(x); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (x - mean(x, na.rm = TRUE)) / s }
has_variation <- function(x) { z <- x[!is.na(x)]; length(unique(z)) >= 2L }

d <- as.data.table(readRDS(need(file.path(prep_dir, "04_ADNI_R_only_longitudinal_structural_frame.rds"))))
d[, RID := as.character(RID)]
d[, mri_date := as.IDate(as.character(mri_date))]
d[, age_mri := num(age_mri)]

# Baseline covariates come from the clinical/PGS baseline record.
pgs_file <- file.path(project_dir, "results", "10A_ADNI_gMAV_PGS_assembly", "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
if (file.exists(pgs_file)) {
  pgs <- fread(file = pgs_file, sep = "\t", showProgress = FALSE)
  pgs[, RID := as.character(RID)]
  pgs[, exam_date__ := as.IDate(as.character(EXAMDATE))]
  setorder(pgs, RID, exam_date__)
  base_cov <- pgs[, .(
    baseline_age = num(age_bl[which(is.finite(num(age_bl)))[1L]]),
    baseline_sex = as.character(sex[which(!is.na(sex) & nzchar(as.character(sex)))[1L]]),
    baseline_education = num(education[which(is.finite(num(education)))[1L]]),
    baseline_diagnosis = as.character(dx_bl[which(!is.na(dx_bl) & nzchar(as.character(dx_bl)))[1L]])
  ), by = RID]
  d <- merge(d, base_cov, by = "RID", all.x = TRUE, sort = FALSE)
} else {
  d[, `:=`(baseline_age = age_mri, baseline_sex = NA_character_, baseline_education = NA_real_, baseline_diagnosis = NA_character_)]
}
d[, `:=`(
  gMAV = num(clinical_gMAV_noAPOE_z_eur),
  gAD = num(clinical_gAD_noAPOE_z_eur),
  baseline_age = num(baseline_age),
  education = num(baseline_education),
  sex = factor(baseline_sex),
  diagnosis = factor(baseline_diagnosis),
  scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_"))
)]
pc_cols <- paste0("clinical_PC", 1:10)
for (v in pc_cols) if (v %in% names(d)) d[, (v) := num(get(v))]
d[, core_complete := finite(ad_vulnerable_index)]

# Deterministic one record per person/date; no duplicate dates are treated as
# independent observations. Prefer complete primary endpoints, then IMAGEUID.
setorderv(d, c("RID", "mri_date", "core_complete", "IMAGEUID"), order = c(1, 1, -1, -1), na.last = TRUE)
d <- d[, .SD[1L], by = .(RID, mri_date)]
setorder(d, RID, mri_date)
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
d[, age_mri := baseline_age + time_years]

# The sensitivity model excludes participants with any automatically flagged primary T1 image.
flagged_file <- file.path(project_dir, "results", "13F_automated_raw_T1_QC_and_processing_manifest", "04_automatic_QC_flagged_images.csv")
flagged_rids <- character()
raw_qc_available <- file.exists(flagged_file)
if (raw_qc_available) {
  flagged <- fread(file = flagged_file, showProgress = FALSE)
  flagged_rids <- unique(as.character(flagged[cohort == "ADNI", RID]))
}
if (!raw_qc_available) message("Raw-image QC file absent; raw-QC exclusion sensitivities are skipped.")
d[, raw_auto_flagged_participant := RID %in% flagged_rids]
d[, `:=`(
  ad_vulnerable_z = zscore(ad_vulnerable_index),
  global_cortical_z = zscore(global_cortical_thickness),
  hippocampus_z = zscore(hip_hippocampus_icv_ratio),
  wmh_z = zscore(wmh_wmh_brain_ratio),
  mPACCdigit_z = zscore(clinical_mPACCdigit),
  mPACCtrailsB_z = zscore(clinical_mPACCtrailsB)
)]
d <- d[finite(ad_vulnerable_z) & finite(gMAV) & finite(gAD)]

# Preserve a secondary cohort with >=2 distinct dates irrespective of spacing.
# The frozen primary cohort below still requires the prespecified >=270-day
# interval; this alternative checks whether conclusions depend on that rule.
d_no270 <- copy(d)

elig <- d[, .(
  n_dates = uniqueN(mri_date),
  second_date_270d = { z <- sort(unique(mri_date)); z <- z[as.integer(z - z[1L]) >= 270L]; if (length(z)) z[1L] else as.IDate(NA) }
), by = RID]
elig[, longitudinal_eligible := n_dates >= 2L & !is.na(second_date_270d)]
d <- d[RID %in% elig[longitudinal_eligible == TRUE, RID]]

covariate_names <- c("baseline_age", "sex", "education", "diagnosis", "scanner", pc_cols)
usable_covariates <- function(x) {
  out <- character()
  for (v in covariate_names) {
    if (!v %in% names(x)) next
    if (is.factor(x[[v]])) { if (nlevels(droplevels(x[[v]])) >= 2L) out <- c(out, v) }
    else if (has_variation(x[[v]])) out <- c(out, v)
  }
  unique(out)
}

fit_mixed <- function(x, outcome, rhs, model_id) {
  vars <- unique(c(outcome, "RID", all.vars(as.formula(paste("~", rhs)))))
  vars <- vars[vars %in% names(x)]
  a <- x[complete.cases(x[, vars, with = FALSE])]
  np <- uniqueN(a$RID); nr <- nrow(a); nrep <- a[, .N, by = RID][N >= 2L, .N]
  if (nr < 200L || np < 100L) return(data.table(model_id = model_id, outcome = outcome, term = NA_character_, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p_value = NA_real_, status = "NOT_RUN_INSUFFICIENT_DATA", n_people = np, n_rows = nr, n_repeated = nrep, random_structure = NA_character_))
  a <- droplevels(a)
  f1 <- as.formula(paste0(outcome, " ~ ", rhs, " + (1 + time_years | RID)"))
  f0 <- as.formula(paste0(outcome, " ~ ", rhs, " + (1 | RID)"))
  fit <- tryCatch(lmer(f1, data = a, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))), error = function(e) e)
  structure <- "random_intercept_and_time_slope"
  if (inherits(fit, "error") || isSingular(fit, tol = 1e-4)) { fit <- tryCatch(lmer(f0, data = a, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))), error = function(e) e); structure <- "random_intercept_fallback" }
  if (inherits(fit, "error")) return(data.table(model_id = model_id, outcome = outcome, term = NA_character_, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p_value = NA_real_, status = "FAILED", n_people = np, n_rows = nr, n_repeated = nrep, random_structure = structure, message = conditionMessage(fit)))
  sm <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(sm, c("Estimate", "Std. Error", "t value"), c("estimate", "std_error", "statistic"), skip_absent = TRUE)
  sm[, `:=`(model_id = model_id, outcome = outcome, status = "SUCCESS", p_value = 2 * pnorm(-abs(statistic)), n_people = np, n_rows = nr, n_repeated = nrep, random_structure = structure)]
  sm[, .(model_id, outcome, term, estimate, std_error, statistic, p_value, status, n_people, n_rows, n_repeated, random_structure)]
}

results <- list(); audit <- list(); k <- 0L
structural <- c("ad_vulnerable_z", "global_cortical_z", "hippocampus_z", "wmh_z")
run_specs <- c("primary", "exclude_raw_auto_flagged", "without_baseline_diagnosis")
for (outcome in structural) {
  if (!outcome %in% names(d)) next
  x <- d[finite(get(outcome))]
  covars <- usable_covariates(x)
  rhs_main <- paste(c("time_years", "gMAV", "gAD", covars), collapse = " + ")
  rhs_time <- paste(c("time_years * gMAV", "time_years * gAD", covars), collapse = " + ")
  for (spec in c("main", "time_interaction")) {
    for (variant in run_specs) {
      if (variant == "exclude_raw_auto_flagged" && !raw_qc_available) next
      xx <- x
      if (variant == "exclude_raw_auto_flagged") xx <- xx[raw_auto_flagged_participant == FALSE]
      covars_variant <- covars
      if (variant == "without_baseline_diagnosis") covars_variant <- setdiff(covars_variant, "diagnosis")
      rhs_main_variant <- paste(c("time_years", "gMAV", "gAD", covars_variant), collapse = " + ")
      rhs_time_variant <- paste(c("time_years * gMAV", "time_years * gAD", covars_variant), collapse = " + ")
      rhs_variant <- if (spec == "main") rhs_main_variant else rhs_time_variant
      id_prefix <- if (variant == "primary") "PRIMARY" else toupper(variant)
      model_id <- paste(id_prefix, outcome, spec, sep = "__")
      k <- k + 1L
      ans <- fit_mixed(xx, outcome, rhs_variant, model_id)
      results[[k]] <- ans; audit[[k]] <- ans[1L, .(model_id, outcome, status, n_people, n_rows, n_repeated, random_structure)]
    }
  }
}

# Spacing-threshold sensitivity for the primary structural endpoint only.
elig_no270 <- d_no270[, .(n_dates = uniqueN(mri_date)), by = RID]
d_no270 <- d_no270[RID %in% elig_no270[n_dates >= 2L, RID]]
if (nrow(d_no270)) {
  outcome <- "ad_vulnerable_z"
  x <- d_no270[finite(get(outcome))]
  covars <- usable_covariates(x)
  rhs <- paste(c("time_years * gMAV", "time_years * gAD", covars), collapse = " + ")
  for (variant in c("NO270", "NO270_EXCLUDE_RAW_AUTO_FLAGGED")) {
    if (grepl("EXCLUDE_RAW_AUTO_FLAGGED|exclude_raw_auto_flagged", variant) && !raw_qc_available) next
    xx <- if (variant == "NO270") x else x[raw_auto_flagged_participant == FALSE]
    k <- k + 1L
    ans <- fit_mixed(xx, outcome, rhs, paste(variant, outcome, "time_interaction", sep = "__"))
    results[[k]] <- ans
    audit[[k]] <- ans[1L, .(model_id, outcome, status, n_people, n_rows, n_repeated, random_structure)]
  }
}

# First-MRI structure predicts subsequent cognitive slope.
baseline_structure <- d[finite(ad_vulnerable_z), .SD[1L], by = RID][, .(RID, baseline_ad_vulnerable_z = ad_vulnerable_z)]
d <- merge(d, baseline_structure, by = "RID", all.x = TRUE, sort = FALSE)
for (outcome in c("mPACCdigit_z", "mPACCtrailsB_z")) {
  if (!outcome %in% names(d)) next
  x <- d[finite(baseline_ad_vulnerable_z) & finite(get(outcome))]
  covars <- usable_covariates(x)
  rhs <- paste(c("time_years * baseline_ad_vulnerable_z", covars), collapse = " + ")
  for (variant in c("primary", "exclude_raw_auto_flagged")) {
    if (grepl("EXCLUDE_RAW_AUTO_FLAGGED|exclude_raw_auto_flagged", variant) && !raw_qc_available) next
    xx <- if (variant == "primary") x else x[raw_auto_flagged_participant == FALSE]
    k <- k + 1L
    ans <- fit_mixed(xx, outcome, rhs, paste(variant, "cognitive_bridge", outcome, sep = "__"))
    results[[k]] <- ans; audit[[k]] <- ans[1L, .(model_id, outcome, status, n_people, n_rows, n_repeated, random_structure)]
  }
}

model_results <- rbindlist(results, fill = TRUE)
model_audit <- rbindlist(audit, fill = TRUE)
model_results[, family := fifelse(grepl("PRIMARY__ad_vulnerable_z__time_interaction", model_id), "primary_time_PGS", fifelse(grepl("global_cortical_z|cognitive_bridge", model_id), "secondary", "sensitivity"))]
model_results[, fdr := NA_real_]
ix <- model_results$family == "primary_time_PGS" & model_results$term %in% c("time_years:gMAV", "time_years:gAD")
if (sum(ix) > 0L) model_results[ix, fdr := p.adjust(p_value, method = "holm")]
ix2 <- model_results$family == "secondary" & grepl(":", model_results$term)
if (sum(ix2) > 0L) model_results[ix2, fdr := p.adjust(p_value, method = "BH")]

person_flow <- d[, .(n_dates = uniqueN(mri_date)), by = RID]
sample_flow <- data.table(
  set = c("PGS_complete_after_dedup", "people_with_at_least_2_dates", "people_with_270_day_spacing", "primary_structure_rows", "primary_structure_people", "hippocampus_rows", "WMH_rows", "mPACCdigit_rows", "mPACCtrailsB_rows"),
  rows = c(nrow(d_no270), NA_integer_, NA_integer_, sum(finite(d$ad_vulnerable_z)), NA_integer_, sum(finite(d$hippocampus_z)), sum(finite(d$wmh_z)), sum(finite(d$mPACCdigit_z)), sum(finite(d$mPACCtrailsB_z))),
  people = c(uniqueN(d_no270$RID), uniqueN(d_no270$RID), uniqueN(d$RID), uniqueN(d[finite(ad_vulnerable_z), RID]), uniqueN(d[finite(ad_vulnerable_z), RID]), uniqueN(d[finite(hippocampus_z), RID]), uniqueN(d[finite(wmh_z), RID]), uniqueN(d[finite(mPACCdigit_z), RID]), uniqueN(d[finite(mPACCtrailsB_z), RID]))
)
interpretation <- c(

  "Primary structure: frozen AD-vulnerable cortical index. Primary genetic terms: time_years:gMAV and time_years:gAD.",
  "mPACCdigit is the primary cognitive bridge; mPACCtrailsB is sensitivity."

)
fwrite(model_results, file.path(out_dir, "01_longitudinal_model_coefficients.csv"), bom = TRUE)
fwrite(model_audit, file.path(out_dir, "02_model_execution_audit.csv"), bom = TRUE)
fwrite(sample_flow, file.path(out_dir, "03_sample_flow.csv"), bom = TRUE)
writeLines(interpretation, file.path(out_dir, "06_interpretation_scope.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
message("13F-R1 longitudinal models completed: ", out_dir)
