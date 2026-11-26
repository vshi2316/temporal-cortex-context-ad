# AIBL risk set and incident-AD models
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

# Reconstruct incident AD from AIBL visit-level diagnoses and baseline CAT12 volume.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

PROJECT_ROOT <- study_path("project")
OUTDIR <- file.path(
  PROJECT_ROOT,
  "results",
  "13J_outcome_harmonization_and_AIBL_incident_AD_validation"
)

CAT12_ENDPOINT <- file.path(
  PROJECT_ROOT,
  "results",
  "13F_CAT12_R_direct_minimal",
  "tables",
  "04_AIBL_CAT12_endpoint.csv"
)
AIBL_CLINICAL <- study_path("aibl_clinical")
AIBL_DX <- study_path("aibl_data", "aibl_pdxconv_01-Jun-2018.csv")
AIBL_MRI_MANIFEST <- study_path("aibl_mri_manifest")
AIBL_MRI_META <- study_path("aibl_data", "aibl_mrimeta_01-Jun-2018.csv")

required_files <- c(
  CAT12_ENDPOINT,
  AIBL_CLINICAL,
  AIBL_DX,
  AIBL_MRI_MANIFEST,
  AIBL_MRI_META
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Required input file(s) missing:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

norm_id <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("^0+", "", y)
  y[y == ""] <- "0"
  y
}

zscore <- function(x) {
  x <- as_num(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

parse_visit_month <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x))
  out[x == "bl"] <- 0
  is_month <- grepl("^m[0-9]+$", x)
  out[is_month] <- as_num(sub("^m", "", x[is_month]))
  out
}

parse_us_date <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "-4", "NA")] <- NA_character_
  as.IDate(x, format = "%m/%d/%Y")
}

write_csv <- function(x, name) {
  fwrite(x, file.path(OUTDIR, name), bom = TRUE)
}

# ----------------------------------------------------------------------------
# 2. Read CAT12, visit-level AIBL diagnosis, and covariates
# ----------------------------------------------------------------------------

cat12 <- fread(CAT12_ENDPOINT, showProgress = FALSE)
clinical <- fread(AIBL_CLINICAL, showProgress = FALSE)
dx <- fread(AIBL_DX, showProgress = FALSE)
mri_manifest <- fread(AIBL_MRI_MANIFEST, showProgress = FALSE)
mri_meta <- fread(AIBL_MRI_META, showProgress = FALSE)

required_cat12 <- c(
  "RID", "image_id", "report_status", "n_target_rois",
  "cat12_temporal_gm_tiv_ratio", "cat12_ad_vulnerable_gm"
)
required_clinical <- c(
  "RID", "Gender", "Age", "APOE4_Carrier", "DX_Baseline",
  "MMSE", "Converted_to_AD", "Conversion_Month", "Max_Followup", "Event"
)
required_dx <- c("RID", "VISCODE", "DXCURREN", "DXAD", "DXOTHDEM")

for (spec in list(
  list(n = names(cat12), r = required_cat12, label = "CAT12 endpoint"),
  list(n = names(clinical), r = required_clinical, label = "AIBL clinical"),
  list(n = names(dx), r = required_dx, label = "AIBL pdxconv")
)) {
  miss <- setdiff(spec$r, spec$n)
  if (length(miss)) {
    stop(spec$label, " is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

cat12[, RID := norm_id(RID)]
cat12[, image_id := trimws(as.character(image_id))]
if (cat12[, anyDuplicated(RID)] > 0L) {
  stop("CAT12 endpoint contains duplicate RID values.", call. = FALSE)
}

clinical[, RID := norm_id(RID)]
clinical <- unique(clinical, by = "RID")

dx[, `:=`(
  RID = norm_id(RID),
  visit_month = parse_visit_month(VISCODE),
  dx_current = as.integer(as_num(DXCURREN)),
  dx_ad = as.integer(as_num(DXAD)),
  dx_other_dementia = as.integer(as_num(DXOTHDEM))
)]
dx <- dx[is.finite(visit_month)]
setorder(dx, RID, visit_month)

# AIBL coding used by the supplied release:
# DXCURREN: 1 = cognitively normal, 2 = MCI, 3 = dementia.
# AD events additionally require DXAD == 1, preventing other dementias from
# being counted as Alzheimer disease.
dx[, `:=`(
  valid_diagnosis = dx_current %in% c(1L, 2L, 3L),
  ad_diagnosis = dx_current == 3L & dx_ad == 1L,
  other_dementia = dx_current == 3L & dx_ad != 1L
)]

baseline_dx <- dx[
  visit_month == 0 & valid_diagnosis,
  .SD[1L],
  by = RID
][, .(
  RID,
  baseline_dx_code = dx_current,
  baseline_ad = ad_diagnosis,
  baseline_other_dementia = other_dementia
)]

followup_state <- dx[visit_month > 0 & valid_diagnosis, {
  ad_months <- visit_month[ad_diagnosis]
  other_months <- visit_month[other_dementia]
  first_ad <- if (length(ad_months)) min(ad_months) else NA_real_
  first_other <- if (length(other_months)) min(other_months) else NA_real_
  last_valid <- max(visit_month)

  # A preceding non-AD dementia ends the AD risk period. Same-visit AD coding
  # takes precedence because it is explicitly coded as AD.
  incident_ad <- is.finite(first_ad) &&
    (!is.finite(first_other) || first_ad <= first_other)
  analysis_month <- if (incident_ad) {
    first_ad
  } else if (is.finite(first_other)) {
    min(first_other, last_valid)
  } else {
    last_valid
  }

  list(
    first_ad_month = first_ad,
    first_other_dementia_month = first_other,
    last_valid_dx_month = last_valid,
    incident_ad = incident_ad,
    analysis_month = analysis_month,
    n_postbaseline_dx = .N
  )
}, by = RID]

# ----------------------------------------------------------------------------
# 3. Audit MRI baseline alignment
# ----------------------------------------------------------------------------

setnames(
  mri_manifest,
  old = c("Subject ID", "Image ID", "Visit", "Study Date", "Age"),
  new = c("RID", "image_id", "manifest_visit", "manifest_study_date", "manifest_mri_age"),
  skip_absent = TRUE
)
mri_manifest[, `:=`(
  RID = norm_id(RID),
  image_id = trimws(as.character(image_id)),
  manifest_study_date = parse_us_date(manifest_study_date),
  manifest_mri_age = as_num(manifest_mri_age)
)]
mri_manifest <- unique(
  mri_manifest[, .(
    RID, image_id, manifest_visit, manifest_study_date, manifest_mri_age
  )],
  by = c("RID", "image_id")
)

mri_meta[, `:=`(
  RID = norm_id(RID),
  meta_exam_date = parse_us_date(EXAMDATE)
)]
mri_meta_bl <- mri_meta[
  tolower(trimws(VISCODE)) == "bl",
  .(meta_exam_date = meta_exam_date[which.max(!is.na(meta_exam_date))]),
  by = RID
]

image_link <- merge(
  cat12[, .(RID, image_id)],
  mri_manifest,
  by = c("RID", "image_id"),
  all.x = TRUE,
  sort = FALSE
)
image_link <- merge(image_link, mri_meta_bl, by = "RID", all.x = TRUE, sort = FALSE)
image_link[, `:=`(
  manifest_baseline = tolower(trimws(manifest_visit)) == "baseline",
  exact_date_gap_days = abs(as.numeric(manifest_study_date - meta_exam_date))
)]

# ----------------------------------------------------------------------------
# 4. Assemble the incident-AD risk set
# ----------------------------------------------------------------------------

covariates <- clinical[, .(
  RID,
  age = as_num(Age),
  sex_code = as.integer(as_num(Gender)),
  apoe4 = as.integer(as_num(APOE4_Carrier)),
  mmse = as_num(MMSE),
  integrated_baseline_dx = as.integer(as_num(DX_Baseline)),
  old_event = as.integer(as_num(Event) == 1),
  old_converted_to_ad = as.integer(as_num(Converted_to_AD) == 1),
  old_conversion_month = as_num(Conversion_Month),
  old_max_followup_month = as_num(Max_Followup)
)]

assembled <- Reduce(
  function(x, y) merge(x, y, by = "RID", all.x = TRUE, sort = FALSE),
  list(
    cat12,
    baseline_dx,
    followup_state,
    covariates,
    image_link[, .(
      RID, image_id, manifest_visit, manifest_study_date,
      manifest_mri_age, meta_exam_date, manifest_baseline,
      exact_date_gap_days
    )]
  )
)

assembled[, `:=`(
  eligible_baseline_dx = baseline_dx_code %in% c(1L, 2L) &
    baseline_ad == FALSE & baseline_other_dementia == FALSE,
  cat12_qc = report_status == "PARSED" &
    as.integer(n_target_rois) == 12L &
    is.finite(as_num(cat12_temporal_gm_tiv_ratio)) &
    as_num(cat12_temporal_gm_tiv_ratio) > 0,
  positive_followup = is.finite(analysis_month) & analysis_month > 0
)]

analysis <- assembled[
  eligible_baseline_dx == TRUE &
    cat12_qc == TRUE &
    positive_followup == TRUE &
    manifest_baseline == TRUE
]

analysis[, `:=`(
  baseline_dx = factor(
    baseline_dx_code,
    levels = c(1L, 2L),
    labels = c("CN", "MCI")
  ),
  sex = factor(
    sex_code,
    levels = c(1L, 2L),
    labels = c("Male", "Female")
  ),
  event = as.integer(incident_ad),
  followup_years = analysis_month / 12,
  age_z = zscore(age),
  mmse_z = zscore(mmse),
  # Re-standardize the frozen CAT12 ROI ratio inside the baseline at-risk
  # population. Higher values continue to mean greater structural vulnerability.
  cat12_vulnerability_z = -zscore(log(as_num(cat12_temporal_gm_tiv_ratio)))
)]

if (!nrow(analysis)) {
  stop("The corrected AIBL incident-AD risk set is empty.", call. = FALSE)
}

event_n <- sum(analysis$event)
if (event_n < 10L) {
  warning("Corrected AIBL analysis has fewer than 10 incident AD events; Cox estimates will not be fitted.")
}

# ----------------------------------------------------------------------------
# 5. Sample-flow, contamination, and time-alignment audits
# ----------------------------------------------------------------------------

sample_flow <- data.table(
  stage = c(
    "CAT12 endpoints",
    "visit-level baseline diagnosis available",
    "excluded prevalent baseline AD",
    "excluded baseline other dementia",
    "baseline CN/MCI",
    "baseline CN/MCI with post-baseline diagnosis",
    "CAT12 and baseline-image QC complete",
    "corrected incident-AD analysis",
    "corrected incident AD events"
  ),
  n = c(
    nrow(cat12),
    sum(!is.na(assembled$baseline_dx_code)),
    sum(assembled$baseline_ad %in% TRUE, na.rm = TRUE),
    sum(assembled$baseline_other_dementia %in% TRUE, na.rm = TRUE),
    sum(assembled$eligible_baseline_dx %in% TRUE, na.rm = TRUE),
    sum(assembled$eligible_baseline_dx %in% TRUE & assembled$positive_followup %in% TRUE, na.rm = TRUE),
    sum(assembled$eligible_baseline_dx %in% TRUE & assembled$positive_followup %in% TRUE & assembled$cat12_qc %in% TRUE & assembled$manifest_baseline %in% TRUE, na.rm = TRUE),
    nrow(analysis),
    event_n
  )
)

event_reclassification <- assembled[, .N, by = .(
  baseline_dx_code,
  old_event,
  corrected_incident_ad = as.integer(incident_ad %in% TRUE)
)][order(baseline_dx_code, old_event, corrected_incident_ad)]

event_by_baseline_dx <- analysis[, .(
  n = .N,
  incident_ad_events = sum(event),
  median_followup_years = median(followup_years),
  maximum_followup_years = max(followup_years)
), by = baseline_dx]

alignment_audit <- data.table(
  metric = c(
    "CAT12 endpoints with baseline-labelled image",
    "Corrected risk-set participants with baseline-labelled image",
    "CAT12 endpoints with exact MRI metadata date available",
    "Exact MRI metadata dates matching image-manifest date",
    "Corrected risk-set incident events",
    "Corrected risk-set incident events with exact MRI date match"
  ),
  n = c(
    sum(image_link$manifest_baseline %in% TRUE, na.rm = TRUE),
    sum(analysis$manifest_baseline %in% TRUE, na.rm = TRUE),
    sum(is.finite(image_link$exact_date_gap_days)),
    sum(image_link$exact_date_gap_days == 0, na.rm = TRUE),
    event_n,
    sum(analysis$event == 1L & analysis$exact_date_gap_days == 0, na.rm = TRUE)
  ),
  interpretation = c(
    "Categorical visit alignment",
    "Primary time-origin requirement",
    "Exact-date sensitivity availability",
    "Exact-date concordance when both dates are present",
    "Corrected event support",
    "Exact-date sensitivity event support"
  )
)

write_csv(sample_flow, "01_AIBL_corrected_sample_flow.csv")
write_csv(event_reclassification, "02_AIBL_old_vs_corrected_event_reclassification.csv")
write_csv(event_by_baseline_dx, "03_AIBL_incident_AD_by_baseline_diagnosis.csv")
write_csv(alignment_audit, "04_AIBL_MRI_time_alignment_audit.csv")

# ----------------------------------------------------------------------------
# 6. Cox models for incident AD
# ----------------------------------------------------------------------------

fit_cox <- function(d, model_id, rhs) {
  vars <- unique(c("followup_years", "event", all.vars(as.formula(paste("~", rhs)))))
  dd <- d[complete.cases(d[, ..vars])]
  ne <- sum(dd$event)

  empty_coef <- data.table(
    model_id = model_id,
    term = NA_character_,
    estimate = NA_real_,
    standard_error = NA_real_,
    hazard_ratio = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = NA_real_,
    n = nrow(dd),
    events = ne,
    status = "INSUFFICIENT_EVENTS"
  )
  empty_diag <- data.table(
    model_id = model_id,
    n = nrow(dd),
    events = ne,
    c_index = NA_real_,
    AIC = NA_real_,
    logLik = NA_real_,
    status = "INSUFFICIENT_EVENTS"
  )

  if (ne < 10L || nrow(dd) < 50L) {
    return(list(coef = empty_coef, diag = empty_diag, ph = data.table()))
  }

  fit <- tryCatch(
    coxph(
      as.formula(paste0("Surv(followup_years, event) ~ ", rhs)),
      data = dd,
      ties = "efron",
      x = TRUE,
      model = TRUE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    empty_coef[, status := paste0("FAILED: ", conditionMessage(fit))]
    empty_diag[, status := paste0("FAILED: ", conditionMessage(fit))]
    return(list(coef = empty_coef, diag = empty_diag, ph = data.table()))
  }

  b <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  z <- b / se
  coef_table <- data.table(
    model_id = model_id,
    term = names(b),
    estimate = as.numeric(b),
    standard_error = as.numeric(se),
    hazard_ratio = exp(as.numeric(b)),
    conf_low = exp(as.numeric(b) - qnorm(0.975) * as.numeric(se)),
    conf_high = exp(as.numeric(b) + qnorm(0.975) * as.numeric(se)),
    p_value = 2 * pnorm(-abs(as.numeric(z))),
    n = nrow(dd),
    events = ne,
    status = "SUCCESS"
  )

  sm <- summary(fit)
  diag_table <- data.table(
    model_id = model_id,
    n = nrow(dd),
    events = ne,
    c_index = unname(sm$concordance[1L]),
    AIC = AIC(fit),
    logLik = as.numeric(logLik(fit)),
    status = "SUCCESS"
  )

  ph_table <- tryCatch({
    ph <- cox.zph(fit)
    out <- as.data.table(ph$table, keep.rownames = "term")
    setnames(out, names(out), make.names(names(out), unique = TRUE))
    out[, model_id := model_id]
    setcolorder(out, c("model_id", setdiff(names(out), "model_id")))
    out
  }, error = function(e) {
    data.table(
      model_id = model_id,
      term = "PH_TEST_FAILED",
      error = conditionMessage(e)
    )
  })

  list(coef = coef_table, diag = diag_table, ph = ph_table)
}

model_specs <- list(
  AIBL_incident_AD_primary =
    "cat12_vulnerability_z + age_z + sex + strata(baseline_dx)",
  AIBL_incident_AD_plus_APOE4 =
    "cat12_vulnerability_z + age_z + sex + apoe4 + strata(baseline_dx)",
  AIBL_incident_AD_plus_MMSE =
    "cat12_vulnerability_z + age_z + sex + mmse_z + strata(baseline_dx)",
  AIBL_incident_AD_unadjusted_stratified =
    "cat12_vulnerability_z + strata(baseline_dx)"
)

model_fits <- lapply(names(model_specs), function(id) {
  fit_cox(analysis, id, model_specs[[id]])
})

coefficients <- rbindlist(lapply(model_fits, `[[`, "coef"), fill = TRUE)
model_diagnostics <- rbindlist(lapply(model_fits, `[[`, "diag"), fill = TRUE)
ph_tests <- rbindlist(lapply(model_fits, `[[`, "ph"), fill = TRUE)

write_csv(coefficients, "05_AIBL_incident_AD_cox_coefficients.csv")
write_csv(model_diagnostics, "06_AIBL_incident_AD_model_diagnostics.csv")
write_csv(ph_tests, "07_AIBL_incident_AD_PH_tests.csv")
write_csv(analysis, "08_AIBL_corrected_incident_AD_analysis_frame.csv")

# ----------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()),file.path(OUTDIR,"99_sessionInfo.txt"))
