# Participant folds and held-out predictions used by the final MCI comparison
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

root <- study_path("project")
population_root <- file.path(root, "results")
out_dir <- file.path(root, "results", paste0("13H46_limited_prediction_audit_",
  "release"))
dir.create(out_dir, recursive = TRUE)
save_csv <- function(x, name) fwrite(x, file.path(out_dir, name), bom = TRUE)

predict_cox_risk <- function(fit, newdata, horizon = 5) {
  bh <- basehaz(fit, centered = FALSE)
  before <- bh$time <= horizon
  h0 <- if (any(before)) max(bh$hazard[before]) else 0
  lp <- predict(fit, newdata = newdata, type = "lp", reference = "zero")
  pmin(1 - 1e-8, pmax(1e-8, 1 - exp(-h0 * exp(lp))))
}

cross_validated_risk <- function(d, base_terms, add_term, seed) {
  set.seed(seed)
  d[, case5 := as.integer(event == 1L & followup_years <= 5)]
  d[, fold := sample(rep(1:5, length.out = .N)),
    by = .(baseline_dx, case5)]
  d[, `:=`(risk_base = NA_real_, risk_plus = NA_real_)]
  formula_for <- function(terms) as.formula(paste(
    "Surv(followup_years, event) ~", paste(terms, collapse = " + ")))
  for (k in 1:5) {
    train <- d[fold != k]
    test <- d[fold == k]
    base <- coxph(formula_for(base_terms), data = train, ties = "efron")
    plus <- coxph(formula_for(c(base_terms, add_term)), data = train,
      ties = "efron")
    if (any(!is.finite(coef(base))) || any(!is.finite(coef(plus))))
      stop("Non-finite Cox coefficient in fold ", k)
    d[fold == k, `:=`(risk_base = predict_cox_risk(base, test),
                       risk_plus = predict_cox_risk(plus, test))]
  }
  if (any(!is.finite(d$risk_base)) || any(!is.finite(d$risk_plus)))
    stop("Incomplete cross-validated predictions")
  d[, case5 := NULL]
  d
}

baseline_file <- file.path(population_root,
  "13F_R1_clinical_closure_confirmatory", "09_baseline_event_frame.csv")
baseline <- fread(baseline_file, showProgress = FALSE)
baseline[, `:=`(linked = as.logical(linked),
  prevalent_dementia = as.logical(prevalent_dementia),
  event = as.integer(as.logical(event)),
  followup_years = as.numeric(followup_years),
  baseline_MMSE = as.numeric(baseline_MMSE),
  baseline_FAQ = as.numeric(baseline_FAQ))]
baseline <- baseline[linked == TRUE & prevalent_dementia == FALSE &
  baseline_dx %in% c("CN", "MCI") & followup_years > 0]
if (nrow(baseline) != 680L || uniqueN(baseline$RID) != 680L)
  stop("Frozen baseline risk set is not 680 unique participants")
baseline_terms <- c("age_index_z", "sex", "education_z",
  "baseline_dx", "baseline_MMSE")
baseline_vars <- c("RID", "followup_years", "event", baseline_terms,
  "structure_vulnerability_z")
baseline <- droplevels(baseline[complete.cases(baseline[, ..baseline_vars])])
baseline[, `:=`(sex = factor(sex), baseline_dx = factor(baseline_dx))]
baseline <- cross_validated_risk(baseline, baseline_terms,
  "structure_vulnerability_z", 20261030L)

landmark_file <- file.path(population_root,
  "13F_R1_clinical_closure_confirmatory", "10_two_MRI_landmark_frame.csv")
landmark <- fread(landmark_file, showProgress = FALSE)
if (nrow(landmark) != 616L || uniqueN(landmark$RID) != 616L)
  stop("Frozen two-MRI frame is not 616 unique participants")
for (field in c("landmark_mri_date", "baseline_dx_date", "event_date",
                "censor_date")) set(landmark, j = field,
  value = as.IDate(landmark[[field]]))
landmark <- landmark[!is.na(baseline_dx_date) &
  !is.na(landmark_mri_date) & !is.na(censor_date)]
landmark[, risk_start := landmark_mri_date]
landmark[baseline_dx_date > risk_start, risk_start := baseline_dx_date]
landmark <- landmark[(is.na(event_date) | event_date > risk_start) &
  censor_date > risk_start]
landmark[, `:=`(
  event = as.integer(!is.na(event_date)),
  followup_years = (ifelse(!is.na(event_date), as.integer(event_date),
    as.integer(censor_date)) - as.integer(risk_start)) / 365.25,
  structure_vulnerability_z = as.numeric(baseline_structure_vulnerability_z),
  baseline_MMSE = as.numeric(baseline_MMSE),
  baseline_FAQ = as.numeric(baseline_FAQ))]
landmark <- landmark[followup_years > 0]
if (nrow(landmark) != 594L || sum(landmark$event) != 162L)
  stop("Strict post-landmark risk set no longer matches 594 / 162")
landmark_terms <- c(baseline_terms, "structure_vulnerability_z")
landmark_vars <- c("RID", "followup_years", "event", landmark_terms,
  "atrophy_rate_z")
landmark <- droplevels(landmark[complete.cases(landmark[, ..landmark_vars])])
landmark[, `:=`(sex = factor(sex), baseline_dx = factor(baseline_dx))]
landmark <- cross_validated_risk(landmark, landmark_terms,
  "atrophy_rate_z", 20261031L)

save_csv(rbindlist(list(
  baseline[, .(analysis = "baseline", RID, fold,
    index_mri_date = as.character(index_mri_date),
    risk_start = as.character(risk_start), baseline_dx,
    baseline_FAQ, followup_years, event, risk_base, risk_plus)],
  landmark[, .(analysis = "two_scan_landmark", RID, fold,
    index_mri_date = as.character(index_mri_date),
    landmark_mri_date = as.character(landmark_mri_date),
    risk_start = as.character(risk_start), baseline_dx,
    baseline_FAQ, followup_years, event, risk_base, risk_plus)]),
  fill = TRUE),
  "01_ADNI_out_of_fold_predictions.csv")

writeLines(capture.output(sessionInfo()),file.path(out_dir,"11_sessionInfo.txt"))
