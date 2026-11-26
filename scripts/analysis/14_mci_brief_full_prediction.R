# 427-person MCI five-year held-out benchmark
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

root <- study_path("project")
in_file <- study_path("project", "results/13F_R1_clinical_closure_confirmatory/09_baseline_event_frame.csv")
fold_file <- study_path("project", "results/13H46_limited_prediction_audit_release/01_ADNI_out_of_fold_predictions.csv")
out_dir <- file.path(root, "results", paste0(
  "13H47_MCI_strong_clinical_MRI_increment_",
  "release"))
dir.create(out_dir, recursive = TRUE)
save_csv <- function(x, name) fwrite(x, file.path(out_dir, name), bom = TRUE)

clinical_core <- c("age_index_z", "sex", "education_z", "baseline_MMSE")
clinical_strong <- c(clinical_core, "baseline_CDRSB", "baseline_FAQ",
  "baseline_mPACC")
model_terms <- list(core = clinical_core,
  core_plus_MRI = c(clinical_core, "structure_vulnerability_z"),
  strong = clinical_strong,
  strong_plus_MRI = c(clinical_strong, "structure_vulnerability_z"))

d <- fread(in_file, showProgress = FALSE)
d[, `:=`(linked = as.logical(linked),
  prevalent_dementia = as.logical(prevalent_dementia),
  event = as.integer(as.logical(event)),
  followup_years = as.numeric(followup_years))]
d <- d[linked == TRUE & prevalent_dementia == FALSE &
  baseline_dx == "MCI" & followup_years > 0]
if (nrow(d) != 428L || uniqueN(d$RID) != 428L)
  stop("Frozen MCI risk set is not 428 unique participants")

folds <- fread(fold_file, select = c("analysis", "RID", "fold"))[
  analysis == "baseline", .(RID, fold)]
if (nrow(folds) != 680L || uniqueN(folds$RID) != 680L)
  stop("13H46 baseline fold assignment changed")
d <- merge(d, folds, by = "RID", all.x = TRUE, sort = FALSE)
vars <- c("RID", "followup_years", "event", "fold",
  clinical_strong, "structure_vulnerability_z")
d <- droplevels(d[complete.cases(d[, ..vars])])
d[, `:=`(sex = factor(sex),
  baseline_MMSE = as.numeric(baseline_MMSE),
  baseline_CDRSB = as.numeric(baseline_CDRSB),
  baseline_FAQ = as.numeric(baseline_FAQ),
  baseline_mPACC = as.numeric(baseline_mPACC))]
if (nrow(d) != 427L || sum(d$event == 1L &
    d$followup_years <= 5) != 138L ||
    !setequal(d$fold, 1:5))
  stop("Complete MCI comparison frame changed")

predict_risk <- function(fit, newdata) {
  bh <- basehaz(fit, centered = FALSE)
  h0 <- if (any(bh$time <= 5)) max(bh$hazard[bh$time <= 5]) else 0
  lp <- predict(fit, newdata = newdata, type = "lp", reference = "zero")
  pmin(1 - 1e-8, pmax(1e-8, 1 - exp(-h0 * exp(lp))))
}

for (name in names(model_terms)) d[, (paste0("risk_", name)) := NA_real_]
for (k in 1:5) {
  train <- d[fold != k]
  test <- d[fold == k]
  for (name in names(model_terms)) {
    form <- as.formula(paste("Surv(followup_years, event) ~",
      paste(model_terms[[name]], collapse = " + ")))
    fit <- coxph(form, data = train, ties = "efron")
    if (any(!is.finite(coef(fit))))
      stop("Non-finite Cox coefficient: ", name, ", fold ", k)
    d[fold == k, (paste0("risk_", name)) := predict_risk(fit, test)]
  }
}
risk_cols <- paste0("risk_", names(model_terms))
if (any(!is.finite(as.matrix(d[, ..risk_cols]))))
  stop("Incomplete out-of-fold predictions")

save_csv(d[, c("RID", "fold", "index_mri_date", "risk_start",
  "followup_years", "event", "baseline_MMSE", "baseline_CDRSB",
  "baseline_FAQ", "baseline_mPACC", risk_cols), with = FALSE],
  "01_MCI_out_of_fold_predictions.csv")

weighted_auc <- function(y, risk, weight) {
  a <- data.table(y, risk, weight)[weight > 0,
    .(case_weight = sum(weight[y == 1L]),
      control_weight = sum(weight[y == 0L])), by = risk]
  setorder(a, risk)
  den <- sum(a$case_weight) * sum(a$control_weight)
  if (den == 0) return(NA_real_)
  sum(a$case_weight * (cumsum(a$control_weight) -
    0.5 * a$control_weight)) / den
}

ipcw <- function(x) {
  fit <- survfit(Surv(followup_years, 1L - event) ~ 1, data = x)
  G <- function(t) {
    j <- findInterval(t, fit$time)
    pmax(0.02, c(1, fit$surv)[j + 1L])
  }
  case <- x$event == 1L & x$followup_years <= 5
  control <- x$followup_years > 5
  w <- numeric(nrow(x))
  w[case] <- 1 / G(pmax(0, x$followup_years[case] - 1e-8))
  w[control] <- 1 / G(5)
  w
}

metrics <- function(y, risk, weight, calibrate = TRUE) {
  risk <- pmin(1 - 1e-8, pmax(1e-8, risk))
  ans <- c(AUC = weighted_auc(y, risk, weight),
    Brier = sum(weight * (y - risk)^2) / length(y),
    observed_risk = sum(weight * y) / length(y),
    predicted_mean = mean(risk),
    calibration_intercept = NA_real_, calibration_slope = NA_real_)
  if (calibrate) {
    z <- data.table(y, lp = qlogis(risk), weight)[weight > 0]
    slope_fit <- suppressWarnings(glm(y ~ lp, data = z,
      weights = weight, family = quasibinomial()))
    intercept_fit <- suppressWarnings(glm(y ~ 1 + offset(lp), data = z,
      weights = weight, family = quasibinomial()))
    ans["calibration_intercept"] <- unname(coef(intercept_fit)[1L])
    ans["calibration_slope"] <- unname(coef(slope_fit)["lp"])
  }
  ans
}

y <- as.integer(d$event == 1L & d$followup_years <= 5)
w <- ipcw(d)
flow <- data.table(initial_MCI = 428L, complete_MCI = nrow(d),
  events_all_followup = sum(d$event), cases_5y = sum(y),
  controls_5y = sum(d$followup_years > 5),
  early_censored = sum(w == 0), max_ipcw = max(w))
save_csv(flow, "02_MCI_sample_flow.csv")

performance <- rbindlist(lapply(names(model_terms), function(name) {
  z <- metrics(y, d[[paste0("risk_", name)]], w)
  data.table(model = name, metric = names(z), value = as.numeric(z))
}))
save_csv(performance, "03_MCI_model_performance.csv")

increment <- function(x, weight, calibrate = FALSE) {
  response <- as.integer(x$event == 1L & x$followup_years <= 5)
  s <- lapply(names(model_terms), function(name)
    metrics(response, x[[paste0("risk_", name)]], weight,
      calibrate = calibrate))
  names(s) <- names(model_terms)
  c(core_AUC = unname(s$core["AUC"]),
    core_Brier = unname(s$core["Brier"]),
    strong_AUC = unname(s$strong["AUC"]),
    strong_Brier = unname(s$strong["Brier"]),
    core_MRI_delta_AUC = unname(s$core_plus_MRI["AUC"] - s$core["AUC"]),
    core_MRI_delta_Brier = unname(s$core_plus_MRI["Brier"] - s$core["Brier"]),
    strong_MRI_delta_AUC = unname(s$strong_plus_MRI["AUC"] - s$strong["AUC"]),
    strong_MRI_delta_Brier = unname(s$strong_plus_MRI["Brier"] - s$strong["Brier"]),
    strong_minus_core_delta_AUC = unname((s$strong_plus_MRI["AUC"] -
      s$strong["AUC"]) - (s$core_plus_MRI["AUC"] - s$core["AUC"])),
    strong_minus_core_delta_Brier = unname((s$strong_plus_MRI["Brier"] -
      s$strong["Brier"]) - (s$core_plus_MRI["Brier"] - s$core["Brier"])))
}
observed <- increment(d, w)
set.seed(20261103L)
boot <- replicate(400L, {
  sample <- d[sample.int(nrow(d), replace = TRUE)]
  increment(sample, ipcw(sample))
})
intervals <- rbindlist(lapply(names(observed), function(name) {
  values <- boot[name, ]
  data.table(contrast = name, estimate = unname(observed[name]),
    ci_low = unname(quantile(values, .025, na.rm = TRUE)),
    ci_high = unname(quantile(values, .975, na.rm = TRUE)),
    bootstrap_replicates = sum(is.finite(values)))
}))
save_csv(intervals, "04_MCI_MRI_increment_by_tier.csv")

thresholds <- seq(.10, .30, .05)
dca <- rbindlist(lapply(names(model_terms), function(name) {
  risk <- d[[paste0("risk_", name)]]
  rbindlist(lapply(thresholds, function(t) {
    positive <- risk >= t
    all_benefit <- sum(w * (y - (1 - y) * t / (1 - t))) / nrow(d)
    data.table(model = name, threshold = t,
      net_benefit = sum(w * positive *
        (y - (1 - y) * t / (1 - t))) / nrow(d),
      treat_all = all_benefit, treat_none = 0)
  }))
}))
save_csv(dca, "05_MCI_decision_curve.csv")

calibration_bins <- rbindlist(lapply(names(model_terms), function(name) {
  risk <- d[[paste0("risk_", name)]]
  bin <- pmin(5L, ceiling(frank(risk, ties.method = "average") /
    length(risk) * 5L))
  data.table(bin, risk, y, w)[, .(n = .N, cases = sum(y),
    predicted_mean = mean(risk), observed_ipcw = sum(w * y) / .N),
    by = bin][, model := name]
}))
save_csv(calibration_bins, "06_MCI_calibration_bins.csv")

status <- c("STATUS=COMPLETE", "HORIZON_YEARS=5",
  "SAMPLE=MCI_ONLY_SAME_427_COMPLETE_PEOPLE",
  "VALIDATION=FIVE_FOLD_OUT_OF_FOLD_SAME_13H46_FOLDS",
  "CI=PAIRED_SUBJECT_BOOTSTRAP_CONDITIONAL_ON_OUT_OF_FOLD_PREDICTIONS",
  "CORE=AGE_SEX_EDUCATION_MMSE",
  "STRONG=CORE_PLUS_CDRSB_FAQ_MPACC",
  "MRI=BASELINE_TEMPORAL_STRUCTURE_VULNERABILITY",
  paste0("MCI_PEOPLE=", nrow(d)),
  paste0("FIVE_YEAR_CASES=", sum(y)),
  paste0("STRONG_DELTA_AUC=", signif(
    unname(observed["strong_MRI_delta_AUC"]), 5)),
  paste0("STRONG_DELTA_BRIER=", signif(
    unname(observed["strong_MRI_delta_Brier"]), 5)),
  paste0("OUTPUT_DIR=", out_dir))
writeLines(status, file.path(out_dir, "00_status.txt"))
writeLines(capture.output(sessionInfo()),
  file.path(out_dir, "07_sessionInfo.txt"))
cat(paste(status, collapse = "\n"), "\n")
