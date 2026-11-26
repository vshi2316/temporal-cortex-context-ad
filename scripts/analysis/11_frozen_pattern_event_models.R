# Frozen 34-region pattern at the second MRI origin
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

project <- study_path("project")
main_results <- study_path("project", "results")
aligned_dir <- study_path("project", "results/13H1_aligned_regional_slope_map_release")
clinical_file <- file.path(main_results, "13F_R1_clinical_closure_confirmatory",
                           "10_two_MRI_landmark_frame.csv")
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
template[, genetic_slope_difference :=
           aligned$estimate[match(region, aligned$region)]]
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
setorderv(fs, c("RID", "mri_date", "IMAGEUID"),
          order = c(1L, 1L, -1L), na.last = TRUE)
fs <- fs[, .SD[1L], by = .(RID, mri_date)]
fs <- merge(fs, keys, by = c("RID", "mri_date"), sort = FALSE)
if (nrow(fs) != nrow(keys) ||
    any(as.character(fs$IMAGEUID) != fs$main_imageuid))
  stop("FreeSurfer scan keys or Image IDs differ from 13K.")
for (i in seq_len(nrow(pairs))) {
  left <- pairs$left_field[i]
  right <- pairs$right_field[i]
  fs[, (regions[i]) := rowMeans(.SD, na.rm = TRUE), .SDcols = c(left, right)]
}
setorder(fs, RID, mri_date)
first_scan <- fs[, .SD[1L], by = RID]
template[, baseline_sd_mm := vapply(region, function(v)
  sd(first_scan[[v]], na.rm = TRUE), numeric(1))]
if (any(!is.finite(template$baseline_sd_mm) | template$baseline_sd_mm <= 0))
  stop("At least one regional baseline thickness SD is invalid.")
fs[, temporal_thickness_mm := rowMeans(.SD, na.rm = FALSE), .SDcols = temporal]
scan <- fs[, c("RID", "mri_date", "temporal_thickness_mm", "IMAGEUID", regions),
           with = FALSE]
setnames(scan, "IMAGEUID", "imageuid")
index <- copy(scan)
setnames(index, c("mri_date", "temporal_thickness_mm", "imageuid", regions),
         c("index_mri_date", "index_thickness_mm", "index_imageuid",
           paste0("index__", regions)))
landmark <- copy(scan)
setnames(landmark, c("mri_date", "temporal_thickness_mm", "imageuid", regions),
         c("landmark_mri_date", "landmark_thickness_mm", "landmark_imageuid",
           paste0("landmark__", regions)))

clinical <- fread(clinical_file)
clinical[, `:=`(RID = as.character(RID),
                index_mri_date = as.IDate(index_mri_date),
                landmark_mri_date = as.IDate(landmark_mri_date),
                baseline_dx_date = as.IDate(baseline_dx_date),
                event_date = as.IDate(event_date),
                censor_date = as.IDate(censor_date))]
if (nrow(clinical) != 616L || uniqueN(clinical$RID) != 616L)
  stop("The frozen two-MRI clinical landmark frame changed.")
x <- merge(clinical, index, by = c("RID", "index_mri_date"), all.x = TRUE,
           sort = FALSE)
x <- merge(x, landmark, by = c("RID", "landmark_mri_date"), all.x = TRUE,
           sort = FALSE)
if (nrow(x) != nrow(clinical)) stop("Clinical-to-MRI linkage is not one-to-one.")
x[, regional_thinning_mm_per_year :=
    (index_thickness_mm - landmark_thickness_mm) /
    (as.numeric(landmark_mri_date - index_mri_date) / 365.25)]
equivalence_error <- max(abs(x$regional_thinning_mm_per_year -
                             num(x$atrophy_rate)), na.rm = TRUE)
if (!is.finite(equivalence_error) || equivalence_error > 1e-8)
  stop("The six-region mean no longer reproduces the original atrophy rate.")
interval <- as.numeric(x$landmark_mri_date - x$index_mri_date) / 365.25
regional_rates <- sapply(seq_along(regions), function(i) {
  (num(x[[paste0("index__", regions[i])]]) -
     num(x[[paste0("landmark__", regions[i])]])) /
    interval / template$baseline_sd_mm[i]
})
x[, spatial_pattern_score := as.vector(regional_rates %*% template$weight)]
x[, risk_start_strict := landmark_mri_date]
x[!is.na(baseline_dx_date) & baseline_dx_date > risk_start_strict,
  risk_start_strict := baseline_dx_date]
x[, `:=`(event_before_risk_start = !is.na(event_date) &
           event_date <= risk_start_strict,
         censor_before_risk_start = !is.na(censor_date) &
           censor_date <= risk_start_strict)]
x[, `:=`(event_strict = as.integer(!is.na(event_date)),
         followup_strict = (fifelse(!is.na(event_date), as.numeric(event_date),
                                   as.numeric(censor_date)) -
                            as.numeric(risk_start_strict)) / 365.25)]
flow <- data.table(stage = "frozen_two_MRI_landmark", people = nrow(x),
                   events = sum(x$event_strict))
add_flow <- function(label, dat) data.table(stage = label, people = nrow(dat),
                                             events = sum(dat$event_strict))
linked <- x[!is.na(index_imageuid) & !is.na(landmark_imageuid)]
flow <- rbind(flow, add_flow("both_FreeSurfer_scans_matched", linked))
linked <- linked[is.finite(spatial_pattern_score)]
flow <- rbind(flow, add_flow("all_34_pattern_regions_complete", linked))
linked <- linked[!is.na(baseline_dx_date) &
                   event_before_risk_start == FALSE &
                   censor_before_risk_start == FALSE &
                   is.finite(followup_strict) & followup_strict > 0]
flow <- rbind(flow, add_flow("strict_post_landmark_risk_start", linked))

linked[, `:=`(spatial_pattern_z = zscore(spatial_pattern_score),
              event = event_strict, followup_years = followup_strict,
              gAD = num(gAD), gMAV = num(gMAV),
              baseline_structure_vulnerability_z = num(baseline_structure_vulnerability_z),
              atrophy_rate_z = num(atrophy_rate_z),
              baseline_mPACC_z = zscore(baseline_mPACC),
              baseline_CDRSB_z = zscore(baseline_CDRSB),
              baseline_FAQ_z = zscore(baseline_FAQ),
              age_index_z = num(age_index_z),
              education_z = num(education_z), apoe4_z = num(apoe4_z),
              sex = factor(sex), scanner = factor(scanner),
              baseline_dx = factor(baseline_dx))]
pc <- paste0("pc", 1:10, "_z")
for (v in pc) linked[, (v) := num(get(v))]
strong <- c("age_index_z", "sex", "education_z", "apoe4_z", "scanner",
            "baseline_dx", pc, "gAD", "gMAV", "baseline_mPACC_z",
            "baseline_CDRSB_z", "baseline_FAQ_z",
            "baseline_structure_vulnerability_z", "atrophy_rate_z")
vars <- c("followup_years", "event", strong, "spatial_pattern_z")
complete_rows <- complete.cases(linked[, vars, with = FALSE])
dat <- droplevels(linked[complete_rows])
flow <- rbind(flow, add_flow("common_complete_model_frame", dat))

out <- file.path(project, "results", paste0("13H2_spatial_pattern_landmark_validation_",
                                           "release"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(flow, file.path(out, "01_sample_flow.csv"), bom = TRUE)
fwrite(template, file.path(out, "01_frozen_34_region_template.csv"), bom = TRUE)
fwrite(x[, .(RID, index_mri_date, landmark_mri_date, index_imageuid,
             landmark_imageuid, index_thickness_mm, landmark_thickness_mm,
             regional_thinning_mm_per_year, spatial_pattern_score,
             baseline_dx_date,
             risk_start_strict, event_date, censor_date,
             event_before_risk_start, censor_before_risk_start,
             event_strict, followup_strict)],
       file.path(out, "02_landmark_date_and_image_audit.csv"), bom = TRUE)
if (nrow(dat) < 100L || sum(dat$event) < 50L)
  stop("Fewer than 100 complete people or 50 later events; inspect 01_sample_flow.csv. Output: ", out)

fit_cox <- function(train, terms, label = "bootstrap", strict = FALSE) {
  form <- as.formula(paste("Surv(followup_years, event) ~",
                           paste(terms, collapse = " + ")))
  fit <- tryCatch(coxph(form, data = train, ties = "efron", x = TRUE,
                        y = TRUE, model = TRUE), error = function(e) e)
  if (inherits(fit, "error")) {
    if (strict) stop(label, " Cox error: ", conditionMessage(fit))
    return(NULL)
  }
  bad <- names(coef(fit))[!is.finite(coef(fit))]
  if (length(bad)) {
    if (strict) stop(label, " has non-finite coefficients: ",
                     paste(bad, collapse = ", "))
    return(NULL)
  }
  fit
}
predict_risk <- function(fit, newdata, horizon = 5) {
  if (is.null(fit)) return(rep(NA_real_, nrow(newdata)))
  bh <- basehaz(fit, centered = FALSE)
  h0 <- bh[bh$time <= horizon, "hazard"]
  h0 <- if (length(h0)) max(h0) else 0
  lp <- tryCatch(predict(fit, newdata = newdata, type = "lp",
                         reference = "zero"), error = function(e) rep(NA_real_, nrow(newdata)))
  pmin(pmax(1 - exp(-h0 * exp(lp)), 1e-8), 1 - 1e-8)
}
metrics <- function(test, risk, horizon = 5) {
  if (any(!is.finite(risk))) return(c(AUC = NA_real_, Brier = NA_real_))
  time <- test$followup_years
  event <- test$event
  case <- event == 1L & time <= horizon
  control <- time > horizon
  if (!any(case) || !any(control)) return(c(AUC = NA_real_, Brier = NA_real_))
  censor_fit <- survfit(Surv(time, 1L - event) ~ 1)
  g <- function(at) {
    j <- findInterval(at, censor_fit$time)
    pmax(c(1, censor_fit$surv)[j + 1L], 1e-6)
  }
  w <- numeric(nrow(test))
  w[case] <- 1 / g(pmax(time[case] - 1e-8, 0))
  w[control] <- 1 / g(horizon)
  cmp <- outer(risk[case], risk[control], `-`)
  pair_w <- outer(w[case], w[control], `*`)
  c(AUC = sum(pair_w * (cmp > 0) + 0.5 * pair_w * (cmp == 0)) /
      sum(pair_w),
    Brier = sum(w * (as.integer(case) - risk)^2) / nrow(test))
}

base_fit <- fit_cox(dat, strong,
                    "strong clinical and genetic benchmark plus temporal mean",
                    strict = TRUE)
aug_fit <- fit_cox(dat, c(strong, "spatial_pattern_z"),
                   "strong benchmark plus temporal mean and spatial pattern", strict = TRUE)
beta <- coef(aug_fit)[["spatial_pattern_z"]]
se <- sqrt(vcov(aug_fit)["spatial_pattern_z", "spatial_pattern_z"])
hr <- data.table(term = "spatial_pattern_z", HR = exp(beta),
                 lower_95 = exp(beta - 1.96 * se),
                 upper_95 = exp(beta + 1.96 * se),
                 p_value = 2 * pnorm(-abs(beta / se)),
                 n = nrow(dat), events = sum(dat$event))
lr <- 2 * (as.numeric(logLik(aug_fit)) - as.numeric(logLik(base_fit)))
hr[, `:=`(lr_chisq = lr, lr_p_value = pchisq(lr, df = 1, lower.tail = FALSE))]
fwrite(hr, file.path(out, "03_regional_increment_cox.csv"), bom = TRUE)
ph <- tryCatch(cox.zph(aug_fit), error = function(e) NULL)
if (!is.null(ph))
  fwrite(as.data.table(ph$table, keep.rownames = "term"),
         file.path(out, "03_proportional_hazards.csv"), bom = TRUE)

apparent_base <- metrics(dat, predict_risk(base_fit, dat))
apparent_aug <- metrics(dat, predict_risk(aug_fit, dat))
apparent_delta <- apparent_aug - apparent_base
set.seed(20261019L)
B <- 500L
boot <- rbindlist(lapply(seq_len(B), function(i) {
  train <- dat[sample.int(nrow(dat), nrow(dat), replace = TRUE)]
  f0 <- fit_cox(train, strong)
  f1 <- fit_cox(train, c(strong, "spatial_pattern_z"))
  if (is.null(f0) || is.null(f1))
    return(data.table(replicate = i, metric = names(apparent_delta),
                      bootstrap_delta = NA_real_, test_delta = NA_real_))
  boot_delta <- metrics(train, predict_risk(f1, train)) -
    metrics(train, predict_risk(f0, train))
  test_delta <- metrics(dat, predict_risk(f1, dat)) -
    metrics(dat, predict_risk(f0, dat))
  data.table(replicate = i, metric = names(boot_delta),
             bootstrap_delta = as.numeric(boot_delta),
             test_delta = as.numeric(test_delta))
}))
boot[, optimism := bootstrap_delta - test_delta]
fwrite(boot, file.path(out, "04_bootstrap_increment.csv"), bom = TRUE)
performance <- boot[, .(successful_bootstraps = sum(is.finite(optimism)),
                        mean_optimism = mean(optimism, na.rm = TRUE)), by = metric]
performance[, `:=`(apparent_base = as.numeric(apparent_base[metric]),
                   apparent_plus_region = as.numeric(apparent_aug[metric]),
                   apparent_delta = as.numeric(apparent_delta[metric]))]
performance[, optimism_corrected_delta := apparent_delta - mean_optimism]
fwrite(performance, file.path(out, "05_five_year_increment.csv"), bom = TRUE)

status <- c(
  paste0("LINKED_PEOPLE=", flow[stage == "both_FreeSurfer_scans_matched", people]),
  paste0("COMPLETE_PEOPLE=", nrow(dat)),
  paste0("LATER_EVENTS=", sum(dat$event)),
  paste0("SPATIAL_PATTERN_HR_PER_SD=", signif(hr$HR, 4)),
  paste0("SPATIAL_PATTERN_LR_P=", signif(hr$lr_p_value, 4)),
  paste0("SIX_REGION_MEAN_DUPLICATE_MAX_ERROR=", signif(equivalence_error, 4)),
  paste0("BOOTSTRAP_SUCCESS_AUC=", performance[metric == "AUC", successful_bootstraps]),
  "NOTE=The added score projects 34 regional thinning rates onto the outcome-blind 13H1 genetic contrast map; the original six-region mean remains in the benchmark.",
  "NOTE=Risk starts at the later of landmark MRI and baseline diagnosis; event dates and censor dates are inherited from the frozen clinical frame.",

  paste0("OUTPUT_DIR=", out)
)
writeLines(status, file.path(out, "00_status.txt"), useBytes = TRUE)
message(paste(status, collapse = "\n"))
