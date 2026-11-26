# Diagnosis-conditioned regional maps and paired participant bootstrap
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

project <- study_path("project")
root <- study_path("project", "results")
reference_dir <- file.path(project,
  "results/13H1_aligned_regional_slope_map_release")
model_dir <- file.path(root,
  "13H44_ADNI_diagnosis_time_specification_audit_release")
structural_file <- file.path(root, "13F_R_only_structural_analysis",
  "04_ADNI_R_only_longitudinal_structural_frame.rds")
pgs_file <- file.path(root, "10A_ADNI_gMAV_PGS_assembly",
  "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
fs_file <- study_path("cohorts", "ADNI/data/MRI/MRI_FreeSurfer.csv")
key_file <- study_path("project", "input/analysis_definitions/scan_keys.csv")
map_file <- file.path(REPO_ROOT, "definitions/bilateral_region_mapping.csv")
regional_file <- file.path(reference_dir,
  "01_aligned_34_region_contrasts.csv")
model_frame_file <- file.path(model_dir, "01_analysis_frame_audit.csv")
inputs <- c(structural_file, pgs_file, fs_file, key_file, map_file,
            regional_file, model_frame_file)
if (any(!file.exists(inputs)))
  stop("Missing frozen input(s): ", paste(inputs[!file.exists(inputs)],
       collapse = "; "))

temporal <- c("entorhinal", "fusiform", "inferiortemporal",
              "middletemporal", "temporalpole", "parahippocampal")
pairs <- fread(map_file)
original_regions <- fread(regional_file)
if (nrow(pairs) != 34L || anyDuplicated(pairs$region) ||
    !setequal(pairs$region, original_regions$region) ||
    !setequal(original_regions[temporal_target == TRUE, region], temporal) ||
    any(original_regions$status != "SUCCESS"))
  stop("The frozen 34-region map or six temporal targets changed.")
rest <- setdiff(pairs$region, temporal)

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
first_num <- function(x) {
  y <- num(x)
  if (any(is.finite(y))) y[which(is.finite(y))[1L]] else NA_real_
}
first_chr <- function(x) {
  y <- trimws(as.character(x))
  y <- y[!is.na(y) & nzchar(y)]
  if (length(y)) y[1L] else NA_character_
}
varies <- function(x) {
  if (is.factor(x)) return(nlevels(droplevels(x)) >= 2L)
  length(unique(x[is.finite(num(x))])) >= 2L
}
zscore <- function(x) {
  x <- num(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

keys <- fread(key_file)[in_main_model == TRUE,
  .(RID = as.character(RID), mri_date = as.IDate(mri_date),
    main_imageuid = as.character(main_imageuid))]
if (nrow(keys) != 4659L || uniqueN(keys$RID) != 695L ||
    anyDuplicated(keys, by = c("RID", "mri_date")))
  stop("The frozen 13K scan keys changed.")

pgs <- fread(pgs_file, sep = "\t", showProgress = FALSE)
pgs[, `:=`(RID = as.character(RID),
           exam_date__ = as.IDate(as.character(EXAMDATE)))]
setorder(pgs, RID, exam_date__)
base <- pgs[, .(baseline_age = first_num(age_bl),
                baseline_sex = first_chr(sex),
                baseline_education = first_num(education),
                baseline_diagnosis = first_chr(dx_bl)), by = RID]
d <- as.data.table(readRDS(structural_file))
d[, `:=`(RID = as.character(RID),
         mri_date = as.IDate(as.character(mri_date)),
         gMAV = num(clinical_gMAV_noAPOE_z_eur),
         gAD = num(clinical_gAD_noAPOE_z_eur))]
d <- merge(d, base, by = "RID", all.x = TRUE, sort = FALSE)
d[, `:=`(baseline_age = num(baseline_age),
         sex = factor(baseline_sex),
         education = num(baseline_education),
         diagnosis = factor(baseline_diagnosis),
         scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_")),
         primary_complete = is.finite(num(ad_vulnerable_index)))]
pc_cols <- paste0("clinical_PC", 1:10)
for (v in pc_cols) d[, (v) := num(get(v))]
setorderv(d, c("RID", "mri_date", "primary_complete", "IMAGEUID"),
          c(1L, 1L, -1L, -1L), na.last = TRUE)
d <- d[, .SD[1L], by = .(RID, mri_date)]
setorder(d, RID, mri_date)
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
d <- merge(d, keys, by = c("RID", "mri_date"), sort = FALSE)
if (nrow(d) != nrow(keys) ||
    any(as.character(d$IMAGEUID) != d$main_imageuid))
  stop("Structural Image IDs do not reproduce the frozen scan frame.")

frozen_model <- fread(model_frame_file)
frozen_model[, `:=`(RID = as.character(RID),
                    mri_date = as.IDate(mri_date))]
frame_check <- merge(
  d[, .(RID, mri_date, diagnosis = as.character(diagnosis), gAD, gMAV)],
  frozen_model[, .(RID, mri_date,
                   frozen_diagnosis = as.character(diagnosis),
                   frozen_gAD = gAD, frozen_gMAV = gMAV)],
  by = c("RID", "mri_date"), sort = FALSE)
if (nrow(frame_check) != 4659L ||
    anyNA(frame_check$diagnosis) ||
    any(frame_check$diagnosis != frame_check$frozen_diagnosis) ||
    any(abs(frame_check$gAD - frame_check$frozen_gAD) > 1e-10) ||
    any(abs(frame_check$gMAV - frame_check$frozen_gMAV) > 1e-10))
  stop("The regional frame does not match the frozen 13H44 model frame.")

fields <- unique(c(pairs$left_field, pairs$right_field))
fs <- fread(fs_file,
  select = c("RID", "EXAMDATE", "OVERALLQC", "IMAGEUID", fields),
  showProgress = FALSE)
fs[, `:=`(RID = as.character(RID),
          mri_date = as.IDate(as.character(EXAMDATE)))]
fs <- fs[!is.na(RID) & !is.na(mri_date)]
fs <- fs[!grepl("fail|exclude|poor|unusable",
                 tolower(trimws(as.character(OVERALLQC))))]
setorderv(fs, c("RID", "mri_date", "IMAGEUID"), c(1L, 1L, -1L),
          na.last = TRUE)
fs <- fs[, .SD[1L], by = .(RID, mri_date)]
fs[, fs_imageuid := as.character(IMAGEUID)]
fs[, c("EXAMDATE", "IMAGEUID", "OVERALLQC") := NULL]
for (v in fields) set(fs, j = v, value = num(fs[[v]]))
a <- merge(d, fs, by = c("RID", "mri_date"), sort = FALSE)
if (nrow(a) != nrow(keys) ||
    any(a$main_imageuid != a$fs_imageuid))
  stop("FreeSurfer Image IDs do not match the frozen scans.")

regional_z <- vapply(seq_len(nrow(pairs)), function(i) {
  cols <- c(pairs$left_field[i], pairs$right_field[i])
  zscore(rowMeans(as.matrix(a[, ..cols]), na.rm = TRUE))
}, numeric(nrow(a)))
colnames(regional_z) <- pairs$region
coverage <- data.table(region = pairs$region,
  temporal_target = pairs$region %in% temporal,
  missing_scans = colSums(!is.finite(regional_z)))
candidate_covars <- c("baseline_age", "sex", "education", "diagnosis",
                      "scanner", pc_cols)
covars <- candidate_covars[vapply(candidate_covars,
                                  function(v) varies(a[[v]]), logical(1))]
if (!"diagnosis" %in% covars)
  stop("Baseline diagnosis has no variation in the regional frame.")
rhs_original <- paste(c("time_years * gMAV", "time_years * gAD", covars),
                      collapse = " + ")
rhs_stage <- paste(rhs_original, "time_years * diagnosis", sep = " + ")
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

fit_contrast <- function(dat, outcome, rhs) {
  formula_for <- function(random) as.formula(paste0(
    outcome, " ~ ", rhs,
    if (random == "slope") " + (1 + time_years | RID)"
    else " + (1 | RID)"))
  fit <- tryCatch(lmer(formula_for("slope"), data = dat,
                       REML = FALSE, control = ctrl), error = function(e) e)
  structure <- "random_intercept_and_time_slope"
  if (inherits(fit, "error") || isSingular(fit, tol = 1e-4)) {
    fit <- tryCatch(lmer(formula_for("intercept"), data = dat,
                         REML = FALSE, control = ctrl), error = function(e) e)
    structure <- "random_intercept_fallback"
  }
  if (inherits(fit, "error"))
    stop("Model failed for ", outcome, ": ", conditionMessage(fit))
  b <- fixef(fit)
  V <- as.matrix(vcov(fit))
  terms <- c("time_years:gAD", "time_years:gMAV")
  if (!all(terms %in% names(b)))
    stop("Genetic slope term missing for ", outcome)
  estimate <- unname(b[terms[1L]] - b[terms[2L]])
  variance <- V[terms[1L], terms[1L]] + V[terms[2L], terms[2L]] -
    2 * V[terms[1L], terms[2L]]
  se <- sqrt(max(variance, 0))
  list(estimate = estimate, se = se,
       lower = estimate - 1.96 * se, upper = estimate + 1.96 * se,
       p = 2 * pnorm(-abs(estimate / se)),
       gAD = unname(b[terms[1L]]), gMAV = unname(b[terms[2L]]),
       random_structure = structure)
}

adjusted <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  dat <- copy(a)
  dat[, outcome_z := regional_z[, i]]
  needed <- c("outcome_z", "RID", "time_years", "gMAV", "gAD", covars)
  dat <- droplevels(dat[complete.cases(dat[, ..needed])])
  ref <- original_regions[region == pairs$region[i]]
  if (nrow(dat) != ref$n_rows || uniqueN(dat$RID) != ref$n_people)
    stop("Regional sample changed for ", pairs$region[i])
  fit <- fit_contrast(dat, "outcome_z", rhs_stage)
  data.table(region = pairs$region[i],
    adjusted_estimate = fit$estimate, adjusted_se = fit$se,
    adjusted_lower_95 = fit$lower, adjusted_upper_95 = fit$upper,
    adjusted_p = fit$p, adjusted_gAD = fit$gAD,
    adjusted_gMAV = fit$gMAV, n_people = uniqueN(dat$RID),
    n_scans = nrow(dat), random_structure = fit$random_structure)
}))
adjusted[, adjusted_q := p.adjust(adjusted_p, method = "BH")]
paired_regions <- merge(
  original_regions[, .(region, temporal_target, original_estimate = estimate,
                       original_se = std_error, original_p = p_value,
                       original_q = q_value, original_n_scans = n_rows,
                       original_random_structure = random_structure)],
  adjusted, by = "region", sort = FALSE)
paired_regions[, change_adjusted_minus_original :=
  adjusted_estimate - original_estimate]
setorder(paired_regions, region)

full <- rowSums(is.finite(regional_z)) == 34L
dat <- copy(a[full])
dat[, `:=`(
  temporal_mean_z = rowMeans(regional_z[full, temporal, drop = FALSE]),
  rest_mean_z = rowMeans(regional_z[full, rest, drop = FALSE]))]
dat[, temporal_minus_rest_z := temporal_mean_z - rest_mean_z]
needed <- c("RID", "time_years", "gMAV", "gAD", covars,
            "temporal_mean_z", "rest_mean_z", "temporal_minus_rest_z")
dat <- droplevels(dat[complete.cases(dat[, ..needed])])
if (uniqueN(dat$RID) != 695L || nrow(dat) != 4658L)
  stop("The complete-region paired frame does not reproduce 13H21.")

outcomes <- c("temporal_mean_z", "rest_mean_z", "temporal_minus_rest_z")
specifications <- c(original = rhs_original, stage_adjusted = rhs_stage)
group_fits <- list()
group_rows <- list()
for (spec in names(specifications)) {
  for (outcome in outcomes) {
    fit <- fit_contrast(dat, outcome, specifications[[spec]])
    group_fits[[paste(spec, outcome, sep = "/")]] <- fit
    group_rows[[length(group_rows) + 1L]] <- data.table(
      specification = spec, outcome = outcome,
      estimate = fit$estimate, se = fit$se,
      lower_wald = fit$lower, upper_wald = fit$upper,
      p_wald = fit$p, n_people = uniqueN(dat$RID), n_scans = nrow(dat),
      random_structure = fit$random_structure)
  }
}
group_rows <- rbindlist(group_rows)
out_dir <- file.path(project, "results", paste0(
  "13H50B_ADNI_stage_adjusted_regional_",
  "release"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(coverage, file.path(out_dir, "01_region_coverage.csv"), bom = TRUE)
fwrite(paired_regions,
       file.path(out_dir, "02_paired_34_region_estimates.csv"), bom = TRUE)
fwrite(group_rows,
       file.path(out_dir, "03_paired_temporal_vs_rest.csv"), bom = TRUE)

bootstrap_form <- function(spec) {
  random <- group_fits[[paste(spec, "temporal_minus_rest_z", sep = "/")]]$random_structure
  as.formula(paste0("temporal_minus_rest_z ~ ", specifications[[spec]],
    if (random == "random_intercept_and_time_slope")
      " + (1 + time_years | boot_id)" else " + (1 | boot_id)"))
}
forms <- lapply(names(specifications), bootstrap_form)
names(forms) <- names(specifications)
boot_fit <- function(form, sample_data) {
  fit <- suppressWarnings(tryCatch(lmer(form, data = sample_data,
    REML = FALSE, control = ctrl), error = function(e) e))
  if (inherits(fit, "error"))
    return(list(estimate = NA_real_, singular = NA,
                status = conditionMessage(fit)))
  b <- fixef(fit)
  terms <- c("time_years:gAD", "time_years:gMAV")
  if (!all(terms %in% names(b)))
    return(list(estimate = NA_real_, singular = NA,
                status = "MISSING_GENETIC_SLOPE"))
  list(estimate = unname(b[terms[1L]] - b[terms[2L]]),
       singular = isSingular(fit, tol = 1e-4), status = "SUCCESS")
}

set.seed(13021L)
ids <- unique(dat$RID)
groups <- split(seq_len(nrow(dat)), dat$RID)
draws <- rbindlist(lapply(seq_len(500L), function(i) {
  selected <- sample(ids, length(ids), replace = TRUE)
  selected_groups <- groups[selected]
  boot <- copy(dat[unlist(selected_groups, use.names = FALSE)])
  boot[, boot_id := factor(rep(seq_along(selected), lengths(selected_groups)))]
  orig <- boot_fit(forms$original, boot)
  adj <- boot_fit(forms$stage_adjusted, boot)
  if (i %% 50L == 0L) message("PAIRED_BOOTSTRAP=", i, "/500")
  data.table(draw = i, original_estimate = orig$estimate,
    adjusted_estimate = adj$estimate,
    change_adjusted_minus_original = adj$estimate - orig$estimate,
    original_singular = orig$singular, adjusted_singular = adj$singular,
    original_status = orig$status, adjusted_status = adj$status)
}))
fwrite(draws, file.path(out_dir, "04_paired_subject_bootstrap.csv"),
       bom = TRUE)
success <- draws[original_status == "SUCCESS" &
                 adjusted_status == "SUCCESS" &
                 is.finite(change_adjusted_minus_original)]
if (nrow(success) < 450L)
  stop("Fewer than 450 paired bootstrap fits succeeded. See ", out_dir)

intervals <- rbindlist(lapply(
  c("original_estimate", "adjusted_estimate",
    "change_adjusted_minus_original"), function(v) {
  ci <- quantile(success[[v]], c(0.025, 0.975), names = FALSE)
  data.table(quantity = v, lower_95 = ci[1L], upper_95 = ci[2L],
             successful_draws = nrow(success))
}))
fwrite(intervals, file.path(out_dir, "05_paired_bootstrap_intervals.csv"),
       bom = TRUE)
status <- c("STATUS=COMPLETE",
  "MODEL_DIAGNOSIS_FIELD=dx_bl",
  "REGIONAL_OUTCOME=bilateral_FreeSurfer_thickness_zscore",
  "MAP_VALUES=bilateral_region_estimates_for_both_hemispheres",
  paste0("PEOPLE=", uniqueN(a$RID)),
  paste0("SCAN_KEYS=", nrow(a)),
  paste0("COMMON_COMPLETE_SCANS=", nrow(dat)),
  paste0("REGIONS=", nrow(paired_regions)),
  paste0("PAIRED_BOOTSTRAP_SUCCESS=", nrow(success), "/500"),
  paste0("OUTPUT_DIR=", out_dir))
writeLines(status, file.path(out_dir, "00_status.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
message(paste(status, collapse = "\n"))
