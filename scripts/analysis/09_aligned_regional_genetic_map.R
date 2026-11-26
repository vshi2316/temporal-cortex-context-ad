# Aligned 34-region gAD and gMAV slopes
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

project <- study_path("project")
main_results <- study_path("project", "results")
fs_file <- study_path("cohorts", "ADNI/data/MRI/MRI_FreeSurfer.csv")
structural_file <- file.path(main_results, "13F_R_only_structural_analysis", "04_ADNI_R_only_longitudinal_structural_frame.rds")
pgs_file <- file.path(main_results, "10A_ADNI_gMAV_PGS_assembly", "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
key_file <- study_path("project", "input/analysis_definitions/scan_keys.csv")
map_file <- file.path(REPO_ROOT, "definitions/bilateral_region_mapping.csv")
required <- c(fs_file, structural_file, pgs_file, key_file, map_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing frozen inputs: ", paste(missing, collapse = "; "))

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
  y <- num(x)
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(y)))
  (y - mean(y, na.rm = TRUE)) / s
}

keys <- fread(key_file)[in_main_model == TRUE,
                        .(RID = as.character(RID), mri_date = as.IDate(mri_date),
                          main_imageuid = as.character(main_imageuid))]
if (nrow(keys) != 4659L || uniqueN(keys$RID) != 695L ||
    anyDuplicated(keys, by = c("RID", "mri_date")))
  stop("The frozen 13K scan keys are not 695 people and 4659 unique scans.")

pgs <- fread(pgs_file, sep = "\t", showProgress = FALSE)
pgs[, `:=`(RID = as.character(RID), exam_date__ = as.IDate(as.character(EXAMDATE)))]
setorder(pgs, RID, exam_date__)
base <- pgs[, .(baseline_age = first_num(age_bl),
                baseline_sex = first_chr(sex),
                baseline_education = first_num(education),
                baseline_diagnosis = first_chr(dx_bl)), by = RID]
d <- as.data.table(readRDS(structural_file))
d[, `:=`(RID = as.character(RID), mri_date = as.IDate(as.character(mri_date)),
         gMAV = num(clinical_gMAV_noAPOE_z_eur),
         gAD = num(clinical_gAD_noAPOE_z_eur))]
d <- merge(d, base, by = "RID", all.x = TRUE, sort = FALSE)
d[, `:=`(baseline_age = num(baseline_age), sex = factor(baseline_sex),
         education = num(baseline_education), diagnosis = factor(baseline_diagnosis),
         scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_")),
         primary_complete = is.finite(num(ad_vulnerable_index)))]
pc_cols <- paste0("clinical_PC", 1:10)
for (v in pc_cols) d[, (v) := num(get(v))]
setorderv(d, c("RID", "mri_date", "primary_complete", "IMAGEUID"),
          order = c(1L, 1L, -1L, -1L), na.last = TRUE)
d <- d[, .SD[1L], by = .(RID, mri_date)]
setorder(d, RID, mri_date)
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
d <- merge(d, keys, by = c("RID", "mri_date"), sort = FALSE)
if (nrow(d) != nrow(keys) || any(as.character(d$IMAGEUID) != d$main_imageuid))
  stop("Structural records do not reproduce the frozen 13K Image IDs.")

pairs <- fread(map_file)
if (nrow(pairs) != 34L || anyDuplicated(pairs$region))
  stop("The frozen bilateral cortical map is not 34 unique regions.")
fields <- unique(c(pairs$left_field, pairs$right_field))
fs <- fread(fs_file, select = c("RID", "EXAMDATE", "OVERALLQC", "IMAGEUID", fields),
            showProgress = FALSE)
fs[, `:=`(RID = as.character(RID), mri_date = as.IDate(as.character(EXAMDATE)))]
fs <- fs[!is.na(RID) & !is.na(mri_date)]
fs <- fs[!grepl("fail|exclude|poor|unusable", tolower(trimws(as.character(OVERALLQC))))]
setorderv(fs, c("RID", "mri_date", "IMAGEUID"),
          order = c(1L, 1L, -1L), na.last = TRUE)
fs <- fs[, .SD[1L], by = .(RID, mri_date)]
fs[, fs_imageuid := as.character(IMAGEUID)]
fs[, c("EXAMDATE", "IMAGEUID", "OVERALLQC") := NULL]
for (v in fields) set(fs, j = v, value = num(fs[[v]]))
a <- merge(d, fs, by = c("RID", "mri_date"), sort = FALSE)
if (nrow(a) != nrow(keys) || any(a$main_imageuid != a$fs_imageuid))
  stop("FreeSurfer records do not match all frozen 13K scan keys and Image IDs.")

candidate_covars <- c("baseline_age", "sex", "education", "diagnosis",
                      "scanner", pc_cols)
covars <- candidate_covars[vapply(candidate_covars,
                                  function(v) varies(a[[v]]), logical(1))]
rhs <- paste(c("time_years * gMAV", "time_years * gAD", covars), collapse = " + ")

fit_region <- function(region, left, right) {
  dat <- copy(a)
  dat[, thickness := rowMeans(.SD, na.rm = TRUE), .SDcols = c(left, right)]
  dat[, outcome_z := zscore(thickness)]
  vars <- c("outcome_z", "RID", "time_years", "gMAV", "gAD", covars)
  dat <- droplevels(dat[complete.cases(dat[, ..vars])])
  result <- data.table(region = region, estimate = NA_real_, std_error = NA_real_,
                       p_value = NA_real_, gAD_estimate = NA_real_,
                       gAD_se = NA_real_, gAD_p = NA_real_,
                       gMAV_estimate = NA_real_, gMAV_se = NA_real_,
                       gMAV_p = NA_real_, n_people = uniqueN(dat$RID),
                       n_rows = nrow(dat),
                       n_repeated = dat[, .N, by = RID][N >= 2L, .N],
                       random_structure = NA_character_, status = "NOT_RUN")
  if (result$n_people < 100L || result$n_rows < 200L || result$n_repeated < 50L)
    return(result)
  slope <- as.formula(paste0("outcome_z ~ ", rhs, " + (1 + time_years | RID)"))
  intercept <- as.formula(paste0("outcome_z ~ ", rhs, " + (1 | RID)"))
  fit <- tryCatch(lmer(slope, data = dat, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 200000))),
                  error = function(e) e)
  structure <- "random_intercept_and_time_slope"
  if (inherits(fit, "error") || isSingular(fit, tol = 1e-4)) {
    fit <- tryCatch(lmer(intercept, data = dat, REML = FALSE,
                         control = lmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 200000))),
                    error = function(e) e)
    structure <- "random_intercept_fallback"
  }
  result[, random_structure := structure]
  if (inherits(fit, "error")) {
    result[, status := paste0("FAILED: ", conditionMessage(fit))]
    return(result)
  }
  b <- fixef(fit)
  V <- as.matrix(vcov(fit))
  if (!all(c("time_years:gAD", "time_years:gMAV") %in% names(b))) {
    result[, status := "MISSING_GENETIC_INTERACTION"]
    return(result)
  }
  beta <- b[["time_years:gAD"]] - b[["time_years:gMAV"]]
  variance <- V["time_years:gAD", "time_years:gAD"] +
    V["time_years:gMAV", "time_years:gMAV"] -
    2 * V["time_years:gAD", "time_years:gMAV"]
  se_hat <- sqrt(max(variance, 0))
  se_ad <- sqrt(V["time_years:gAD", "time_years:gAD"])
  se_mav <- sqrt(V["time_years:gMAV", "time_years:gMAV"])
  result[, `:=`(estimate = beta, std_error = se_hat,
                p_value = 2 * pnorm(-abs(beta / se_hat)),
                gAD_estimate = b[["time_years:gAD"]], gAD_se = se_ad,
                gAD_p = 2 * pnorm(-abs(b[["time_years:gAD"]] / se_ad)),
                gMAV_estimate = b[["time_years:gMAV"]], gMAV_se = se_mav,
                gMAV_p = 2 * pnorm(-abs(b[["time_years:gMAV"]] / se_mav)),
                status = "SUCCESS")]
  result
}

results <- rbindlist(lapply(seq_len(nrow(pairs)), function(i)
  fit_region(pairs$region[i], pairs$left_field[i], pairs$right_field[i])))
results[, `:=`(lower_95 = estimate - 1.96 * std_error,
               upper_95 = estimate + 1.96 * std_error,
               q_value = p.adjust(p_value, method = "BH"),
               gAD_q = p.adjust(gAD_p, method = "BH"),
               gMAV_q = p.adjust(gMAV_p, method = "BH"),
               temporal_target = region %in% c("entorhinal", "fusiform",
                   "inferiortemporal", "middletemporal", "temporalpole",
                   "parahippocampal"))]
out <- file.path(project, "results", paste0("13H1_aligned_regional_slope_map_",
                                           "release"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(results, file.path(out, "01_aligned_34_region_contrasts.csv"), bom = TRUE)
status <- c(
  "FROZEN_MAIN_FRAME_MATCHED=TRUE",
  paste0("PEOPLE=", uniqueN(a$RID)),
  paste0("SCAN_KEYS=", nrow(a)),
  paste0("IMAGE_ID_MISMATCHES=", sum(a$main_imageuid != a$fs_imageuid)),
  paste0("SUCCESS_REGIONS=", sum(results$status == "SUCCESS")),
  paste0("TEMPORAL_REGIONS=", sum(results$temporal_target)),
  "SEPARATE_gAD_gMAV_SLOPES=TRUE",
  paste0("gAD_FDR_REGIONS=", sum(results$gAD_q < 0.05, na.rm = TRUE)),
  paste0("gMAV_FDR_REGIONS=", sum(results$gMAV_q < 0.05, na.rm = TRUE)),
  "NOTE=Regional FreeSurfer thickness uses the fixed participant and scan keys.",

  paste0("OUTPUT_DIR=", out)
)
writeLines(status, file.path(out, "00_status.txt"), useBytes = TRUE)
message(paste(status, collapse = "\n"))
