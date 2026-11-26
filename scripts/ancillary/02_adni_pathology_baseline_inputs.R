# Baseline-linked A/T/V/N preparation for CAT12 inclusion
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

#!/usr/bin/env Rscript

# Step 10B: Prespecified ADNI pathology and clinical-expression models.
# Endpoint-specific samples are used; no outcome is used to tune PGS parameters.

options(stringsAsFactors = FALSE, scipen = 999, width = 180, warn = 1)

PROJECT_DIR <- study_path("project")
PGS_FILE <- file.path(PROJECT_DIR, "results", "10A_ADNI_gMAV_PGS_assembly",
                      "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
ADNI_ROOT <- study_path("cohorts", "ADNI/data")
CLINICAL_FILE <- file.path(ADNI_ROOT, "Clinical/ADNIMERGE_03Jan2026.csv")
FILES <- c(
  amyloid = file.path(ADNI_ROOT, "PET/PET_Amyloid_Berkeley.csv"),
  tau = file.path(ADNI_ROOT, "PET/PET_Tau_Berkeley.csv"),
  wmh = file.path(ADNI_ROOT, "MRI/UCD_WMH_03Jan2026.csv"),
  hippocampus = file.path(ADNI_ROOT, "MRI/MRI_Hippocampus_ASHS.csv")
)
OUTDIR <- file.path(PROJECT_DIR, "results", "10B_ADNI_ATVN_PGS_clinical_expression")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.", call. = FALSE)
if (!requireNamespace("lme4", quietly = TRUE)) stop("lme4 is required for the prespecified mixed models.", call. = FALSE)
need <- c(PGS_FILE, CLINICAL_FILE, FILES)
if (any(!file.exists(need))) stop("Missing input(s):\n", paste(need[!file.exists(need)], collapse = "\n"), call. = FALSE)

LOG <- file.path(OUTDIR, "10B_ADNI_ATVN_PGS_clinical_expression_models.log")
if (file.exists(LOG)) unlink(LOG, force = TRUE)
log_line <- function(level = "INFO", ...) {
  z <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(..., collapse = ""))
  cat(z, "\n"); cat(z, "\n", file = LOG, append = TRUE)
}
as_date <- function(x) as.Date(as.character(x))
zscore <- function(x) {
  x <- as.numeric(x); s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
tidy_lm <- function(fit, endpoint, model_id, n_participants, n_rows) {
  x <- as.data.frame(summary(fit)$coefficients); x$term <- rownames(x); rownames(x) <- NULL
  names(x)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  x$endpoint <- endpoint; x$model_id <- model_id; x$n_participants <- n_participants; x$n_rows <- n_rows
  x[, c("endpoint", "model_id", "term", "estimate", "std_error", "statistic", "p_value", "n_participants", "n_rows")]
}
tidy_lmer <- function(fit, endpoint, outcome, model_id, n_participants, n_rows, random_structure) {
  x <- as.data.frame(summary(fit)$coefficients); x$term <- rownames(x); rownames(x) <- NULL
  names(x)[1:3] <- c("estimate", "std_error", "statistic")
  x$p_value <- 2 * stats::pnorm(-abs(x$statistic))
  x$endpoint <- endpoint; x$outcome <- outcome; x$model_id <- model_id
  x$n_participants <- n_participants; x$n_rows <- n_rows; x$random_structure <- random_structure
  x[, c("endpoint", "outcome", "model_id", "term", "estimate", "std_error", "statistic", "p_value", "n_participants", "n_rows", "random_structure")]
}
fit_mixed <- function(formula_slope, formula_intercept, data) {
  out <- tryCatch(lme4::lmer(formula_slope, data = data, REML = FALSE,
                             control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))),
                  error = function(e) e)
  if (inherits(out, "error") || lme4::isSingular(out, tol = 1e-4)) {
    out2 <- lme4::lmer(formula_intercept, data = data, REML = FALSE,
                       control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))
    return(list(fit = out2, random_structure = "random_intercept_fallback"))
  }
  list(fit = out, random_structure = "random_intercept_and_time_slope")
}

log_line("INFO", "Step 10B started. Primary analysis is European-ancestry ADNI participants with PC1-PC10 adjustment.")
long <- data.table::fread(PGS_FILE, sep = "\t", showProgress = TRUE)
req <- c("PTID", "RID", "EXAMDATE", "time_years", "age_bl", "sex", "education", "dx_bl",
         "primary_european_ancestry", "gMAV_noAPOE_z_eur", "gAD_noAPOE_z_eur", paste0("PC", 1:10),
         "mPACCdigit", "mPACCtrailsB")
if (length(setdiff(req, names(long)))) stop("10A data lacks: ", paste(setdiff(req, names(long)), collapse = ", "), call. = FALSE)
# 10A's longitudinal export intentionally prioritizes visit-level clinical fields.
# Read the fixed baseline APOE4 covariate directly from ADNIMERGE.
clinical_apoe <- data.table::fread(CLINICAL_FILE, select = c("PTID", "VISCODE", "EXAMDATE", "APOE4"), showProgress = TRUE)
clinical_apoe[, `:=`(PTID = as.character(PTID), EXAMDATE = as_date(EXAMDATE), APOE4 = suppressWarnings(as.numeric(APOE4)))]
clinical_apoe <- clinical_apoe[is.finite(APOE4)]
clinical_apoe[, bl_priority := ifelse(tolower(as.character(VISCODE)) %in% c("bl", "sc"), 0L, 1L)]
data.table::setorder(clinical_apoe, PTID, bl_priority, EXAMDATE)
clinical_apoe <- clinical_apoe[, .SD[1L], by = PTID][, .(PTID, APOE4)]
long <- merge(long, clinical_apoe, by = "PTID", all.x = TRUE)
long[, `:=`(PTID = as.character(PTID), EXAMDATE = as_date(EXAMDATE), baseline_date = as.Date(NA))]
long <- long[primary_european_ancestry %in% TRUE & is.finite(gMAV_noAPOE_z_eur) &
               is.finite(gAD_noAPOE_z_eur) & is.finite(APOE4)]
long[, sex := factor(sex)]
long[, dx_bl := factor(dx_bl)]

# One baseline record per participant supplies fixed covariates and the anchor date.
base <- long[order(PTID, time_years, EXAMDATE), .SD[1L], by = PTID]
base <- base[is.finite(age_bl) & is.finite(education) & is.finite(APOE4)]
base[, baseline_date := EXAMDATE]

# Predeclared pathology construction. A/T/V/N are continuous standardized endpoints;
# dates are paired to the nearest baseline assessment within 730 days.
am <- data.table::fread(FILES[["amyloid"]], showProgress = TRUE)
am <- am[, .(PTID = as.character(PTID), endpoint_date = as_date(SCANDATE), value = suppressWarnings(as.numeric(CENTILOIDS)))]
ta <- data.table::fread(FILES[["tau"]], showProgress = TRUE)
ta <- ta[, .(PTID = as.character(PTID), endpoint_date = as_date(SCANDATE), value = suppressWarnings(as.numeric(META_TEMPORAL_SUVR)))]
wm <- data.table::fread(FILES[["wmh"]], showProgress = TRUE)
wm[, total_wmh := suppressWarnings(as.numeric(TOTAL_WMH))]
wm[, total_brain := suppressWarnings(as.numeric(TOTAL_BRAIN))]
wm <- wm[is.finite(total_wmh) & is.finite(total_brain) & total_wmh >= 0 & total_brain > 0,
         .(PTID = as.character(PTID), endpoint_date = as_date(EXAMDATE), value = log1p(total_wmh / total_brain))]
hp <- data.table::fread(FILES[["hippocampus"]], showProgress = TRUE)
hp[, hippo := suppressWarnings(as.numeric(LEFT_HIPP_VOL)) + suppressWarnings(as.numeric(RIGHT_HIPP_VOL))]
hp[, icv := suppressWarnings(as.numeric(ICV))]
hp <- hp[is.finite(hippo) & is.finite(icv) & hippo > 0 & icv > 0,
         .(PTID = as.character(PTID), endpoint_date = as_date(EXAMDATE), value = hippo / icv)]
endpoint_raw <- list(A = am, T = ta, V = wm, N = hp)

pair_baseline <- function(ep, base) {
  z <- merge(ep[is.finite(value) & !is.na(endpoint_date)], base, by = "PTID", all = FALSE, allow.cartesian = TRUE)
  z[, date_gap_days := abs(as.numeric(endpoint_date - baseline_date))]
  z <- z[is.finite(date_gap_days) & date_gap_days <= 730]
  data.table::setorder(z, PTID, date_gap_days, endpoint_date)
  z <- z[, .SD[1L], by = PTID]
  z[, endpoint_z := zscore(value)]
  z
}

paired_sets <- lapply(endpoint_raw, pair_baseline, base = base)
for (nm in names(paired_sets)) data.table::fwrite(paired_sets[[nm]], file.path(OUTDIR, paste0("05_paired_", nm, "_baseline_data.tsv")), sep = "\t", na = "NA")
