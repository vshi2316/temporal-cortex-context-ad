# Endpoint-indexed A/T/V/N preparation for CAT12 models
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

#!/usr/bin/env Rscript

# Step 10C: Endpoint-anchored ADNI A/T/V/N clinical-expression analysis.
# Each participant's earliest available endpoint measurement is the fixed index.
# Cognition is followed from that index forward; no cognitive result selects the index.

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
OUTDIR <- file.path(PROJECT_DIR, "results", "10C_ADNI_endpoint_anchored_ATVN")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.", call. = FALSE)
if (!requireNamespace("lme4", quietly = TRUE)) stop("lme4 is required.", call. = FALSE)
need <- c(PGS_FILE, CLINICAL_FILE, FILES)
if (any(!file.exists(need))) stop("Missing input(s):\n", paste(need[!file.exists(need)], collapse = "\n"), call. = FALSE)

LOG <- file.path(OUTDIR, "10C_ADNI_endpoint_anchored_ATVN_clinical_expression.log")
if (file.exists(LOG)) unlink(LOG, force = TRUE)
log_line <- function(level = "INFO", ...) {
  z <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(..., collapse = ""))
  cat(z, "\n"); cat(z, "\n", file = LOG, append = TRUE)
}
as_date <- function(x) suppressWarnings(as.Date(as.character(x)))
zscore <- function(x) {
  x <- as.numeric(x); s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}
tidy_lm <- function(fit, endpoint, model_id, n_participants) {
  x <- as.data.frame(summary(fit)$coefficients); x$term <- rownames(x); rownames(x) <- NULL
  names(x)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  x$endpoint <- endpoint; x$model_id <- model_id; x$n_participants <- n_participants
  x[, c("endpoint", "model_id", "term", "estimate", "std_error", "statistic", "p_value", "n_participants")]
}
tidy_lmer <- function(fit, endpoint, outcome, random_structure, n_id, n_row) {
  x <- as.data.frame(summary(fit)$coefficients); x$term <- rownames(x); rownames(x) <- NULL
  names(x)[1:3] <- c("estimate", "std_error", "statistic")
  x$p_value <- 2 * stats::pnorm(-abs(x$statistic))
  x$endpoint <- endpoint; x$outcome <- outcome; x$random_structure <- random_structure
  x$n_participants <- n_id; x$n_rows <- n_row
  x[, c("endpoint", "outcome", "term", "estimate", "std_error", "statistic", "p_value", "n_participants", "n_rows", "random_structure")]
}
fit_mixed <- function(f_slope, f_int, d) {
  x <- tryCatch(lme4::lmer(f_slope, data = d, REML = FALSE,
                           control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))),
                error = function(e) e)
  if (inherits(x, "error") || lme4::isSingular(x, tol = 1e-4)) {
    y <- lme4::lmer(f_int, data = d, REML = FALSE,
                    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))
    return(list(fit = y, random_structure = "random_intercept_fallback"))
  }
  list(fit = x, random_structure = "random_intercept_and_time_slope")
}

log_line("INFO", "Step 10C started. Endpoint index = earliest available endpoint measurement per participant.")
long <- data.table::fread(PGS_FILE, sep = "\t", showProgress = TRUE)
need_long <- c("PTID", "RID", "EXAMDATE", "time_years", "age_bl", "sex", "education", "dx_bl",
               "primary_european_ancestry", "gMAV_noAPOE_z_eur", "gAD_noAPOE_z_eur", paste0("PC", 1:10),
               "mPACCdigit", "mPACCtrailsB")
if (length(setdiff(need_long, names(long)))) stop("10A file lacks: ", paste(setdiff(need_long, names(long)), collapse = ", "), call. = FALSE)
long[, `:=`(PTID = as.character(PTID), EXAMDATE = as_date(EXAMDATE))]

# Recover baseline APOE4 from the original clinical source; this avoids dependence
# on column suffixes generated while producing the 10A longitudinal export.
apoe <- data.table::fread(CLINICAL_FILE, select = c("PTID", "VISCODE", "EXAMDATE", "APOE4"), showProgress = FALSE)
apoe[, `:=`(PTID = as.character(PTID), EXAMDATE = as_date(EXAMDATE), APOE4 = suppressWarnings(as.numeric(APOE4)))]
apoe <- apoe[is.finite(APOE4)]
apoe[, bl_priority := ifelse(tolower(as.character(VISCODE)) %in% c("bl", "sc"), 0L, 1L)]
data.table::setorder(apoe, PTID, bl_priority, EXAMDATE)
apoe <- apoe[, .SD[1L], by = PTID][, .(PTID, APOE4)]
long <- merge(long, apoe, by = "PTID", all.x = TRUE)
long <- long[primary_european_ancestry %in% TRUE & is.finite(gMAV_noAPOE_z_eur) &
               is.finite(gAD_noAPOE_z_eur) & is.finite(APOE4) & !is.na(EXAMDATE)]
long[, `:=`(sex = factor(sex), dx_bl = factor(dx_bl))]

# Endpoint definitions are frozen before model fitting.
am <- data.table::fread(FILES[["amyloid"]], showProgress = TRUE)
am <- am[, .(PTID = as.character(PTID), endpoint_date = as_date(SCANDATE), endpoint_value = suppressWarnings(as.numeric(CENTILOIDS)))]
ta <- data.table::fread(FILES[["tau"]], showProgress = TRUE)
ta <- ta[, .(PTID = as.character(PTID), endpoint_date = as_date(SCANDATE), endpoint_value = suppressWarnings(as.numeric(META_TEMPORAL_SUVR)))]
wm <- data.table::fread(FILES[["wmh"]], showProgress = TRUE)
wm[, endpoint_value := log1p(suppressWarnings(as.numeric(TOTAL_WMH)) / suppressWarnings(as.numeric(TOTAL_BRAIN)))]
wm <- wm[, .(PTID = as.character(PTID), endpoint_date = as_date(EXAMDATE), endpoint_value)]
hp <- data.table::fread(FILES[["hippocampus"]], showProgress = TRUE)
hp[, endpoint_value := (suppressWarnings(as.numeric(LEFT_HIPP_VOL)) + suppressWarnings(as.numeric(RIGHT_HIPP_VOL))) /
     suppressWarnings(as.numeric(ICV))]
hp <- hp[, .(PTID = as.character(PTID), endpoint_date = as_date(EXAMDATE), endpoint_value)]
endpoint_raw <- list(A = am, T = ta, V = wm, N = hp)

# Fixed participant-level covariates are taken once, before pathology linkage.
person <- long[order(PTID, EXAMDATE), .SD[1L], by = PTID]
person <- person[, c("PTID", "age_bl", "sex", "education", "dx_bl", "APOE4", "gMAV_noAPOE_z_eur", "gAD_noAPOE_z_eur", paste0("PC", 1:10)), with = FALSE]

index_sets <- list()
for (ep_name in names(endpoint_raw)) {
  log_line("INFO", "Constructing endpoint-anchored index for ", ep_name)
  ep <- endpoint_raw[[ep_name]][is.finite(endpoint_value) & !is.na(endpoint_date)]
  # Do not select on cognition: use the earliest valid endpoint record per participant.
  data.table::setorder(ep, PTID, endpoint_date)
  ep <- ep[, .SD[1L], by = PTID]
  idx <- merge(ep, person, by = "PTID", all = FALSE)
  idx[, endpoint_z := zscore(endpoint_value)]
  idx <- idx[is.finite(endpoint_z)]
  index_sets[[ep_name]] <- idx
}
for (nm in names(index_sets)) data.table::fwrite(index_sets[[nm]], file.path(OUTDIR, paste0("05_", nm, "_endpoint_index_data.tsv")), sep = "\t", na = "NA")
