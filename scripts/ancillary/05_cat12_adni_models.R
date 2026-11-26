# CAT12 volume, vascular/pathology and cognitive models
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

# CAT12 temporal volume models in the 402-participant sample.

options(stringsAsFactors = FALSE, scipen = 999, width = 180, warn = 1)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.", call. = FALSE)
}
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required for longitudinal mixed models.", call. = FALSE)
}
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

PROJECT_ROOT <- study_path("project")
CAT12_ROOT <- file.path(PROJECT_ROOT, "results", "13F_ADNI_CAT12_bridge")
CAT12_REPORT <- file.path(CAT12_ROOT, "tables", "01_CAT12_report_metrics.csv")
CAT12_ROI <- file.path(CAT12_ROOT, "tables", "02_CAT12_roi_values.csv")
CAT12_JOBS <- file.path(CAT12_ROOT, "tables", "03_CAT12_job_manifest.csv")
BRIDGE_MANIFEST <- study_path("cat12_config", "02_ADNI_CAT12_bridge_baseline_manifest.csv")
PGS_FILE <- file.path(
  PROJECT_ROOT, "results", "10A_ADNI_gMAV_PGS_assembly",
  "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"
)
ATVN_ROOT <- file.path(PROJECT_ROOT, "results", "10C_ADNI_endpoint_anchored_ATVN")
OUTDIR <- file.path(
  PROJECT_ROOT, "results", "13F_ADNI_CAT12_mechanistic_bridge_analysis"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

INPUTS <- c(CAT12_REPORT, CAT12_ROI, CAT12_JOBS, BRIDGE_MANIFEST, PGS_FILE)
if (any(!file.exists(INPUTS))) {
  stop("Missing required input(s):\n",
       paste(INPUTS[!file.exists(INPUTS)], collapse = "\n"), call. = FALSE)
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
as_date <- function(x) as.IDate(as.character(x))
norm_id <- function(x) {
  z <- sub("^.*?([0-9]+)$", "\\1", as.character(x))
  z <- sub("^0+", "", z)
  z[z == ""] <- "0"
  z
}
safe_z <- function(x) {
  x <- as_num(x)
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(m) || !is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - m) / s
}
required_columns <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) {
    stop(label, " is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
}
write_empty <- function(path, cols) {
  z <- as.data.table(setNames(
    replicate(length(cols), character(), simplify = FALSE), cols
  ))
  fwrite(z, path, bom = TRUE)
}
tidy_lm <- function(fit, model_id, outcome, endpoint = NA_character_,
                    n_people, n_rows) {
  sm <- as.data.frame(summary(fit)$coefficients)
  sm$term <- rownames(sm)
  rownames(sm) <- NULL
  names(sm)[seq_len(4L)] <- c("estimate", "std_error", "statistic", "p_value")
  out <- as.data.table(sm)
  out[, c("endpoint", "outcome", "model_id", "n_people", "n_rows") :=
        list(endpoint, outcome, model_id, n_people, n_rows)]
  crit <- stats::qt(0.975, df = max(1, fit$df.residual))
  out[, c("conf_low", "conf_high") :=
        list(estimate - crit * std_error, estimate + crit * std_error)]
  out[, .(endpoint, outcome, model_id, term, estimate, std_error, statistic,
          p_value, conf_low, conf_high, n_people, n_rows)]
}
tidy_lmer <- function(fit, model_id, outcome, n_people, n_rows,
                      random_structure) {
  sm <- as.data.frame(summary(fit)$coefficients)
  sm$term <- rownames(sm)
  rownames(sm) <- NULL
  names(sm)[seq_len(3L)] <- c("estimate", "std_error", "statistic")
  out <- as.data.table(sm)
  out[, c("endpoint", "outcome", "model_id", "p_value", "n_people",
          "n_rows", "random_structure") :=
        list(NA_character_, outcome, model_id,
             2 * stats::pnorm(-abs(statistic)), n_people, n_rows,
             random_structure)]
  out[, .(endpoint, outcome, model_id, term, estimate, std_error, statistic,
          p_value, n_people, n_rows, random_structure)]
}
fit_mixed <- function(formula_slope, formula_intercept, data) {
  fit1 <- tryCatch(
    lme4::lmer(
      formula_slope, data = data, REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
      )
    ),
    error = function(e) e
  )
  if (!inherits(fit1, "error") && !lme4::isSingular(fit1, tol = 1e-4)) {
    return(list(fit = fit1,
                random_structure = "random_intercept_and_time_slope"))
  }
  fit2 <- tryCatch(
    lme4::lmer(
      formula_intercept, data = data, REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
      )
    ),
    error = function(e) e
  )
  if (inherits(fit2, "error")) {
    msg <- if (inherits(fit1, "error")) fit1$message else fit2$message
    stop("Both mixed-model specifications failed: ", msg, call. = FALSE)
  }
  list(fit = fit2, random_structure = "random_intercept_fallback")
}

# ---------------------------------------------------------------------------
# 1. Frozen bridge manifest and 10A longitudinal PGS/cognition
# ---------------------------------------------------------------------------
bridge_manifest <- fread(BRIDGE_MANIFEST, showProgress = FALSE)
if (!"PTID" %in% names(bridge_manifest) &&
    "participant_id" %in% names(bridge_manifest)) {
  bridge_manifest[, PTID := as.character(participant_id)]
}
required_columns(
  bridge_manifest,
  c("PTID", "RID", "cohort", "role", "processing_image_code",
    "image_id", "analysis_timepoint"),
  "Bridge manifest"
)
bridge_manifest[, c("PTID", "RID", "cohort", "role") :=
  list(as.character(PTID), norm_id(RID), as.character(cohort), as.character(role))]
bridge_manifest <- bridge_manifest[cohort == "ADNI" & role == "PRIMARY"]
if (!nrow(bridge_manifest) || anyDuplicated(bridge_manifest$PTID)) {
  stop("Bridge manifest is empty or contains duplicate PTIDs.", call. = FALSE)
}
if (nrow(bridge_manifest) != 402L) {
  stop("Frozen ADNI bridge manifest must contain 402 participants; found ",
       nrow(bridge_manifest), call. = FALSE)
}

pgs <- fread(PGS_FILE, sep = "\t", showProgress = FALSE)
required_columns(
  pgs,
  c("PTID", "EXAMDATE", "time_years", "age_bl", "sex", "education",
    "dx_bl", "primary_european_ancestry",
    "gMAV_noAPOE_z_eur", "gAD_noAPOE_z_eur",
    paste0("PC", 1:10), "MMSE", "mPACCdigit", "mPACCtrailsB",
    "ADAS13", "CDRSB"),
  "10A longitudinal PGS export"
)
pgs[, c("PTID", "EXAMDATE", "time_years", "age_bl", "education",
        "primary_european_ancestry", "gMAV_noAPOE_z_eur",
        "gAD_noAPOE_z_eur") :=
  list(as.character(PTID), as_date(EXAMDATE), as_num(time_years),
       as_num(age_bl), as_num(education),
       toupper(as.character(primary_european_ancestry)) == "TRUE",
       as_num(gMAV_noAPOE_z_eur), as_num(gAD_noAPOE_z_eur))]
for (cc in c("MMSE", "mPACCdigit", "mPACCtrailsB", "ADAS13", "CDRSB",
             paste0("PC", 1:10))) {
  pgs[, (cc) := as_num(get(cc))]
}
pgs <- pgs[
  primary_european_ancestry == TRUE &
    is.finite(gMAV_noAPOE_z_eur) & is.finite(gAD_noAPOE_z_eur)
]
setorder(pgs, PTID, time_years, EXAMDATE)
pgs_base <- pgs[, .SD[1L], by = PTID]

# ---------------------------------------------------------------------------
# 2. CAT12 structural phenotype
# ---------------------------------------------------------------------------
report <- fread(CAT12_REPORT, showProgress = FALSE)
roi <- fread(CAT12_ROI, showProgress = FALSE)
jobs <- fread(CAT12_JOBS, showProgress = FALSE)
required_columns(
  report,
  c("job_id", "RID", "image_id", "tiv_ml", "gm_ml", "wm_ml",
    "iqr", "siqr", "ncr", "icr", "wmh_abs_ml"),
  "CAT12 report"
)
required_columns(roi, c("job_id", "atlas", "roi_name", "Vgm_ml"),
                  "CAT12 ROI table")
required_columns(jobs, c("job_id", "status"), "CAT12 job manifest")
if (nrow(jobs) != 402L ||
    any(!jobs$status %in% c("SUCCESS", "ALREADY_COMPLETE"))) {
  stop("CAT12 job manifest does not show 402 complete CAT12 jobs.", call. = FALSE)
}
if ("participant_id" %in% names(report)) {
  report[, PTID := as.character(participant_id)]
}
if ("participant_id" %in% names(roi)) {
  roi[, PTID := as.character(participant_id)]
}
required_columns(report, "PTID", "CAT12 report")
required_columns(roi, "PTID", "CAT12 ROI table")

report[, c("PTID", "RID", "tiv_ml", "gm_ml", "wm_ml", "iqr", "siqr",
           "ncr", "icr", "wmh_abs_ml") :=
  list(as.character(PTID), norm_id(RID), as_num(tiv_ml), as_num(gm_ml),
       as_num(wm_ml), as_num(iqr), as_num(siqr), as_num(ncr),
       as_num(icr), as_num(wmh_abs_ml))]
if (anyDuplicated(report$PTID)) {
  stop("CAT12 report contains duplicate PTIDs.", call. = FALSE)
}
roi[, c("PTID", "Vgm_ml") := list(as.character(PTID), as_num(Vgm_ml))]

targets <- c(
  "Right Ent entorhinal area", "Left Ent entorhinal area",
  "Right FuG fusiform gyrus", "Left FuG fusiform gyrus",
  "Right ITG inferior temporal gyrus", "Left ITG inferior temporal gyrus",
  "Right MTG middle temporal gyrus", "Left MTG middle temporal gyrus",
  "Right TMP temporal pole", "Left TMP temporal pole",
  "Right PHG parahippocampal gyrus", "Left PHG parahippocampal gyrus"
)
roi_target <- roi[
  atlas == "neuromorphometrics" &
    roi_name %in% targets &
    is.finite(Vgm_ml)
]
temporal <- roi_target[, .(
  cat12_temporal_gm_ml = sum(Vgm_ml, na.rm = TRUE),
  n_target_rois = uniqueN(roi_name)
), by = PTID]
if (nrow(temporal) != 402L || any(temporal$n_target_rois != length(targets))) {
  stop("CAT12 temporal ROI extraction is incomplete; expected 12 ROIs for all 402 participants.",
       call. = FALSE)
}
cat12 <- merge(
  report[, .(PTID, RID, image_id, analysis_timepoint, report_xml,
             iqr, siqr, ncr, icr, csf_ml, gm_ml, wm_ml, wmh_ml,
             tiv_ml, wmh_abs_ml)],
  temporal, by = "PTID", all = FALSE
)
cat12[, c("cat12_global_gm_tiv_ratio", "cat12_wm_tiv_ratio",
          "cat12_temporal_gm_tiv_ratio") :=
  list(gm_ml / tiv_ml, wm_ml / tiv_ml,
       cat12_temporal_gm_ml / tiv_ml)]
cat12[, cat12_ad_vulnerable_gm :=
  -safe_z(log(cat12_temporal_gm_tiv_ratio))]
cat12 <- cat12[
  is.finite(cat12_temporal_gm_tiv_ratio) &
    cat12_temporal_gm_tiv_ratio > 0 &
    is.finite(cat12_ad_vulnerable_gm)
]
if (!nrow(cat12)) stop("No valid CAT12 structural phenotype.", call. = FALSE)

bridge_meta <- bridge_manifest[, .(
  PTID,
  bridge_RID = norm_id(RID),
  t1_anchor_date = if ("baseline_date" %in% names(bridge_manifest))
    as_date(baseline_date) else as.IDate(NA),
  APOE4 = if ("APOE4" %in% names(bridge_manifest)) as_num(APOE4) else NA_real_,
  has_ATVN = if ("has_ATVN" %in% names(bridge_manifest)) as.logical(has_ATVN) else NA,
  eligible_bridge = if ("eligible_bridge" %in% names(bridge_manifest)) as.logical(eligible_bridge) else NA
)]
base <- merge(cat12, bridge_meta, by = "PTID", all = FALSE)
base <- merge(
  base,
  pgs_base[, .(
    PTID, age_bl, sex, education, dx_bl, primary_european_ancestry,
    gMAV_noAPOE_z_eur, gAD_noAPOE_z_eur,
    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10
  )],
  by = "PTID", all = FALSE
)
base[, c("sex", "dx_bl", "APOE4", "age_bl", "education") :=
  list(factor(sex), factor(dx_bl), as_num(APOE4),
       as_num(age_bl), as_num(education))]
if (anyDuplicated(base$PTID)) {
  stop("CAT12/PGS bridge contains duplicate PTIDs.", call. = FALSE)
}
if (nrow(base) < 100L) {
  stop("Fewer than 100 complete CAT12/PGS bridge participants.", call. = FALSE)
}
if (nrow(base) != 402L) {
  stop("Complete CAT12/PGS bridge must contain 402 participants; found ",
       nrow(base), call. = FALSE)
}

# ---------------------------------------------------------------------------
# 3. Endpoint-anchored A/T/V/N from 10C
# ---------------------------------------------------------------------------
read_endpoint <- function(ep) {
  f <- file.path(ATVN_ROOT, paste0("05_", ep, "_endpoint_index_data.tsv"))
  if (!file.exists(f) || file.info(f)$size <= 400) return(NULL)
  d <- fread(f, sep = "\t", showProgress = FALSE)
  if (!nrow(d)) return(NULL)
  required_columns(
    d, c("PTID", "endpoint_date", "endpoint_value"),
    paste0("10C ", ep, " endpoint file")
  )
  d[, c("PTID", "endpoint_date", "endpoint_value") :=
    list(as.character(PTID), as_date(endpoint_date), as_num(endpoint_value))]
  if ("endpoint_z" %in% names(d)) {
    d[, endpoint_z := as_num(endpoint_z)]
  } else {
    d[, endpoint_z := safe_z(endpoint_value)]
  }
  d <- d[is.finite(endpoint_value) & is.finite(endpoint_z)]
  setorder(d, PTID, endpoint_date)
  d <- d[, .SD[1L], by = PTID]
  d[, endpoint := ep]
  d
}
endpoint_names <- c("A", "T", "V", "N")
endpoint_raw <- setNames(lapply(endpoint_names, read_endpoint), endpoint_names)
endpoint_flow <- rbindlist(lapply(endpoint_names, function(ep) {
  d <- endpoint_raw[[ep]]
  if (is.null(d)) {
    return(data.table(endpoint = ep, source = "10C_endpoint_index",
                      n_source_people = 0L, n_bridge_people = 0L,
                      status = "NOT_AVAILABLE"))
  }
  z <- merge(d[, .(PTID, endpoint_date)],
             base[, .(PTID, t1_anchor_date)], by = "PTID", all = FALSE)
  z[, absolute_t1_endpoint_gap_days := abs(as.numeric(endpoint_date - t1_anchor_date))]
  data.table(
    endpoint = ep,
    source = "10C_endpoint_index",
    n_source_people = uniqueN(d$PTID),
    n_bridge_people = uniqueN(z$PTID),
    median_absolute_t1_endpoint_gap_days = stats::median(z$absolute_t1_endpoint_gap_days, na.rm = TRUE),
    n_within_730_days = sum(z$absolute_t1_endpoint_gap_days <= 730, na.rm = TRUE),
    status = "AVAILABLE"
  )
}), fill = TRUE)
fwrite(endpoint_flow, file.path(OUTDIR, "01_endpoint_bridge_flow.csv"),
       bom = TRUE)

# ---------------------------------------------------------------------------
# 4. Sample audit and CAT12 structural association with PGS
# ---------------------------------------------------------------------------
bridge_audit <- data.table(
  gate = c(
    "CAT12_report_people", "CAT12_target_ROI_people", "CAT12_PGS_people",
    "complete_baseline_bridge_people", "complete_mPACCdigit_followup_people",
    "complete_mPACCtrailsB_followup_people", "complete_MMSE_followup_people"
  ),
  n_people = c(
    uniqueN(report$PTID), uniqueN(temporal$PTID), uniqueN(base$PTID),
    uniqueN(base$PTID),
    uniqueN(pgs[PTID %in% base$PTID & is.finite(mPACCdigit), PTID]),
    uniqueN(pgs[PTID %in% base$PTID & is.finite(mPACCtrailsB), PTID]),
    uniqueN(pgs[PTID %in% base$PTID & is.finite(MMSE), PTID])
  )
)
fwrite(bridge_audit, file.path(OUTDIR, "02_bridge_sample_audit.csv"),
       bom = TRUE)
fwrite(base, file.path(OUTDIR, "03_CAT12_PGS_baseline_bridge_frame.csv"),
       bom = TRUE)

pc_terms <- paste0("PC", 1:10, collapse = " + ")
base_formula <- stats::as.formula(paste(
  "cat12_ad_vulnerable_gm ~ gMAV_noAPOE_z_eur + gAD_noAPOE_z_eur + APOE4 +",
  "age_bl + sex + education + dx_bl +", pc_terms
))
base_vars <- all.vars(base_formula)
d_base <- base[complete.cases(base[, ..base_vars])]
baseline_audit <- data.table(
  model_id = "CAT12_structural_on_PGS", status = "NOT_RUN",
  n_people = uniqueN(d_base$PTID), n_rows = nrow(d_base)
)
baseline_out <- data.table()
if (nrow(d_base) >= 100L) {
  fit <- stats::lm(base_formula, data = d_base)
  baseline_out <- tidy_lm(
    fit, model_id = "CAT12_structural_on_PGS",
    outcome = "cat12_ad_vulnerable_gm",
    n_people = uniqueN(d_base$PTID), n_rows = nrow(d_base)
  )
  baseline_out[, fdr_p_value := stats::p.adjust(p_value, method = "BH")]
  baseline_audit[, status := "RUN"]
}
if (nrow(baseline_out)) fwrite(
  baseline_out, file.path(OUTDIR, "04_CAT12_PGS_baseline_models.csv"),
  bom = TRUE
) else write_empty(
  file.path(OUTDIR, "04_CAT12_PGS_baseline_models.csv"),
  c("endpoint", "outcome", "model_id", "term", "estimate", "std_error",
    "statistic", "p_value", "conf_low", "conf_high", "n_people",
    "n_rows", "fdr_p_value")
)
fwrite(baseline_audit,
       file.path(OUTDIR, "05_baseline_model_execution_audit.csv"), bom = TRUE)

# ---------------------------------------------------------------------------
# 5. CAT12 structural association with endpoint-specific A/T/V/N
# ---------------------------------------------------------------------------
atvn_formula_text <- paste(
  "endpoint_z ~ cat12_ad_vulnerable_gm + gMAV_noAPOE_z_eur +",
  "gAD_noAPOE_z_eur + APOE4 + age_bl + sex + education + dx_bl +", pc_terms
)
atvn_results <- list()
atvn_audit <- list()
for (ep in endpoint_names) {
  d0 <- endpoint_raw[[ep]]
  if (is.null(d0)) {
    atvn_audit[[ep]] <- data.table(
      endpoint = ep, model_id = "CAT12_on_endpoint",
      status = "NOT_AVAILABLE", n_people = 0L, n_rows = 0L
    )
    next
  }
  d_ep <- merge(
    d0[, .(PTID, endpoint, endpoint_date, endpoint_value, endpoint_z)],
    base, by = "PTID", all = FALSE
  )
  d_ep[, absolute_t1_endpoint_gap_days :=
    abs(as.numeric(endpoint_date - t1_anchor_date))]
  # A/T/V/N bridge models are restricted to the same prespecified 730-day
  # linkage window used in the original endpoint-pairing branch. This avoids
  # describing later tau or hippocampal measures as concurrent CAT12 effects.
  d_ep <- d_ep[is.finite(absolute_t1_endpoint_gap_days) &
                 absolute_t1_endpoint_gap_days <= 730]
  f_ep <- stats::as.formula(atvn_formula_text)
  ep_vars <- all.vars(f_ep)
  d_ep <- d_ep[complete.cases(d_ep[, ..ep_vars])]
  n_people <- uniqueN(d_ep$PTID)
  if (n_people < 100L) {
    atvn_audit[[ep]] <- data.table(
      endpoint = ep, model_id = "CAT12_on_endpoint",
      status = "NOT_RUN_N_LT_100", n_people = n_people, n_rows = nrow(d_ep)
    )
    next
  }
  fit <- stats::lm(f_ep, data = d_ep)
  atvn_results[[ep]] <- tidy_lm(
    fit, model_id = "CAT12_on_endpoint", outcome = "endpoint_z",
    endpoint = ep, n_people = n_people, n_rows = nrow(d_ep)
  )
  atvn_audit[[ep]] <- data.table(
    endpoint = ep, model_id = "CAT12_on_endpoint",
    status = "RUN", n_people = n_people, n_rows = nrow(d_ep)
  )
}
atvn_out <- if (length(atvn_results)) rbindlist(atvn_results, fill = TRUE) else data.table()
if (nrow(atvn_out)) {
  atvn_out[, fdr_p_value := stats::p.adjust(p_value, method = "BH")]
  atvn_out[, cat12_term := term == "cat12_ad_vulnerable_gm"]
  atvn_out[, cat12_term_fdr_p := NA_real_]
  ix <- atvn_out$cat12_term == TRUE & is.finite(atvn_out$p_value)
  if (any(ix)) {
    atvn_out[ix, cat12_term_fdr_p :=
      stats::p.adjust(p_value, method = "BH")]
  }
  fwrite(atvn_out, file.path(OUTDIR, "06_CAT12_ATVN_baseline_models.csv"),
         bom = TRUE)
} else write_empty(
  file.path(OUTDIR, "06_CAT12_ATVN_baseline_models.csv"),
  c("endpoint", "outcome", "model_id", "term", "estimate", "std_error",
    "statistic", "p_value", "conf_low", "conf_high", "n_people",
    "n_rows", "fdr_p_value", "cat12_term", "cat12_term_fdr_p")
)
fwrite(rbindlist(atvn_audit, fill = TRUE),
       file.path(OUTDIR, "07_ATVN_model_execution_audit.csv"), bom = TRUE)

# ---------------------------------------------------------------------------
# 6. Longitudinal cognitive expression of CAT12 and PGS
# ---------------------------------------------------------------------------
cog_outcomes <- c("mPACCdigit", "mPACCtrailsB", "MMSE")
long <- merge(
  pgs[PTID %in% base$PTID],
  base[, .(PTID, cat12_ad_vulnerable_gm, APOE4, age_bl, sex, education,
           dx_bl, gMAV_noAPOE_z_eur, gAD_noAPOE_z_eur,
           PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)],
  by = "PTID", all = FALSE, suffixes = c("", "_cat12")
)
long[, time_years := as_num(time_years)]
long <- long[is.finite(time_years)]
cog_results <- list()
cog_audit <- list()
for (cog in cog_outcomes) {
  d <- copy(long)
  d[, cognition_raw := as_num(get(cog))]
  d <- d[is.finite(cognition_raw)]
  d[, cognition_z := safe_z(cognition_raw)]
  d <- d[is.finite(cognition_z)]
  d[, sex := factor(sex)]
  d[, dx_bl := factor(dx_bl)]
  n_people <- uniqueN(d$PTID)
  n_rows <- nrow(d)
  if (n_people < 100L || n_rows < 250L) {
    cog_audit[[cog]] <- data.table(
      outcome = cog, model_id = "CAT12_PGS_cognition_expression",
      status = "NOT_RUN_INSUFFICIENT_LONGITUDINAL_DATA",
      n_people = n_people, n_rows = n_rows
    )
    next
  }
  fixed <- paste(
    "cognition_z ~ time_years * (cat12_ad_vulnerable_gm +",
    "gMAV_noAPOE_z_eur + gAD_noAPOE_z_eur) + APOE4 + age_bl + sex +",
    "education + dx_bl +", pc_terms
  )
  f_slope <- stats::as.formula(paste(fixed, "+ (1 + time_years | PTID)"))
  f_intercept <- stats::as.formula(paste(fixed, "+ (1 | PTID)"))
  fit <- fit_mixed(f_slope, f_intercept, d)
  cog_results[[cog]] <- tidy_lmer(
    fit$fit, model_id = "CAT12_PGS_cognition_expression",
    outcome = cog, n_people = n_people, n_rows = n_rows,
    random_structure = fit$random_structure
  )
  cog_audit[[cog]] <- data.table(
    outcome = cog, model_id = "CAT12_PGS_cognition_expression",
    status = "RUN", n_people = n_people, n_rows = n_rows,
    random_structure = fit$random_structure
  )
}
cog_out <- if (length(cog_results)) rbindlist(cog_results, fill = TRUE) else data.table()
if (nrow(cog_out)) {
  cog_out[, interaction_term := grepl("^time_years:", term)]
  cog_out[, cat12_time_term := term == "time_years:cat12_ad_vulnerable_gm"]
  cog_out[, fdr_interaction_p := NA_real_]
  ix <- cog_out$interaction_term == TRUE & is.finite(cog_out$p_value)
  if (any(ix)) {
    cog_out[ix, fdr_interaction_p :=
      stats::p.adjust(p_value, method = "BH")]
  }
  fwrite(cog_out,
         file.path(OUTDIR, "08_CAT12_cognition_longitudinal_models.csv"),
         bom = TRUE)
} else write_empty(
  file.path(OUTDIR, "08_CAT12_cognition_longitudinal_models.csv"),
  c("endpoint", "outcome", "model_id", "term", "estimate", "std_error",
    "statistic", "p_value", "n_people", "n_rows", "random_structure",
    "interaction_term", "cat12_time_term", "fdr_interaction_p")
)
fwrite(rbindlist(cog_audit, fill = TRUE),
       file.path(OUTDIR, "09_cognition_model_execution_audit.csv"), bom = TRUE)

# ---------------------------------------------------------------------------
# 7. Provenance and interpretation
# ---------------------------------------------------------------------------
fwrite(data.table(
  component = c(
    "CAT12_primary_structural_phenotype", "CAT12_target_regions",
    "PGS_primary", "PGS_AD_adjustment", "APOE_adjustment",
    "ancestry_adjustment", "ATVN_source", "cognition_primary",
    "cognition_sensitivity", "causal_scope"
  ),
  specification = c(
    "bilateral temporal GM/TIV; larger score means greater atrophy vulnerability",
    paste(targets, collapse = "; "),
    "gMAV-noAPOE European-ancestry score",
    "gAD-noAPOE European-ancestry score",
    "APOE4 carrier indicator from the frozen bridge manifest",
    "PC1-PC10, age, sex, education and baseline diagnosis",
    "10C earliest endpoint-anchored A/T/V/N records",
    "mPACCdigit longitudinal mixed model",
    "mPACCtrailsB and MMSE longitudinal mixed models",
    "Covariate-adjusted structural and longitudinal cognitive associations"
  )
), file.path(OUTDIR, "10_frozen_bridge_specification.csv"), bom = TRUE)

writeLines(c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Complete CAT12/PGS baseline bridge participants: ", uniqueN(base$PTID)),
  "No CAT12 process was launched by this script.",
  "A/T/V/N measurements use the endpoint-anchored records from ancillary/03_adni_pathology_index_inputs.R.",
  "The primary structural bridge coefficient is cat12_ad_vulnerable_gm.",
  "The primary clinical-expression coefficient is time_years:cat12_ad_vulnerable_gm in mPACCdigit."
), file.path(OUTDIR, "11_INTERPRETATION.txt"), useBytes = TRUE)

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "12_sessionInfo.txt"), useBytes = TRUE)
cat("ADNI CAT12 models completed: ", OUTDIR, "\n", sep = "")
