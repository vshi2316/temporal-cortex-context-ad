# First and second MRI event risk sets
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

# Construct clinical event risk sets from regional MRI derivatives and ADNIMERGE.

project_dir <- study_path("project")
project_dir_native <- study_path("project")
result_root <- file.path(project_dir, "results")
prep_file <- file.path(result_root, "13F_R_only_structural_analysis", "04_ADNI_R_only_longitudinal_structural_frame.rds")
pgs_file <- file.path(result_root, "10A_ADNI_gMAV_PGS_assembly", "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
clinical_candidates <- c(
  study_path("adni_data", "Clinical/ADNIMERGE_03Jan2026.csv"),
  study_path("cohorts", "ADNI/data/Clinical/ADNIMERGE_03Jan2026.csv")
)
clinical_file <- clinical_candidates[file.exists(clinical_candidates)][1L]
out_dir <- file.path(project_dir_native, "results", "13F_R1_clinical_event_bridge")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  path
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
finite <- function(x) is.finite(num(x))
norm_rid <- function(x) {
  x <- sub("^.*?([0-9]+)$", "\\1", as.character(x))
  x <- sub("^0+", "", x)
  x[x == ""] <- "0"
  x
}
zscore <- function(x) {
  x <- num(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
first_nonmissing <- function(x, numeric = FALSE) {
  if (numeric) {
    z <- num(x)
    return(if (any(is.finite(z))) z[which(is.finite(z))[1L]] else NA_real_)
  }
  z <- trimws(as.character(x)); z <- z[!is.na(z) & nzchar(z)]
  if (length(z)) z[1L] else NA_character_
}
parse_date <- function(x) {
  x <- as.character(x)
  out <- as.IDate(rep(NA_character_, length(x)))
  for (fmt in c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y")) {
    ix <- is.na(out) & !is.na(x) & nzchar(x)
    if (!any(ix)) break
    out[ix] <- as.IDate(x[ix], format = fmt)
  }
  out
}
clean_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("CN", "NORMAL", "COGNITIVELY NORMAL")] <- "CN"
  out[x %in% c("MCI", "EMCI", "LMCI", "SMC")] <- "MCI"
  out[x %in% c("AD", "DEMENTIA", "ALZHEIMER DISEASE", "ALZHEIMER'S DISEASE")] <- "Dementia"
  out
}

extract_cox <- function(fit, model_id, n_input, n_events) {
  if (is.null(fit)) {
    return(data.table(model_id = model_id, term = NA_character_, estimate = NA_real_,
                      hazard_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
                      p_value = NA_real_, n = n_input, events = n_events,
                      status = "FAILED"))
  }
  b <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  keep <- is.finite(b) & is.finite(se)
  z <- b[keep] / se[keep]
  data.table(
    model_id = model_id, term = names(b)[keep], estimate = as.numeric(b[keep]),
    hazard_ratio = exp(as.numeric(b[keep])),
    conf_low = exp(as.numeric(b[keep]) - qnorm(0.975) * se[keep]),
    conf_high = exp(as.numeric(b[keep]) + qnorm(0.975) * se[keep]),
    p_value = 2 * pnorm(-abs(z)), n = n_input, events = n_events,
    status = "SUCCESS"
  )
}

fit_cox_model <- function(dat, rhs, model_id) {
  vars <- unique(c("followup_years", "event", all.vars(as.formula(paste("~", rhs)))))
  vars <- vars[vars %in% names(dat)]
  a <- dat[complete.cases(dat[, vars, with = FALSE])]
  a <- a[is.finite(followup_years) & followup_years > 0 & event %in% c(FALSE, TRUE)]
  ne <- sum(a$event)
  if (nrow(a) < 100L || ne < 20L) {
    return(list(
      fit = NULL,
      coefficients = extract_cox(NULL, model_id, nrow(a), ne),
      diagnostics = data.table(model_id = model_id, status = "NOT_RUN_INSUFFICIENT_DATA",
                               n = nrow(a), events = ne, c_index = NA_real_,
                               ph_global_p = NA_real_, ph_status = NA_character_),
      ph_table = data.table()
    ))
  }
  f <- as.formula(paste0("Surv(followup_years, event) ~ ", rhs))
  fit <- tryCatch(
    coxph(f, data = a, ties = "efron", x = TRUE, y = TRUE, model = TRUE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      coefficients = data.table(model_id = model_id, term = NA_character_, estimate = NA_real_,
                                hazard_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
                                p_value = NA_real_, n = nrow(a), events = ne,
                                status = "FAILED", message = conditionMessage(fit)),
      diagnostics = data.table(model_id = model_id, status = "FAILED", n = nrow(a), events = ne,
                               c_index = NA_real_, ph_global_p = NA_real_, ph_status = conditionMessage(fit)),
      ph_table = data.table()
    ))
  }
  sm <- summary(fit)
  c_index <- if (!is.null(sm$concordance)) unname(sm$concordance[1L]) else NA_real_
  ph <- tryCatch(cox.zph(fit), error = function(e) e)
  ph_global <- if (inherits(ph, "cox.zph")) as.numeric(ph$table[nrow(ph$table), "p"])[1L] else NA_real_
  ph_status <- if (inherits(ph, "cox.zph")) "TESTED" else if (inherits(ph, "error")) conditionMessage(ph) else NA_character_
  ph_table <- if (inherits(ph, "cox.zph")) {
    out <- as.data.table(ph$table, keep.rownames = "term")
    out[, model_id := model_id]
    setcolorder(out, c("model_id", "term", setdiff(names(out), c("model_id", "term"))))
    out
  } else data.table()
  list(
    fit = fit,
    coefficients = extract_cox(fit, model_id, nrow(a), ne),
    diagnostics = data.table(model_id = model_id, status = "SUCCESS", n = nrow(a), events = ne,
                             c_index = c_index, ph_global_p = ph_global, ph_status = ph_status,
                             AIC = AIC(fit), logLik = as.numeric(logLik(fit))),
    ph_table = ph_table
  )
}

# -----------------------------------------------------------------------------
# Build a deterministic one-record-per-RID/date regional MRI frame.
# -----------------------------------------------------------------------------
d <- as.data.table(readRDS(need(prep_file)))
d[, RID := norm_rid(RID)]
d[, mri_date := parse_date(mri_date)]
d[, primary_complete := finite(ad_vulnerable_index)]
d[, IMAGEUID_num := num(IMAGEUID)]
setorderv(d, c("RID", "mri_date", "primary_complete", "IMAGEUID_num"),
          order = c(1, 1, -1, -1), na.last = TRUE)
d <- d[, .SD[1L], by = .(RID, mri_date)]
d <- d[finite(ad_vulnerable_index) & finite(clinical_gMAV_noAPOE_z_eur) &
         finite(clinical_gAD_noAPOE_z_eur) & !is.na(mri_date)]
setorder(d, RID, mri_date)
d[, `:=`(
  gMAV = num(clinical_gMAV_noAPOE_z_eur),
  gAD = num(clinical_gAD_noAPOE_z_eur),
  structure_z = zscore(ad_vulnerable_index),
  global_structure_z = zscore(global_cortical_thickness)
)]
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]

pgs <- fread(file = need(pgs_file), sep = "\t", showProgress = FALSE)
pgs[, RID := norm_rid(RID)]
pgs[, pgs_exam_date := parse_date(EXAMDATE)]
setorder(pgs, RID, pgs_exam_date)
base_cov <- pgs[, .(
  pgs_baseline_age = first_nonmissing(age_bl, TRUE),
  pgs_baseline_sex = first_nonmissing(sex),
  pgs_baseline_education = first_nonmissing(education, TRUE),
  pgs_baseline_dx = first_nonmissing(dx_bl)
), by = RID]
d <- merge(d, base_cov, by = "RID", all.x = TRUE, sort = FALSE)
d[, `:=`(
  baseline_age = num(pgs_baseline_age),
  education = num(pgs_baseline_education),
  sex = factor(pgs_baseline_sex),
  pgs_baseline_dx = trimws(as.character(pgs_baseline_dx)),
  scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_"))
)]
for (v in paste0("clinical_PC", 1:10)) if (v %in% names(d)) d[, (v) := num(get(v))]

# -----------------------------------------------------------------------------
# Read complete ADNI clinical follow-up and normalize diagnosis labels.
# -----------------------------------------------------------------------------
clinical_cols <- c("RID", "EXAMDATE", "DX", "DX_bl", "APOE4", "AGE", "PTGENDER", "PTEDUCAT",
                   "mPACCdigit", "mPACCtrailsB", "ADAS13", "CDRSB", "MMSE", "MOCA")
clinical <- fread(file = need(clinical_file), select = clinical_cols, showProgress = FALSE)
clinical[, RID := norm_rid(RID)]
clinical[, exam_date := parse_date(EXAMDATE)]
clinical[, dx_clean := trimws(as.character(DX))]
clinical <- clinical[!is.na(exam_date)]
clinical <- unique(clinical, by = c("RID", "exam_date", "dx_clean"))
setkey(clinical, RID)

# The event is the observed ADNIMERGE Dementia classification.
event_label <- "Dementia"

make_event_table <- function(index_table, start_date_col = "index_mri_date", label = "baseline") {
  out <- vector("list", nrow(index_table))
  for (i in seq_len(nrow(index_table))) {
    rid <- index_table$RID[i]
    start <- index_table[[start_date_col]][i]
    z <- clinical[.(rid)]
    z <- z[!is.na(exam_date)]
    if (!nrow(z) || is.na(start)) {
      out[[i]] <- data.table(RID = rid, baseline_dx_clinical = NA_character_,
                             prevalent_dementia = NA, event = NA, event_date = as.IDate(NA),
                             censor_date = as.IDate(NA), followup_years = NA_real_,
                             clinical_visits = 0L, bridge_label = label)
      next
    }
    pre <- z[exam_date <= start & !is.na(dx_clean) & nzchar(dx_clean)]
  if (nrow(pre)) {
      setorder(pre, exam_date)
      baseline_dx <- clean_stage(pre$dx_clean[nrow(pre)])
    } else {
      baseline_dx <- NA_character_
    }
    # A 90-day buffer treats near-index dementia ascertainment as prevalent.
    prev <- any(z$exam_date <= start + 90 & z$dx_clean == event_label, na.rm = TRUE)
    post <- z[exam_date > start + 90]
    ev <- post[dx_clean == event_label]
    event_date <- if (nrow(ev)) min(ev$exam_date) else as.IDate(NA)
    censor_date <- max(z$exam_date, na.rm = TRUE)
    end_date <- if (!is.na(event_date)) event_date else censor_date
    fu <- as.numeric(end_date - start) / 365.25
    out[[i]] <- data.table(
      RID = rid, baseline_dx_clinical = baseline_dx,
      prevalent_dementia = prev, event = !is.na(event_date), event_date = event_date,
      censor_date = censor_date, followup_years = fu,
      clinical_visits = nrow(z), bridge_label = label
    )
  }
  rbindlist(out, fill = TRUE)
}

# -----------------------------------------------------------------------------
# Baseline MRI -> incident dementia analysis.
# -----------------------------------------------------------------------------
base <- d[, .SD[1L], by = RID]
setnames(base, "mri_date", "index_mri_date")
event_base <- make_event_table(base)
analysis <- merge(base, event_base, by = "RID", all.x = TRUE, sort = FALSE)

# Match APOE4 and age to the closest dated clinical assessment to the index MRI.
analysis[, `:=`(
  age_index = NA_real_, apoe4 = NA_real_,
  baseline_dx = clean_stage(ifelse(!is.na(baseline_dx_clinical), baseline_dx_clinical, pgs_baseline_dx)),
  pc1 = num(clinical_PC1), pc2 = num(clinical_PC2), pc3 = num(clinical_PC3), pc4 = num(clinical_PC4),
  pc5 = num(clinical_PC5), pc6 = num(clinical_PC6), pc7 = num(clinical_PC7), pc8 = num(clinical_PC8),
  pc9 = num(clinical_PC9), pc10 = num(clinical_PC10)
)]
for (i in seq_len(nrow(analysis))) {
  z <- clinical[.(analysis$RID[i])]
  z <- z[!is.na(exam_date)]
  if (!nrow(z)) next
  z[, distance := abs(as.numeric(exam_date - analysis$index_mri_date[i]))]
  z <- z[order(distance, exam_date)]
  analysis$age_index[i] <- num(z$AGE[1L])
  analysis$apoe4[i] <- num(z$APOE4[1L])
}
analysis[, `:=`(
  sex = droplevels(factor(pgs_baseline_sex)),
  baseline_dx = droplevels(factor(baseline_dx)),
  scanner = droplevels(scanner),
  apoe4_z = zscore(apoe4), age_index_z = zscore(age_index),
  education_z = zscore(education)
)]
analysis[baseline_dx == "Dementia", prevalent_dementia := TRUE]

# No prevalent Dementia at or before the index MRI, and a positive follow-up.
analysis_incident <- analysis[prevalent_dementia == FALSE & is.finite(followup_years) & followup_years > 0]

usable <- function(dat, vars) {
  vars[vapply(vars, function(v) {
    if (!v %in% names(dat)) return(FALSE)
    x <- dat[[v]]
    if (is.factor(x)) nlevels(droplevels(x)) >= 2L else length(unique(x[is.finite(num(x))])) >= 2L
  }, logical(1))]
}
pc_terms <- paste0("pc", 1:10)
covariates <- c("age_index_z", "sex", "education_z", "apoe4_z", "scanner", pc_terms)
cov_full <- usable(analysis_incident, covariates)
rhs_structure <- paste(c("structure_z", "gMAV", "gAD", cov_full, "strata(baseline_dx)"), collapse = " + ")
rhs_global <- paste(c("global_structure_z", "gMAV", "gAD", cov_full, "strata(baseline_dx)"), collapse = " + ")
rhs_reduced <- paste(c("structure_z", "gMAV", "gAD", usable(analysis_incident, setdiff(covariates, pc_terms)),
                       "strata(baseline_dx)"), collapse = " + ")

fits <- list()
fits[["baseline_structure_full"]] <- fit_cox_model(analysis_incident, rhs_structure, "baseline_structure_full")
fits[["baseline_global_sensitivity"]] <- fit_cox_model(analysis_incident, rhs_global, "baseline_global_sensitivity")
fits[["baseline_structure_reduced"]] <- fit_cox_model(analysis_incident, rhs_reduced, "baseline_structure_reduced")

# -----------------------------------------------------------------------------
# First-to-second-MRI structural change -> later dementia landmark analysis.
# The event clock begins at the first MRI at least 270 days after the index MRI.
# -----------------------------------------------------------------------------
setorder(d, RID, mri_date)
pair <- d[, {
  dates <- sort(unique(mri_date))
  eligible_second <- if (length(dates) >= 2L) dates[as.integer(dates - dates[1L]) >= 270L] else as.IDate(character())
  if (!length(eligible_second)) {
    NULL
  } else {
    first_date <- dates[1L]
    second_date <- eligible_second[1L]
    first_row <- .SD[mri_date == first_date][1L]
    second_row <- .SD[mri_date == second_date][1L]
    interval <- as.numeric(second_date - first_date) / 365.25
    data.table(
      index_mri_date = first_date,
      landmark_mri_date = second_date,
      interval_years = interval,
      baseline_structure_z = first_row$structure_z,
      landmark_structure_z = second_row$structure_z,
      annual_change = (num(second_row$ad_vulnerable_index) - num(first_row$ad_vulnerable_index)) / interval,
      gMAV = first_row$gMAV, gAD = first_row$gAD,
      education = first_row$education, sex = as.character(first_row$sex),
      scanner = as.character(first_row$scanner),
      pc1 = first_row$clinical_PC1, pc2 = first_row$clinical_PC2,
      pc3 = first_row$clinical_PC3, pc4 = first_row$clinical_PC4,
      pc5 = first_row$clinical_PC5, pc6 = first_row$clinical_PC6,
      pc7 = first_row$clinical_PC7, pc8 = first_row$clinical_PC8,
      pc9 = first_row$clinical_PC9, pc10 = first_row$clinical_PC10
    )
  }
}, by = RID]

if (nrow(pair)) {
  pair <- merge(pair, analysis[, .(RID, age_index, apoe4)], by = "RID", all.x = TRUE, sort = FALSE)
  pair[, `:=`(
    atrophy_rate_z = zscore(-annual_change),
    age_landmark_z = zscore(age_index + interval_years),
    education_z = zscore(education), apoe4_z = zscore(apoe4),
    sex = factor(sex), scanner = factor(scanner)
  )]
  landmark_event <- make_event_table(pair, start_date_col = "landmark_mri_date", label = "landmark")
  pair <- merge(pair, landmark_event[, .(
    RID, landmark_dx = baseline_dx_clinical,
    prevalent_dementia, event, event_date, censor_date,
    followup_years, clinical_visits
  )], by = "RID", all.x = TRUE, sort = FALSE)
  pair[, landmark_dx := clean_stage(landmark_dx)]
  # Fall back to the original clinical baseline stage only when no dated
  # diagnosis exists before the landmark; do not use a post-landmark label.
  pair <- merge(pair, analysis[, .(RID, baseline_dx)], by = "RID", all.x = TRUE, sort = FALSE)
  pair[is.na(landmark_dx), landmark_dx := as.character(baseline_dx)]
  pair[, landmark_dx := droplevels(factor(landmark_dx))]
  pair_incident <- pair[prevalent_dementia == FALSE & is.finite(followup_years) & followup_years > 0]

  landmark_covars <- c("age_landmark_z", "sex", "education_z", "apoe4_z", "scanner", pc_terms)
  cov_landmark <- usable(pair_incident, landmark_covars)
  rhs_landmark <- paste(c("atrophy_rate_z", "baseline_structure_z", "gMAV", "gAD", cov_landmark,
                          "strata(landmark_dx)"), collapse = " + ")
  rhs_landmark_reduced <- paste(c("atrophy_rate_z", "baseline_structure_z", "gMAV", "gAD",
                                  usable(pair_incident, setdiff(landmark_covars, c(pc_terms, "scanner"))),
                                  "strata(landmark_dx)"), collapse = " + ")
  fits[["landmark_atrophy_full"]] <- fit_cox_model(pair_incident, rhs_landmark, "landmark_atrophy_full")
  fits[["landmark_atrophy_reduced"]] <- fit_cox_model(pair_incident, rhs_landmark_reduced, "landmark_atrophy_reduced")
} else {
  pair_incident <- data.table()
}

# -----------------------------------------------------------------------------
# Export model estimates and risk-set counts.
# -----------------------------------------------------------------------------
coefficients <- rbindlist(lapply(fits, `[[`, "coefficients"), fill = TRUE)
diagnostics <- rbindlist(lapply(fits, `[[`, "diagnostics"), fill = TRUE)
ph_results <- rbindlist(lapply(fits, `[[`, "ph_table"), fill = TRUE)

sample_flow <- rbindlist(list(
  data.table(stage = "deduplicated_PGS_complete_structural_frame",
             people = uniqueN(d$RID), rows = nrow(d)),
  data.table(stage = "baseline_MRI_index",
             people = uniqueN(base$RID), rows = nrow(base)),
  data.table(stage = "clinical_table_linked",
             people = sum(!is.na(analysis$clinical_visits) & analysis$clinical_visits > 0), rows = nrow(analysis)),
  data.table(stage = "excluded_prevalent_dementia",
             people = sum(analysis$prevalent_dementia %in% TRUE, na.rm = TRUE), rows = NA_integer_),
  data.table(stage = "baseline_incident_dementia_analysis",
             people = uniqueN(analysis_incident$RID), rows = nrow(analysis_incident)),
  data.table(stage = "baseline_incident_dementia_events",
             people = sum(analysis_incident$event %in% TRUE, na.rm = TRUE), rows = NA_integer_),
  data.table(stage = "two_MRI_landmark_pairs",
             people = uniqueN(pair$RID), rows = nrow(pair)),
  data.table(stage = "landmark_incident_dementia_analysis",
             people = uniqueN(pair_incident$RID), rows = nrow(pair_incident)),
  data.table(stage = "landmark_incident_dementia_events",
             people = sum(pair_incident$event %in% TRUE, na.rm = TRUE), rows = NA_integer_)
), fill = TRUE)

diagnosis_distribution <- analysis[, .N, by = .(baseline_dx_clinical, prevalent_dementia)]
setorder(diagnosis_distribution, baseline_dx_clinical, prevalent_dementia)
event_label_distribution <- clinical[, .N, by = dx_clean][order(-N)]

fwrite(coefficients, file.path(out_dir, "01_clinical_event_cox_coefficients.csv"), bom = TRUE)
fwrite(diagnostics, file.path(out_dir, "02_clinical_event_model_diagnostics.csv"), bom = TRUE)
fwrite(ph_results, file.path(out_dir, "02B_proportional_hazards_tests.csv"), bom = TRUE)
fwrite(sample_flow, file.path(out_dir, "03_clinical_event_sample_flow.csv"), bom = TRUE)
fwrite(diagnosis_distribution, file.path(out_dir, "04_baseline_diagnosis_distribution.csv"), bom = TRUE)
fwrite(event_label_distribution, file.path(out_dir, "05_clinical_dx_label_distribution.csv"), bom = TRUE)
fwrite(analysis, file.path(out_dir, "06_baseline_event_analysis_frame.csv"), bom = TRUE)
if (nrow(pair)) fwrite(pair, file.path(out_dir, "07_landmark_event_analysis_frame.csv"), bom = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
message("ADNI event risk sets completed: ", out_dir)
