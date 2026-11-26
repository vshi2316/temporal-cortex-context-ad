# Event timing, clinical adjustment and sensitivity models
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

# The event is the first detected ADNI Dementia classification.
# Baseline diagnosis is matched within 90 days of the first MRI.
# Risk starts at the later MRI or diagnosis date in baseline CN/MCI participants.
# Valid diagnostic visits define follow-up. Cox models use diagnosis strata
# and site-clustered variance.

project_dir <- study_path("project")
project_native <- study_path("project")
result_root <- file.path(project_dir, "results")
prep_file <- file.path(result_root, "13F_R_only_structural_analysis", "04_ADNI_R_only_longitudinal_structural_frame.rds")
pgs_file <- file.path(result_root, "10A_ADNI_gMAV_PGS_assembly", "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
clinical_candidates <- c(
  study_path("adni_data", "Clinical/ADNIMERGE_03Jan2026.csv"),
  study_path("cohorts", "ADNI/data/Clinical/ADNIMERGE_03Jan2026.csv")
)
clinical_file <- clinical_candidates[file.exists(clinical_candidates)][1L]
out_dir <- file.path(project_native, "results", "13F_R1_clinical_closure_confirmatory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (is.na(path) || !file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  path
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
finite <- function(x) is.finite(num(x))
norm_rid <- function(x) {
  x <- sub("^.*?([0-9]+)$", "\\1", as.character(x))
  x <- sub("^0+", "", x); x[x == ""] <- "0"; x
}
parse_date <- function(x) {
  x <- as.character(x); out <- as.IDate(rep(NA_character_, length(x)))
  for (fmt in c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y")) {
    ix <- is.na(out) & !is.na(x) & nzchar(x)
    if (!any(ix)) break
    out[ix] <- as.IDate(x[ix], format = fmt)
  }
  out
}
clean_stage <- function(x) {
  x <- toupper(trimws(as.character(x))); out <- rep(NA_character_, length(x))
  out[x %in% c("CN", "NORMAL", "COGNITIVELY NORMAL")] <- "CN"
  out[x %in% c("MCI", "EMCI", "LMCI", "SMC")] <- "MCI"
  out[x %in% c("AD", "DEMENTIA", "ALZHEIMER DISEASE", "ALZHEIMER'S DISEASE")] <- "Dementia"
  out
}
zscore <- function(x) {
  x <- num(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
first_nonmissing <- function(x, numeric = FALSE) {
  if (numeric) {
    z <- num(x); return(if (any(is.finite(z))) z[which(is.finite(z))[1L]] else NA_real_)
  }
  z <- trimws(as.character(x)); z <- z[!is.na(z) & nzchar(z)]
  if (length(z)) z[1L] else NA_character_
}

# Robust extraction with site-clustered standard errors.
extract_cox <- function(fit, model_id, n, events) {
  if (is.null(fit)) return(data.table(model_id = model_id, term = NA_character_, estimate = NA_real_,
                                      hazard_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
                                      p_value = NA_real_, n = n, events = events, status = "FAILED"))
  b <- coef(fit); se <- sqrt(diag(vcov(fit))); keep <- is.finite(b) & is.finite(se)
  zz <- b[keep] / se[keep]
  data.table(model_id = model_id, term = names(b)[keep], estimate = as.numeric(b[keep]),
             hazard_ratio = exp(as.numeric(b[keep])),
             conf_low = exp(as.numeric(b[keep]) - qnorm(.975) * se[keep]),
             conf_high = exp(as.numeric(b[keep]) + qnorm(.975) * se[keep]),
             p_value = 2 * pnorm(-abs(zz)), n = n, events = events, status = "SUCCESS")
}

fit_cox <- function(dat, rhs, model_id) {
  vars <- unique(c("followup_years", "event", all.vars(as.formula(paste("~", rhs)))))
  vars <- vars[vars %in% names(dat)]
  a <- dat[complete.cases(dat[, vars, with = FALSE])]
  a <- a[is.finite(followup_years) & followup_years > 0 & event %in% c(FALSE, TRUE)]
  ne <- sum(a$event)
  if (nrow(a) < 100L || ne < 20L) {
    return(list(fit = NULL, coefficients = extract_cox(NULL, model_id, nrow(a), ne),
                diagnostics = data.table(model_id = model_id, status = "NOT_RUN_INSUFFICIENT_DATA",
                                         n = nrow(a), events = ne, c_index = NA_real_,
                                         ph_global_p = NA_real_, ph_status = NA_character_),
                ph_table = data.table()))
  }
  f <- as.formula(paste0("Surv(followup_years, event) ~ ", rhs))
  fit <- tryCatch(coxph(f, data = a, ties = "efron", x = TRUE, y = TRUE,
                        model = TRUE), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(fit = NULL,
                coefficients = data.table(model_id = model_id, term = NA_character_, estimate = NA_real_,
                                          hazard_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
                                          p_value = NA_real_, n = nrow(a), events = ne, status = "FAILED",
                                          message = conditionMessage(fit)),
                diagnostics = data.table(model_id = model_id, status = "FAILED", n = nrow(a), events = ne,
                                         c_index = NA_real_, ph_global_p = NA_real_, ph_status = conditionMessage(fit)),
                ph_table = data.table()))
  }
  sm <- summary(fit); c_index <- if (!is.null(sm$concordance)) unname(sm$concordance[1L]) else NA_real_
  ph <- tryCatch(cox.zph(fit), error = function(e) e)
  ph_table <- if (inherits(ph, "cox.zph")) {
    q <- as.data.table(ph$table, keep.rownames = "term"); q[, model_id := model_id]
    setcolorder(q, c("model_id", "term", setdiff(names(q), c("model_id", "term")))); q
  } else data.table()
  global_p <- if (nrow(ph_table)) ph_table[term == "GLOBAL", p][1L] else NA_real_
  list(fit = fit, coefficients = extract_cox(fit, model_id, nrow(a), ne),
       diagnostics = data.table(model_id = model_id, status = "SUCCESS", n = nrow(a), events = ne,
                                c_index = c_index, ph_global_p = global_p,
                                ph_status = if (nrow(ph_table)) "TESTED" else "NOT_TESTED",
                                AIC = AIC(fit), logLik = as.numeric(logLik(fit))),
       ph_table = ph_table)
}

# -----------------------------------------------------------------------------
# MRI and PGS frame.
# -----------------------------------------------------------------------------
d <- as.data.table(readRDS(need(prep_file)))
d[, `:=`(RID = norm_rid(RID), mri_date = parse_date(mri_date),
         primary_complete = finite(ad_vulnerable_index), IMAGEUID_num = num(IMAGEUID))]
setorderv(d, c("RID", "mri_date", "primary_complete", "IMAGEUID_num"),
          c(1, 1, -1, -1), na.last = TRUE)
d <- d[, .SD[1L], by = .(RID, mri_date)]
d <- d[finite(ad_vulnerable_index) & finite(global_cortical_thickness) &
         finite(clinical_gMAV_noAPOE_z_eur) & finite(clinical_gAD_noAPOE_z_eur) & !is.na(mri_date)]
setorder(d, RID, mri_date)
d[, `:=`(
  gMAV = num(clinical_gMAV_noAPOE_z_eur),
  gAD = num(clinical_gAD_noAPOE_z_eur),
  structure_vulnerability_z = -zscore(ad_vulnerable_index),
  global_vulnerability_z = -zscore(global_cortical_thickness),
  scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_"))
)]
for (v in paste0("clinical_PC", 1:10)) d[, (sub("clinical_", "", v)) := num(get(v))]

pgs <- fread(need(pgs_file), sep = "\t", showProgress = FALSE)
pgs[, `:=`(RID = norm_rid(RID), pgs_date = parse_date(EXAMDATE))]
setorder(pgs, RID, pgs_date)
pgs_base <- pgs[, .(
  pgs_age = first_nonmissing(age_bl, TRUE),
  pgs_sex = first_nonmissing(sex),
  pgs_education = first_nonmissing(education, TRUE),
  pgs_dx = clean_stage(first_nonmissing(dx_bl))
), by = RID]
d <- merge(d, pgs_base, by = "RID", all.x = TRUE, sort = FALSE)
setorder(d, RID, mri_date)
base <- d[, .SD[1L], by = RID]
setnames(base, "mri_date", "index_mri_date")

# -----------------------------------------------------------------------------
# Complete valid-DX clinical sequence.
# -----------------------------------------------------------------------------
clinical_cols <- c("RID", "SITE", "EXAMDATE", "DX", "DX_bl", "APOE4", "AGE", "Years_bl",
                   "PTGENDER", "PTEDUCAT", "mPACCdigit", "mPACCtrailsB", "ADAS13",
                   "CDRSB", "FAQ", "MMSE", "MOCA")
clinical <- fread(need(clinical_file), select = clinical_cols, showProgress = FALSE)
clinical[, `:=`(RID = norm_rid(RID), exam_date = parse_date(EXAMDATE), dx_stage = clean_stage(DX))]
clinical <- clinical[!is.na(exam_date) & !is.na(dx_stage)]
clinical[, completeness__ := rowSums(!is.na(.SD)), .SDcols = c("APOE4", "AGE", "PTGENDER", "PTEDUCAT",
                                                                "mPACCdigit", "CDRSB", "FAQ", "MMSE")]
setorderv(clinical, c("RID", "exam_date", "completeness__"), c(1, 1, -1))
clinical <- clinical[, .SD[1L], by = .(RID, exam_date)]
setkey(clinical, RID)

safe_idate_max <- function(a, b) as.IDate(max(as.numeric(c(a, b)), na.rm = TRUE), origin = "1970-01-01")
safe_idate_mid <- function(a, b) as.IDate(round(mean(as.numeric(c(a, b)))), origin = "1970-01-01")

link_clinical_event <- function(index_table, start_col, window_days = 90L, label = "baseline") {
  out <- vector("list", nrow(index_table))
  for (i in seq_len(nrow(index_table))) {
    rid <- index_table$RID[i]; mri_date <- index_table[[start_col]][i]
    z <- clinical[.(rid)]
    if (!nrow(z) || is.na(mri_date)) {
      out[[i]] <- data.table(RID = rid, linked = FALSE, label = label)
      next
    }
    near <- copy(z[abs(as.numeric(exam_date - mri_date)) <= window_days])
    if (!nrow(near)) {
      out[[i]] <- data.table(RID = rid, linked = FALSE, label = label)
      next
    }
    near[, `:=`(gap_abs = abs(as.numeric(exam_date - mri_date)),
                after_mri = exam_date < mri_date)]
    setorder(near, gap_abs, after_mri, exam_date)
    b <- near[1L]
    risk_start <- safe_idate_max(mri_date, b$exam_date)
    prevalent <- b$dx_stage == "Dementia" || any(z$exam_date <= risk_start & z$dx_stage == "Dementia")
    post <- copy(z[exam_date > risk_start])
    setorder(post, exam_date)
    event_rows <- which(post$dx_stage == "Dementia")
    first_event_i <- if (length(event_rows)) event_rows[1L] else NA_integer_
    event_date <- if (is.finite(first_event_i)) post$exam_date[first_event_i] else as.IDate(NA)
    next_dx <- if (is.finite(first_event_i) && first_event_i < nrow(post)) post$dx_stage[first_event_i + 1L] else NA_character_
    confirmed_next <- !is.na(next_dx) && next_dx == "Dementia"
    confirmed_ever <- if (is.finite(first_event_i) && first_event_i < nrow(post)) {
      any(post$dx_stage[(first_event_i + 1L):nrow(post)] == "Dementia")
    } else FALSE
    confirmed_next_date <- if (confirmed_next) event_date else as.IDate(NA)
    confirmed_ever_date <- if (confirmed_ever) event_date else as.IDate(NA)
    last_valid_date <- max(z$exam_date)
    previous_non_dementia <- if (!is.na(event_date)) {
      q <- z[exam_date < event_date & dx_stage != "Dementia", exam_date]
      if (length(q)) max(q) else risk_start
    } else as.IDate(NA)
    midpoint_date <- if (!is.na(event_date)) safe_idate_mid(previous_non_dementia, event_date) else as.IDate(NA)
    matched_age <- num(b$AGE) + ifelse(is.finite(num(b$Years_bl)), num(b$Years_bl), 0)
    out[[i]] <- data.table(
      RID = rid, linked = TRUE, label = label, window_days = window_days,
      mri_index_date = mri_date, baseline_dx_date = b$exam_date,
      baseline_dx_gap_days = as.numeric(b$exam_date - mri_date), risk_start = risk_start,
      baseline_dx = b$dx_stage, prevalent_dementia = prevalent,
      event_date = event_date, confirmed_next_date = confirmed_next_date,
      confirmed_ever_date = confirmed_ever_date, next_dx_after_event = next_dx,
      censor_date = last_valid_date, previous_non_dementia_date = previous_non_dementia,
      midpoint_event_date = midpoint_date, n_valid_dx = nrow(z),
      age_at_dx = matched_age, sex_clinical = as.character(b$PTGENDER),
      education_clinical = num(b$PTEDUCAT), apoe4 = num(b$APOE4),
      site = as.character(b$SITE), baseline_mPACC = num(b$mPACCdigit),
      baseline_CDRSB = num(b$CDRSB), baseline_FAQ = num(b$FAQ),
      baseline_MMSE = num(b$MMSE), baseline_ADAS13 = num(b$ADAS13)
    )
  }
  rbindlist(out, fill = TRUE)
}

# -----------------------------------------------------------------------------
# Index cohort and event-time variables.
# -----------------------------------------------------------------------------
event90 <- link_clinical_event(base, "index_mri_date", 90L, "baseline")
analysis <- merge(base, event90, by = "RID", all.x = TRUE, sort = FALSE)
analysis[, `:=`(
  baseline_dx = droplevels(factor(clean_stage(baseline_dx))),
  sex = droplevels(factor(ifelse(!is.na(sex_clinical) & nzchar(sex_clinical), sex_clinical, as.character(pgs_sex)))),
  education = ifelse(is.finite(education_clinical), education_clinical, num(pgs_education)),
  age_index = age_at_dx,
  scanner = droplevels(factor(paste(FIELD_STRENGTH, FSVER, sep = "_"))),
  site = droplevels(factor(site)),
  structure_vulnerability_z = num(structure_vulnerability_z),
  global_vulnerability_z = num(global_vulnerability_z)
)]
for (v in paste0("PC", 1:10)) if (!v %in% names(analysis)) analysis[, (v) := num(get(paste0("clinical_", v)))]
analysis[, `:=`(
  gMAV = num(gMAV), gAD = num(gAD), apoe4 = num(apoe4),
  age_index_z = zscore(age_index), education_z = zscore(education), apoe4_z = zscore(apoe4),
  baseline_mPACC_z = zscore(baseline_mPACC), baseline_CDRSB_z = zscore(baseline_CDRSB),
  baseline_FAQ_z = zscore(baseline_FAQ), baseline_MMSE_z = zscore(baseline_MMSE),
  baseline_ADAS13_z = zscore(baseline_ADAS13),
  pc1_z = zscore(PC1), pc2_z = zscore(PC2), pc3_z = zscore(PC3), pc4_z = zscore(PC4),
  pc5_z = zscore(PC5), pc6_z = zscore(PC6), pc7_z = zscore(PC7), pc8_z = zscore(PC8),
  pc9_z = zscore(PC9), pc10_z = zscore(PC10)
)]
analysis[, `:=`(
  followup_years = (fifelse(!is.na(event_date), as.numeric(event_date), as.numeric(censor_date)) - as.numeric(risk_start)) / 365.25,
  event = !is.na(event_date),
  followup_confirmed_years = (fifelse(!is.na(confirmed_ever_date), as.numeric(confirmed_ever_date), as.numeric(censor_date)) - as.numeric(risk_start)) / 365.25,
  event_confirmed = !is.na(confirmed_ever_date),
  followup_next_confirmed_years = (fifelse(!is.na(confirmed_next_date), as.numeric(confirmed_next_date), as.numeric(censor_date)) - as.numeric(risk_start)) / 365.25,
  event_next_confirmed = !is.na(confirmed_next_date),
  followup_midpoint_years = (fifelse(!is.na(midpoint_event_date), as.numeric(midpoint_event_date), as.numeric(censor_date)) - as.numeric(risk_start)) / 365.25,
  event_midpoint = !is.na(midpoint_event_date)
)]
analysis <- analysis[linked == TRUE & baseline_dx %in% c("CN", "MCI") & prevalent_dementia == FALSE]

pc_terms <- paste0("pc", 1:10, "_z")
base_covars <- c("age_index_z", "sex", "education_z", "apoe4_z", "scanner", pc_terms)
clinical_terms <- c("baseline_mPACC_z", "baseline_CDRSB_z", "baseline_FAQ_z")
rhs_suffix <- function(extra = character(), stage = "baseline_dx") {
  paste(c(extra, base_covars, paste0("strata(", stage, ")"), "cluster(site)"), collapse = " + ")
}

# Apparent C-index comparisons use a common complete-case sample.
all_core_vars <- unique(c("followup_years", "event", "structure_vulnerability_z", "gMAV", "gAD", base_covars, "baseline_dx", "site"))
common_core <- analysis[complete.cases(analysis[, all_core_vars, with = FALSE])]
all_benchmark_vars <- unique(c(all_core_vars, clinical_terms))
common_benchmark <- analysis[complete.cases(analysis[, all_benchmark_vars, with = FALSE])]

model_specs <- list(
  clinical_core = list(data = common_core, rhs = rhs_suffix(), endpoint = "first_detection"),
  clinical_plus_PGS = list(data = common_core, rhs = rhs_suffix(c("gMAV", "gAD")), endpoint = "first_detection"),
  clinical_plus_structure = list(data = common_core, rhs = rhs_suffix("structure_vulnerability_z"), endpoint = "first_detection"),
  clinical_full = list(data = common_core, rhs = rhs_suffix(c("structure_vulnerability_z", "gMAV", "gAD")), endpoint = "first_detection"),
  benchmark_core = list(data = common_benchmark, rhs = rhs_suffix(clinical_terms), endpoint = "first_detection"),
  benchmark_plus_PGS = list(data = common_benchmark, rhs = rhs_suffix(c(clinical_terms, "gMAV", "gAD")), endpoint = "first_detection"),
  benchmark_plus_structure = list(data = common_benchmark, rhs = rhs_suffix(c(clinical_terms, "structure_vulnerability_z")), endpoint = "first_detection"),
  benchmark_full = list(data = common_benchmark, rhs = rhs_suffix(c(clinical_terms, "structure_vulnerability_z", "gMAV", "gAD")), endpoint = "first_detection")
)

fits <- list()
for (id in names(model_specs)) fits[[id]] <- fit_cox(model_specs[[id]]$data, model_specs[[id]]$rhs, id)

# -----------------------------------------------------------------------------
# Prespecified robustness: confirmed DX, midpoint timing, lag, horizon, and
# +/-60/120-day MRI-DX linkage windows.
# -----------------------------------------------------------------------------
fit_endpoint_variant <- function(dat, time_col, event_col, rhs, id) {
  x <- copy(dat)
  x[, `:=`(followup_years = num(get(time_col)), event = as.logical(get(event_col)))]
  fit_cox(x, rhs, id)
}
full_rhs <- rhs_suffix(c("structure_vulnerability_z", "gMAV", "gAD"))
fits[["confirmed_ever_sensitivity"]] <- fit_endpoint_variant(common_core, "followup_confirmed_years", "event_confirmed", full_rhs, "confirmed_ever_sensitivity")
fits[["confirmed_next_sensitivity"]] <- fit_endpoint_variant(common_core, "followup_next_confirmed_years", "event_next_confirmed", full_rhs, "confirmed_next_sensitivity")
fits[["midpoint_sensitivity"]] <- fit_endpoint_variant(common_core, "followup_midpoint_years", "event_midpoint", full_rhs, "midpoint_sensitivity")

lag1 <- common_core[followup_years > 1]
lag1[, followup_years := followup_years - 1]
fits[["exclude_first_1y"]] <- fit_cox(lag1, full_rhs, "exclude_first_1y")
lag2 <- common_core[followup_years > 2]
lag2[, followup_years := followup_years - 2]
fits[["exclude_first_2y"]] <- fit_cox(lag2, full_rhs, "exclude_first_2y")

horizon5 <- copy(common_core)
horizon5[, `:=`(event = event & followup_years <= 5, followup_years = pmin(followup_years, 5))]
fits[["five_year_horizon"]] <- fit_cox(horizon5, full_rhs, "five_year_horizon")

rebuild_window <- function(window_days) {
  ev <- link_clinical_event(base, "index_mri_date", window_days, paste0("window_", window_days))
  x <- merge(base, ev, by = "RID", all.x = FALSE, sort = FALSE)
  x[, `:=`(
    baseline_dx = droplevels(factor(clean_stage(baseline_dx))),
    sex = factor(ifelse(!is.na(sex_clinical) & nzchar(sex_clinical), sex_clinical, as.character(pgs_sex))),
    education = ifelse(is.finite(education_clinical), education_clinical, num(pgs_education)),
    age_index = age_at_dx, scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_")),
    site = factor(site), gMAV = num(gMAV), gAD = num(gAD), apoe4 = num(apoe4),
    followup_years = (fifelse(!is.na(event_date), as.numeric(event_date), as.numeric(censor_date)) - as.numeric(risk_start)) / 365.25,
    event = !is.na(event_date)
  )]
  for (v in paste0("PC", 1:10)) x[, (v) := num(get(paste0("clinical_", v)))]
  x[, `:=`(age_index_z = zscore(age_index), education_z = zscore(education), apoe4_z = zscore(apoe4),
            pc1_z = zscore(PC1), pc2_z = zscore(PC2), pc3_z = zscore(PC3), pc4_z = zscore(PC4),
            pc5_z = zscore(PC5), pc6_z = zscore(PC6), pc7_z = zscore(PC7), pc8_z = zscore(PC8),
            pc9_z = zscore(PC9), pc10_z = zscore(PC10))]
  x <- x[linked == TRUE & baseline_dx %in% c("CN", "MCI") & prevalent_dementia == FALSE]
  vars <- unique(c("followup_years", "event", "structure_vulnerability_z", "gMAV", "gAD", base_covars, "baseline_dx", "site"))
  x[complete.cases(x[, vars, with = FALSE])]
}
for (w in c(60L, 120L)) {
  xw <- rebuild_window(w)
  fits[[paste0("window_", w, "d")]] <- fit_cox(xw, full_rhs, paste0("window_", w, "d"))
}

# -----------------------------------------------------------------------------
# Two-MRI landmark: atrophy rate -> subsequent first Dementia detection.
# -----------------------------------------------------------------------------
pair <- d[, {
  dates <- sort(unique(mri_date)); eligible <- if (length(dates) >= 2L) dates[as.integer(dates - dates[1L]) >= 270L] else as.IDate(character())
  if (!length(eligible)) NULL else {
    d0 <- dates[1L]; d1 <- eligible[1L]; a <- .SD[mri_date == d0][1L]; b <- .SD[mri_date == d1][1L]
    interval <- as.numeric(d1 - d0) / 365.25
    data.table(index_mri_date = d0, landmark_mri_date = d1, interval_years = interval,
               baseline_structure_vulnerability_z = a$structure_vulnerability_z,
               atrophy_rate = -(num(b$ad_vulnerable_index) - num(a$ad_vulnerable_index)) / interval,
               gMAV = a$gMAV, gAD = a$gAD, scanner = as.character(a$scanner),
               PC1 = a$PC1, PC2 = a$PC2, PC3 = a$PC3, PC4 = a$PC4, PC5 = a$PC5,
               PC6 = a$PC6, PC7 = a$PC7, PC8 = a$PC8, PC9 = a$PC9, PC10 = a$PC10,
               pgs_sex = a$pgs_sex, pgs_education = a$pgs_education)
  }
}, by = RID]
pair_event <- link_clinical_event(pair, "landmark_mri_date", 90L, "two_MRI_landmark")
pair <- merge(pair, pair_event, by = "RID", all.x = FALSE, sort = FALSE)
pair[, `:=`(
  atrophy_rate_z = zscore(atrophy_rate), baseline_dx = factor(clean_stage(baseline_dx)),
  sex = factor(ifelse(!is.na(sex_clinical) & nzchar(sex_clinical), sex_clinical, pgs_sex)),
  education = ifelse(is.finite(education_clinical), education_clinical, num(pgs_education)),
  age_index_z = zscore(age_at_dx), education_z = zscore(ifelse(is.finite(education_clinical), education_clinical, num(pgs_education))),
  apoe4_z = zscore(apoe4), scanner = factor(scanner), site = factor(site),
  followup_years = (fifelse(!is.na(event_date), as.numeric(event_date), as.numeric(censor_date)) - as.numeric(risk_start)) / 365.25,
  event = !is.na(event_date)
)]
for (v in paste0("PC", 1:10)) pair[, (paste0(tolower(v), "_z")) := zscore(get(v))]
pair <- pair[linked == TRUE & baseline_dx %in% c("CN", "MCI") & prevalent_dementia == FALSE]
landmark_vars <- unique(c("followup_years", "event", "atrophy_rate_z", "baseline_structure_vulnerability_z",
                          "gMAV", "gAD", base_covars, "baseline_dx", "site"))
pair_common <- pair[complete.cases(pair[, landmark_vars, with = FALSE])]
landmark_rhs <- rhs_suffix(c("atrophy_rate_z", "baseline_structure_vulnerability_z", "gMAV", "gAD"))
fits[["two_MRI_landmark"]] <- fit_cox(pair_common, landmark_rhs, "two_MRI_landmark")

# -----------------------------------------------------------------------------
# Synthesis, nested-model contrasts, and reproducible outputs.
# -----------------------------------------------------------------------------
coefficients <- rbindlist(lapply(fits, `[[`, "coefficients"), fill = TRUE)
diagnostics <- rbindlist(lapply(fits, `[[`, "diagnostics"), fill = TRUE)
ph_results <- rbindlist(lapply(fits, `[[`, "ph_table"), fill = TRUE)

lrt_pair <- function(smaller, larger, label) {
  a <- fits[[smaller]]$fit; b <- fits[[larger]]$fit
  na <- fits[[smaller]]$diagnostics$n[1L]
  nb <- fits[[larger]]$diagnostics$n[1L]
  if (is.null(a) || is.null(b) || !is.finite(na) || !is.finite(nb) || na != nb) {
    return(data.table(comparison = label, smaller = smaller, larger = larger,
                      n = NA_integer_, events = NA_integer_, lr_chisq = NA_real_, df = NA_integer_, p_value = NA_real_))
  }
  lr <- 2 * (as.numeric(logLik(b)) - as.numeric(logLik(a)))
  df <- length(coef(b)) - length(coef(a))
  data.table(comparison = label, smaller = smaller, larger = larger,
             n = nb, events = b$nevent, lr_chisq = lr, df = df,
             p_value = pchisq(lr, df = df, lower.tail = FALSE))
}
lrt_results <- rbindlist(list(
  lrt_pair("clinical_core", "clinical_plus_PGS", "PGS_increment_over_clinical_core"),
  lrt_pair("clinical_core", "clinical_plus_structure", "structure_increment_over_clinical_core"),
  lrt_pair("clinical_core", "clinical_full", "structure_and_PGS_increment_over_clinical_core"),
  lrt_pair("benchmark_core", "benchmark_plus_PGS", "PGS_increment_over_strong_benchmark"),
  lrt_pair("benchmark_core", "benchmark_plus_structure", "structure_increment_over_strong_benchmark"),
  lrt_pair("benchmark_core", "benchmark_full", "structure_and_PGS_increment_over_strong_benchmark")
))

performance <- diagnostics[, .(model_id, n, events, c_index, AIC, logLik, ph_global_p, status)]
core_c <- performance[model_id == "clinical_core", c_index]
benchmark_c <- performance[model_id == "benchmark_core", c_index]
performance[, delta_c_vs_tier_core := fifelse(grepl("^clinical", model_id) & length(core_c), c_index - core_c,
                                               fifelse(grepl("^benchmark", model_id) & length(benchmark_c), c_index - benchmark_c, NA_real_))]

sample_flow <- rbindlist(list(
  data.table(stage = "PGS_complete_first_MRI", people = nrow(base), events = NA_integer_),
  data.table(stage = "MRI_DX_linked_within_90d", people = sum(event90$linked %in% TRUE), events = NA_integer_),
  data.table(stage = "eligible_CN_MCI_after_prevalent_exclusion", people = nrow(analysis), events = sum(analysis$event)),
  data.table(stage = "common_core_model", people = nrow(common_core), events = sum(common_core$event)),
  data.table(stage = "strong_clinical_benchmark", people = nrow(common_benchmark), events = sum(common_benchmark$event)),
  data.table(stage = "confirmed_ever_later_Dementia", people = nrow(common_core), events = sum(common_core$event_confirmed)),
  data.table(stage = "confirmed_next_Dementia", people = nrow(common_core), events = sum(common_core$event_next_confirmed)),
  data.table(stage = "two_MRI_landmark", people = nrow(pair_common), events = sum(pair_common$event))
))

event_confirmation <- common_core[event == TRUE, .N, by = .(
  next_dx_after_event,
  ever_later_Dementia = event_confirmed,
  next_Dementia_confirmed = event_next_confirmed
)][order(-N)]
interval_days <- as.numeric(common_core[event == TRUE, event_date - previous_non_dementia_date])
interval_summary <- data.table(
  n_events = length(interval_days), median_days = median(interval_days, na.rm = TRUE),
  q1_days = quantile(interval_days, .25, na.rm = TRUE),
  q3_days = quantile(interval_days, .75, na.rm = TRUE),
  min_days = min(interval_days, na.rm = TRUE), max_days = max(interval_days, na.rm = TRUE)
)

fwrite(coefficients, file.path(out_dir, "01_cox_coefficients.csv"), bom = TRUE)
fwrite(diagnostics, file.path(out_dir, "02_model_diagnostics.csv"), bom = TRUE)
fwrite(ph_results, file.path(out_dir, "03_proportional_hazards_tests.csv"), bom = TRUE)
fwrite(performance, file.path(out_dir, "04_model_performance_apparent.csv"), bom = TRUE)
fwrite(lrt_results, file.path(out_dir, "05_nested_model_likelihood_ratio_tests.csv"), bom = TRUE)
fwrite(sample_flow, file.path(out_dir, "06_sample_flow.csv"), bom = TRUE)
fwrite(event_confirmation, file.path(out_dir, "07_event_confirmation_audit.csv"), bom = TRUE)
fwrite(interval_summary, file.path(out_dir, "08_event_interval_summary.csv"), bom = TRUE)
fwrite(analysis, file.path(out_dir, "09_baseline_event_frame.csv"), bom = TRUE)
fwrite(pair, file.path(out_dir, "10_two_MRI_landmark_frame.csv"), bom = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
message("ADNI event models completed: ", out_dir)
