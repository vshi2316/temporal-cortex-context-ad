# Endpoint-specific pathology and memory associations
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
suppressPackageStartupMessages({ library(data.table); library(lme4) })

pick_root <- function(ascii, native) if (file.exists(ascii)) ascii else native
project_dir <- pick_root(study_path("project"), study_path("project"))
result_root <- file.path(project_dir, "results")
repair_dir <- file.path(result_root, "13F_HABS_linkage_repair")
out_dir <- file.path(result_root, "13F_HABS_repaired_external_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
need <- function(path) { if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE); path }
num <- function(x) { z <- suppressWarnings(as.numeric(as.character(x))); z[z <= -9000] <- NA_real_; z }
finite <- function(x) is.finite(num(x))
zscore <- function(x) { x <- num(x); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (x - mean(x, na.rm = TRUE)) / s }
normalize_id <- function(x) { z <- trimws(as.character(x)); n <- suppressWarnings(as.numeric(z)); z[is.finite(n)] <- format(n[is.finite(n)], scientific = FALSE, trim = TRUE, digits = 15); z <- sub("\\.0+$", "", z); z[z %in% c("", "NA", "N/A", ".")] <- NA_character_; z }
month_num <- function(x) { y <- suppressWarnings(as.integer(x)); bad <- !is.finite(y); if (any(bad)) y[bad] <- match(tolower(substr(trimws(as.character(x[bad])), 1, 3)), tolower(month.abb)); y }
make_date <- function(y, m) { yy <- suppressWarnings(as.integer(y)); mm <- month_num(m); out <- as.IDate(rep(NA_character_, length(yy))); ok <- is.finite(yy) & is.finite(mm) & mm >= 1 & mm <= 12; out[ok] <- as.IDate(sprintf("%04d-%02d-15", yy[ok], mm[ok])); out }

z <- readRDS(need(file.path(repair_dir, "02_repaired_HABSHD_linked_frames.rds")))
h <- as.data.table(z$linked)
h[, participant_id := normalize_id(participant_id)]
h[, mri_date := as.IDate(as.character(mri_date))]
h[, scanner_factor := factor(scanner)]
setorder(h, participant_id, mri_date)
h <- h[!is.na(participant_id) & !is.na(mri_date) & finite(ad_vulnerable_index)]
h <- h[, .SD[1L], by = .(participant_id, mri_date)]
setorder(h, participant_id, mri_date)
h[, first_mri_date := min(mri_date), by = participant_id]
h[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
h[, `:=`(ad_vulnerable_z = zscore(ad_vulnerable_index), global_cortical_z = zscore(global_cortical_thickness), hippocampus_z = zscore(hip_hippocampus_icv_ratio), brain_z = zscore(brain_brain_icv_ratio))]

# Baseline demographic covariates, selected independently of pathology and
# cognition from the earliest available interview/residence record.
demo <- fread(file = file.path(pick_root(study_path("habs_release"), study_path("cohorts", "0MRI")), "interviewdemographics_12Aug2026.csv"), showProgress = FALSE)
demo[, participant_id := normalize_id(MED_ID)]
demo[, visit_date__ := make_date(VISIT_YEAR, VISIT_MONTH)]
setorder(demo, participant_id, visit_date__)
demo <- demo[!is.na(participant_id), .SD[1L], by = participant_id]
demo[, sex := factor(as.character(ID_GENDER))]
demo[, race := fifelse(num(ID_RACE_BLACK) == 1, "Black", fifelse(num(ID_RACE_ASIAN) == 1, "Asian", fifelse(num(ID_RACE_WHITE) == 1, "White", "Other")))]
demo <- demo[, .(participant_id, sex, race)]
h <- merge(h, demo, by = "participant_id", all.x = TRUE, sort = FALSE)

# Same-visit cognitive scores define the cross-sectional structure/cognition association.
zs <- fread(file = file.path(pick_root(study_path("habs_release"), study_path("cohorts", "0MRI")), "zscores_12Aug2026.csv"), showProgress = FALSE)
zs[, participant_id := normalize_id(MED_ID)]
zs[, cog_date := make_date(VISIT_YEAR, VISIT_MONTH)]
mem_cols <- intersect(c("LM1_AB_ZSCORE", "LM2_AB_ZSCORE", "SEVLT_T1235_ZSCORE", "SEVLT_DR_ZSCORE"), names(zs))
for (v in mem_cols) set(zs, j = v, value = num(zs[[v]]))
zs[, memory_composite := rowMeans(.SD, na.rm = TRUE), .SDcols = mem_cols]
zs[!is.finite(memory_composite), memory_composite := NA_real_]
zs <- zs[is.finite(memory_composite), .(participant_id, cog_date, memory_composite)]

nearest <- function(index, candidate, date_col, value_col, prefix, max_days = 365L) {
  a <- copy(index); b <- copy(candidate)
  a[, row_id__ := .I]; setnames(b, date_col, "candidate_date__")
  q <- merge(a, b, by = "participant_id", allow.cartesian = TRUE, all.x = TRUE)
  q[, gap_days__ := abs(as.integer(mri_date - candidate_date__))]
  setorder(q, row_id__, gap_days__); q <- q[, head(.SD, 1L), by = row_id__]
  q[!is.finite(gap_days__) | gap_days__ > max_days, (value_col) := NA_real_]
  setnames(q, c("candidate_date__", value_col), c(paste0(prefix, "date"), paste0(prefix, value_col)))
  q[, c("row_id__", "gap_days__") := NULL]
  q
}
h <- nearest(h, zs, "cog_date", "memory_composite", "cog_", 365L)
h[, memory_z := zscore(cog_memory_composite)]

# Cross-sectional models use one deterministic MRI record per participant to
# avoid treating repeated scans as independent observations. Repeated scans
# are reserved for the mixed-effects sensitivity model below.
index_one <- function(x, value_col) {
  x <- x[finite(get(value_col))]
  if (!nrow(x)) return(x)
  setorder(x, participant_id, mri_date)
  x[, .SD[1L], by = participant_id]
}

extract_lm <- function(fit, family, outcome, exposure, n) {
  sm <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(sm, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), c("estimate", "std_error", "statistic", "p_value"), skip_absent = TRUE)
  sm[, `:=`(family = family, outcome = outcome, exposure = exposure, n = n)]
  sm[, .(family, outcome, exposure, term, estimate, std_error, statistic, p_value, n)]
}
fit_lm_one <- function(data, outcome, exposure, family) {
  vars <- c(outcome, exposure, "age", "sex", "race", "scanner_factor")
  vars <- vars[vars %in% names(data)]
  a <- data[complete.cases(data[, vars, with = FALSE])]
  if (nrow(a) < 200L || !any(finite(a[[outcome]])) || !any(finite(a[[exposure]]))) return(NULL)
  a <- droplevels(a)
  covars <- c("age", "sex", "race", "scanner_factor"); covars <- covars[covars %in% names(a)]
  covars <- covars[vapply(a[, covars, with = FALSE], function(x) if (is.factor(x)) nlevels(x) >= 2 else length(unique(x[is.finite(x)])) >= 2, logical(1))]
  f <- lm(as.formula(paste(outcome, "~", paste(c(exposure, covars), collapse = " + "))), data = a)
  extract_lm(f, family, outcome, exposure, nrow(a))
}

h[, age := num(age)]
path_results <- rbindlist(list(
  fit_lm_one(index_one(h, "amyloid_amyloid_value"), "ad_vulnerable_z", "amyloid_amyloid_value", "amyloid_to_fixed_structure"),
  fit_lm_one(index_one(h, "tau_tau_value"), "ad_vulnerable_z", "tau_tau_value", "tau_to_fixed_structure"),
  fit_lm_one(index_one(h, "amyloid_amyloid_value"), "global_cortical_z", "amyloid_amyloid_value", "amyloid_to_global_structure"),
  fit_lm_one(index_one(h, "tau_tau_value"), "global_cortical_z", "tau_tau_value", "tau_to_global_structure")
), fill = TRUE)
if (nrow(path_results)) path_results[, fdr := p.adjust(p_value, method = "BH")]

cog_results <- fit_lm_one(index_one(h, "memory_z"), "memory_z", "ad_vulnerable_z", "fixed_structure_to_memory")
if (!is.null(cog_results) && nrow(cog_results)) cog_results[, fdr := p.adjust(p_value, method = "BH")]

# Longitudinal external sensitivity: earliest amyloid value per participant is
# attached to all subsequent MRI rows; interpretation remains associative.
amy_base <- h[finite(amyloid_amyloid_value) & !is.na(amyloid_date), .SD[which.min(amyloid_date)], by = participant_id]
amy_base <- amy_base[, .(participant_id, amyloid_baseline = amyloid_amyloid_value, amyloid_baseline_date = amyloid_date)]
hl <- merge(h, amy_base, by = "participant_id", all = FALSE, sort = FALSE)
hl[, amyloid_baseline_z := zscore(amyloid_baseline)]
hl[, time_since_amyloid := as.numeric(mri_date - amyloid_baseline_date) / 365.25]
hl <- hl[finite(ad_vulnerable_z) & finite(amyloid_baseline_z) & finite(time_since_amyloid) & time_since_amyloid >= 0 & time_since_amyloid <= 5]
longitudinal_fit <- NULL
longitudinal_result <- data.table()
if (uniqueN(hl$participant_id) >= 200L) {
  hl <- droplevels(hl)
  covars <- c("age", "sex", "race", "scanner_factor"); covars <- covars[covars %in% names(hl)]
  covars <- covars[vapply(hl[, covars, with = FALSE], function(x) if (is.factor(x)) nlevels(x) >= 2 else length(unique(x[is.finite(x)])) >= 2, logical(1))]
  rhs <- paste(c("time_since_amyloid * amyloid_baseline_z", covars), collapse = " + ")
  f1 <- as.formula(paste0("ad_vulnerable_z ~ ", rhs, " + (1 + time_since_amyloid | participant_id)"))
  f0 <- as.formula(paste0("ad_vulnerable_z ~ ", rhs, " + (1 | participant_id)"))
  longitudinal_fit <- tryCatch(lmer(f1, data = hl, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))), error = function(e) e)
  random_structure <- "random_intercept_and_time_slope"
  if (inherits(longitudinal_fit, "error") || isSingular(longitudinal_fit, tol = 1e-4)) { longitudinal_fit <- tryCatch(lmer(f0, data = hl, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))), error = function(e) e); random_structure <- "random_intercept_fallback" }
  if (!inherits(longitudinal_fit, "error")) {
    longitudinal_result <- as.data.table(summary(longitudinal_fit)$coefficients, keep.rownames = "term")
    setnames(longitudinal_result, c("Estimate", "Std. Error", "t value"), c("estimate", "std_error", "statistic"), skip_absent = TRUE)
    longitudinal_result[, `:=`(family = "longitudinal_amyloid_structure", outcome = "ad_vulnerable_z", p_value = 2 * pnorm(-abs(statistic)), n_people = uniqueN(hl$participant_id), n_rows = nrow(hl), random_structure = random_structure)]
  }
}

flow <- data.table(
  set = c("cortical_repaired_rows", "cortical_repaired_people", "amyloid_linked_people_548d", "tau_linked_people_548d", "hippocampus_linked_people_548d", "cognition_linked_rows", "longitudinal_amyloid_people"),
  value = c(nrow(h), uniqueN(h$participant_id), uniqueN(h[finite(amyloid_amyloid_value) & !is.na(amyloid_date), participant_id]), uniqueN(h[finite(tau_tau_value) & !is.na(tau_date), participant_id]), uniqueN(h[finite(hip_hippocampus_icv_ratio) & !is.na(hip_date), participant_id]), sum(finite(h$memory_z)), uniqueN(hl$participant_id))
)
scope <- c(
  "HABS-HD models use a shared identifier/date rule across cortical measures and the linked pathology, volume and cognition files.",
  "Regions and directions are frozen from 12B; no region or weight was re-selected using outcomes.",
  "Cross-sectional models are external biological associations adjusted for available age, sex, race and scanner.",
  "The longitudinal sensitivity model links earliest amyloid to repeated MRI measurements."
)
fwrite(path_results, file.path(out_dir, "01_repaired_pathology_structure_models.csv"), bom = TRUE)
fwrite(cog_results, file.path(out_dir, "02_repaired_structure_cognition_models.csv"), bom = TRUE)
fwrite(longitudinal_result, file.path(out_dir, "03_repaired_longitudinal_amyloid_structure_model.csv"), bom = TRUE)
fwrite(flow, file.path(out_dir, "04_repaired_external_flow.csv"), bom = TRUE)
writeLines(scope, file.path(out_dir, "05_interpretation_scope.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
message("HABS repaired external models completed: ", out_dir)
