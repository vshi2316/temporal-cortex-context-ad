# Direct genetic contrast and 679-person bootstrap optimism benchmark
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(survival)
})

# Direct genetic time contrast and bootstrap-corrected five-year performance.

project_dir <- study_path("project")
result_root <- file.path(project_dir, "results")
out_dir <- file.path(result_root, "13K_core_claim_completion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

B <- getOption("AD_13K_BOOTSTRAP_B", 500L)
seed <- getOption("AD_13K_SEED", 20260914L)
horizon <- 5

need <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  path
}

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

zscore <- function(x) {
  x <- num(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

has_variation <- function(x) {
  if (is.factor(x)) return(nlevels(droplevels(x)) >= 2L)
  length(unique(x[is.finite(num(x))])) >= 2L
}

write_csv <- function(x, name) {
  fwrite(as.data.table(x), file.path(out_dir, name), bom = TRUE, na = "NA")
}

# -----------------------------------------------------------------------------
# 1. Direct gAD-versus-gMAV longitudinal contrast
# -----------------------------------------------------------------------------

structural_file <- file.path(
  result_root, "13F_R_only_structural_analysis",
  "04_ADNI_R_only_longitudinal_structural_frame.rds"
)
pgs_file <- file.path(
  result_root, "10A_ADNI_gMAV_PGS_assembly",
  "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"
)

d <- as.data.table(readRDS(need(structural_file)))
pgs <- fread(need(pgs_file), sep = "\t", showProgress = FALSE)

d[, `:=`(
  RID = as.character(RID),
  mri_date = as.IDate(as.character(mri_date)),
  gMAV = num(clinical_gMAV_noAPOE_z_eur),
  gAD = num(clinical_gAD_noAPOE_z_eur)
)]

pgs[, `:=`(
  RID = as.character(RID),
  exam_date__ = as.IDate(as.character(EXAMDATE))
)]
setorder(pgs, RID, exam_date__)

first_num <- function(x) {
  z <- num(x)
  if (any(is.finite(z))) z[which(is.finite(z))[1L]] else NA_real_
}
first_chr <- function(x) {
  z <- trimws(as.character(x))
  z <- z[!is.na(z) & nzchar(z)]
  if (length(z)) z[1L] else NA_character_
}

base_cov <- pgs[, .(
  baseline_age = first_num(age_bl),
  baseline_sex = first_chr(sex),
  baseline_education = first_num(education),
  baseline_diagnosis = first_chr(dx_bl)
), by = RID]

d <- merge(d, base_cov, by = "RID", all.x = TRUE, sort = FALSE)
d[, `:=`(
  baseline_age = num(baseline_age),
  sex = factor(baseline_sex),
  education = num(baseline_education),
  diagnosis = factor(baseline_diagnosis),
  scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_")),
  primary_complete = is.finite(num(ad_vulnerable_index))
)]

pc_cols <- paste0("clinical_PC", 1:10)
for (v in pc_cols) d[, (v) := num(get(v))]

setorderv(
  d,
  c("RID", "mri_date", "primary_complete", "IMAGEUID"),
  order = c(1, 1, -1, -1), na.last = TRUE
)
d <- d[, .SD[1L], by = .(RID, mri_date)]
setorder(d, RID, mri_date)
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
d[, ad_vulnerable_z := zscore(ad_vulnerable_index)]
d <- d[
  is.finite(ad_vulnerable_z) & is.finite(gMAV) & is.finite(gAD) &
    !is.na(mri_date)
]

elig <- d[, {
  dates <- sort(unique(mri_date))
  later <- dates[as.integer(dates - dates[1L]) >= 270L]
  list(n_dates = length(dates), has_270d_second = length(later) > 0L)
}, by = RID]
d <- d[RID %in% elig[n_dates >= 2L & has_270d_second == TRUE, RID]]

candidate_covars <- c(
  "baseline_age", "sex", "education", "diagnosis", "scanner", pc_cols
)
covars <- candidate_covars[vapply(candidate_covars, function(v) {
  v %in% names(d) && has_variation(d[[v]])
}, logical(1))]

model_vars <- unique(c(
  "ad_vulnerable_z", "RID", "time_years", "gMAV", "gAD", covars
))
model_data <- droplevels(d[complete.cases(d[, ..model_vars])])

fixed_rhs <- paste(
  c("time_years * gMAV", "time_years * gAD", covars),
  collapse = " + "
)
f_slope <- as.formula(paste0(
  "ad_vulnerable_z ~ ", fixed_rhs, " + (1 + time_years | RID)"
))
f_intercept <- as.formula(paste0(
  "ad_vulnerable_z ~ ", fixed_rhs, " + (1 | RID)"
))

fit_long <- tryCatch(
  lmer(
    f_slope, data = model_data, REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
    )
  ),
  error = function(e) e
)
random_structure <- "random_intercept_and_time_slope"
if (inherits(fit_long, "error") || isSingular(fit_long, tol = 1e-4)) {
  fit_long <- lmer(
    f_intercept, data = model_data, REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
    )
  )
  random_structure <- "random_intercept_fallback"
}

b <- fixef(fit_long)
V <- vcov(fit_long)
required_terms <- c("time_years:gAD", "time_years:gMAV")
if (!all(required_terms %in% names(b))) {
  stop("The fitted longitudinal model lacks the two required interaction terms.",
       call. = FALSE)
}

L <- setNames(rep(0, length(b)), names(b))
L["time_years:gAD"] <- 1
L["time_years:gMAV"] <- -1
contrast_est <- as.numeric(sum(L * b))
contrast_se <- sqrt(as.numeric(t(L) %*% V %*% L))
contrast_z <- contrast_est / contrast_se
contrast_p <- 2 * pnorm(-abs(contrast_z))

contrast <- data.table(
  contrast = "time_x_gAD_minus_time_x_gMAV",
  estimate_difference_SD_per_year = contrast_est,
  standard_error = contrast_se,
  conf_low = contrast_est - 1.96 * contrast_se,
  conf_high = contrast_est + 1.96 * contrast_se,
  z_value = contrast_z,
  p_value = contrast_p,
  n_people = uniqueN(model_data$RID),
  n_MRI_rows = nrow(model_data),
  random_structure = random_structure
)
write_csv(contrast, "01_gAD_vs_gMAV_direct_slope_contrast.csv")

key_long <- data.table(
  term = required_terms,
  estimate = b[required_terms],
  standard_error = sqrt(diag(V))[required_terms]
)
key_long[, `:=`(
  conf_low = estimate - 1.96 * standard_error,
  conf_high = estimate + 1.96 * standard_error,
  p_value = 2 * pnorm(-abs(estimate / standard_error))
)]
write_csv(key_long, "02_longitudinal_interaction_coefficients.csv")

time_grid <- seq(0, 5, by = 0.05)
base_time <- unname(b["time_years"])
trajectory <- rbindlist(list(
  rbindlist(lapply(c(-1, 0, 1), function(z) data.table(
    panel = "gAD-noAPOE",
    score = z,
    time_years = time_grid,
    predicted_change_z = time_grid *
      (base_time + z * b["time_years:gAD"])
  ))),
  rbindlist(lapply(c(-1, 0, 1), function(z) data.table(
    panel = "gMAV-noAPOE",
    score = z,
    time_years = time_grid,
    predicted_change_z = time_grid *
      (base_time + z * b["time_years:gMAV"])
  )))
))
write_csv(trajectory, "03_standardized_predicted_trajectories.csv")

# -----------------------------------------------------------------------------
# 2. Five-year ADNI incremental performance with optimism correction
# -----------------------------------------------------------------------------

event_file <- file.path(
  result_root, "13F_R1_clinical_closure_confirmatory",
  "09_baseline_event_frame.csv"
)
clinical <- fread(need(event_file), showProgress = FALSE)

clinical[, `:=`(
  linked = as.logical(linked),
  prevalent_dementia = as.logical(prevalent_dementia),
  event = as.integer(as.logical(event)),
  followup_years = num(followup_years),
  baseline_dx = factor(baseline_dx),
  sex = factor(sex),
  scanner = factor(scanner),
  structure_vulnerability_z = num(structure_vulnerability_z)
)]
clinical <- clinical[
  linked == TRUE & prevalent_dementia == FALSE &
    baseline_dx %in% c("CN", "MCI") & followup_years > 0
]

pc_terms <- paste0("pc", 1:10, "_z")
core_terms <- c(
  "age_index_z", "sex", "education_z", "apoe4_z", "scanner",
  "baseline_dx", pc_terms
)
benchmark_terms <- c(
  core_terms, "baseline_mPACC_z", "baseline_CDRSB_z", "baseline_FAQ_z"
)
all_prediction_vars <- unique(c(
  "followup_years", "event", "structure_vulnerability_z", benchmark_terms
))
pred_data <- droplevels(clinical[complete.cases(clinical[, ..all_prediction_vars])])

prediction_specs <- list(
  clinical_core = core_terms,
  clinical_plus_structure = c(core_terms, "structure_vulnerability_z"),
  strong_benchmark = benchmark_terms,
  strong_benchmark_plus_structure =
    c(benchmark_terms, "structure_vulnerability_z")
)

fit_prediction_model <- function(dat, terms) {
  f <- as.formula(paste0(
    "Surv(followup_years, event) ~ ", paste(terms, collapse = " + ")
  ))
  tryCatch(
    coxph(f, data = dat, ties = "efron", x = TRUE, y = TRUE, model = TRUE),
    error = function(e) NULL
  )
}

predict_risk <- function(fit, newdata, horizon_years) {
  if (is.null(fit) || any(!is.finite(coef(fit)))) {
    return(rep(NA_real_, nrow(newdata)))
  }
  bh <- basehaz(fit, centered = FALSE)
  h0 <- bh[bh$time <= horizon_years, "hazard"]
  h0 <- if (length(h0)) max(h0) else 0
  lp <- predict(fit, newdata = newdata, type = "lp", reference = "zero")
  pmin(pmax(1 - exp(-h0 * exp(lp)), 1e-8), 1 - 1e-8)
}

km_censor_survival <- function(time, event, at) {
  fit <- survfit(Surv(time, 1L - event) ~ 1)
  out <- summary(fit, times = pmax(at, 0), extend = TRUE)$surv
  pmax(out, 1e-6)
}

fixed_horizon_metrics <- function(dat, risk, horizon_years) {
  if (length(risk) != nrow(dat) || any(!is.finite(risk))) {
    return(c(AUC = NA_real_, Brier = NA_real_,
             calibration_intercept = NA_real_, calibration_slope = NA_real_))
  }
  time <- dat$followup_years
  event <- dat$event
  is_case <- event == 1L & time <= horizon_years
  is_control <- time > horizon_years
  y <- as.integer(is_case)
  w <- rep(0, nrow(dat))
  if (any(is_case)) {
    w[is_case] <- 1 / km_censor_survival(
      time, event, pmax(time[is_case] - 1e-8, 0)
    )
  }
  if (any(is_control)) {
    g_h <- km_censor_survival(time, event, horizon_years)[1L]
    w[is_control] <- 1 / g_h
  }

  case_risk <- risk[is_case]
  control_risk <- risk[is_control]
  case_w <- w[is_case]
  control_w <- w[is_control]
  if (!length(case_risk) || !length(control_risk)) {
    auc <- NA_real_
  } else {
    cmp <- outer(case_risk, control_risk, `-`)
    pair_w <- outer(case_w, control_w, `*`)
    auc <- sum(pair_w * (cmp > 0) + 0.5 * pair_w * (cmp == 0)) /
      sum(pair_w)
  }

  brier <- sum(w * (y - risk)^2) / nrow(dat)
  known <- w > 0 & is.finite(w)
  logit_risk <- qlogis(risk)

  cal_slope_fit <- tryCatch(
    suppressWarnings(glm(
      y ~ logit_risk, family = binomial(), weights = w,
      subset = known
    )),
    error = function(e) NULL
  )
  cal_int_fit <- tryCatch(
    suppressWarnings(glm(
      y ~ 1 + offset(logit_risk), family = binomial(), weights = w,
      subset = known
    )),
    error = function(e) NULL
  )

  c(
    AUC = auc,
    Brier = brier,
    calibration_intercept = if (is.null(cal_int_fit)) NA_real_ else
      unname(coef(cal_int_fit)[1L]),
    calibration_slope = if (is.null(cal_slope_fit)) NA_real_ else
      unname(coef(cal_slope_fit)["logit_risk"])
  )
}

apparent_fits <- lapply(prediction_specs, function(x) {
  fit_prediction_model(pred_data, x)
})
apparent_metrics <- rbindlist(lapply(names(apparent_fits), function(id) {
  risk <- predict_risk(apparent_fits[[id]], pred_data, horizon)
  m <- fixed_horizon_metrics(pred_data, risk, horizon)
  data.table(
    model_id = id, horizon_years = horizon,
    n = nrow(pred_data), events_by_horizon = sum(
      pred_data$event == 1L & pred_data$followup_years <= horizon
    ),
    metric = names(m), apparent = as.numeric(m)
  )
}))

set.seed(seed)
boot_rows <- vector("list", B * length(prediction_specs))
k <- 0L
for (bb in seq_len(B)) {
  index <- sample.int(nrow(pred_data), nrow(pred_data), replace = TRUE)
  boot_data <- pred_data[index]
  for (id in names(prediction_specs)) {
    k <- k + 1L
    fit_b <- fit_prediction_model(boot_data, prediction_specs[[id]])
    if (is.null(fit_b)) {
      boot_rows[[k]] <- data.table(
        replicate = bb, model_id = id, metric = c(
          "AUC", "Brier", "calibration_intercept", "calibration_slope"
        ), bootstrap_apparent = NA_real_, original_test = NA_real_
      )
      next
    }
    m_boot <- fixed_horizon_metrics(
      boot_data, predict_risk(fit_b, boot_data, horizon), horizon
    )
    m_test <- fixed_horizon_metrics(
      pred_data, predict_risk(fit_b, pred_data, horizon), horizon
    )
    boot_rows[[k]] <- data.table(
      replicate = bb, model_id = id, metric = names(m_boot),
      bootstrap_apparent = as.numeric(m_boot),
      original_test = as.numeric(m_test)
    )
  }
  if (bb %% 50L == 0L) {
    message("13K bootstrap completed: ", bb, "/", B)
  }
}
bootstrap <- rbindlist(boot_rows, fill = TRUE)
bootstrap[, optimism := bootstrap_apparent - original_test]
write_csv(bootstrap, "04_ADNI_5y_bootstrap_replicates.csv")

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

optimism <- bootstrap[, .(
  successful_bootstraps = sum(is.finite(optimism)),
  mean_optimism = safe_mean(optimism)
), by = .(model_id, metric)]
performance <- merge(
  apparent_metrics, optimism,
  by = c("model_id", "metric"), all.x = TRUE, sort = FALSE
)
performance[, optimism_corrected := apparent - mean_optimism]
write_csv(performance, "05_ADNI_5y_optimism_corrected_performance.csv")

pairs <- data.table(
  comparison = c("clinical_tier", "strong_benchmark_tier"),
  base_model = c("clinical_core", "strong_benchmark"),
  structure_model = c(
    "clinical_plus_structure", "strong_benchmark_plus_structure"
  )
)

delta_rows <- list()
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  for (metric_name in c("AUC", "Brier")) {
    app_base <- performance[
      model_id == p$base_model & metric == metric_name, apparent
    ]
    app_struct <- performance[
      model_id == p$structure_model & metric == metric_name, apparent
    ]
    cor_base <- performance[
      model_id == p$base_model & metric == metric_name, optimism_corrected
    ]
    cor_struct <- performance[
      model_id == p$structure_model & metric == metric_name, optimism_corrected
    ]

    x <- bootstrap[
      model_id == p$base_model & metric == metric_name,
      .(replicate, base_test = original_test)
    ]
    y <- bootstrap[
      model_id == p$structure_model & metric == metric_name,
      .(replicate, structure_test = original_test)
    ]
    dd <- merge(x, y, by = "replicate")
    dd[, delta := structure_test - base_test]
    q <- quantile(dd$delta[is.finite(dd$delta)], c(0.025, 0.975), na.rm = TRUE)

    delta_rows[[length(delta_rows) + 1L]] <- data.table(
      comparison = p$comparison,
      metric = metric_name,
      apparent_delta = app_struct - app_base,
      optimism_corrected_delta = cor_struct - cor_base,
      bootstrap_test_delta_ci_low = q[1L],
      bootstrap_test_delta_ci_high = q[2L],
      successful_bootstraps = sum(is.finite(dd$delta))
    )
  }
}
increment <- rbindlist(delta_rows)
write_csv(increment, "06_ADNI_5y_structural_increment.csv")

# -----------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()),file.path(out_dir,"99_sessionInfo.txt"))
