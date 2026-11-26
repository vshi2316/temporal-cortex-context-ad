# HABS-HD within-cohort prediction and anatomical benchmark
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

local({
  options(stringsAsFactors = FALSE, scipen = 999)
  suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(survival)
  })

  root <- study_path("project", "results")
  release <- study_path("habs_release")
  mri_file <- file.path(root, "12B_ADNI_HABSHD_structural_MRI_feature_freeze",
                        "01_HABSHD_frozen_cortical_features.rds")
  date_file <- file.path(release,
    "HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx")
  out <- file.path(root, paste0("13H58_HABSHD_MCI_clinical_context_transport_",
                                 "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, name) fwrite(x, file.path(out, name), bom = TRUE)
  num <- function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    x[x <= -8888] <- NA_real_
    x
  }
  id <- function(x) {
    x <- trimws(as.character(x))
    sub("\\.0+$", "", x)
  }
  visit <- function(x) trimws(as.character(x))
  date <- function(year, month) {
    m <- match(tolower(substr(trimws(as.character(month)), 1L, 3L)),
               tolower(month.abb))
    numeric_month <- suppressWarnings(as.integer(as.character(month)))
    m[is.na(m)] <- numeric_month[is.na(m)]
    as.IDate(sprintf("%04d-%02d-15", as.integer(year), m))
  }
  z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  read_visit <- function(name, fields) {
    x <- fread(file.path(release, name), select = c("MED_ID", "VISIT_ID",
      "VISIT_MONTH", "VISIT_YEAR", fields), showProgress = FALSE)
    x[, `:=`(participant_id = id(MED_ID), visit_id = visit(VISIT_ID),
              visit_date = date(VISIT_YEAR, VISIT_MONTH))]
    x
  }

  mri <- as.data.table(readRDS(mri_file))
  mri[, `:=`(participant_id = id(participant_id),
              mri_timepoint = as.character(mri_timepoint))]
  sheets <- excel_sheets(date_file)
  sheets <- sheets[grepl("inclusion", sheets, ignore.case = TRUE) &
                     !grepl("excl|pending", sheets, ignore.case = TRUE)]
  dates <- rbindlist(lapply(sheets, function(s) {
    x <- as.data.table(read_excel(date_file, sheet = s, col_types = "text"))
    setnames(x, gsub("[^A-Za-z0-9]+", "_", names(x)))
    x[, .(participant_id = id(Med_ID),
          mri_timepoint = as.character(Time_point_MRI),
          mri_date = date(MRI_year, MRI_month))]
  }))
  dates <- unique(dates)
  mri <- merge(mri, dates, by = c("participant_id", "mri_timepoint"),
               all.x = TRUE, sort = FALSE)
  mri <- mri[!is.na(mri_date) & is.finite(num(ad_vulnerable_index)) &
               is.finite(num(global_cortical_thickness))]
  setorder(mri, participant_id, mri_date)
  mri <- mri[, .SD[1L], by = .(participant_id, mri_date)]

  dx <- read_visit("dx_12Aug2026.csv",
                   c("CDX_COG", "CDX_QUESTIONABLEVALIDITY"))
  dx[, `:=`(stage = num(CDX_COG), questionable = num(CDX_QUESTIONABLEVALIDITY))]
  dx <- dx[stage %in% 0:2 & (is.na(questionable) | questionable != 1) &
             !is.na(visit_date)]
  setorder(dx, participant_id, visit_date, visit_id)
  dx <- dx[, .SD[1L], by = .(participant_id, visit_date)]
  dx <- dx[, .(participant_id, visit_id, dx_date = visit_date, stage)]

  aligned <- merge(mri, dx, by = "participant_id", allow.cartesian = TRUE)
  aligned[, gap_days := abs(as.integer(mri_date - dx_date))]
  aligned <- aligned[gap_days <= 90L]
  setorder(aligned, participant_id, mri_date, gap_days, dx_date)
  aligned <- aligned[, .SD[1L], by = .(participant_id, mri_date)]
  aligned <- aligned[stage == 1L]
  setorder(aligned, participant_id, mri_date)
  aligned <- aligned[, .SD[1L], by = participant_id]
  aligned[, risk_start := pmax(mri_date, dx_date)]

  mmse <- read_visit("mmse_12Aug2026.csv", "MMSE_TOTAL")
  mmse <- mmse[, .(participant_id, visit_id, mmse = num(MMSE_TOTAL))]
  cdr <- read_visit("consensuscdr_12Aug2026.csv", "CDR_SUM")
  cdr <- cdr[, .(participant_id, visit_id, cdr_sb = num(CDR_SUM))]
  cog <- read_visit("zscores_12Aug2026.csv",
    c("DIGIT_SYMBOL_SUBSTITUTION_ZSCORE", "LM2_AB_ZSCORE",
      "SEVLT_DR_ZSCORE"))
  cog[, memory := rowMeans(cbind(num(LM2_AB_ZSCORE),
                                  num(SEVLT_DR_ZSCORE)), na.rm = TRUE)]
  cog[!is.finite(memory), memory := NA_real_]
  cog <- cog[, .(participant_id, visit_id,
                 digit_symbol = num(DIGIT_SYMBOL_SUBSTITUTION_ZSCORE), memory)]
  dem <- read_visit("interviewdemographics_12Aug2026.csv", "ID_GENDER")
  dem <- dem[, .(participant_id, visit_id, sex = num(ID_GENDER))]
  edu <- read_visit("residence_12Aug2026.csv", "ID_EDUCATION")
  edu <- edu[, .(participant_id, visit_id, education = num(ID_EDUCATION))]
  for (x in list(mmse, cdr, cog, dem, edu))
    aligned <- merge(aligned,
                     x[, .SD[1L], by = .(participant_id, visit_id)],
                     by = c("participant_id", "visit_id"),
                     all.x = TRUE, sort = FALSE)
  if (anyDuplicated(aligned$participant_id))
    stop("Baseline matching produced duplicate people")

  event <- dx[, {
    this_id <- .BY$participant_id
    origin <- aligned[participant_id == this_id, risk_start]
    if (!length(origin)) NULL else {
      at_risk <- !any(stage == 2L & dx_date <= origin)
      later <- .SD[dx_date > origin]
      event_date <- if (any(later$stage == 2L))
        min(later[stage == 2L, dx_date]) else as.IDate(NA)
      censor_date <- if (nrow(later)) max(later$dx_date) else as.IDate(NA)
      data.table(at_risk = at_risk, event_date = event_date,
                 censor_date = censor_date, later_visits = nrow(later))
    }
  }, by = participant_id]
  aligned <- merge(aligned, event, by = "participant_id", all.x = TRUE)
  aligned <- aligned[at_risk == TRUE & later_visits > 0L]
  aligned[, `:=`(
    followup_years = (as.numeric(fifelse(!is.na(event_date), event_date,
                                         censor_date)) - as.numeric(risk_start)) / 365.25,
    event = as.integer(!is.na(event_date)),
    six_z = -z(num(ad_vulnerable_index)),
    global_z = -z(num(global_cortical_thickness)),
    sex = factor(sex, levels = 0:1, labels = c("Male", "Female"))
  )]
  aligned <- aligned[is.finite(followup_years) & followup_years > 0]
  setorder(aligned, participant_id)

  brief_vars <- c("age", "sex", "education", "mmse", "six_z", "global_z")
  extended_vars <- c(brief_vars, "cdr_sb", "digit_symbol", "memory")
  brief <- droplevels(aligned[complete.cases(aligned[, ..brief_vars])])
  ext <- droplevels(aligned[complete.cases(aligned[, ..extended_vars])])
  flow <- data.table(
    step = c("MRI_people", "MCI_MRI_with_later_diagnosis", "brief_complete",
             "extended_complete"),
    people = c(uniqueN(mri$participant_id), nrow(aligned), nrow(brief), nrow(ext)),
    events = c(NA_integer_, sum(aligned$event), sum(brief$event), sum(ext$event)),
    five_year_cases = c(NA_integer_,
      sum(aligned$event == 1L & aligned$followup_years <= 5),
      sum(brief$event == 1L & brief$followup_years <= 5),
      sum(ext$event == 1L & ext$followup_years <= 5)),
    five_year_controls = c(NA_integer_,
      sum(aligned$followup_years > 5), sum(brief$followup_years > 5),
      sum(ext$followup_years > 5))
  )
  write_csv(flow, "00_sample_flow.csv")
  write_csv(ext[, .(participant_id, mri_date, dx_date, risk_start,
                    gap_days, followup_years, event, age, sex, education,
                    mmse, cdr_sb, digit_symbol, memory, six_z, global_z)],
            "01_extended_analysis_frame.csv")

  if (nrow(ext) >= 100L && sum(ext$event == 1L &
       ext$followup_years <= 5) >= 20L && sum(ext$followup_years > 5) >= 20L) {
    d <- copy(ext)
    set.seed(20261118L)
    d[, fold := sample(rep(1:5, length.out = .N))]
    base <- c("age", "sex", "education", "mmse")
    strong <- c(base, "cdr_sb", "digit_symbol", "memory")
    models <- list(brief = base, brief_six = c(base, "six_z"),
      brief_global = c(base, "global_z"),
      brief_both = c(base, "six_z", "global_z"),
      extended = strong, extended_six = c(strong, "six_z"),
      extended_global = c(strong, "global_z"),
      extended_both = c(strong, "six_z", "global_z"))
    predict_5y <- function(fit, x) {
      h <- basehaz(fit, centered = FALSE)
      h5 <- if (any(h$time <= 5)) max(h$hazard[h$time <= 5]) else 0
      lp <- predict(fit, x, type = "lp", reference = "zero")
      pmin(1 - 1e-8, pmax(1e-8, 1 - exp(-h5 * exp(lp))))
    }
    for (name in names(models)) d[, (paste0("risk_", name)) := NA_real_]
    for (k in 1:5) {
      train <- d[fold != k]
      test <- d[fold == k]
      for (name in names(models)) {
        form <- as.formula(paste("Surv(followup_years, event) ~",
          paste(models[[name]], collapse = " + ")))
        fit <- coxph(form, data = train, ties = "efron")
        if (any(!is.finite(coef(fit))))
          stop("Cox coefficient was not estimable: ", name, " fold ", k)
        d[fold == k, (paste0("risk_", name)) := predict_5y(fit, test)]
      }
    }
    write_csv(d[, c("participant_id", "fold", "followup_years", "event",
      paste0("risk_", names(models))), with = FALSE],
      "02_out_of_fold_predictions.csv")

    ipcw <- function(x) {
      fit <- survfit(Surv(followup_years, 1L - event) ~ 1, data = x)
      G <- function(t) {
        j <- findInterval(t, fit$time)
        pmax(.02, c(1, fit$surv)[j + 1L])
      }
      case <- x$event == 1L & x$followup_years <= 5
      control <- x$followup_years > 5
      w <- numeric(nrow(x))
      w[case] <- 1 / G(pmax(0, x$followup_years[case] - 1e-8))
      w[control] <- 1 / G(5)
      w
    }
    weighted_auc <- function(y, risk, w) {
      a <- data.table(y, risk, w)[w > 0,
        .(case = sum(w[y == 1L]), control = sum(w[y == 0L])), by = risk]
      setorder(a, risk)
      den <- sum(a$case) * sum(a$control)
      if (den == 0) return(NA_real_)
      sum(a$case * (cumsum(a$control) - .5 * a$control)) / den
    }
    score <- function(x) {
      w <- ipcw(x)
      y <- as.integer(x$event == 1L & x$followup_years <= 5)
      unlist(lapply(names(models), function(name) {
        r <- x[[paste0("risk_", name)]]
        c(AUC = weighted_auc(y, r, w),
          Brier = sum(w * (y - r)^2) / nrow(x))
      }), use.names = FALSE)
    }
    score_names <- as.vector(t(outer(names(models), c("AUC", "Brier"),
      paste, sep = ".")))
    point <- setNames(score(d), score_names)
    set.seed(20261119L)
    boot <- replicate(400L, {
      sample <- d[sample.int(nrow(d), replace = TRUE)]
      score(sample)
    })
    intervals <- data.table(model = rep(names(models), each = 2L),
      metric = rep(c("AUC", "Brier"), length(models)),
      estimate = as.numeric(point),
      lower = apply(boot, 1L, quantile, .025, na.rm = TRUE),
      upper = apply(boot, 1L, quantile, .975, na.rm = TRUE))
    write_csv(intervals, "03_oof_performance.csv")

    pairs <- list(brief_six = c("brief_six", "brief"),
      brief_global = c("brief_global", "brief"),
      brief_six_beyond_global = c("brief_both", "brief_global"),
      extended_six = c("extended_six", "extended"),
      extended_global = c("extended_global", "extended"),
      extended_six_beyond_global = c("extended_both", "extended_global"))
    contrast <- rbindlist(lapply(names(pairs), function(name) {
      a <- pairs[[name]][1L]
      b <- pairs[[name]][2L]
      rbindlist(lapply(c("AUC", "Brier"), function(metric) {
        ia <- match(paste(a, metric, sep = "."), score_names)
        ib <- match(paste(b, metric, sep = "."), score_names)
        v <- boot[ia, ] - boot[ib, ]
        data.table(comparison = name, metric = metric,
          estimate = unname(point[ia] - point[ib]),
          lower = unname(quantile(v, .025, na.rm = TRUE)),
          upper = unname(quantile(v, .975, na.rm = TRUE)))
      }))
    }))
    depth <- rbindlist(lapply(c("six", "global"), function(target) {
      rbindlist(lapply(c("AUC", "Brier"), function(metric) {
        i_ext <- match(paste0("extended_", target, ".", metric), score_names)
        i_ext_base <- match(paste0("extended.", metric), score_names)
        i_brief <- match(paste0("brief_", target, ".", metric), score_names)
        i_brief_base <- match(paste0("brief.", metric), score_names)
        v <- (boot[i_ext, ] - boot[i_ext_base, ]) -
             (boot[i_brief, ] - boot[i_brief_base, ])
        data.table(comparison = paste0(target, "_extended_minus_brief_increment"),
          metric = metric,
          estimate = unname((point[i_ext] - point[i_ext_base]) -
                              (point[i_brief] - point[i_brief_base])),
          lower = unname(quantile(v, .025, na.rm = TRUE)),
          upper = unname(quantile(v, .975, na.rm = TRUE)))
      }))
    }))
    contrast <- rbind(contrast, depth)
    write_csv(contrast, "04_paired_MRI_increment.csv")

    w <- ipcw(d)
    y <- as.integer(d$event == 1L & d$followup_years <= 5)
    dca <- rbindlist(lapply(names(models), function(name) {
      rbindlist(lapply(seq(.10, .30, .05), function(t) {
        chosen <- d[[paste0("risk_", name)]] >= t
        data.table(model = name, threshold = t,
          net_benefit = sum(w * chosen *
            (y - (1 - y) * t / (1 - t))) / nrow(d),
          treat_all = sum(w * (y - (1 - y) * t / (1 - t))) / nrow(d),
          treat_none = 0)
      }))
    }))
    write_csv(dca, "05_decision_curve.csv")
    status <- "FIVE_YEAR_PAIRED_COMPARISON_COMPLETED"
  } else {
    status <- "FIVE_YEAR_PAIRED_COMPARISON_NOT_ESTIMABLE_FROM_CURRENT_MCI_FOLLOWUP"
  }
  writeLines(c(
    paste0("STATUS=", status),
    "ENDPOINT=first later final-consensus dementia classification (CDX_COG=2)",
    "ORIGIN=first MCI-matched MRI within 90 days; month-level dates use the 15th",
    "MRI=same frozen six-region and global mean cortical thickness constructs",
    "BRIEF=age+sex+education+MMSE",
    "EXTENDED=brief+CDR-SB+digit-symbol z+episodic-memory z",
    "Clinical assessment tiers use the available HABS-HD measures specified in this script.",
    "Five-year performance intervals resample fixed out-of-fold predictions",
    paste0("OUTPUT=", out)
  ), file.path(out, "06_status_and_scope.txt"))
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  cat("13H58 complete. Results: ", out, "\n", sep = "")
})
