# Direct EMCI/LMCI event interaction and increment comparison
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

local({
  options(stringsAsFactors = FALSE, scipen = 999)
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })

  root <- study_path("project", "results")
  event_file <- file.path(root, "13F_R1_clinical_closure_confirmatory",
                          "09_baseline_event_frame.csv")
  oof_file <- file.path(root,
    "13H57_six_region_vs_global_clinical_value_release",
    "01_paired_oof_predictions.csv")
  pgs_file <- file.path(root, "10A_ADNI_gMAV_PGS_assembly",
                        "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
  out <- file.path(root, paste0("13H60_ADNI_MCI_entry_stage_event_increment_",
                                "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(x, name) fwrite(x, file.path(out, name), bom = TRUE)
  norm_rid <- function(x) {
    x <- sub("^.*?([0-9]+)$", "\\1", as.character(x))
    x <- sub("^0+", "", x)
    x[x == ""] <- "0"
    x
  }

  oof <- fread(oof_file)
  oof[, RID := norm_rid(RID)]
  if (nrow(oof) != 427L || uniqueN(oof$RID) != 427L ||
      sum(oof$event == 1L & oof$followup_years <= 5) != 138L)
    stop("The frozen 13H57 MCI prediction frame was not reproduced")
  event <- fread(event_file, select = c("RID", "index_mri_date",
    "age_index_z", "sex", "education_z", "baseline_MMSE",
    "baseline_CDRSB", "baseline_FAQ", "baseline_mPACC", "site"))
  event[, RID := norm_rid(RID)]
  if (anyDuplicated(event$RID)) stop("Duplicate RID in event frame")
  d <- merge(oof, event, by = "RID", all.x = TRUE, sort = FALSE)
  if (nrow(d) != 427L || anyNA(d$index_mri_date))
    stop("The event frame did not match all 427 prediction participants")

  pgs <- fread(pgs_file, sep = "\t",
               select = c("RID", "EXAMDATE", "dx_bl"), showProgress = FALSE)
  pgs[, `:=`(RID = norm_rid(RID), exam_date = as.IDate(EXAMDATE))]
  setorder(pgs, RID, exam_date)
  entry <- pgs[, {
    z <- toupper(trimws(as.character(dx_bl)))
    z <- z[!is.na(z) & nzchar(z)]
    data.table(entry_stage = if (length(z)) z[1L] else NA_character_,
               entry_exam_date = if (any(!is.na(exam_date)))
                 min(exam_date, na.rm = TRUE) else as.IDate(NA))
  }, by = RID]
  d <- merge(d, entry, by = "RID", all.x = TRUE, sort = FALSE)
  d[, entry_to_MRI_days := as.numeric(as.IDate(index_mri_date) -
                                       entry_exam_date)]
  d[, stage_group := fifelse(entry_stage %in% c("EMCI", "LMCI"),
                             entry_stage, "other_or_missing")]
  d[, five_year_status := fifelse(event == 1L & followup_years <= 5,
    "case", fifelse(followup_years > 5, "observed_non_case",
                    "censored_before_5y"))]
  flow <- d[, .(people = .N,
                cases_5y = sum(five_year_status == "case"),
                observed_non_cases_5y = sum(five_year_status == "observed_non_case"),
                censored_before_5y = sum(five_year_status == "censored_before_5y"),
                median_entry_to_MRI_days = if (all(is.na(entry_to_MRI_days)))
                  NA_real_ else median(entry_to_MRI_days, na.rm = TRUE),
                MRI_over_1y_after_entry = sum(entry_to_MRI_days > 365,
                                               na.rm = TRUE)),
            by = stage_group]
  save_csv(flow, "00_stage_flow.csv")

  d <- d[stage_group %in% c("EMCI", "LMCI")]
  d[, `:=`(entry_stage = factor(stage_group, levels = c("EMCI", "LMCI")),
            sex = factor(sex), site = factor(site),
            time5 = pmin(followup_years, 5),
            event5 = as.integer(event == 1L & followup_years <= 5),
            six_by_LMCI = structure_vulnerability_z *
              as.integer(stage_group == "LMCI"))]
  needed <- c("time5", "event5", "structure_vulnerability_z", "six_by_LMCI",
              "age_index_z", "sex", "education_z", "baseline_MMSE",
              "baseline_CDRSB", "baseline_FAQ", "baseline_mPACC", "site")
  if (nrow(d) == 0L || !setequal(as.character(unique(d$entry_stage)),
                                c("EMCI", "LMCI")) ||
      any(!complete.cases(d[, ..needed])))
    stop("EMCI/LMCI overlap or complete predictors are unavailable; inspect 00_stage_flow.csv")
  if (any(d[, .(cases = sum(event5), controls = sum(followup_years > 5)),
            by = entry_stage][, cases == 0L | controls == 0L]))
    stop("At least one stage lacks five-year cases or observed non-cases")

  cox_result <- function(tier, vars) {
    rhs <- paste(c("structure_vulnerability_z", "six_by_LMCI", vars,
                   "strata(entry_stage)", "cluster(site)"), collapse = " + ")
    fit <- coxph(as.formula(paste("Surv(time5, event5) ~", rhs)),
                 data = d, ties = "efron")
    b <- coef(fit)
    V <- vcov(fit)
    targets <- list(EMCI = "structure_vulnerability_z",
                    LMCI = c("structure_vulnerability_z", "six_by_LMCI"),
                    LMCI_vs_EMCI = "six_by_LMCI")
    rbindlist(lapply(names(targets), function(target) {
      L <- setNames(rep(0, length(b)), names(b))
      L[targets[[target]]] <- 1
      estimate <- sum(L * b)
      se <- sqrt(as.numeric(t(L) %*% V %*% L))
      data.table(tier = tier, contrast = target, log_HR = estimate,
        HR = exp(estimate), lower = exp(estimate - 1.96 * se),
        upper = exp(estimate + 1.96 * se),
        p = 2 * pnorm(-abs(estimate / se)),
        people = nrow(d), events_5y = sum(d$event5))
    }))
  }
  brief <- c("age_index_z", "sex", "education_z", "baseline_MMSE")
  extended <- c(brief, "baseline_CDRSB", "baseline_FAQ", "baseline_mPACC")
  cox <- rbind(cox_result("brief", brief),
               cox_result("extended", extended))
  save_csv(cox, "01_stage_MRI_event_interaction.csv")

  ipcw <- function(x) {
    fit <- survfit(Surv(followup_years, 1L - event) ~ 1, data = x)
    G <- function(t) {
      j <- findInterval(t, fit$time)
      pmax(0.02, c(1, fit$surv)[j + 1L])
    }
    y <- as.integer(x$event == 1L & x$followup_years <= 5)
    w <- numeric(nrow(x))
    w[y == 1L] <- 1 / G(pmax(0, x$followup_years[y == 1L] - 1e-8))
    w[x$followup_years > 5] <- 1 / G(5)
    list(y = y, w = w, G5 = G(5))
  }
  weighted_auc <- function(y, risk, w) {
    x <- data.table(y, risk, w)[w > 0,
      .(case = sum(w[y == 1L]), control = sum(w[y == 0L])), by = risk]
    setorder(x, risk)
    den <- sum(x$case) * sum(x$control)
    if (den == 0) return(NA_real_)
    sum(x$case * (cumsum(x$control) - 0.5 * x$control)) / den
  }
  increment <- function(x) {
    a <- ipcw(x)
    unlist(lapply(c("brief", "extended"), function(tier) {
      clinical <- x[[paste0("risk_", tier)]]
      with_six <- x[[paste0("risk_", tier, "_six")]]
      c(AUC = weighted_auc(a$y, with_six, a$w) -
                weighted_auc(a$y, clinical, a$w),
        Brier = sum(a$w * ((a$y - with_six)^2 -
                           (a$y - clinical)^2)) / nrow(x))
    }), use.names = FALSE)
  }
  score_names <- as.vector(t(outer(c("brief", "extended"),
                                   c("AUC", "Brier"), paste, sep = ".")))
  emci <- d[entry_stage == "EMCI"]
  lmci <- d[entry_stage == "LMCI"]
  performance <- rbindlist(lapply(list(EMCI = emci, LMCI = lmci), function(x) {
    a <- ipcw(x)
    rbindlist(lapply(c("brief", "extended"), function(tier) {
      rbindlist(lapply(c("clinical", "with_six"), function(model) {
        col <- paste0("risk_", tier, if (model == "with_six") "_six" else "")
        risk <- x[[col]]
        data.table(tier = tier, model = model,
          AUC = weighted_auc(a$y, risk, a$w),
          Brier = sum(a$w * (a$y - risk)^2) / nrow(x))
      }))
    }))
  }), idcol = "entry_stage")
  save_csv(performance, "02_stage_model_performance.csv")
  point_e <- setNames(increment(emci), score_names)
  point_l <- setNames(increment(lmci), score_names)
  set.seed(20261121L)
  boot <- replicate(400L, {
    e <- emci[sample.int(nrow(emci), replace = TRUE)]
    l <- lmci[sample.int(nrow(lmci), replace = TRUE)]
    c(increment(e), increment(l))
  })
  ci <- function(x) unname(quantile(x, c(.025, .975), na.rm = TRUE))
  estimates <- rbindlist(lapply(seq_along(score_names), function(i) {
    e <- ci(boot[i, ])
    l <- ci(boot[i + length(score_names), ])
    data.table(entry_stage = c("EMCI", "LMCI"),
               comparison = score_names[i],
               estimate = c(point_e[i], point_l[i]),
               lower = c(e[1L], l[1L]), upper = c(e[2L], l[2L]),
               valid_resamples = c(sum(is.finite(boot[i, ])),
                                   sum(is.finite(boot[i + length(score_names), ]))))
  }))
  save_csv(estimates, "03_stage_specific_MRI_increment.csv")
  heterogeneity <- rbindlist(lapply(seq_along(score_names), function(i) {
    delta <- boot[i + length(score_names), ] - boot[i, ]
    interval <- ci(delta)
    data.table(comparison = score_names[i],
      LMCI_minus_EMCI = unname(point_l[i] - point_e[i]),
      lower = interval[1L], upper = interval[2L],
      valid_resamples = sum(is.finite(delta)))
  }))
  save_csv(heterogeneity, "04_stage_difference_in_MRI_increment.csv")
  weights <- rbindlist(lapply(list(EMCI = emci, LMCI = lmci), function(x) {
    a <- ipcw(x)
    data.table(people = nrow(x), cases_5y = sum(a$y),
               observed_non_cases_5y = sum(x$followup_years > 5),
               censor_survival_G5 = a$G5, max_IPCW = max(a$w))
  }), idcol = "entry_stage")
  save_csv(weights, "05_stage_censoring.csv")
  writeLines(c(
    "Entry-stage EMCI/LMCI is the first nonmissing dx_bl, matching the 13H55 genetic-stage definition.",
    "The event origin and fixed out-of-fold risks come from the 13H57 427-person MCI frame.",
    "Cox models stop at five years, stratify baseline hazard by entry stage and use site-clustered variance.",
    "LMCI_vs_EMCI is the direct ratio of MRI hazard ratios per one ADNI SD; it is the interaction test.",
    "Stage-specific MRI increments reuse the fixed folds and out-of-fold risks.",
    "AUC/Brier stage-difference intervals resample participants within each stage and condition on fixed predictions.",
    "Negative Brier increments indicate improvement. Sample-flow outputs include entry-to-MRI gaps and overlap counts.",
    paste0("OUTPUT=", out)
  ), file.path(out, "06_scope.txt"))
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  cat("13H60 complete. Results: ", out, "\n", sep = "")
})
