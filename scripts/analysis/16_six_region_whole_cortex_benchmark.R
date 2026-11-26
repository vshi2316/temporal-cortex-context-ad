# Six-region and whole-cortex MRI additions in matched ADNI folds
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
  fold_file <- file.path(root,
    "13H47_MCI_strong_clinical_MRI_increment_release",
    "01_MCI_out_of_fold_predictions.csv")
  out <- file.path(root, paste0("13H57_six_region_vs_global_clinical_value_",
                                 "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(x, name) fwrite(x, file.path(out, name), bom = TRUE)

  d <- fread(event_file, showProgress = FALSE)
  d[, `:=`(linked = as.logical(linked),
            prevalent_dementia = as.logical(prevalent_dementia),
            event = as.integer(as.logical(event)),
            followup_years = as.numeric(followup_years))]
  d <- d[linked == TRUE & prevalent_dementia == FALSE &
           baseline_dx == "MCI" & followup_years > 0]
  folds <- unique(fread(fold_file, select = c("RID", "fold")))
  d <- merge(d, folds, by = "RID", all.x = TRUE, sort = FALSE)

  brief <- c("age_index_z", "sex", "education_z", "baseline_MMSE")
  extended <- c(brief, "baseline_CDRSB", "baseline_FAQ", "baseline_mPACC")
  six <- "structure_vulnerability_z"
  global <- "global_vulnerability_z"
  needed <- c("RID", "fold", "followup_years", "event", extended, six, global)
  d <- droplevels(d[complete.cases(d[, ..needed])])
  d[, sex := factor(sex)]
  if (nrow(d) != 427L || uniqueN(d$RID) != 427L ||
      sum(d$event == 1L & d$followup_years <= 5) != 138L ||
      !setequal(d$fold, 1:5))
    stop("The paired 427-person / 138-case MCI frame was not reproduced")

  models <- list(
    brief = brief,
    brief_global = c(brief, global),
    brief_six = c(brief, six),
    brief_both = c(brief, global, six),
    extended = extended,
    extended_global = c(extended, global),
    extended_six = c(extended, six),
    extended_both = c(extended, global, six)
  )
  predict_5y <- function(fit, x) {
    h <- basehaz(fit, centered = FALSE)
    h5 <- if (any(h$time <= 5)) max(h$hazard[h$time <= 5]) else 0
    lp <- predict(fit, newdata = x, type = "lp", reference = "zero")
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
  risk_cols <- paste0("risk_", names(models))
  if (any(!is.finite(as.matrix(d[, ..risk_cols]))))
    stop("Out-of-fold predictions are incomplete")
  save_csv(d[, c("RID", "fold", "followup_years", "event", six, global,
                 risk_cols), with = FALSE], "01_paired_oof_predictions.csv")

  ipcw <- function(x) {
    fit <- survfit(Surv(followup_years, 1L - event) ~ 1, data = x)
    G <- function(t) {
      j <- findInterval(t, fit$time)
      pmax(0.02, c(1, fit$surv)[j + 1L])
    }
    case <- x$event == 1L & x$followup_years <= 5
    control <- x$followup_years > 5
    w <- numeric(nrow(x))
    w[case] <- 1 / G(pmax(0, x$followup_years[case] - 1e-8))
    w[control] <- 1 / G(5)
    w
  }
  weighted_auc <- function(y, risk, weight) {
    a <- data.table(y, risk, weight)[weight > 0,
      .(case_weight = sum(weight[y == 1L]),
        control_weight = sum(weight[y == 0L])), by = risk]
    setorder(a, risk)
    den <- sum(a$case_weight) * sum(a$control_weight)
    if (den == 0) return(NA_real_)
    sum(a$case_weight *
          (cumsum(a$control_weight) - 0.5 * a$control_weight)) / den
  }
  scores <- function(x, w) {
    y <- as.integer(x$event == 1L & x$followup_years <= 5)
    unlist(lapply(names(models), function(name) {
      risk <- x[[paste0("risk_", name)]]
      c(AUC = weighted_auc(y, risk, w),
        Brier = sum(w * (y - risk)^2) / nrow(x))
    }), use.names = FALSE)
  }
  score_names <- as.vector(t(outer(names(models), c("AUC", "Brier"), paste,
                                   sep = ".")))
  observed <- setNames(scores(d, ipcw(d)), score_names)
  performance <- data.table(model = rep(names(models), each = 2L),
                            metric = rep(c("AUC", "Brier"), length(models)),
                            value = as.numeric(observed))
  save_csv(performance, "02_oof_performance.csv")

  pairs <- list(
    brief_six_increment = c("brief_six", "brief"),
    brief_global_increment = c("brief_global", "brief"),
    brief_six_beyond_global = c("brief_both", "brief_global"),
    brief_global_beyond_six = c("brief_both", "brief_six"),
    extended_six_increment = c("extended_six", "extended"),
    extended_global_increment = c("extended_global", "extended"),
    extended_six_beyond_global = c("extended_both", "extended_global"),
    extended_global_beyond_six = c("extended_both", "extended_six")
  )
  contrasts <- function(x, w) {
    s <- setNames(scores(x, w), score_names)
    v <- unlist(lapply(names(pairs), function(name) {
      a <- pairs[[name]][1L]
      b <- pairs[[name]][2L]
      c(AUC = s[paste0(a, ".AUC")] - s[paste0(b, ".AUC")],
        Brier = s[paste0(a, ".Brier")] - s[paste0(b, ".Brier")])
    }), use.names = FALSE)
    setNames(v, as.vector(t(outer(names(pairs), c("AUC", "Brier"), paste,
                                 sep = "."))))
  }
  point <- contrasts(d, ipcw(d))
  set.seed(20261118L)
  boot <- replicate(400L, {
    sample <- d[sample.int(nrow(d), replace = TRUE)]
    contrasts(sample, ipcw(sample))
  })
  summary <- data.table(contrast = names(point), estimate = as.numeric(point),
                        lower = apply(boot, 1L, quantile, .025, na.rm = TRUE),
                        upper = apply(boot, 1L, quantile, .975, na.rm = TRUE))
  save_csv(summary, "03_paired_MRI_contrasts.csv")

  y <- as.integer(d$event == 1L & d$followup_years <= 5)
  w <- ipcw(d)
  dca <- rbindlist(lapply(names(models), function(name) {
    rbindlist(lapply(seq(.10, .30, .05), function(t) {
      chosen <- d[[paste0("risk_", name)]] >= t
      data.table(model = name, threshold = t,
        net_benefit = sum(w * chosen * (y - (1 - y) * t / (1 - t))) / nrow(d))
    }))
  }))
  save_csv(dca, "04_decision_curve.csv")
  save_csv(data.table(people = nrow(d), cases_5y = sum(y),
                      controls_5y = sum(d$followup_years > 5),
                      early_censored = sum(w == 0)), "00_sample_flow.csv")
  writeLines(c(
    "Eight models use the same 427 MCI participants, the same five folds and the same five-year outcome.",
    "Six-region temporal vulnerability is compared with global cortical vulnerability in the same baseline FreeSurfer frame.",
    "Paired intervals resample fixed out-of-fold risks and are conditional on the fitted folds.",
    "Paired six-region/global comparisons use person-aligned AUC and Brier estimates with the corresponding net-benefit curves."
  ), file.path(out, "05_readme.txt"))
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  cat("13H57 complete. Results: ", out, "\n", sep = "")
})
