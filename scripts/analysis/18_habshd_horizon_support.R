# Follow-up support, censoring and paired horizon-specific evaluation
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

local({
  options(stringsAsFactors = FALSE, scipen = 999)
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })

  root <- study_path("project", "results")
  input <- file.path(root,
    "13H58_HABSHD_MCI_clinical_context_transport_release")
  out <- file.path(root, paste0("13H59_HABSHD_MCI_horizon_assessment_",
                                "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(x, name) fwrite(x, file.path(out, name), bom = TRUE)

  d <- fread(file.path(input, "01_extended_analysis_frame.csv"))
  folds <- fread(file.path(input, "02_out_of_fold_predictions.csv"),
                 select = c("participant_id", "fold"))
  d <- merge(d, unique(folds), by = "participant_id", all.x = TRUE,
             sort = FALSE)
  d[, sex := factor(sex, levels = c("Male", "Female"))]
  if (nrow(d) != 489L || uniqueN(d$participant_id) != 489L ||
      anyNA(d$fold) || !setequal(d$fold, 1:5))
    stop("The frozen HABS-HD 489-person frame/folds were not reproduced")

  brief <- c("age", "sex", "education", "mmse")
  extended <- c(brief, "cdr_sb", "digit_symbol", "memory")
  models <- list(brief = brief, brief_six = c(brief, "six_z"),
    brief_global = c(brief, "global_z"),
    brief_both = c(brief, "six_z", "global_z"),
    extended = extended, extended_six = c(extended, "six_z"),
    extended_global = c(extended, "global_z"),
    extended_both = c(extended, "six_z", "global_z"))
  horizons <- c(2, 3, 5)
  predict_at <- function(fit, newdata, horizon) {
    h <- basehaz(fit, centered = FALSE)
    ht <- if (any(h$time <= horizon))
      max(h$hazard[h$time <= horizon]) else 0
    lp <- predict(fit, newdata = newdata, type = "lp", reference = "zero")
    pmin(1 - 1e-8, pmax(1e-8, 1 - exp(-ht * exp(lp))))
  }
  for (h in horizons)
    for (name in names(models))
      d[, (paste0("risk_", name, "_", h, "y")) := NA_real_]
  for (k in 1:5) {
    train <- d[fold != k]
    test <- d[fold == k]
    for (name in names(models)) {
      form <- as.formula(paste("Surv(followup_years, event) ~",
                               paste(models[[name]], collapse = " + ")))
      fit <- coxph(form, data = train, ties = "efron")
      if (any(!is.finite(coef(fit))))
        stop("Cox coefficient not estimable: ", name, " fold ", k)
      for (h in horizons)
        d[fold == k, (paste0("risk_", name, "_", h, "y")) :=
            predict_at(fit, test, h)]
    }
  }
  risk_cols <- unlist(lapply(horizons, function(h)
    paste0("risk_", names(models), "_", h, "y")), use.names = FALSE)
  if (any(!is.finite(as.matrix(d[, ..risk_cols]))))
    stop("Out-of-fold predictions are incomplete")
  save_csv(d[, c("participant_id", "fold", "followup_years", "event",
                 risk_cols), with = FALSE], "01_oof_predictions.csv")

  ipcw <- function(x, h) {
    censor <- survfit(Surv(followup_years, 1L - event) ~ 1, data = x)
    G <- function(t) {
      j <- findInterval(t, censor$time)
      c(1, censor$surv)[j + 1L]
    }
    case <- x$event == 1L & x$followup_years <= h
    control <- x$followup_years > h
    w <- numeric(nrow(x))
    w[case] <- 1 / G(pmax(0, x$followup_years[case] - 1e-8))
    w[control] <- 1 / G(h)
    list(y = as.integer(case), w = w, case = case, control = control,
         G_h = G(h))
  }
  ess <- function(w) if (!length(w) || sum(w^2) == 0) NA_real_ else
    sum(w)^2 / sum(w^2)
  weighted_auc <- function(y, risk, w) {
    x <- data.table(y, risk, w)[w > 0,
      .(case = sum(w[y == 1L]), control = sum(w[y == 0L])), by = risk]
    setorder(x, risk)
    den <- sum(x$case) * sum(x$control)
    if (den == 0) return(NA_real_)
    sum(x$case * (cumsum(x$control) - .5 * x$control)) / den
  }
  score <- function(x, h) {
    ip <- ipcw(x, h)
    unlist(lapply(names(models), function(name) {
      r <- x[[paste0("risk_", name, "_", h, "y")]]
      c(AUC = weighted_auc(ip$y, r, ip$w),
        Brier = sum(ip$w * (ip$y - r)^2) / nrow(x))
    }), use.names = FALSE)
  }
  score_names <- as.vector(t(outer(names(models), c("AUC", "Brier"),
                                   paste, sep = ".")))
  pairs <- list(brief_six = c("brief_six", "brief"),
    brief_six_beyond_global = c("brief_both", "brief_global"),
    extended_six = c("extended_six", "extended"),
    extended_six_beyond_global = c("extended_both", "extended_global"))

  flow <- list(); perf <- list(); increments <- list(); calibration <- list()
  set.seed(20261120L)
  for (h in horizons) {
    ip <- ipcw(d, h)
    observed <- 1 - summary(survfit(Surv(followup_years, event) ~ 1,
                                     data = d), times = h, extend = TRUE)$surv
    flow[[as.character(h)]] <- data.table(
      horizon_years = h, people = nrow(d), cases = sum(ip$case),
      observed_non_cases = sum(ip$control),
      censored_before_horizon = nrow(d) - sum(ip$case) - sum(ip$control),
      censor_survival_G = ip$G_h,
      case_effective_n = ess(ip$w[ip$case]),
      non_case_effective_n = ess(ip$w[ip$control]),
      max_ipcw = max(ip$w), km_event_risk = observed)
    point <- setNames(score(d, h), score_names)
    boot <- replicate(400L, score(d[sample.int(nrow(d), replace = TRUE)], h))
    perf[[as.character(h)]] <- data.table(
      horizon_years = h, model = rep(names(models), each = 2L),
      metric = rep(c("AUC", "Brier"), length(models)),
      estimate = as.numeric(point),
      lower = apply(boot, 1L, quantile, .025, na.rm = TRUE),
      upper = apply(boot, 1L, quantile, .975, na.rm = TRUE))
    increments[[as.character(h)]] <- rbindlist(lapply(names(pairs), function(p) {
      rbindlist(lapply(c("AUC", "Brier"), function(metric) {
        a <- match(paste(pairs[[p]][1L], metric, sep = "."), score_names)
        b <- match(paste(pairs[[p]][2L], metric, sep = "."), score_names)
        delta <- boot[a, ] - boot[b, ]
        data.table(horizon_years = h, comparison = p, metric = metric,
          estimate = unname(point[a] - point[b]),
          lower = unname(quantile(delta, .025, na.rm = TRUE)),
          upper = unname(quantile(delta, .975, na.rm = TRUE)))
      }))
    }))
    depth <- rbindlist(lapply(c("AUC", "Brier"), function(metric) {
      i <- match(paste(c("extended_six", "extended", "brief_six", "brief"),
                       metric, sep = "."), score_names)
      delta <- (boot[i[1L], ] - boot[i[2L], ]) -
               (boot[i[3L], ] - boot[i[4L], ])
      data.table(horizon_years = h,
        comparison = "six_extended_minus_brief_increment", metric = metric,
        estimate = unname((point[i[1L]] - point[i[2L]]) -
                            (point[i[3L]] - point[i[4L]])),
        lower = unname(quantile(delta, .025, na.rm = TRUE)),
        upper = unname(quantile(delta, .975, na.rm = TRUE)))
    }))
    increments[[as.character(h)]] <- rbind(increments[[as.character(h)]], depth)
    calibration[[as.character(h)]] <- rbindlist(lapply(names(models), function(name) {
      pred <- mean(d[[paste0("risk_", name, "_", h, "y")]])
      data.table(horizon_years = h, model = name,
                 mean_predicted_risk = pred, km_observed_risk = observed,
                 observed_to_expected = observed / pred)
    }))
  }
  save_csv(rbindlist(flow), "00_horizon_information.csv")
  save_csv(rbindlist(perf), "02_oof_performance.csv")
  save_csv(rbindlist(increments), "03_paired_MRI_increment.csv")
  save_csv(rbindlist(calibration), "04_average_calibration.csv")

  old <- fread(file.path(input, "03_oof_performance.csv"))
  now <- rbindlist(perf)[horizon_years == 5]
  check <- merge(now[, .(model, metric, new_estimate = estimate)],
                 old[, .(model, metric, old_estimate = estimate)],
                 by = c("model", "metric"))
  if (nrow(check) != nrow(now) ||
      any(abs(check$new_estimate - check$old_estimate) > 1e-6))
    stop("Five-year results do not reproduce 13H58; inspect inputs")
  writeLines(c(
    "HABS-HD uses the fixed 489-participant frame and five folds.",
    "Two- and three-year windows assess follow-up suitability; five years reproduces 13H58.",
    "Models are fitted and evaluated within HABS-HD.",
    "IPCW uses marginal censoring KM, consistent with 13H58; inspect G(t), effective N and maximum weight.",
    "Bootstrap resamples fixed out-of-fold predictions conditional on the fitted models.",
    "Brier: lower is better. AUC: higher is better. Average calibration is descriptive.",
    paste0("OUTPUT=", out)
  ), file.path(out, "05_scope.txt"))
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  cat("13H59 complete. Results: ", out, "\n", sep = "")
})
