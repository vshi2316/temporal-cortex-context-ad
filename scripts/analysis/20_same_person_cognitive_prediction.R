# Ordered same-person held-out continuous cognitive prediction
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

local({
  options(stringsAsFactors = FALSE, scipen = 999)
  suppressPackageStartupMessages(library(data.table))

  root <- study_path("project", "results")
  h55 <- file.path(root,
    "13H55_transition_ordered_structural_chain_release")
  triple_file <- file.path(h55, "04_nonoverlapping_interval_triples.csv")
  frame_file <- file.path(h55, "01_frozen_frame.csv")
  clinical_file <- file.path(root, "10A_ADNI_gMAV_PGS_assembly",
    "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
  inputs <- c(triple_file, frame_file, clinical_file)
  if (any(!file.exists(inputs)))
    stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = "; "))

  out <- file.path(root, paste0("13H61_ADNI_same_person_cognitive_increment_",
    "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(x, name)
    fwrite(as.data.table(x), file.path(out, name), bom = TRUE)

  triples <- fread(triple_file, showProgress = FALSE)
  triples[, RID := as.character(RID)]
  for (v in c("t0", "t1", "c0", "c1", "c2"))
    triples[, (v) := as.IDate(get(v))]
  n_forward <- uniqueN(triples[forward_order == TRUE, RID])
  setorder(triples, RID, t0, t1)
  d <- triples[forward_order == TRUE & stage_t1 %in% c("CN", "MCI"),
               .SD[1L], by = RID]

  frame <- fread(frame_file,
    select = c("RID", "mri_date", "global_z"), showProgress = FALSE)
  frame[, `:=`(RID = as.character(RID), mri_date = as.IDate(mri_date))]
  if (anyDuplicated(frame[, .(RID, mri_date)]))
    stop("13H55 MRI frame has duplicate RID and MRI date")
  global_at_t1 <- frame[, .(RID, t1 = mri_date, global1 = global_z)]
  d <- merge(d, global_at_t1, by = c("RID", "t1"), all.x = TRUE,
             sort = FALSE)

  cog <- fread(clinical_file, sep = "\t",
    select = c("RID", "EXAMDATE", "mPACCdigit"), showProgress = FALSE)
  cog[, `:=`(RID = as.character(RID),
             cog_date = as.IDate(EXAMDATE),
             mPACC = suppressWarnings(as.numeric(mPACCdigit)))]
  cog <- cog[is.finite(mPACC) & !is.na(cog_date)]
  setorder(cog, RID, cog_date)
  cog <- cog[, .SD[1L], by = .(RID, cog_date)]
  d <- merge(d, cog[, .(RID, c0 = cog_date, m0 = mPACC)],
             by = c("RID", "c0"), all.x = TRUE, sort = FALSE)
  d <- merge(d, cog[, .(RID, c1 = cog_date, m1 = mPACC)],
             by = c("RID", "c1"), all.x = TRUE, sort = FALSE)
  d <- merge(d, cog[, .(RID, c2 = cog_date, m2 = mPACC)],
             by = c("RID", "c2"), all.x = TRUE, sort = FALSE)
  if (anyNA(d[, .(m0, m1, m2, global1)]))
    stop("Frozen 13H55 triples could not be matched to original mPACC or global thickness")
  d[, `:=`(prior_mPACC_rate = (m1 - m0) / dtc01,
           future_mPACC_rate = (m2 - m1) / dtc12,
           sex = factor(sex), stage_t1 = factor(stage_t1))]

  numeric_vars <- c("future_mPACC_rate", "prior_mPACC_rate", "m1",
    "loss01", "global_loss01", "thickness1", "global1", "age1",
    "education", "gAD", "gMAV", "dt01", "dtc01")
  n_first <- nrow(d)
  d <- d[complete.cases(d[, c(numeric_vars, "sex", "stage_t1"),
                          with = FALSE])]
  for (v in numeric_vars) d <- d[is.finite(get(v))]
  d <- droplevels(d)
  if (nrow(d) < 100L || any(d[, .N, by = .(stage_t1, sex)]$N < 5L))
    stop("Insufficient complete first triples for stratified five-fold evaluation")
  if (anyDuplicated(d$RID)) stop("More than one interval per participant")
  save_csv(data.table(
    step = c("13H55_forward_ordered_people", "first_forward_triple_people",
             "complete_model_people"),
    people = c(n_forward, n_first, nrow(d)),
    intervals = c(nrow(triples[forward_order == TRUE]), n_first, nrow(d))),
    "00_sample_flow.csv")

  # All predictors were available by the c1 cognitive visit. The future visit
  # interval dtc12 and the later MRI are never supplied to either model.
  base_formula <- future_mPACC_rate ~ m1 + prior_mPACC_rate + thickness1 +
    global1 + global_loss01 + age1 + sex + education + stage_t1 +
    gAD + gMAV + dt01 + dtc01
  plus_formula <- update(base_formula, . ~ . + loss01)
  set.seed(20261120)
  d[, fold := sample(rep(1:5, length.out = .N)), by = .(stage_t1, sex)]
  d[, `:=`(pred_base = NA_real_, pred_plus = NA_real_)]
  for (k in 1:5) {
    train <- d[fold != k]
    test <- d[fold == k]
    fit_base <- lm(base_formula, data = train)
    fit_plus <- lm(plus_formula, data = train)
    d[fold == k, `:=`(pred_base = as.numeric(predict(fit_base, test)),
                      pred_plus = as.numeric(predict(fit_plus, test)))]
  }
  if (any(!is.finite(d$pred_base)) || any(!is.finite(d$pred_plus)))
    stop("Out-of-fold predictions are incomplete")
  save_csv(d[, .(RID, fold, stage_t1, sex, t0, t1, c1, c2,
                 dt01, dtc01, dtc12, m1, prior_mPACC_rate,
                 future_mPACC_rate, loss01, global_loss01,
                 pred_base, pred_plus)], "01_paired_oof_predictions.csv")

  metrics <- function(x) {
    y <- x$future_mPACC_rate
    e0 <- y - x$pred_base
    e1 <- y - x$pred_plus
    sst <- sum((y - mean(y))^2)
    c(rmse_base = sqrt(mean(e0^2)), rmse_plus = sqrt(mean(e1^2)),
      delta_rmse = sqrt(mean(e0^2)) - sqrt(mean(e1^2)),
      mae_base = mean(abs(e0)), mae_plus = mean(abs(e1)),
      delta_mae = mean(abs(e0)) - mean(abs(e1)),
      r2_base = 1 - sum(e0^2) / sst,
      r2_plus = 1 - sum(e1^2) / sst,
      delta_r2 = (sum(e0^2) - sum(e1^2)) / sst)
  }
  point <- metrics(d)
  save_csv(data.table(people = nrow(d),
    outcome = "annual mPACCdigit change after c1, original units/year",
    t1_to_c1_days_median = median(as.numeric(d$c1 - d$t1)),
    later_visit_years_median = median(d$dtc12),
    metric = names(point), estimate = unname(point)),
    "02_oof_performance.csv")

# Paired participant bootstrap conditional on fixed held-out predictions.
  set.seed(20261121)
  boot <- replicate(1000L,
    metrics(d[sample.int(nrow(d), nrow(d), replace = TRUE)]))
  increments <- c("delta_rmse", "delta_mae", "delta_r2")
  ci <- t(apply(boot[increments, , drop = FALSE], 1L, quantile,
                probs = c(0.025, 0.975), na.rm = TRUE))
  save_csv(data.table(metric = increments, estimate = point[increments],
    lower_95 = ci[, 1L], upper_95 = ci[, 2L],
    positive_means = "earlier six-region thinning improves held-out prediction"),
    "03_paired_increment_bootstrap.csv")

  full_fit <- lm(plus_formula, data = d)
  beta <- unname(coef(full_fit)["loss01"])
  se <- unname(sqrt(diag(vcov(full_fit)))["loss01"])
  scale_loss <- sd(d$loss01)
  save_csv(data.table(
    term = "earlier annual six-region thickness loss per 1 sample SD",
    outcome = "subsequent mPACCdigit units/year",
    estimate = beta * scale_loss,
    lower_95 = (beta - qt(0.975, df.residual(full_fit)) * se) * scale_loss,
    upper_95 = (beta + qt(0.975, df.residual(full_fit)) * se) * scale_loss,
    people = nrow(d)), "04_first_triple_association.csv")

  writeLines(c(
    "Estimand: Among ADNI participants with a first time-ordered, nonoverlapping MRI and cognition triple, does earlier six-region thinning improve prediction of the next annual mPACCdigit change?",
    "Index: c1, the cognitive visit on or after the second MRI t1. Outcome: (mPACC at c2 - mPACC at c1) / years from c1 to c2.",
    "Comparator: current mPACC, earlier mPACC change, current six-region and global thickness, earlier global thinning, age, sex, education, current diagnosis, gAD, gMAV, and observed pre-index intervals.",
    "Added feature: earlier six-region thickness loss, (thickness at t0 - thickness at t1) / years from t0 to t1.",
    "Validation: first eligible triple per person; deterministic five-fold participant-level evaluation stratified by current diagnosis and sex; same participants and folds in both models.",
    "Positive delta_R2, delta_RMSE and delta_MAE indicate improvement. Intervals bootstrap participants on the fixed out-of-fold predictions."

  ), file.path(out, "05_design_and_interpretation.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(out, "99_sessionInfo.txt"))
  cat("13H61 completed. Results: ", out, "\n", sep = "")
})
