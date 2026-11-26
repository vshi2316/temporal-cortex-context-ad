# Clinical component decomposition, reciprocal hippocampal benchmark and CSF linkage
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

local({
  options(stringsAsFactors = FALSE, scipen = 999)
  packages <- c("data.table", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })

  root <- study_path("project", "results")
  adni <- study_path("cohorts", "ADNI/data")
  inputs <- c(
    event = file.path(root, "13F_R1_clinical_closure_confirmatory", "09_baseline_event_frame.csv"),
    frozen_oof = file.path(root, "13H57_six_region_vs_global_clinical_value_release", "01_paired_oof_predictions.csv"),
    genetic_frame = file.path(root, "13H55_transition_ordered_structural_chain_release", "01_frozen_frame.csv"),
    clinical = file.path(adni, "Clinical", "ADNIMERGE_03Jan2026.csv"),
    amyloid = file.path(adni, "PET", "PET_Amyloid_Berkeley.csv"),
    tau_pet = file.path(adni, "PET", "PET_Tau_Berkeley.csv"),
    elecsys = file.path(adni, "Biomarkers", "CSF_Ab_Tau_Elecsys.csv"),
    alzbio3 = file.path(adni, "Biomarkers", "CSF_Ab_Tau_Alzbio3.csv")
  )
  required <- c("event", "frozen_oof", "genetic_frame", "clinical")
  absent <- required[!file.exists(inputs[required])]
  if (length(absent)) stop("Missing input: ", paste(inputs[absent], collapse = "; "))
  out <- file.path(root, paste0("13H62_ADNI_components_hippocampus_AT_audit_",
                                "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(x, name) fwrite(as.data.table(x), file.path(out, name), bom = TRUE)
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  norm_rid <- function(x) {
    z <- trimws(as.character(x))
    z <- sub("[.]0+$", "", z)
    z <- sub("^.*?([0-9]+)$", "\\1", z)
    z <- sub("^0+", "", z)
    z[z == ""] <- "0"
    z[is.na(x) | !nzchar(trimws(as.character(x)))] <- NA_character_
    z
  }
  norm_image <- function(x) {
    z <- toupper(trimws(as.character(x)))
    z <- sub("^I", "", z)
    z <- sub("[.]0+$", "", z)
    z[is.na(x) | z %in% c("", "NA", "NAN", "0")] <- NA_character_
    z
  }
  assert_columns <- function(x, cols, label) {
    absent <- setdiff(cols, names(x))
    if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
  }
  bool_event <- function(x) {
    z <- toupper(trimws(as.character(x)))
    answer <- rep(NA_integer_, length(z))
    answer[z %in% c("TRUE", "1")] <- 1L
    answer[z %in% c("FALSE", "0")] <- 0L
    answer
  }
  interval <- function(x) {
    z <- x[is.finite(x)]
    if (length(z) < 2L) return(c(NA_real_, NA_real_))
    unname(quantile(z, c(.025, .975), names = FALSE))
  }

  # Outcome definitions and folds are inherited from 13H57.
  horizon <- 5
  bootstrap_n <- 1000L
  bootstrap_seed <- 20261120L
  manifest <- data.table(input = names(inputs), path = unname(inputs),
                        exists = file.exists(inputs))
  manifest[, bytes := file.info(path)$size]
  manifest[, modified := as.character(file.info(path)$mtime)]
  manifest[, md5 := NA_character_]
  manifest[exists == TRUE, md5 := unname(tools::md5sum(path))]
  save_csv(manifest, "00_input_manifest.csv")

  frozen <- fread(inputs["frozen_oof"], showProgress = FALSE)
  assert_columns(frozen, c("RID", "fold", "followup_years", "event"), "13H57 predictions")
  frozen[, `:=`(RID = norm_rid(RID), event = bool_event(event), fold = as.integer(fold))]
  if (nrow(frozen) != 427L || uniqueN(frozen$RID) != 427L ||
      anyNA(frozen$event) || sum(frozen$event == 1L & frozen$followup_years <= horizon) != 138L ||
      !setequal(frozen$fold, 1:5)) stop("The frozen 427-person / 138-case frame was not reproduced")
  event <- fread(inputs["event"], showProgress = FALSE)
  event[, RID := norm_rid(RID)]
  assert_columns(event, c("RID", "index_mri_date", "IMAGEUID", "risk_start", "baseline_dx_date",
    "age_index", "sex", "education", "baseline_MMSE", "baseline_mPACC", "baseline_CDRSB",
    "baseline_FAQ", "ad_vulnerable_index", "global_cortical_thickness", "followup_years", "event"),
    "Clinical event frame")
  if (anyDuplicated(event$RID)) stop("Duplicate RID in the clinical event frame")
  event[, `:=`(source_followup = num(followup_years), source_event = bool_event(event))]
  event[, c("event", "followup_years") := NULL]
  d <- merge(frozen[, .(RID, fold, followup_years, event)], event, by = "RID", all.x = TRUE, sort = FALSE)
  if (nrow(d) != 427L || anyNA(d$index_mri_date) ||
      any(d$event != d$source_event) || any(abs(d$followup_years - d$source_followup) > 1e-9))
    stop("Event outcomes or prediction participants changed since 13H57")
  d[, `:=`(index_mri_date = as.IDate(index_mri_date), risk_start = as.IDate(risk_start),
            baseline_dx_date = as.IDate(baseline_dx_date),
            age_raw = num(age_index), education_raw = num(education),
            six_raw = -num(ad_vulnerable_index), global_raw = -num(global_cortical_thickness),
            sex = factor(sex))]
  for (v in c("baseline_MMSE", "baseline_mPACC", "baseline_CDRSB", "baseline_FAQ"))
    d[, (v) := num(get(v))]
  numeric_core <- c("followup_years", "age_raw", "education_raw", "baseline_MMSE",
                    "baseline_mPACC", "baseline_CDRSB", "baseline_FAQ", "six_raw", "global_raw")
  if (any(!is.finite(as.matrix(d[, ..numeric_core]))) || anyNA(d$sex) ||
      anyNA(d$risk_start) || anyNA(d$baseline_dx_date) ||
      any(d$baseline_dx_date > d$risk_start) || any(d$index_mri_date > d$risk_start))
    stop("Predictors, dates, or frozen risk origin are inconsistent")
  setorder(d, RID)
  save_csv(d[, .(RID, fold, index_mri_date, baseline_dx_date, risk_start,
    clinical_to_origin_days = as.integer(risk_start - baseline_dx_date),
    MRI_to_origin_days = as.integer(risk_start - index_mri_date))], "01_predictor_timing.csv")

  brief <- c("age_raw", "sex", "education_raw", "baseline_MMSE")
  tiers <- list(brief = brief,
                cognitive = c(brief, "baseline_mPACC"),
                functional = c(brief, "baseline_CDRSB", "baseline_FAQ"),
                full = c(brief, "baseline_mPACC", "baseline_CDRSB", "baseline_FAQ"))
  # The same training-fold transformation is applied to its held-out fold.
  # Linear rescaling does not alter Cox predictions relative to 13H57.
  fit_oof <- function(dat, models, label) {
    x <- copy(dat)
    cols <- unique(unlist(models, use.names = FALSE))
    assert_columns(x, c("RID", "fold", "followup_years", "event", cols), label)
    if (anyNA(x[, ..cols])) stop(label, ": incomplete predictors")
    for (name in names(models)) x[, (paste0("risk_", name)) := NA_real_]
    log <- list(); scaling <- list(); coefficients <- list(); i <- 0L
    for (k in sort(unique(x$fold))) {
      train <- copy(x[fold != k]); test <- copy(x[fold == k])
      if (!nrow(test) || sum(train$event) == 0L) stop(label, ": invalid fold ", k)
      for (v in cols) {
        if (is.numeric(train[[v]])) {
          m <- mean(train[[v]]); s <- sd(train[[v]])
          if (!is.finite(s) || s == 0) stop(label, ": constant training predictor ", v, " in fold ", k)
          train[, (v) := (get(v) - m) / s]
          test[, (v) := (get(v) - m) / s]
          scaling[[length(scaling) + 1L]] <- data.table(fold = k, variable = v, mean = m, sd = s)
        } else if (!all(unique(as.character(test[[v]])) %in% unique(as.character(train[[v]])))) {
          stop(label, ": held-out factor level absent from training: ", v)
        }
      }
      for (name in names(models)) {
        warning_text <- character()
        form <- as.formula(paste("Surv(followup_years, event) ~", paste(models[[name]], collapse = " + ")))
        fit <- withCallingHandlers(
          coxph(form, data = train, ties = "efron", na.action = na.fail,
                control = coxph.control(iter.max = 100)),
          warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") })
        if (any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit))) ||
            any(grepl("converg|infinite", warning_text, ignore.case = TRUE)))
          stop(label, ": nonestimable model ", name, " fold ", k, "; ", paste(warning_text, collapse = " / "))
        # Centered hazard and sample-referenced LP avoid exp() of a large intercept.
        h <- basehaz(fit, centered = TRUE)
        h5 <- if (any(h$time <= horizon)) max(h$hazard[h$time <= horizon]) else 0
        lp <- as.numeric(predict(fit, newdata = test, type = "lp", reference = "sample"))
        risk <- -expm1(-h5 * exp(lp))
        if (any(!is.finite(risk))) stop(label, ": invalid prediction ", name, " fold ", k)
        x[fold == k, (paste0("risk_", name)) := pmin(1 - 1e-8, pmax(1e-8, risk))]
        i <- i + 1L
        log[[i]] <- data.table(module = label, fold = k, model = name,
          training_people = nrow(train), training_events = sum(train$event),
          test_people = nrow(test), parameters = length(coef(fit)),
          warning = paste(warning_text, collapse = " / "))
        coefficients[[i]] <- data.table(fold = k, model = name, term = names(coef(fit)),
                                        coefficient = as.numeric(coef(fit)))
      }
    }
    risk_cols <- paste0("risk_", names(models))
    if (any(!is.finite(as.matrix(x[, ..risk_cols])))) stop(label, ": out-of-fold predictions incomplete")
    list(data = x, fits = rbindlist(log), scaling = unique(rbindlist(scaling)),
         coefficients = rbindlist(coefficients))
  }

  ipcw <- function(x) {
    s <- survfit(Surv(followup_years, 1L - event) ~ 1, data = x)
    G <- function(t, left = FALSE) {
      j <- findInterval(t, s$time, left.open = left)
      c(1, as.numeric(s$surv))[j + 1L]
    }
    y <- as.integer(x$event == 1L & x$followup_years <= horizon)
    control <- x$followup_years > horizon
    raw_G <- rep(NA_real_, nrow(x)); w <- numeric(nrow(x))
    raw_G[y == 1L] <- G(x$followup_years[y == 1L], left = TRUE)
    raw_G[control] <- G(horizon)
    eligible <- y == 1L | control
    w[eligible] <- 1 / pmax(.02, raw_G[eligible])
    list(y = y, w = w, G5 = G(horizon), raw_G = raw_G,
         truncated = sum(raw_G < .02, na.rm = TRUE))
  }
  weighted_auc <- function(y, risk, weight) {
    a <- data.table(y = y, risk = risk, weight = weight)[weight > 0,
      .(case_weight = sum(weight[y == 1L]), control_weight = sum(weight[y == 0L])), by = risk]
    setorder(a, risk)
    den <- sum(a$case_weight) * sum(a$control_weight)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    sum(a$case_weight * (cumsum(a$control_weight) - .5 * a$control_weight)) / den
  }
  metric_vector <- function(x, model_names) {
    a <- ipcw(x)
    values <- unlist(lapply(model_names, function(name) {
      risk <- x[[paste0("risk_", name)]]
      c(AUC = weighted_auc(a$y, risk, a$w), Brier = sum(a$w * (a$y - risk)^2) / nrow(x))
    }), use.names = FALSE)
    setNames(values, as.vector(t(outer(model_names, c("AUC", "Brier"), paste, sep = "."))))
  }
  paired_performance <- function(x, model_names, contrast_matrix, prefix) {
    point <- metric_vector(x, model_names)
    if (!identical(colnames(contrast_matrix), names(point))) stop("Contrast metric columns do not align")
    set.seed(bootstrap_seed)
    boot <- vapply(seq_len(bootstrap_n), function(b) {
      z <- x[sample.int(nrow(x), replace = TRUE)]
      metric_vector(z, model_names)
    }, numeric(length(point)))
    rownames(boot) <- names(point)
    ci <- t(apply(boot, 1L, interval))
    performance <- data.table(model = rep(model_names, each = 2L),
      metric = rep(c("AUC", "Brier"), length(model_names)),
      estimate = as.numeric(point), lower = ci[, 1L], upper = ci[, 2L],
      valid_resamples = rowSums(is.finite(boot)), people = nrow(x))
    delta <- as.numeric(contrast_matrix %*% point)
    delta_boot <- t(vapply(seq_len(nrow(contrast_matrix)), function(i) {
      use <- which(contrast_matrix[i, ] != 0)
      colSums(boot[use, , drop = FALSE] * contrast_matrix[i, use])
    }, numeric(ncol(boot))))
    dci <- t(apply(delta_boot, 1L, interval))
    contrasts <- data.table(comparison = rownames(contrast_matrix), estimate = delta,
      lower = dci[, 1L], upper = dci[, 2L], valid_resamples = rowSums(is.finite(delta_boot)),
      people = nrow(x), uncertainty = "paired bootstrap of fixed held-out predictions")
    save_csv(performance, paste0(prefix, "_performance.csv"))
    save_csv(contrasts, paste0(prefix, "_contrasts.csv"))
    saveRDS(list(point = point, bootstrap = boot, contrast_matrix = contrast_matrix),
            file.path(out, paste0(prefix, "_paired_bootstrap.rds")))
    list(performance = performance, contrasts = contrasts)
  }
  contrast_rows <- function(model_names, definitions) {
    cols <- as.vector(t(outer(model_names, c("AUC", "Brier"), paste, sep = ".")))
    mat <- matrix(0, nrow = length(definitions) * 2L, ncol = length(cols),
                  dimnames = list(NULL, cols))
    labels <- character(nrow(mat)); j <- 0L
    for (name in names(definitions)) for (metric in c("AUC", "Brier")) {
      j <- j + 1L; labels[j] <- paste(name, metric, sep = ".")
      z <- definitions[[name]]
      mat[j, paste(names(z), metric, sep = ".")] <- as.numeric(z)
    }
    rownames(mat) <- labels
    mat
  }
  sample_summary <- function(x, label) {
    a <- ipcw(x)
    data.table(sample = label, people = nrow(x), events_any_followup = sum(x$event),
      cases_5y = sum(a$y), observed_non_cases_5y = sum(x$followup_years > horizon),
      censored_before_or_at_5y = sum(a$w == 0), G5 = a$G5,
      max_weight = max(a$w), weights_using_floor = a$truncated)
  }

  models <- list()
  definitions <- list()
  for (tier in names(tiers)) {
    models[[tier]] <- tiers[[tier]]
    models[[paste0(tier, "_six")]] <- c(tiers[[tier]], "six_raw")
    definitions[[paste0(tier, "_MRI_increment")]] <- setNames(c(1, -1), c(paste0(tier, "_six"), tier))
  }
  depth <- list(cognitive_minus_brief = c("cognitive", "brief"),
    functional_minus_brief = c("functional", "brief"),
    full_minus_cognitive = c("full", "cognitive"),
    full_minus_functional = c("full", "functional"), full_minus_brief = c("full", "brief"))
  for (name in names(depth)) {
    z <- depth[[name]]
    definitions[[name]] <- setNames(c(1, -1, -1, 1),
      c(paste0(z[1], "_six"), z[1], paste0(z[2], "_six"), z[2]))
  }
  cat("13H62 part 1/3: assessment components, fixed 427 participants\n")
  result <- fit_oof(d, models, "assessment_components")
  oof <- result$data
  save_csv(oof[, c("RID", "fold", "followup_years", "event", paste0("risk_", names(models))), with = FALSE],
           "10_components_oof_predictions.csv")
  save_csv(result$fits, "11_components_fold_diagnostics.csv")
  save_csv(result$scaling, "12_components_training_scaling.csv")
  save_csv(result$coefficients, "13_components_fold_coefficients.csv")
  save_csv(sample_summary(d, "frozen_427_MCI"), "02_sample_information.csv")
  specs <- rbindlist(lapply(names(models), function(name) data.table(model = name, predictors = paste(models[[name]], collapse = " + "))))
  save_csv(specs, "15_components_model_specifications.csv")

  # Reproduce frozen brief/full risks to catch an unintended change in the estimand.
  check_names <- c(brief = "brief", brief_six = "brief_six", full = "extended", full_six = "extended_six")
  assert_columns(frozen, paste0("risk_", unname(check_names)), "Frozen brief/full predictions")
  merged <- merge(oof, frozen, by = "RID", suffixes = c("_new", "_frozen"), sort = FALSE)
  concordance <- rbindlist(lapply(names(check_names), function(name) {
    old <- paste0("risk_", check_names[[name]])
    new <- paste0("risk_", name)
    if (new %in% names(frozen)) new <- paste0(new, "_new")
    if (old %in% paste0("risk_", names(models))) old <- paste0(old, "_frozen")
    delta <- merged[[new]] - merged[[old]]
    data.table(new_model = name, frozen_model = check_names[[name]],
      max_absolute_risk_difference = max(abs(delta)), mean_absolute_risk_difference = mean(abs(delta)))
  }))
  save_csv(concordance, "16_frozen_prediction_reproduction.csv")
  if (any(!is.finite(concordance$max_absolute_risk_difference)) ||
      any(concordance$max_absolute_risk_difference > 1e-5))
    stop("Frozen brief/full predictions changed. Inspect 16_frozen_prediction_reproduction.csv before interpreting the new contrasts")
  component <- paired_performance(oof, names(models), contrast_rows(names(models), definitions), "14_components")

  # Hippocampal comparator: FreeSurfer volume and ICV from the same MRI.
  cat("13H62 part 2/3: same-MRI FreeSurfer hippocampal comparator\n")
  clinical_header <- names(fread(inputs["clinical"], nrows = 0L))
  hip_cols <- intersect(c("RID", "EXAMDATE", "IMAGEUID", "Hippocampus", "ICV", "FSVERSION", "FLDSTRENG"), clinical_header)
  clinical <- fread(inputs["clinical"], select = hip_cols, showProgress = FALSE)
  assert_columns(clinical, c("RID", "EXAMDATE", "IMAGEUID", "Hippocampus", "ICV"), "ADNIMERGE hippocampal source")
  clinical[, `:=`(RID = norm_rid(RID), record_date = as.IDate(EXAMDATE),
    source_image = norm_image(IMAGEUID), hippocampus_mm3 = num(Hippocampus), icv_mm3 = num(ICV))]
  anchors <- d[, .(RID, index_mri_date, risk_start, index_image = norm_image(IMAGEUID))]
  candidates <- merge(anchors, clinical, by = "RID", allow.cartesian = TRUE)
  candidates[, match_type := fifelse(!is.na(index_image) & !is.na(source_image) & index_image == source_image,
    "same_IMAGEUID", fifelse((is.na(index_image) | is.na(source_image)) &
      !is.na(record_date) & record_date == index_mri_date & record_date <= risk_start,
      "same_date_image_ID_missing", NA_character_))]
  candidates <- candidates[!is.na(match_type) & is.finite(hippocampus_mm3) & hippocampus_mm3 > 0 &
    is.finite(icv_mm3) & icv_mm3 > 0]
  candidates[, `:=`(hippocampus_icv_ratio = hippocampus_mm3 / icv_mm3,
    match_rank = fifelse(match_type == "same_IMAGEUID", 1L, 2L),
    record_gap_days = abs(as.integer(record_date - index_mri_date)))]
  setorderv(candidates, c("RID", "match_rank", "record_gap_days"), na.last = TRUE)
  candidates <- candidates[, .SD[match_rank == min(match_rank)], by = RID]
  save_csv(candidates, "20_hippocampus_source_matches.csv")
  if (nrow(candidates)) {
    # Different reprocessings of the same scan are exposed instead of averaged.
    hip <- candidates[, .(distinct_measurements = uniqueN(paste(hippocampus_mm3, icv_mm3)),
      hippocampus_mm3 = first(hippocampus_mm3), icv_mm3 = first(icv_mm3),
      hippocampus_icv_ratio = first(hippocampus_icv_ratio),
      match_type = first(match_type), source_record_date = first(record_date),
      source_rows = .N), by = RID]
  } else hip <- data.table(RID = character(), distinct_measurements = integer(),
    hippocampus_mm3 = numeric(), icv_mm3 = numeric(), hippocampus_icv_ratio = numeric(),
    match_type = character(), source_record_date = as.IDate(character()), source_rows = integer())
  hip <- merge(d[, .(RID, fold)], hip, by = "RID", all.x = TRUE, sort = FALSE)
  hip[, status := fifelse(is.na(distinct_measurements), "no_same_MRI_FreeSurfer_volume",
    fifelse(distinct_measurements > 1L, "different_measurements_same_MRI", "available"))]
  save_csv(hip, "21_hippocampus_availability.csv")
  hd <- merge(d, hip[status == "available", .(RID, hippocampus_icv_ratio, hippocampus_mm3, icv_mm3)],
              by = "RID", sort = FALSE)
  hd[, hippocampus_raw := -log(hippocampus_icv_ratio)]
  save_csv(data.table(criterion = c("frozen_MCI", "available_unique_same_MRI", "unavailable", "ambiguous_reprocessing"),
    people = c(nrow(d), nrow(hd), sum(hip$status == "no_same_MRI_FreeSurfer_volume"),
               sum(hip$status == "different_measurements_same_MRI"))), "22_hippocampus_sample_flow.csv")
  hip_status <- "not_estimable_in_matched_sample"
  hip_result <- NULL
  if (nrow(hd) > 0L && setequal(hd$fold, 1:5) && length(unique(hd$sex)) > 1L &&
      sum(hd$event == 1L & hd$followup_years <= horizon) > 0L && any(hd$followup_years > horizon)) {
    anatomy <- list(clinical = character(), six = "six_raw", hippocampus = "hippocampus_raw",
      global = "global_raw", six_hippocampus = c("six_raw", "hippocampus_raw"),
      six_global = c("six_raw", "global_raw"), all_three = c("six_raw", "hippocampus_raw", "global_raw"))
    hm <- list(); hp <- list()
    for (tier in c("brief", "full")) {
      for (name in names(anatomy)) hm[[paste(tier, name, sep = "_")]] <- c(tiers[[tier]], anatomy[[name]])
      relations <- list(six_increment = c("six", "clinical"),
        hippocampus_increment = c("hippocampus", "clinical"), global_increment = c("global", "clinical"),
        six_beyond_hippocampus = c("six_hippocampus", "hippocampus"),
        hippocampus_beyond_six = c("six_hippocampus", "six"),
        six_beyond_global = c("six_global", "global"),
        six_beyond_hippocampus_and_global = c("all_three", "hippocampus_global"))
      hm[[paste0(tier, "_hippocampus_global")]] <- c(tiers[[tier]], "hippocampus_raw", "global_raw")
      for (name in names(relations)) hp[[paste(tier, name, sep = "_")]] <-
        setNames(c(1, -1), paste(tier, relations[[name]], sep = "_"))
    }
    hip_result <- tryCatch({
      hf <- fit_oof(hd, hm, "hippocampal_benchmark")
      save_csv(hf$data[, c("RID", "fold", "followup_years", "event", "hippocampus_icv_ratio",
        paste0("risk_", names(hm))), with = FALSE], "23_hippocampus_oof_predictions.csv")
      save_csv(hf$fits, "24_hippocampus_fold_diagnostics.csv")
      save_csv(hf$scaling, "25_hippocampus_training_scaling.csv")
      save_csv(hf$coefficients, "26_hippocampus_fold_coefficients.csv")
      paired_performance(hf$data, names(hm), contrast_rows(names(hm), hp), "27_hippocampus")
    }, error = function(e) {
      writeLines(conditionMessage(e), file.path(out, "28_hippocampus_model_issue.txt"))
      NULL
    })
    if (!is.null(hip_result)) hip_status <- "completed_on_common_same_MRI_subset"
    save_csv(sample_summary(hd, "common_same_MRI_hippocampus_subset"), "29_hippocampus_sample_information.csv")
    save_csv(data.table(model = names(hm), predictors = vapply(hm, paste, character(1), collapse = " + ")),
             "29_hippocampus_model_specifications.csv")
  }
  save_csv(data.table(module = "hippocampal_benchmark", status = hip_status,
    frozen_people = nrow(d), matched_people = nrow(hd)), "29_hippocampus_status.csv")

# Summarize A/T measurement coverage.
  cat("13H62 part 3/3: A/T coverage within the original 695-person genetic frame\n")
  frame <- fread(inputs["genetic_frame"], showProgress = FALSE)
  assert_columns(frame, c("RID", "mri_date", "diagnosis"), "13H55 genetic frame")
  frame[, `:=`(RID = norm_rid(RID), mri_date = as.IDate(mri_date))]
  if (nrow(frame) != 4659L || uniqueN(frame$RID) != 695L ||
      anyDuplicated(frame[, .(RID, mri_date)])) stop("The original 695-person / 4659-scan frame changed")
  setorder(frame, RID, mri_date)
  index <- frame[, .(first_mri = first(mri_date), last_mri = last(mri_date),
    MRI_scans = .N, MRI_followup_years = as.numeric(last(mri_date) - first(mri_date)) / 365.25,
    diagnosis = first(as.character(diagnosis))), by = RID]

  prepare_marker <- function(path, date_column, value_column, label, tracer_column = NULL,
                             tracer_filter = NULL, platform_column = NULL) {
    if (!file.exists(path)) return(NULL)
    cols <- names(fread(path, nrows = 0L))
    required <- c("RID", date_column, value_column)
    assert_columns(as.data.table(setNames(rep(list(logical()), length(cols)), cols)), required, label)
    take <- intersect(c(required, tracer_column, platform_column, "qc_flag", "TRACER_SUVR_WARNING", "COMMENT"), cols)
    z <- fread(path, select = take, showProgress = FALSE)
    z[, RID := norm_rid(RID)]
    z <- z[RID %in% index$RID]
    z[, `:=`(marker_date = as.IDate(get(date_column)), raw_value = as.character(get(value_column)))]
    z[, value := num(raw_value)]
    z[, platform := if (!is.null(platform_column) && platform_column %in% names(z))
      as.character(get(platform_column)) else label]
    if (!is.null(tracer_column)) {
      assert_columns(z, tracer_column, label)
      z[, tracer := as.character(get(tracer_column))]
    } else z[, tracer := NA_character_]
    if (!is.null(tracer_filter)) z <- z[tracer %in% tracer_filter]
    z[, `:=`(marker = label,
      qc_flag_raw = if ("qc_flag" %in% names(z)) as.character(qc_flag) else NA_character_,
      warning_raw = if ("TRACER_SUVR_WARNING" %in% names(z)) as.character(TRACER_SUVR_WARNING) else NA_character_,
      comment_raw = if ("COMMENT" %in% names(z)) as.character(COMMENT) else NA_character_,
      bound_censored = grepl("^[[:space:]]*[<>]", raw_value))]
    z[, numeric_available := is.finite(value) & !is.na(marker_date) &
      (label == "A_PET_Centiloids" | value > 0)]
    z[, .(RID, marker, marker_date, value, raw_value, numeric_available, bound_censored,
           platform, tracer, qc_flag_raw, warning_raw, comment_raw)]
  }
  marker_data <- list(
    prepare_marker(inputs["amyloid"], "SCANDATE", "CENTILOIDS", "A_PET_Centiloids", "TRACER"),
    prepare_marker(inputs["tau_pet"], "SCANDATE", "META_TEMPORAL_SUVR", "T_PET_FTP", "TRACER", "FTP"),
    prepare_marker(inputs["elecsys"], "EXAMDATE", "ABETA42", "A_CSF_Elecsys_Abeta42", platform_column = "BATCH"),
    prepare_marker(inputs["alzbio3"], "DRAWDTE", "ABETA", "A_CSF_Alzbio3_Abeta", platform_column = "BATCH"),
    prepare_marker(inputs["elecsys"], "EXAMDATE", "PTAU", "T_CSF_Elecsys_pTau", platform_column = "BATCH"),
    prepare_marker(inputs["alzbio3"], "DRAWDTE", "PTAU", "T_CSF_Alzbio3_pTau", platform_column = "BATCH"))
  markers <- rbindlist(marker_data, use.names = TRUE, fill = TRUE)
  if (nrow(markers)) {
    save_csv(markers, "30_AT_source_records_within_695.csv")
    save_csv(markers[, .(records = .N, people = uniqueN(RID), numeric_records = sum(numeric_available),
      bounded_records = sum(bound_censored, na.rm = TRUE)),
      by = .(marker, platform, tracer, qc_flag_raw, warning_raw)], "31_AT_platform_QC_inventory.csv")
    nearest_marker <- function(marker_name, window_rule) {
      z <- markers[marker == marker_name & numeric_available == TRUE]
      z <- merge(index, z, by = "RID", allow.cartesian = TRUE)
      z[, gap_days := as.integer(marker_date - first_mri)]
      if (window_rule == "prior_365d") z <- z[gap_days >= -365L & gap_days <= 0L]
      if (window_rule == "near_90d") z <- z[abs(gap_days) <= 90L]
      z[, distance := abs(gap_days)]
      if (!nrow(z)) return(data.table(RID = character(), marker = character(), rule = character(),
        marker_date = as.IDate(character()), gap_days = integer(), value = numeric(),
        same_date_distinct_values = integer(), platform = character(),
        qc_flag_raw = character(), warning_raw = character()))
      z <- z[, .SD[distance == min(distance)], by = RID]
      # Prefer the earlier date if two dates are equally distant from MRI.
      z <- z[, .SD[marker_date == min(marker_date)], by = RID]
      z[, .(marker = marker_name, rule = window_rule, marker_date = first(marker_date),
        gap_days = first(gap_days), value = if (uniqueN(value) == 1L) first(value) else NA_real_,
        same_date_distinct_values = uniqueN(value), platform = paste(sort(unique(platform)), collapse = "|"),
        qc_flag_raw = paste(sort(unique(qc_flag_raw)), collapse = "|"),
        warning_raw = paste(sort(unique(warning_raw)), collapse = "|")), by = RID]
    }
    matches <- rbindlist(lapply(unique(markers$marker), function(m) rbindlist(lapply(
      c("prior_365d", "near_90d"), function(rule) nearest_marker(m, rule)))))
    save_csv(matches, "32_AT_baseline_matches.csv")
    cov <- copy(index)
    labels <- unique(markers$marker)
    for (window_rule in c("prior_365d", "near_90d")) for (m in labels) {
      name <- paste(m, window_rule, sep = "_")
      ids <- matches[marker == m & rule == window_rule & is.finite(value), RID]
      cov[, (name) := RID %in% ids]
    }
    save_csv(cov, "33_AT_person_coverage.csv")
    overlap <- list(); j <- 0L
    joint_pairs <- list(c("A_PET_Centiloids", "T_PET_FTP"),
      c("A_PET_Centiloids", "T_CSF_Elecsys_pTau"), c("A_PET_Centiloids", "T_CSF_Alzbio3_pTau"),
      c("A_CSF_Elecsys_Abeta42", "T_CSF_Elecsys_pTau"),
      c("A_CSF_Alzbio3_Abeta", "T_CSF_Alzbio3_pTau"))
    for (window_rule in c("prior_365d", "near_90d")) for (sources in joint_pairs) {
      a_marker <- sources[1L]; t_marker <- sources[2L]
      a <- paste(a_marker, window_rule, sep = "_"); b <- paste(t_marker, window_rule, sep = "_")
      if (!all(c(a, b) %in% names(cov))) next
      z <- cov[get(a) & get(b)]
      if (nrow(z)) {
        pair <- merge(matches[marker == a_marker & rule == window_rule,
          .(RID, A_date = marker_date)], matches[marker == t_marker & rule == window_rule,
          .(RID, T_date = marker_date)], by = "RID")
        pair <- pair[RID %in% z$RID]
        pair[, AT_gap_days := abs(as.integer(A_date - T_date))]
        z <- merge(z, pair[, .(RID, AT_gap_days)], by = "RID", all.x = TRUE)
      } else z[, AT_gap_days := integer()]
      by_stage <- z[, .(people = .N, MRI_scans = sum(MRI_scans),
        median_MRI_followup_years = as.numeric(median(MRI_followup_years)),
        median_AT_gap_days = as.numeric(median(AT_gap_days)), pairs_within_90_days = sum(AT_gap_days <= 90)), by = diagnosis]
      total <- data.table(diagnosis = "ALL", people = nrow(z), MRI_scans = sum(z$MRI_scans),
        median_MRI_followup_years = if (nrow(z)) median(z$MRI_followup_years) else NA_real_,
        median_AT_gap_days = if (nrow(z)) median(z$AT_gap_days) else NA_real_,
        pairs_within_90_days = sum(z$AT_gap_days <= 90))
      j <- j + 1L
      overlap[[j]] <- cbind(data.table(A_source = a_marker, T_source = t_marker,
        matching_rule = window_rule), rbind(total, by_stage, fill = TRUE))
    }
    save_csv(rbindlist(overlap, fill = TRUE), "34_AT_joint_coverage_by_entry_stage.csv")
    save_csv(matches[, .(people = .N, post_MRI_measurements = sum(gap_days > 0),
      ambiguous_same_date_measurements = sum(same_date_distinct_values > 1L),
      median_signed_gap_days = as.numeric(median(gap_days))), by = .(marker, rule)], "35_AT_timing_summary.csv")
  } else writeLines("No usable marker source records were available. Core clinical results remain complete.",
                     file.path(out, "35_AT_no_source_records.txt"))

  writeLines(c(
    "13H62 reuses the frozen 427-person / 138-case ADNI MCI frame and the original five folds.",
    "The brief screen contains age, sex, education and MMSE. Objective expansion adds mPACCdigit. Functional expansion adds CDR-SB and FAQ. Full expansion includes both.",
    "CDR-SB covers cognitive and functional domains. The functional block represents function and clinical severity ratings, while mPACCdigit represents the objective cognitive-test expansion.",
    "Models use baseline clinical and imaging values. Measurement timing is exported separately.",
    "Training-fold means and SDs transform continuous predictors. Held-out risks must reproduce the frozen brief/full predictions within 0.00001.",
    "Performance uses the original five-year IPCW definitions. Brier uses the full sample denominator. Positive AUC increments and negative Brier changes indicate improvement.",
    "Component contrasts are paired differences in MRI performance increments.",
    "Intervals use 1,000 paired participant resamples conditional on the fitted folds.",
    "Hippocampal volume and ICV come from ADNIMERGE FreeSurfer. Match MRI IMAGEUID, or date when IMAGEUID is missing. Exclude conflicting reprocessings.",
    "Hippocampal vulnerability is minus log(hippocampal volume / ICV), standardized within training folds.",
    "Hippocampal comparisons use a common matched subset with the original fold assignments.",
    "A/T outputs are coverage audits only. Prior-365-day matches are distinguished from +/-90-day matches; future measurements, platforms, warning fields and raw QC labels remain explicit.",
    "Amyloid PET, Elecsys Abeta42 and Alzbio3 Abeta are inventoried separately. Joint coverage pairs PET amyloid with each tau source, and CSF amyloid with pTau from the same assay family.",
    "FTP tau PET and the two CSF assay families are inventoried separately using their continuous values."

  ), file.path(out, "90_methods_and_interpretation.txt"))
  save_csv(data.table(module = c("assessment_components", "hippocampal_benchmark", "AT_coverage"),
    status = c("completed", hip_status, if (nrow(markers)) "coverage_exported_no_stage_model" else "no_marker_source_records")),
    "98_module_status.csv")
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  writeLines(c("13H62_COMPLETE", paste0("Results: ", out), paste0("Hippocampus: ", hip_status)),
             file.path(out, "13H62_COMPLETED.txt"))
  cat("\n13H62 complete. Results: ", out, "\n", sep = "")
})
