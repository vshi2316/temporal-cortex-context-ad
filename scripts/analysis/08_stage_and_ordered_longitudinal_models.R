# Entry-stage heterogeneity and forward/reverse longitudinal comparisons
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

local({
  options(stringsAsFactors = FALSE, scipen = 999)
  packages <- c("data.table", "lme4", "survival", "sandwich")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({
    library(data.table)
    library(lme4)
    library(survival)
  })

  root <- study_path("project", "results")
  structural_file <- file.path(root, "13F_R_only_structural_analysis",
                               "04_ADNI_R_only_longitudinal_structural_frame.rds")
  pgs_file <- file.path(root, "10A_ADNI_gMAV_PGS_assembly",
                        "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv")
  event_file <- file.path(root, "13F_R1_clinical_closure_confirmatory",
                          "09_baseline_event_frame.csv")
  inputs <- c(structural_file, pgs_file, event_file)
  if (any(!file.exists(inputs)))
    stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = "; "))

  out <- file.path(root, paste0("13H55_transition_ordered_structural_chain_",
                                "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, name) fwrite(as.data.table(x), file.path(out, name), bom = TRUE)
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  first_num <- function(x) {
    z <- num(x)
    if (any(is.finite(z))) z[which(is.finite(z))[1L]] else NA_real_
  }
  first_text <- function(x) {
    z <- trimws(as.character(x))
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z)) z[1L] else NA_character_
  }
  zscore <- function(x) {
    z <- num(x)
    (z - mean(z, na.rm = TRUE)) / sd(z, na.rm = TRUE)
  }
  clean_stage <- function(x) {
    z <- toupper(trimws(as.character(x)))
    ans <- rep(NA_character_, length(z))
    ans[z %in% c("CN", "NORMAL", "COGNITIVELY NORMAL")] <- "CN"
    ans[z %in% c("MCI", "EMCI", "LMCI", "SMC")] <- "MCI"
    ans[z %in% c("AD", "DEMENTIA", "ALZHEIMER DISEASE", "ALZHEIMER'S DISEASE")] <- "Dementia"
    ans
  }

  pgs <- fread(pgs_file, sep = "\t", showProgress = FALSE)
  pgs[, `:=`(RID = as.character(RID), exam_date__ = as.IDate(EXAMDATE))]
  setorder(pgs, RID, exam_date__)
  baseline <- pgs[, .(
    baseline_age = first_num(age_bl), sex = first_text(sex),
    education = first_num(education), diagnosis = first_text(dx_bl)
  ), by = RID]

  d <- as.data.table(readRDS(structural_file))
  d[, `:=`(RID = as.character(RID), mri_date = as.IDate(as.character(mri_date)),
            gAD = num(clinical_gAD_noAPOE_z_eur),
            gMAV = num(clinical_gMAV_noAPOE_z_eur),
            primary_complete = is.finite(num(ad_vulnerable_index)),
            IMAGEUID = as.character(IMAGEUID))]
  setorderv(d, c("RID", "mri_date", "primary_complete", "IMAGEUID"),
           c(1L, 1L, -1L, -1L), na.last = TRUE)
  d <- d[, .SD[1L], by = .(RID, mri_date)]
  setorder(d, RID, mri_date)
  d[, time_years := as.numeric(mri_date - first(mri_date)) / 365.25, by = RID]
  d[, `:=`(thickness_z = zscore(ad_vulnerable_index),
            global_z = zscore(global_cortical_thickness),
            scanner = factor(paste(FIELD_STRENGTH, FSVER, sep = "_")))]
  d <- merge(d, baseline, by = "RID", all.x = TRUE, sort = FALSE)
  d[, `:=`(sex = factor(sex), diagnosis = factor(diagnosis))]
  d <- d[is.finite(thickness_z) & is.finite(gAD) & is.finite(gMAV) & !is.na(mri_date)]
  eligible <- d[, .(has_second = any(as.integer(mri_date - min(mri_date)) >= 270L)), by = RID]
  d <- d[RID %in% eligible[has_second == TRUE, RID]]
  pc <- paste0("clinical_PC", 1:10)
  d[, (pc) := lapply(.SD, num), .SDcols = pc]
  fixed_covars <- c("baseline_age", "sex", "education", "diagnosis", "scanner", pc)
  d <- droplevels(d[complete.cases(d[, c("RID", "time_years", "gAD", "gMAV",
                                       "thickness_z", fixed_covars), with = FALSE])])
  if (uniqueN(d$RID) != 695L || nrow(d) != 4659L)
    stop("The original 695-person / 4659-scan genetic frame was not reproduced")
  d[, diagnosis := relevel(diagnosis, ref = "CN")]
  cog <- pgs[is.finite(num(mPACCdigit)) & !is.na(exam_date__),
             .(RID, cog_date = exam_date__, cognition_raw = num(mPACCdigit),
               cog_dx = as.character(DX))]
  setorder(cog, RID, cog_date)
  cog <- cog[, .SD[1L], by = .(RID, cog_date)]
  matched <- merge(d[, .(RID, mri_date)], cog, by = "RID",
                   all.x = TRUE, allow.cartesian = TRUE)
  matched[, cog_gap_days := abs(as.integer(mri_date - cog_date))]
  matched <- matched[is.finite(cog_gap_days) & cog_gap_days <= 90L]
  setorder(matched, RID, mri_date, cog_gap_days, cog_date)
  matched <- matched[, .SD[1L], by = .(RID, mri_date)]
  d <- merge(d, matched[, .(RID, mri_date, cog_date, cog_gap_days,
                            cognition_raw, cog_dx)],
             by = c("RID", "mri_date"), all.x = TRUE, sort = FALSE)
  d[, cognition_z := zscore(cognition_raw)]
  write_csv(d[, .(RID, mri_date, diagnosis, time_years, thickness_z, global_z,
                  cog_date, cog_gap_days, cognition_z, gAD, gMAV)],
            "01_frozen_frame.csv")

  covars <- c("baseline_age", "sex", "education", "scanner", pc)
  base_rhs <- paste(c("time_years * gAD", "time_years * gMAV",
                      "time_years * diagnosis", "gAD:diagnosis",
                      "gMAV:diagnosis", covars), collapse = " + ")
  full_rhs <- paste(base_rhs, "+ time_years:gAD:diagnosis + time_years:gMAV:diagnosis")
  control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
  fit_pair <- function(dat, label) {
    dat <- droplevels(copy(dat))
    dat[, diagnosis := relevel(diagnosis, ref = "CN")]
    fit_with <- function(random) {
      suffix <- if (random == "slope") "(1 + time_years | RID)" else "(1 | RID)"
      a <- lmer(as.formula(paste("thickness_z ~", base_rhs, "+", suffix)),
                data = dat, REML = FALSE, control = control)
      b <- lmer(as.formula(paste("thickness_z ~", full_rhs, "+", suffix)),
                data = dat, REML = FALSE, control = control)
      list(base = a, full = b)
    }
    fits <- tryCatch(fit_with("slope"), error = function(e) NULL)
    random <- "intercept_and_slope"
    if (is.null(fits) || isSingular(fits$base, tol = 1e-4) ||
        isSingular(fits$full, tol = 1e-4)) {
      fits <- fit_with("intercept")
      random <- "intercept_only"
    }
    comparison <- as.data.table(anova(fits$base, fits$full), keep.rownames = "model")
    b <- fixef(fits$full)
    V <- as.matrix(vcov(fits$full))
    stage <- levels(dat$diagnosis)
    contrasts <- rbindlist(lapply(stage, function(s) {
      L <- setNames(rep(0, length(b)), names(b))
      L["time_years:gAD"] <- 1
      L["time_years:gMAV"] <- -1
      if (s != "CN") {
        a <- paste0("time_years:gAD:diagnosis", s)
        m <- paste0("time_years:gMAV:diagnosis", s)
        if (!all(c(a, m) %in% names(b))) stop("Stage contrast term missing: ", s)
        L[a] <- 1
        L[m] <- -1
      }
      estimate <- sum(L * b)
      se <- sqrt(as.numeric(t(L) %*% V %*% L))
      data.table(sample = label, stage = s, estimate = estimate,
                 se = se, lower = estimate - 1.96 * se,
                 upper = estimate + 1.96 * se,
                 p = 2 * pnorm(-abs(estimate / se)),
                 people = uniqueN(dat[diagnosis == s, RID]),
                 scans = nrow(dat[diagnosis == s]))
    }))
    nonref <- setdiff(stage, "CN")
    L <- matrix(0, nrow = length(nonref), ncol = length(b),
                dimnames = list(nonref, names(b)))
    for (s in nonref) {
      L[s, paste0("time_years:gAD:diagnosis", s)] <- 1
      L[s, paste0("time_years:gMAV:diagnosis", s)] <- -1
    }
    delta <- as.numeric(L %*% b)
    W <- as.numeric(t(delta) %*% solve(L %*% V %*% t(L), delta))
    audit <- data.table(sample = label, people = uniqueN(dat$RID),
                        scans = nrow(dat), random_structure = random,
                        joint_score_stage_LR = comparison$Chisq[2L],
                        joint_score_stage_df = comparison$Df[2L],
                        joint_score_stage_p = comparison$`Pr(>Chisq)`[2L],
                        contrast_heterogeneity_chisq = W,
                        contrast_heterogeneity_df = length(nonref),
                        contrast_heterogeneity_p = pchisq(W, length(nonref),
                                                          lower.tail = FALSE))
    list(audit = audit, contrasts = contrasts, fit = fits$full)
  }

  full_stage <- fit_pair(d, "original_695")

  event <- fread(event_file, select = c("RID", "event_date", "censor_date", "site"),
                 showProgress = FALSE)
  event[, `:=`(RID = as.character(RID), event_date = as.IDate(event_date),
              censor_date = as.IDate(censor_date), site = as.character(site))]
  if (anyDuplicated(event$RID)) stop("Baseline event frame has duplicate participants")
  pre <- merge(d, event, by = "RID", all = FALSE, sort = FALSE)
  pre <- pre[(is.na(event_date) |
                (mri_date < event_date & cog_date < event_date)) &
               !is.na(censor_date) & mri_date <= censor_date &
               cog_date <= censor_date & is.finite(cognition_z)]
  setorder(pre, RID, mri_date)

  triples <- pre[, {
    if (.N < 3L) NULL else {
      i0 <- seq_len(.N - 2L)
      i1 <- i0 + 1L
      i2 <- i0 + 2L
      dt01 <- as.numeric(mri_date[i1] - mri_date[i0]) / 365.25
      dt12 <- as.numeric(mri_date[i2] - mri_date[i1]) / 365.25
      dtc01 <- as.numeric(cog_date[i1] - cog_date[i0]) / 365.25
      dtc12 <- as.numeric(cog_date[i2] - cog_date[i1]) / 365.25
      data.table(
        t0 = mri_date[i0], t1 = mri_date[i1], t2 = mri_date[i2],
        c0 = cog_date[i0], c1 = cog_date[i1], c2 = cog_date[i2],
        dt01 = dt01, dt12 = dt12, dtc01 = dtc01, dtc12 = dtc12,
        loss01 = (thickness_z[i0] - thickness_z[i1]) / dt01,
        global_loss01 = (global_z[i0] - global_z[i1]) / dt01,
        loss12 = (thickness_z[i1] - thickness_z[i2]) / dt12,
        cog01 = (cognition_z[i1] - cognition_z[i0]) / dtc01,
        cog12 = (cognition_z[i2] - cognition_z[i1]) / dtc12,
        thickness1 = thickness_z[i1], thickness2 = thickness_z[i2],
        cognition1 = cognition_z[i1], cognition2 = cognition_z[i2],
        age1 = baseline_age[i1] + time_years[i1],
        age2 = baseline_age[i2] + time_years[i2],
        sex = sex[i1], education = education[i1], diagnosis = diagnosis[i1],
        stage_t1 = clean_stage(cog_dx[i1]),
        stage_t2 = clean_stage(cog_dx[i2]),
        scanner = scanner[i1], gAD = gAD[i1], gMAV = gMAV[i1],
        event_date = event_date[i2], censor_date = censor_date[i2], site = site[i2]
      )
    }
  }, by = RID]
  triples <- triples[dt01 >= 180 / 365.25 & dt01 <= 3 &
                       dt12 >= 180 / 365.25 & dt12 <= 3 &
                       dtc01 >= 180 / 365.25 & dtc01 <= 3 &
                       dtc12 >= 180 / 365.25 & dtc12 <= 3 &
                       stage_t1 %in% c("CN", "MCI")]
  if (!nrow(triples)) stop("No nonoverlapping MRI-to-cognition intervals were available")
  triples[, `:=`(forward_order = t1 <= c1, reverse_order = c1 <= t1)]
  if (!any(triples$forward_order))
    stop("No triple preserves MRI-before-later-cognition order")
  for (v in c("loss01", "global_loss01", "loss12", "cog01", "cog12"))
    triples[, (paste0(v, "_z")) := zscore(get(v))]
  write_csv(triples, "04_nonoverlapping_interval_triples.csv")

  predementia_stage <- fit_pair(
    d[!as.character(diagnosis) %in% c("AD", "Dementia")],
    "entry_CN_EMCI_LMCI"
  )
  stage_audit <- rbind(full_stage$audit, predementia_stage$audit)
  stage_contrasts <- rbind(full_stage$contrasts, predementia_stage$contrasts)
  overlap_ids <- unique(triples[forward_order == TRUE, RID])
  if (length(overlap_ids) >= 100L &&
      "CN" %in% as.character(d[RID %in% overlap_ids, diagnosis]) &&
      uniqueN(d[RID %in% overlap_ids, diagnosis]) >= 2L) {
    overlap_stage <- fit_pair(d[RID %in% overlap_ids], "lagged_overlap")
    stage_audit <- rbind(stage_audit, overlap_stage$audit)
    stage_contrasts <- rbind(stage_contrasts, overlap_stage$contrasts)
  }
  write_csv(stage_audit, "02_stage_interaction_tests.csv")
  write_csv(stage_contrasts, "03_stage_specific_genetic_contrasts.csv")

  flow <- data.table(
    step = c("frozen_genetic_MRI", "clinical_event_overlap_with_cognition",
             "eligible_visit_triples", "forward_ordered_triples",
             "reverse_ordered_triples"),
    people = c(uniqueN(d$RID), uniqueN(pre$RID), uniqueN(triples$RID),
               uniqueN(triples[forward_order == TRUE, RID]),
               uniqueN(triples[reverse_order == TRUE, RID])),
    records = c(nrow(d), nrow(pre), nrow(triples),
                nrow(triples[forward_order == TRUE]),
                nrow(triples[reverse_order == TRUE])),
    events = c(NA_integer_, uniqueN(pre[!is.na(event_date), RID]),
               uniqueN(triples[!is.na(event_date), RID]),
               uniqueN(triples[forward_order == TRUE & !is.na(event_date), RID]),
               uniqueN(triples[reverse_order == TRUE & !is.na(event_date), RID]))
  )

  fit_clustered_lm <- function(x, formula, term, label) {
    vars <- unique(c("RID", all.vars(formula)))
    a <- x[complete.cases(x[, ..vars])]
    a <- droplevels(a)
    if (nrow(a) < 100L || uniqueN(a$RID) < 50L)
      return(data.table(model = label, term = term, status = "INSUFFICIENT_DATA",
                        people = uniqueN(a$RID), intervals = nrow(a)))
    fit <- lm(formula, data = a)
    V <- sandwich::vcovCL(fit, cluster = a$RID, type = "HC1")
    estimate <- unname(coef(fit)[term])
    se <- unname(sqrt(diag(V))[term])
    if (!is.finite(estimate) || !is.finite(se))
      stop(label, " target coefficient was not estimable")
    p <- 2 * pt(-abs(estimate / se), df = uniqueN(a$RID) - 1L)
    data.table(model = label, term = term, status = "FIT",
               estimate = estimate, se = se,
               lower = estimate - 1.96 * se, upper = estimate + 1.96 * se,
               p = p, people = uniqueN(a$RID), intervals = nrow(a))
  }

  adjustment <- paste(c("thickness1", "cognition1", "age1", "sex",
                        "education", "stage_t1", "scanner", "gAD", "gMAV",
                        "dt01", "dt12"), collapse = " + ")
  forward <- fit_clustered_lm(
    triples[forward_order == TRUE],
    as.formula(paste("cog12_z ~ loss01_z + cog01_z + dtc01 + dtc12 +",
                     adjustment)),
    "loss01_z", "earlier_MRI_loss_to_later_cognitive_change"
  )
  forward_global <- fit_clustered_lm(
    triples[forward_order == TRUE],
    as.formula(paste("cog12_z ~ loss01_z + global_loss01_z + cog01_z +",
                     "dtc01 + dtc12 +", adjustment)),
    "loss01_z", "earlier_six_region_loss_beyond_global_loss"
  )
  reverse <- fit_clustered_lm(
    triples[reverse_order == TRUE],
    as.formula(paste("loss12_z ~ cog01_z + loss01_z + dtc01 + dtc12 +",
                     adjustment)),
    "cog01_z", "earlier_cognitive_change_to_later_MRI_loss"
  )
  write_csv(rbind(forward, forward_global, reverse, fill = TRUE),
            "05_lead_lag_models.csv")

  setorder(triples, RID, t2)
  first_triple <- triples[forward_order == TRUE, .SD[1L], by = RID]
  first_triple[, risk_start := pmax(t2, c2)]
  first_triple <- first_triple[
    (is.na(event_date) | event_date > risk_start) &
      censor_date > risk_start &
      stage_t2 %in% c("CN", "MCI")
  ]
  first_triple[, `:=`(
    followup_years = (fifelse(!is.na(event_date), as.numeric(event_date),
                              as.numeric(censor_date)) - as.numeric(risk_start)) / 365.25,
    event = !is.na(event_date), stage_t2 = factor(stage_t2),
    site = factor(site)
  )]
  event_vars <- c("followup_years", "event", "loss01_z", "cog12_z",
                  "thickness2", "cognition2", "age2", "sex", "education",
                  "gAD", "gMAV", "stage_t2", "site")
  first_triple <- droplevels(first_triple[
    complete.cases(first_triple[, ..event_vars]) & followup_years > 0
  ])
  write_csv(first_triple, "06_first_triple_event_frame.csv")
  flow <- rbind(flow, data.table(step = "first_triple_future_event",
                                people = uniqueN(first_triple$RID),
                                records = nrow(first_triple),
                                events = sum(first_triple$event)))
  write_csv(flow, "00_sample_flow.csv")

  if (nrow(first_triple) >= 100L && sum(first_triple$event) >= 20L) {
    event_rhs <- paste(
      "thickness2 + cognition1 + age2 + sex + education + gAD + gMAV",
      "+ strata(stage_t2) + cluster(site)"
    )
    specs <- list(
      earlier_structure = paste("loss01_z +", event_rhs),
      structure_then_cognition = paste("loss01_z + cog12_z +", event_rhs),
      current_state = paste("loss01_z + cog12_z + cognition2 +", event_rhs)
    )
    fits <- lapply(specs, function(rhs) {
      coxph(as.formula(paste("Surv(followup_years, event) ~", rhs)),
            data = first_triple, ties = "efron")
    })
    co <- rbindlist(lapply(names(fits), function(model) {
      x <- as.data.table(summary(fits[[model]])$coefficients,
                         keep.rownames = "term")
      x[, `:=`(model = model, people = nrow(first_triple),
               events = sum(first_triple$event))]
      x
    }), fill = TRUE)
    write_csv(co, "07_subsequent_dementia_cox.csv")
    saveRDS(fits, file.path(out, "07_subsequent_dementia_cox_fits.rds"))
  } else {
    write_csv(data.table(status = "INSUFFICIENT_FUTURE_EVENTS",
                         people = nrow(first_triple),
                         events = sum(first_triple$event)),
              "07_subsequent_dementia_cox.csv")
  }

  writeLines(c(

    "Stage interaction: direct gAD-minus-gMAV slope differences across entry diagnoses.",
    "Lead-lag: clinical exam dates must be within 90 days of MRI and ordered. Forward MRI interval ends before later cognitive interval begins; the reverse uses its own ordered subset.",
    "Global cortex is a comparator for the six-region thickness change.",
    "Event: first future dementia classification after both third-visit MRI and cognitive exam in the first eligible forward triple.",
    "Cox tiers add intervening cognitive change, then cognitive level at the third MRI.",

    paste0("Output directory: ", out)
  ), file.path(out, "08_readme.txt"))
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  cat("13H55 complete. Results: ", out, "\n", sep = "")
})
