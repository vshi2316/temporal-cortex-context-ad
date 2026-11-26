# CSF joint modification and matched-sample genetic trajectories
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

# The outcome and genetic scores retain the original manuscript definitions.
local({
  options(stringsAsFactors = FALSE, scipen = 999)
  packages <- c("data.table", "lme4")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install R packages: ", paste(absent, collapse = ", "))
  suppressPackageStartupMessages({library(data.table); library(lme4)})

  root <- study_path("project", "results")
  prior <- file.path(root, "13H62_ADNI_components_hippocampus_AT_audit_release")
  inputs <- c(
    frozen = file.path(root, "13H55_transition_ordered_structural_chain_release", "01_frozen_frame.csv"),
    matches = file.path(prior, "32_AT_baseline_matches.csv"),
    completed = file.path(prior, "13H62_COMPLETED.txt"),
    structure = file.path(root, "13F_R_only_structural_analysis", "04_ADNI_R_only_longitudinal_structural_frame.rds"),
    pgs = file.path(root, "10A_ADNI_gMAV_PGS_assembly", "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"))
  if (any(!file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = "; "))
  out <- file.path(root, paste0("13H63_ADNI_pathology_conditioned_genetic_trajectories_", "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(x, name) fwrite(as.data.table(x), file.path(out, name), bom = TRUE)
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  first_num <- function(x) {z <- num(x); if (any(is.finite(z))) z[which(is.finite(z))[1L]] else NA_real_}
  first_text <- function(x) {z <- trimws(as.character(x)); z <- z[!is.na(z) & nzchar(z)]; if (length(z)) z[1L] else NA_character_}
  norm_rid <- function(x) {
    z <- sub("[.]0+$", "", trimws(as.character(x)))
    z <- sub("^.*?([0-9]+)$", "\\1", z); z <- sub("^0+", "", z)
    z[z == ""] <- "0"; z[is.na(x) | !nzchar(trimws(as.character(x)))] <- NA_character_; z
  }
  zscore <- function(x) {z <- num(x); (z - mean(z, na.rm = TRUE)) / sd(z, na.rm = TRUE)}
  need_cols <- function(x, cols, label) {
    miss <- setdiff(cols, names(x)); if (length(miss)) stop(label, " missing: ", paste(miss, collapse = ", "))
  }
  bootstrap_n <- 500L
  bootstrap_seed <- 20261120L
  primary_batch <- "UPENNBIOMK9"
  min_followup_days <- 270L
  pc <- paste0("clinical_PC", 1:10)
  control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

  # Exact term lookup is invariant to the order R prints an interaction.
  term_index <- function(coef_names, tokens) {
    same <- vapply(strsplit(coef_names, ":", fixed = TRUE), function(x) setequal(x, tokens), logical(1))
    if (sum(same) != 1L) stop("Interaction term absent or ambiguous: ", paste(tokens, collapse = ":"))
    which(same)
  }
  contrast_vector <- function(fit, stage = "CN", amyloid = 0, ptau = 0,
                              moderator = NULL, include_stage = TRUE, include_pathology = TRUE) {
    cn <- names(fixef(fit)); L <- setNames(numeric(length(cn)), cn)
    add <- function(tokens, value) {j <- term_index(cn, tokens); L[j] <<- L[j] + value}
    if (!is.null(moderator)) {
      add(c("t", "gAD", moderator), 1); add(c("t", "gMAV", moderator), -1)
    } else {
      add(c("t", "gAD"), 1); add(c("t", "gMAV"), -1)
      if (stage != "CN" && include_stage) {
        add(c("t", "gAD", paste0("diagnosis", stage)), 1)
        add(c("t", "gMAV", paste0("diagnosis", stage)), -1)
      }
      if (include_pathology) {
        add(c("t", "gAD", "amyloid_z"), amyloid); add(c("t", "gMAV", "amyloid_z"), -amyloid)
        add(c("t", "gAD", "ptau_z"), ptau); add(c("t", "gMAV", "ptau_z"), -ptau)
      }
    }
    L
  }
  estimate_L <- function(fit, L) {
    b <- fixef(fit); V <- as.matrix(vcov(fit)); value <- sum(L * b)
    se <- sqrt(max(0, as.numeric(t(L) %*% V %*% L)))
    p <- if (se > 0) 2 * pnorm(-abs(value / se)) else NA_real_
    data.table(estimate = value, se = se, lower = value - 1.96 * se, upper = value + 1.96 * se, p_wald = p)
  }
  joint_wald <- function(fit, L) {
    V <- L %*% as.matrix(vcov(fit)) %*% t(L); value <- as.numeric(L %*% fixef(fit))
    if (qr(V)$rank < nrow(L)) stop("Joint contrast covariance is not full rank")
    W <- max(0, as.numeric(t(value) %*% solve(V, value)))
    data.table(chisq = W, df = nrow(L), p = pchisq(W, nrow(L), lower.tail = FALSE))
  }
  fit_once <- function(dat, rhs, random = "intercept_slope", allow_nuisance_drop = FALSE) {
    random_term <- if (random == "intercept_slope") "(1 + t | RID)" else "(1 | RID)"
    form <- as.formula(paste("thickness_z ~", rhs, "+", random_term))
    mm <- model.matrix(as.formula(paste("~", rhs)), dat)
    if (!allow_nuisance_drop && qr(mm)$rank < ncol(mm)) stop("Fixed model matrix is rank deficient")
    warnings <- character()
    fit <- withCallingHandlers(suppressMessages(lmer(form, data = dat, REML = FALSE, control = control, na.action = na.fail)),
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
    messages <- unlist(fit@optinfo$conv$lme4$messages)
    nonboundary <- messages[!grepl("boundary.*singular", messages, ignore.case = TRUE)]
    optcode <- unlist(fit@optinfo$conv$opt)
    if (length(nonboundary) || (length(optcode) && any(optcode != 0)) ||
        any(!is.finite(fixef(fit))) || any(!is.finite(as.matrix(vcov(fit)))))
      stop("Mixed model convergence issue: ", paste(c(warnings, nonboundary), collapse = " / "))
    list(fit = fit, warning = paste(unique(warnings), collapse = " / "), singular = isSingular(fit, tol = 1e-4),
      dropped_terms = names(attr(getME(fit, "X"), "col.dropped")))
  }
  model_rhs <- function(covars) {
    B <- paste(c("t * gAD", "t * gMAV", "t * diagnosis", "gAD:diagnosis", "gMAV:diagnosis", covars), collapse = " + ")
    stage <- "t:gAD:diagnosis + t:gMAV:diagnosis"
    lower <- "(t + gAD + gMAV) * amyloid_z + (t + gAD + gMAV) * ptau_z"
    path <- "t:gAD:amyloid_z + t:gMAV:amyloid_z + t:gAD:ptau_z + t:gMAV:ptau_z"
    list(clinical = paste(B, stage, sep = " + "),
      pathology_lower = paste(B, lower, sep = " + "),
      pathology = paste(B, lower, path, sep = " + "),
      joint_lower = paste(B, stage, lower, sep = " + "),
      joint = paste(B, stage, lower, path, sep = " + "))
  }
  fit_set <- function(dat, covars, label) {
    dat <- droplevels(copy(dat)); dat[, diagnosis := relevel(factor(diagnosis), "CN")]
    if (!all(c("CN", "EMCI", "LMCI") %in% levels(dat$diagnosis))) stop(label, ": entry stages missing")
    covars <- covars[vapply(covars, function(v) uniqueN(dat[[v]]) > 1L, logical(1))]
    batch_terms <- if ("assay_batch" %in% names(dat) && uniqueN(dat$assay_batch) > 1L)
      c("assay_batch", "t:assay_batch") else character()
    rhs <- model_rhs(c(covars, batch_terms))
    fit_all <- function(random) lapply(rhs, function(r) fit_once(dat, r, random))
    failure <- ""; random <- "intercept_slope"
    fits <- tryCatch(fit_all(random), error = function(e) {failure <<- conditionMessage(e); NULL})
    if (is.null(fits) || any(vapply(fits, function(x) x$singular, logical(1)))) {
      if (!nzchar(failure)) failure <- "At least one model had a singular intercept-slope covariance"
      random <- "intercept_only"; fits <- fit_all(random)
    }
    diagnostics <- rbindlist(lapply(names(fits), function(m) {
      f <- fits[[m]]$fit
      data.table(sample = label, model = m, people = uniqueN(dat$RID), scans = nrow(dat),
        fixed_parameters = length(fixef(f)), random_structure = random, singular = fits[[m]]$singular,
        warning = fits[[m]]$warning, shared_fallback_reason = failure, AIC = AIC(f), logLik = as.numeric(logLik(f)))
    }))
    moderators <- rbindlist(lapply(c("amyloid_z", "ptau_z"), function(v) {
      z <- estimate_L(fits$joint$fit, contrast_vector(fits$joint$fit, moderator = v))
      z[, `:=`(sample = label, moderator = v, people = uniqueN(dat$RID), scans = nrow(dat))]; z
    }))
    moderators[, p_BH_two_moderators := p.adjust(p_wald, "BH")]
    L <- rbind(contrast_vector(fits$joint$fit, moderator = "amyloid_z"),
               contrast_vector(fits$joint$fit, moderator = "ptau_z"))
    joint <- joint_wald(fits$joint$fit, L)
    joint[, `:=`(sample = label, test = "gAD_minus_gMAV_time_contrast_modified_by_A_and_pTau", people = uniqueN(dat$RID), scans = nrow(dat))]
    stages <- rbindlist(lapply(c("clinical", "joint_lower", "joint"), function(m) {
      rbindlist(lapply(levels(dat$diagnosis), function(s) {
        z <- estimate_L(fits[[m]]$fit, contrast_vector(fits[[m]]$fit, stage = s,
          include_stage = TRUE, include_pathology = m == "joint"))
        z[, `:=`(sample = label, model = m, stage = s, amyloid_z = 0, ptau_z = 0,
          people_in_stage = uniqueN(dat[diagnosis == s, RID]))]; z
      }))
    }))
    stage_tests <- rbindlist(lapply(c("clinical", "joint_lower", "joint"), function(m) {
      f <- fits[[m]]$fit
      L <- do.call(rbind, lapply(setdiff(levels(dat$diagnosis), "CN"), function(s)
        contrast_vector(f, s, include_pathology = m == "joint") -
          contrast_vector(f, "CN", include_pathology = m == "joint")))
      z <- joint_wald(f, L); z[, `:=`(sample = label, model = m, test = "direct_stage_heterogeneity_of_genetic_time_contrast")]; z
    }))
    lrt <- rbindlist(lapply(list(c("joint_lower", "joint"), c("pathology", "joint")), function(pair) {
      a <- anova(fits[[pair[1]]]$fit, fits[[pair[2]]]$fit, refit = FALSE)
      data.table(sample = label, reduced = pair[1], full = pair[2], chisq = a$Chisq[2],
        df = a$Df[2], p = a$`Pr(>Chisq)`[2], interpretation = if (pair[1] == "joint_lower")
          "All four genetic-time-pathology coefficients jointly" else "All genetic-time-clinical-stage coefficients conditional on pathology")
    }))
    list(data = dat, fits = fits, rhs = rhs, random = random, diagnostics = diagnostics,
      moderators = moderators, joint_test = joint, stage_contrasts = stages, stage_tests = stage_tests, lrt = lrt)
  }
  participant_bootstrap <- function(result, B, seed) {
    x <- result$data; groups <- split(seq_len(nrow(x)), as.character(x$RID)); ids <- names(groups)
    names_coef <- names(fixef(result$fits$joint$fit))
    L <- rbind(contrast_vector(result$fits$joint$fit, moderator = "amyloid_z"),
               contrast_vector(result$fits$joint$fit, moderator = "ptau_z"))
    set.seed(seed)
    values <- matrix(NA_real_, B, 2L, dimnames = list(NULL, c("amyloid_z", "ptau_z")))
    log <- vector("list", B)
    for (b in seq_len(B)) {
      draw <- sample(ids, length(ids), replace = TRUE)
      z <- rbindlist(lapply(seq_along(draw), function(i) {
        a <- copy(x[groups[[draw[i]]]]); a[, RID := paste0("boot_", i)]; a
      }))
      failure <- ""; singular <- NA; dropped <- ""
      f <- tryCatch(fit_once(z, result$rhs$joint, result$random, allow_nuisance_drop = TRUE),
        error = function(e) {failure <<- conditionMessage(e); NULL})
      if (!is.null(f)) {
        needed <- names_coef[colSums(abs(L)) > 0]
        beta <- fixef(f$fit)
        if (!all(needed %in% names(beta))) failure <- "A target genetic-time-pathology coefficient was not estimable" else {
          values[b, ] <- as.numeric(L[, names(beta), drop = FALSE] %*% beta); singular <- f$singular
          dropped <- paste(f$dropped_terms, collapse = "|")
        }
      }
      log[[b]] <- data.table(resample = b, valid = all(is.finite(values[b, ])), singular = singular,
        dropped_nuisance_terms = dropped, error = failure)
      if (b %% 25L == 0L || b == B) cat("13H63 participant bootstrap ", b, "/", B, "\n", sep = "")
    }
    ci <- rbindlist(lapply(colnames(values), function(v) {
      a <- values[, v]; a <- a[is.finite(a)]
      ok <- length(a) >= ceiling(.8 * B)
      q <- if (ok) unname(quantile(a, c(.025, .975), names = FALSE)) else c(NA_real_, NA_real_)
      data.table(moderator = v, lower_bootstrap = q[1], upper_bootstrap = q[2], valid_resamples = length(a),
        requested_resamples = B, interval_status = if (ok) "available" else "too_many_failed_resamples")
    }))
    list(values = as.data.table(values), log = rbindlist(log), intervals = ci)
  }

  manifest <- data.table(input = names(inputs), path = unname(inputs), bytes = file.info(inputs)$size,
    md5 = unname(tools::md5sum(inputs)))
  save_csv(manifest, "00_input_manifest.csv")
  writeLines(c(
    "Question: does baseline soluble AD pathology modify the original gAD-minus-gMAV six-region cortical-thickness time contrast, conditional on entry clinical stage?",
    "Primary: originally nondemented participants, same-date Elecsys Abeta42 and pTau, UPENNBIOMK9 only, both measured within 90 days of the original first MRI.",
    "Use only MRIs on or after both pathology measurements. At least two dates spanning 270 days are required; the first retained MRI defines time zero.",
    "Amyloid burden is -log(Abeta42). pTau is log(pTau). Primary-participant means and SDs define their scales.",
    "Five nested models use identical rows and a shared random-effects structure. The primary test is a two-df direct contrast of the amyloid and pTau modifiers of gAD-minus-gMAV time effects.",
    "The two modifier Wald P values receive BH adjustment as a two-test family. Other model/stage contrasts are labelled secondary. Bootstrap intervals refit the joint model by participant.",
    "Two fixed sensitivities use the same batch across all entry stages, and all matched Elecsys batches among entry-nondemented participants with batch and batch-time terms.",
    "Pooled-batch sensitivity models include batch and batch-time terms."
  ), file.path(out, "00_locked_analysis_plan.txt"))

  cat("13H63 1/4: reproduce the original genetic MRI frame\n")
  frozen <- fread(inputs["frozen"], showProgress = FALSE)
  need_cols(frozen, c("RID", "mri_date", "diagnosis", "time_years", "thickness_z", "gAD", "gMAV"), "13H55")
  frozen[, `:=`(RID = norm_rid(RID), mri_date = as.IDate(mri_date))]
  if (nrow(frozen) != 4659L || uniqueN(frozen$RID) != 695L || anyDuplicated(frozen[, .(RID, mri_date)]))
    stop("Original 695-person / 4659-scan frame changed")
  pgs <- fread(inputs["pgs"], select = c("RID", "EXAMDATE", "age_bl", "sex", "education", "dx_bl"), showProgress = FALSE)
  pgs[, `:=`(RID = norm_rid(RID), exam_date__ = as.IDate(EXAMDATE))]
  setorder(pgs, RID, exam_date__)
  baseline <- pgs[, .(baseline_age = first_num(age_bl), sex = first_text(sex),
                      education = first_num(education), original_diagnosis = first_text(dx_bl)), by = RID]
  raw <- as.data.table(readRDS(inputs["structure"]))
  need_cols(raw, c("RID", "mri_date", "ad_vulnerable_index", "IMAGEUID", "FIELD_STRENGTH", "FSVER",
    "clinical_gAD_noAPOE_z_eur", "clinical_gMAV_noAPOE_z_eur", pc), "Original structural source")
  raw[, `:=`(RID = norm_rid(RID), mri_date = as.IDate(as.character(mri_date)), IMAGEUID = as.character(IMAGEUID),
    primary_complete = is.finite(num(ad_vulnerable_index)))]
  setorderv(raw, c("RID", "mri_date", "primary_complete", "IMAGEUID"), c(1L, 1L, -1L, -1L), na.last = TRUE)
  raw <- raw[, .SD[1L], by = .(RID, mri_date)]
  raw[, `:=`(source_thickness_z = zscore(ad_vulnerable_index),
    source_gAD = num(clinical_gAD_noAPOE_z_eur), source_gMAV = num(clinical_gMAV_noAPOE_z_eur),
    scanner = paste(FIELD_STRENGTH, FSVER, sep = "_"))]
  fields <- c("RID", "mri_date", "source_thickness_z", "source_gAD", "source_gMAV", "scanner", pc)
  d <- merge(frozen, raw[, ..fields], by = c("RID", "mri_date"), all.x = TRUE, sort = FALSE)
  d <- merge(d, baseline, by = "RID", all.x = TRUE, sort = FALSE)
  if (nrow(d) != 4659L || anyNA(d$source_thickness_z) || anyNA(d$baseline_age)) stop("Frozen rows could not be linked")
  check <- rbindlist(lapply(c(thickness_z = "source_thickness_z", gAD = "source_gAD", gMAV = "source_gMAV"), function(v) {
    old <- names(c(thickness_z = "source_thickness_z", gAD = "source_gAD", gMAV = "source_gMAV"))[match(v, c("source_thickness_z", "source_gAD", "source_gMAV"))]
    data.table(variable = old, max_absolute_difference = max(abs(num(d[[old]]) - num(d[[v]]))))
  }))
  save_csv(check, "01_frozen_value_reproduction.csv")
  if (any(!is.finite(check$max_absolute_difference)) || any(check$max_absolute_difference > 1e-9) ||
      any(as.character(d$diagnosis) != d$original_diagnosis)) stop("Original genetic or thickness values changed")
  for (v in c("gAD", "gMAV", "thickness_z", "time_years", pc)) d[, (v) := num(get(v))]
  setorder(d, RID, mri_date)

  cat("13H63 2/4: freeze same-date, same-batch pathology and post-measurement MRI\n")
  matched <- fread(inputs["matches"], showProgress = FALSE)
  need_cols(matched, c("RID", "marker", "rule", "marker_date", "gap_days", "value", "platform", "same_date_distinct_values"), "13H62 matches")
  matched[, `:=`(RID = norm_rid(RID), marker_date = as.IDate(marker_date), value = num(value))]
  a <- matched[marker == "A_CSF_Elecsys_Abeta42" & rule == "near_90d",
    .(RID, A_date = marker_date, A_gap = num(gap_days), abeta42 = value, A_batch = platform, A_distinct = same_date_distinct_values)]
  p <- matched[marker == "T_CSF_Elecsys_pTau" & rule == "near_90d",
    .(RID, T_date = marker_date, T_gap = num(gap_days), ptau = value, T_batch = platform, T_distinct = same_date_distinct_values)]
  if (anyDuplicated(a$RID) || anyDuplicated(p$RID)) stop("Duplicate pathology matching records")
  # One entry diagnosis is frozen for each person.
  people <- d[, .(diagnosis = as.character(first(diagnosis)), first_mri = first(mri_date)), by = RID]
  people <- merge(merge(people, a, by = "RID", all.x = TRUE), p, by = "RID", all.x = TRUE)
  people[, status := "eligible"]
  people[!is.finite(abeta42) | abeta42 <= 0 | !is.finite(ptau) | ptau <= 0 |
    is.na(A_date) | is.na(T_date) | A_distinct != 1L | T_distinct != 1L, status := "missing_or_ambiguous_numeric_pathology"]
  people[status == "eligible" & (abs(as.integer(A_date - first_mri)) > 90L |
    abs(as.integer(T_date - first_mri)) > 90L), status := "outside_locked_90_day_window"]
  people[status == "eligible" & A_date != T_date, status := "pathology_dates_differ"]
  people[status == "eligible" & (is.na(A_batch) | is.na(T_batch) | A_batch != T_batch |
    grepl("|", A_batch, fixed = TRUE)), status := "assay_batch_not_unique_or_not_shared"]
  people[, `:=`(assay_batch = A_batch, pathology_date = A_date)]
  save_csv(people, "02_pathology_matching_by_person.csv")
  save_csv(people[, .(people = .N), by = .(status, assay_batch, diagnosis)], "03_pathology_matching_flow.csv")
  scan <- merge(d, people[status == "eligible", .(RID, abeta42, ptau, assay_batch, pathology_date)], by = "RID")
  scan <- scan[mri_date >= pathology_date]
  setorder(scan, RID, mri_date)
  scan[, `:=`(landmark_date = first(mri_date), t = as.numeric(mri_date - first(mri_date)) / 365.25,
    age_landmark = baseline_age + first(time_years)), by = RID]
  scan[, (pc) := lapply(.SD, num), .SDcols = pc]
  required <- c("RID", "mri_date", "t", "thickness_z", "gAD", "gMAV", "diagnosis", "sex",
    "education", "age_landmark", "scanner", pc, "abeta42", "ptau", "assay_batch")
  complete <- complete.cases(scan[, ..required])
  save_csv(data.table(step = c("matched_post_pathology_MRI_rows", "rows_with_complete_original_covariates"),
    scans = c(nrow(scan), sum(complete)), people = c(uniqueN(scan$RID), uniqueN(scan[complete, RID]))), "04_post_pathology_covariate_flow.csv")
  scan <- scan[complete]
  eligible <- scan[, .(scans = .N, followup_days = as.integer(max(mri_date) - min(mri_date))), by = RID]
  save_csv(eligible, "05_post_pathology_followup_by_person.csv")
  scan <- scan[RID %in% eligible[scans >= 2L & followup_days >= min_followup_days, RID]]
  # Reset time only after complete-row filtering.
  setorder(scan, RID, mri_date)
  scan[, `:=`(landmark_date = first(mri_date), t = as.numeric(mri_date - first(mri_date)) / 365.25,
    age_landmark = baseline_age + first(time_years)), by = RID]
  scan[, `:=`(sex = factor(sex), diagnosis = factor(diagnosis), scanner = factor(scanner), assay_batch = factor(assay_batch))]
  primary <- droplevels(scan[as.character(assay_batch) == primary_batch & as.character(diagnosis) %in% c("CN", "EMCI", "LMCI")])
  if (!nrow(primary)) stop("No eligible participants in the locked primary batch")
  reference <- unique(primary[, c("RID", "abeta42", "ptau", "age_landmark", "education", pc), with = FALSE])
  if (anyDuplicated(reference$RID)) stop("Person-level genetic covariates differ across MRI visits")
  reference[, `:=`(amyloid_log = -log(abeta42), ptau_log = log(ptau))]
  scale_vars <- c("amyloid_log", "ptau_log", "age_landmark", "education", pc)
  scaling <- rbindlist(lapply(scale_vars, function(v) data.table(variable = v,
    mean = mean(reference[[v]]), sd = sd(reference[[v]]), reference_people = nrow(reference))))
  if (any(!is.finite(scaling$sd)) || any(scaling$sd <= 0)) stop("A primary reference variable has no variation")
  scan[, `:=`(amyloid_log = -log(abeta42), ptau_log = log(ptau))]
  for (v in scale_vars) {
    row <- scaling[variable == v]; target <- if (v == "amyloid_log") "amyloid_z" else if (v == "ptau_log") "ptau_z" else paste0(v, "_z")
    scan[, (target) := (get(v) - row$mean) / row$sd]
  }
  save_csv(scaling, "06_primary_participant_scaling.csv")
  sample_definitions <- list(primary_batch9_entry_nondementia = scan[as.character(assay_batch) == primary_batch & as.character(diagnosis) %in% c("CN", "EMCI", "LMCI")],
    sensitivity_batch9_all_entry_stages = scan[as.character(assay_batch) == primary_batch],
    sensitivity_all_Elecsys_batches_entry_nondementia = scan[as.character(diagnosis) %in% c("CN", "EMCI", "LMCI")])
  flow <- rbindlist(lapply(names(sample_definitions), function(label) {
    x <- sample_definitions[[label]]
    person <- x[, .(stage = as.character(first(diagnosis)), followup = max(t), scans = .N,
      landmark_date = first(landmark_date)), by = RID]
    rbind(data.table(sample = label, stage = "ALL", people = nrow(person), scans = nrow(x),
      median_followup_years = median(person$followup)),
      person[, .(sample = label, people = .N, scans = sum(scans), median_followup_years = as.numeric(median(followup))), by = stage], fill = TRUE)
  }))
  save_csv(flow, "07_analysis_sample_flow.csv")
  save_csv(scan[, c("RID", "mri_date", "landmark_date", "pathology_date", "t", "diagnosis", "assay_batch",
    "thickness_z", "gAD", "gMAV", "abeta42", "ptau", "amyloid_z", "ptau_z"), with = FALSE], "08_analysis_measurements.csv")
  save_csv(unique(scan[, .(RID, diagnosis, assay_batch, abeta42, ptau, amyloid_z, ptau_z)])[, .(
    people = .N, amyloid_ptau_correlation = cor(amyloid_z, ptau_z),
    min_abeta42 = min(abeta42), max_abeta42 = max(abeta42), min_ptau = min(ptau), max_ptau = max(ptau)),
    by = .(diagnosis, assay_batch)], "09_pathology_distributions.csv")

  cat("13H63 3/4: locked clinical-stage and pathology-conditioned models\n")
  covars <- c("age_landmark_z", "sex", "education_z", "scanner", paste0(pc, "_z"))
  primary_name <- names(sample_definitions)[1L]
  primary_result <- fit_set(sample_definitions[[primary_name]], covars, primary_name)
  results <- list(); results[[primary_name]] <- primary_result
  issues <- list()
  for (label in names(sample_definitions)[-1L]) {
    try_result <- tryCatch(fit_set(sample_definitions[[label]], covars, label),
      error = function(e) {issues[[label]] <<- data.table(sample = label, issue = conditionMessage(e)); NULL})
    if (!is.null(try_result)) results[[label]] <- try_result
  }
  for (field in c("diagnostics", "moderators", "joint_test", "stage_contrasts", "stage_tests", "lrt"))
    save_csv(rbindlist(lapply(results, `[[`, field), fill = TRUE), paste0("10_", field, ".csv"))
  if (length(issues)) save_csv(rbindlist(issues), "11_sensitivity_model_issues.csv")
  coefficients <- rbindlist(lapply(names(results), function(label) {
    rbindlist(lapply(names(results[[label]]$fits), function(m) {
      f <- results[[label]]$fits[[m]]$fit
      data.table(sample = label, model = m, term = names(fixef(f)), estimate = as.numeric(fixef(f)), se = sqrt(diag(vcov(f))))
    }))
  }))
  save_csv(coefficients, "12_fixed_coefficients.csv")
  specs <- rbindlist(lapply(names(results), function(label) data.table(sample = label,
    model = names(results[[label]]$rhs), fixed_rhs = unlist(results[[label]]$rhs), random_structure = results[[label]]$random)))
  save_csv(specs, "13_model_specifications.csv")
  saveRDS(results, file.path(out, "14_fitted_model_sets.rds"))

  cat("13H63 4/4: participant bootstrap of the primary direct modifiers\n")
  boot <- participant_bootstrap(primary_result, bootstrap_n, bootstrap_seed)
  save_csv(boot$values, "20_primary_participant_bootstrap_estimates.csv")
  save_csv(boot$log, "21_primary_participant_bootstrap_status.csv")
  main <- merge(primary_result$moderators, boot$intervals, by = "moderator", all.x = TRUE, sort = FALSE)
  save_csv(main, "22_primary_pathology_modifiers.csv")
  save_csv(primary_result$joint_test, "23_primary_two_df_joint_test.csv")

  writeLines(c(
    "The outcome is the original six-region thickness z-score; negative genetic contrasts indicate greater thinning per original gAD SD than per original gMAV SD.",
    "A negative amyloid modifier means the contrast becomes more negative as CSF Abeta42 decreases. A negative pTau modifier means it becomes more negative as pTau increases.",
    "The direct two-df contrast test answers whether the two genetic time effects differ in their pathology dependence. The four-df likelihood-ratio test answers a broader question about either genetic score.",
    "The model includes additive amyloid and pTau moderators.",
    "Clinical-only and joint-model stage contrasts use identical rows. Joint contrasts set both pathology variables to the primary sample mean.",
    "Random intercept and slope are tried first. If any of the five nested models is singular or fails, all five are refitted with a common random intercept. The common choice is frozen during participant bootstrap.",
    "Primary inference uses the two-df joint Wald test, two BH-adjusted moderator tests, and 500 participant-resampling intervals for those modifiers. Bootstrap scaling is frozen to the original primary reference population.",
    "At least 80 percent of bootstrap fits must succeed for percentile intervals. Resamples can drop absent nuisance-factor coefficients, but must estimate every target genetic-time-pathology term. Dropped terms, failures and singular fits are recorded.",
    "CSF pTau is analyzed as a continuous soluble tau marker."
  ), file.path(out, "90_methods_and_interpretation.txt"))
  save_csv(data.table(module = c("primary_models", "primary_bootstrap", "sensitivities"),
    status = c("completed", if (all(main$interval_status == "available")) "completed" else "review_bootstrap_failures",
      if (length(issues)) "some_unavailable_see_issues" else "completed")), "98_module_status.csv")
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  writeLines(c("13H63_COMPLETE", paste0("Results: ", out),
    paste0("Primary participants: ", uniqueN(primary_result$data$RID)),
    paste0("Primary MRI rows: ", nrow(primary_result$data))), file.path(out, "13H63_COMPLETED.txt"))
  cat("\n13H63 complete. Results: ", out, "\n", sep = "")
})
