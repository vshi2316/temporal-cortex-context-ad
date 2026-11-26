# Continuous pathology models and donor eligibility
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

# Read the donor scores and linked pathology metadata.
local({
  for (pkg in c("data.table", "sandwich"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)

  prior <- study_path("project", "results/13H34_GSE297004_MTG_neuron_bridge_release")
  metadata_file <- study_path("mtg_metadata")
  out <- file.path(study_path("project", "results"),
    paste0("13H35_GSE297004_pathology_selection_", "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  status_file <- file.path(out, "00_status.txt")
  status <- function(value, note = NULL) {
    x <- c(paste0("STATUS=", value), if (!is.null(note)) paste0("NOTE=", note),
      paste0("OUTPUT_DIR=", out))
    writeLines(x, status_file)
    cat(paste(x, collapse = "\n"), "\n")
  }
  put <- function(x, filename) data.table::fwrite(x, file.path(out, filename), na = "NA", bom = TRUE)
  if (!"STATUS=GSE297004_FROZEN_NEURON_BRIDGE_COMPLETE" %in%
      readLines(file.path(prior, "00_status.txt"), warn = FALSE)) stop("13H34 did not complete")
  status("RUNNING")
  donor <- data.table::fread(file.path(prior, "05_donor_metadata_audit.csv"))
  captured <- data.table::fread(file.path(prior, "04_donor_subtype_capture_audit.csv"))
  paired <- data.table::fread(file.path(prior, "11_paired_donor_contrasts.csv"))
  targets <- c("RORB_LCN15_vulnerable", "RORB_LINC01202_resilient")
  total <- captured[, .(total_neurons = sum(n_nuclei)), by = .(donor, diagnosis)]
  target <- captured[target %in% targets, .(n_nuclei = sum(n_nuclei)), by = .(donor, diagnosis, target)]
  wide <- data.table::dcast(target, donor + diagnosis ~ target, value.var = "n_nuclei", fill = 0)
  selection <- merge(donor, total, by = c("donor", "diagnosis"))
  selection <- merge(selection, wide, by = c("donor", "diagnosis"), all.x = TRUE)
  for (name in targets) {
    if (!name %in% names(selection)) selection[, (name) := 0L]
    data.table::set(selection, which(is.na(selection[[name]])), name, 0L)
  }
  selection[, included_13H34 := donor %in% paired$donor]
  for (cut in c(5L, 10L, 20L))
    selection[, (paste0("both_at_least_", cut)) :=
      get(targets[1]) >= cut & get(targets[2]) >= cut]
  put(selection, "01_all_donor_selection.csv")
  flow <- selection[, .(donors = .N,
    both_5 = sum(both_at_least_5), both_10 = sum(both_at_least_10),
    both_20 = sum(both_at_least_20), median_total_neurons = as.numeric(median(total_neurons)),
    female = sum(toupper(sex) == "F", na.rm = TRUE)),
    by = .(diagnosis, included_13H34)]
  put(flow, "02_selection_summary.csv")

  meta <- readRDS(metadata_file)
  if (!is.data.frame(meta)) meta <- as.data.frame(meta)
  required <- c("individual", "amyloid_beta", "pTau", "braak")
  if (!all(required %in% names(meta))) stop("Author pathology metadata fields have changed")
  metadata <- data.table::as.data.table(meta[, required])
  unique_values <- metadata[, lapply(.SD, function(x) data.table::uniqueN(x[!is.na(x)])),
    by = individual, .SDcols = required[-1]]
  put(unique_values, "03_pathology_values_per_donor.csv")
  path_cols <- required[-1]
  if (any(as.matrix(unique_values[, ..path_cols]) > 1L)) {
    status("PATHOLOGY_UNIT_REVIEW_REQUIRED", "A pathology value varies within donor; no donor-level model was fitted")
    stop("Inspect 03_pathology_values_per_donor.csv and the author definitions. Output: ", out)
  }
  pathology <- unique(metadata, by = "individual")
  data.table::setnames(pathology, "individual", "donor")
  put(pathology, "04_author_donor_pathology.csv")
  analysis <- merge(paired, pathology, by = "donor")
  if (nrow(analysis) != nrow(paired) || anyDuplicated(analysis$donor))
    stop("Pathology-to-paired-donor linkage is incomplete")
  analysis[, AD01 := as.integer(diagnosis == "AD")]
  analysis[, TREM2_risk := as.integer(trem2 %in% c("R47H", "R62H"))]
  put(analysis, "05_paired_donors_with_pathology.csv")

  fit_one <- function(marker, adjusted) {
    d <- as.data.frame(analysis)
    d$marker <- suppressWarnings(as.numeric(d[[marker]]))
    d <- d[complete.cases(d[, c("paired_direct", "marker", "age", "sex", "TREM2_risk", "AD01")]), ]
    if (nrow(d) < 25L || sd(d$marker) == 0) return(NULL)
    d$marker_z <- as.numeric(scale(d$marker))
    formula <- if (adjusted)
      paired_direct ~ marker_z + AD01 + age + sex + TREM2_risk else
      paired_direct ~ marker_z
    model <- tryCatch(lm(formula, data = d), error = function(e) NULL)
    if (is.null(model) || model$rank != length(coef(model)) || df.residual(model) < 12) return(NULL)
    se <- sqrt(sandwich::vcovHC(model, type = "HC3")["marker_z", "marker_z"])
    estimate <- unname(coef(model)["marker_z"])
    if (!is.finite(estimate) || !is.finite(se) || se <= 0) return(NULL)
    data.frame(marker = marker, model = if (adjusted) "diagnosis_and_covariates" else "unadjusted",
      n_donors = nobs(model), estimate_per_SD = estimate, se_HC3 = se,
      ci_low = estimate - qt(.975, df.residual(model)) * se,
      ci_high = estimate + qt(.975, df.residual(model)) * se,
      p = 2 * pt(-abs(estimate / se), df.residual(model)))
  }
  tests <- data.table::rbindlist(list(
    fit_one("amyloid_beta", FALSE), fit_one("amyloid_beta", TRUE),
    fit_one("pTau", FALSE), fit_one("pTau", TRUE)), fill = TRUE)
  if (nrow(tests)) {
    tests[, q_within_model := p.adjust(p, "BH"), by = model]
    put(tests, "06_continuous_pathology_models.csv")
  }

  writeLines(c(
    "Secondary donor-level pathology and selection models use the fixed donor program score.",
    "Captured nucleus counts summarize the sequenced sample.",
    "The adjusted marker coefficient conditions on diagnosis, age, sex and TREM2; two markers are BH-adjusted within each model."),
    file.path(out, "07_interpretation.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  status("GSE297004_PATHOLOGY_SELECTION_AUDIT_COMPLETE",
    paste0("PAIRED_DONORS=", nrow(analysis), "; PATHOLOGY_MODELS=", nrow(tests)))
})
