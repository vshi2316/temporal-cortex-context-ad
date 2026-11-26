# ADNI repeated-MRI eligibility and derivative frame
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)

# Assemble repeated ADNI MRI derivatives and linked HABS-HD features.
suppressPackageStartupMessages(library(data.table))

pick_root <- function(ascii, native) if (file.exists(ascii)) ascii else native
project_dir <- pick_root(study_path("project"), study_path("project"))
adni_root <- pick_root(study_path("cohorts"), study_path("cohorts"))
out_dir <- file.path(project_dir, "results", "13F_R_only_structural_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  path
}
read_csv <- function(path) fread(file = need(path), showProgress = FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
is_true <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
finite <- function(x) is.finite(num(x))

linked_path <- file.path(project_dir, "results", "12C_ADNI_common_structural_indices_multimodal_linkage", "03_ADNI_MRI_PGS_ATVN_linked.rds")
habs_path <- file.path(project_dir, "results", "12E_HABSHD_external_structural_MRI_validation", "06_HABSHD_endpoint_specific_validation_frames.rds")
linked <- as.data.table(readRDS(need(linked_path)))
habs_frames <- readRDS(need(habs_path))

# -----------------------------------------------------------------------------
# 1. ADNI repeated regional MRI frame.
# -----------------------------------------------------------------------------
d <- copy(linked)
d[, RID := num(RID)]
d[, mri_date := as.IDate(as.character(mri_date))]
d[, qc_missing := is.na(qc_missing) | is_true(qc_missing)]
d[, qc_explicit_fail := grepl("fail|exclude|poor|unusable", tolower(trimws(as.character(OVERALLQC))))]
d[, pgs_complete := finite(clinical_gMAV_z_eur) & finite(clinical_gAD_noAPOE_z_eur) & is_true(clinical_primary_european_ancestry)]
d[, structural_core_complete := finite(global_cortical_thickness) & finite(ad_vulnerable_index) &
      finite(frontoposterior_gradient) & finite(spatial_dispersion) & finite(mean_absolute_asymmetry)]
d[, hippocampus_complete := finite(hip_hippocampus_icv_ratio)]
d[, wmh_complete := finite(wmh_wmh_brain_ratio)]
d[, structural_extended_complete := structural_core_complete & hippocampus_complete & wmh_complete]
d <- d[pgs_complete & !is.na(RID) & !is.na(mri_date)]

# Time starts at the first eligible MRI date. Retain all repeated scans.
d[, first_mri_date := min(mri_date), by = RID]
d[, time_years := as.numeric(mri_date - first_mri_date) / 365.25]
d[, age_mri := num(clinical_age_bl) + time_years]
d[, scanner_factor := factor(paste(FIELD_STRENGTH, FSVER, sep = "_"))]

keep <- intersect(c(
  "RID", "PTID", "VISCODE", "VISCODE2", "IMAGEUID", "mri_date", "first_mri_date",
  "time_years", "age_mri", "FIELD_STRENGTH", "FSVER", "OVERALLQC", "qc_missing",
  "qc_explicit_fail", "structural_core_complete", "structural_extended_complete",
  "global_cortical_thickness", "ad_vulnerable_index", "frontal_index", "posterior_index",
  "frontoposterior_gradient", "spatial_dispersion", "mean_absolute_asymmetry",
  "hip_hippocampus_icv_ratio", "wmh_wmh_brain_ratio", "clinical_gMAV_z_eur",
  "clinical_gMAV_noAPOE_z_eur", "clinical_gAD_noAPOE_z_eur", "clinical_APOE4",
  paste0("clinical_PC", 1:10), "clinical_MMSE", "clinical_ADAS13", "clinical_CDRSB",
  "clinical_FAQ", "clinical_MOCA", "clinical_mPACCdigit", "clinical_mPACCtrailsB",
  "clinical_DX", "amyloid_CENTILOIDS", "tau_META_TEMPORAL_SUVR"
), names(d))
adni <- d[, ..keep]
setorder(adni, RID, mri_date, IMAGEUID)

# Explicit analysis flags. The lenient set retains missing QC labels but
# excludes explicit failures; the strict set additionally requires a nonmissing
# QC label. Both sets require the pre-existing structural core to be finite.
adni[, analysis_lenient := !qc_explicit_fail & structural_core_complete]
adni[, analysis_strict := analysis_lenient & !qc_missing]

# -----------------------------------------------------------------------------
# 2. HABS-HD frame is retained as a separate external layer; no refitting or
#    region re-selection is performed here.
# -----------------------------------------------------------------------------
habs <- if (is.data.frame(habs_frames$index_memory)) as.data.table(habs_frames$index_memory) else data.table()
if (nrow(habs)) {
  habs[, external_role := "HABS-HD_frozen_external_structural_frame"]
}

# -----------------------------------------------------------------------------
# 3. Frozen raw-image manifest audit, kept separate from derivative analysis.
# -----------------------------------------------------------------------------
person_flow <- adni[, .(
  scans = .N,
  strict_scans = sum(analysis_strict),
  lenient_scans = sum(analysis_lenient),
  first_date = min(mri_date),
  last_date = max(mri_date),
  span_years = as.numeric(max(mri_date) - min(mri_date)) / 365.25
), by = RID]

sample_flow <- rbindlist(list(
  data.table(layer = "ADNI_derivative", set = "all_PGS_linked", rows = nrow(adni), people = uniqueN(adni$RID)),
  data.table(layer = "ADNI_derivative", set = "lenient_structural_core", rows = sum(adni$analysis_lenient), people = uniqueN(adni[analysis_lenient == TRUE, RID])),
  data.table(layer = "ADNI_derivative", set = "strict_structural_core", rows = sum(adni$analysis_strict), people = uniqueN(adni[analysis_strict == TRUE, RID])),
  data.table(layer = "ADNI_derivative", set = "lenient_people_with_2_scans", rows = NA_integer_, people = sum(person_flow$lenient_scans >= 2L)),
  data.table(layer = "ADNI_derivative", set = "strict_people_with_2_scans", rows = NA_integer_, people = sum(person_flow$strict_scans >= 2L)),
  data.table(layer = "ADNI_derivative", set = "lenient_people_with_3_scans", rows = NA_integer_, people = sum(person_flow$lenient_scans >= 3L)),
  data.table(layer = "ADNI_derivative", set = "strict_people_with_3_scans", rows = NA_integer_, people = sum(person_flow$strict_scans >= 3L)),
  data.table(layer = "HABS-HD", set = "frozen_index_memory_frame", rows = nrow(habs), people = uniqueN(habs$participant_id))
), fill = TRUE)

qc_sensitivity <- rbindlist(list(
  adni[, .(set = "all_PGS_linked", rows = .N, people = uniqueN(RID), missing_qc = sum(qc_missing), explicit_fail = sum(qc_explicit_fail), core_complete = sum(structural_core_complete))],
  adni[analysis_lenient == TRUE, .(set = "lenient_structural_core", rows = .N, people = uniqueN(RID), missing_qc = sum(qc_missing), explicit_fail = sum(qc_explicit_fail), core_complete = sum(structural_core_complete))],
  adni[analysis_strict == TRUE, .(set = "strict_structural_core", rows = .N, people = uniqueN(RID), missing_qc = sum(qc_missing), explicit_fail = sum(qc_explicit_fail), core_complete = sum(structural_core_complete))]
), fill = TRUE)

fwrite(sample_flow, file.path(out_dir, "02_derivative_sample_flow.csv"), bom = TRUE)
fwrite(qc_sensitivity, file.path(out_dir, "03_qc_sensitivity.csv"), bom = TRUE)
saveRDS(adni, file.path(out_dir, "04_ADNI_R_only_longitudinal_structural_frame.rds"), compress = "xz")
saveRDS(habs, file.path(out_dir, "05_HABSHD_R_only_external_structural_frame.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
message("R-only structural preparation completed: ", out_dir)
