# Raw T1 acquisition manifest
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(data.table))

project_dir <- study_path("project")
ad_root <- study_path("cohorts")
out_dir <- file.path(project_dir, "results", "13B_longitudinal_MRI_acquisition_cohort_freeze")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  message(sprintf("[%s] [INFO] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}
need <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  path
}
write_csv <- function(x, name) {
  x <- as.data.table(x)
  fwrite(x, file.path(out_dir, name), bom = TRUE)
  log_msg("Wrote ", name, " (", nrow(x), " rows)")
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
date0 <- function(x) suppressWarnings(as.IDate(as.character(x)))
norm_id <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("^0+([0-9]+)$", "\\1", y)
  y[y %in% c("", "NA", "NaN", ".", "-4")] <- NA_character_
  y
}
is_true <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}
present <- function(x) !is.na(x) & trimws(as.character(x)) != "" &
  !trimws(as.character(x)) %in% c("NA", "NaN", ".", "-4")
write_batches <- function(ids, prefix, batch_size = 100L) {
  ids <- unique(as.character(ids[present(ids)]))
  if (!length(ids)) return(0L)
  groups <- split(ids, ceiling(seq_along(ids) / batch_size))
  for (i in seq_along(groups)) {
    writeLines(paste(groups[[i]], collapse = ","),
               file.path(out_dir, sprintf("%s_batch_%02d.txt", prefix, i)),
               useBytes = TRUE)
  }
  length(groups)
}

log_msg("Step 13B started: freeze image-acquisition cohorts and temporal eligibility only.")
log_msg("Build the acquisition manifest from the specified cohort tables.")

# -----------------------------------------------------------------------------
# 1. Frozen acquisition and temporal rules.
# -----------------------------------------------------------------------------
rules <- data.table(
  rule = c(
    "ADNI genetic eligibility", "ADNI MRI change window", "ADNI MRI spacing",
    "ADNI future clinical follow-up", "ADNI scan selection", "ADNI download identifier rule",
    "HABS-HD role", "AIBL MRI rule", "AIBL follow-up rule", "external validation rule",
    "pilot rule", "outcome-blinding rule"
  ),
  frozen_value = c(
    "Primary European-ancestry indicator plus finite gMAV-noAPOE, gAD-noAPOE and PC1-PC10",
    "MRI records within 120 days are treated as one acquisition visit; the first visit and the earliest later visit >=270 days afterward define the change window; later visits are retained for sensitivity analysis",
    "At least 270 days between the representative dates of the first and second MRI visit clusters",
    "At least one valid cognition or diagnosis assessment >=180 days after the second MRI",
    "One FreeSurfer-linked scan per person/date; prefer 3T, non-failing QC and the highest image identifier only to identify the source visit",
    "Locate original T1 IDs by matching PTID and acquisition metadata to the FreeSurfer source visit.",
    "Core analyses use regional MRI derivatives. CAT12 supplies the supplementary volumetric analysis.",
    "One baseline T1 MPRAGE per participant, selected without using future outcomes",
    "At least one non-baseline diagnosis, MMSE or CDR assessment",
    "Use the fixed anatomical definitions in each cohort; HABS-HD prognostic models are fitted within that cohort",
    "Process a 50-person ADNI plus 50-person AIBL pilot before the full image set",
    "No association with gMAV, gAD, cognition or progression is estimated before image QC and feature extraction are frozen"
  )
)
write_csv(rules, "00_frozen_acquisition_and_temporal_rules.csv")

# -----------------------------------------------------------------------------
# 2. ADNI: identify a true longitudinal discovery set before raw-image download.
# -----------------------------------------------------------------------------
adni_fs_file <- need(file.path(ad_root, "ADNI", "data", "MRI", "MRI_FreeSurfer.csv"))
adni_pgs_file <- need(file.path(project_dir, "results", "10A_ADNI_gMAV_PGS_assembly",
                                "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"))

pgs <- fread(adni_pgs_file, showProgress = FALSE)
pgs[, RID := num(RID)]
pgs[, clinical_date := date0(EXAMDATE)]
pc_cols <- paste0("PC", 1:10)
required_pgs <- c("RID", "primary_european_ancestry", "gMAV_noAPOE_z_eur",
                  "gAD_noAPOE_z_eur", pc_cols)
if (!all(required_pgs %in% names(pgs))) {
  stop("ADNI PGS file lacks: ", paste(setdiff(required_pgs, names(pgs)), collapse = ", "), call. = FALSE)
}

pgs[, genetic_complete := is_true(primary_european_ancestry) &
      is.finite(num(gMAV_noAPOE_z_eur)) & is.finite(num(gAD_noAPOE_z_eur))]
for (v in pc_cols) pgs[, genetic_complete := genetic_complete & is.finite(num(get(v)))]
eligible_rids <- unique(pgs[genetic_complete == TRUE & is.finite(RID), RID])

fs <- fread(adni_fs_file, showProgress = FALSE)
fs[, RID := num(RID)]
fs[, mri_date := date0(EXAMDATE)]
fs[, imageuid_num := num(IMAGEUID)]
qc_text <- tolower(trimws(as.character(fs$OVERALLQC)))
fs[, qc_explicit_fail := grepl("fail|exclude|poor|unusable", qc_text)]
fs[, qc_positive := grepl("pass|good|acceptable", qc_text)]
fs[, is_3T := grepl("3", as.character(FIELD_STRENGTH))]
fs <- fs[RID %in% eligible_rids & is.finite(RID) & !is.na(mri_date) &
           present(IMAGEUID) & qc_explicit_fail == FALSE]

# ADNI FreeSurfer tables may contain several scans or processing records from
# the same clinical acquisition period on dates separated by days or weeks.
# Collapse records into visit clusters before defining longitudinal change.
setorder(fs, RID, mri_date, -is_3T, -qc_positive, -imageuid_num, na.last = TRUE)
fs[, visit_cluster := {
  d <- as.integer(mri_date)
  cluster <- integer(.N)
  current <- 1L
  anchor <- d[1]
  cluster[1] <- current
  if (.N >= 2L) {
    for (j in 2:.N) {
      if (is.finite(d[j]) && is.finite(anchor) && (d[j] - anchor) > 120L) {
        current <- current + 1L
        anchor <- d[j]
      }
      cluster[j] <- current
    }
  }
  cluster
}, by = RID]

# Select one representative processed result per visit cluster. The selected
# IMAGEUID identifies the processed source visit.
setorder(fs, RID, visit_cluster, -is_3T, -qc_positive, mri_date, -imageuid_num, na.last = TRUE)
fs_visit <- fs[, .SD[1], by = .(RID, visit_cluster)]
setorder(fs_visit, RID, mri_date)

adni_span <- fs_visit[, {
  dates <- sort(unique(mri_date))
  first_date <- dates[1]
  eligible_second <- dates[as.integer(dates - first_date) >= 270L]
  if (!length(eligible_second)) {
    list(n_distinct_MRI_dates = length(dates), first_MRI_date = first_date,
         second_MRI_date = as.IDate(NA), last_MRI_date = dates[length(dates)],
         first_to_second_days = NA_integer_, total_MRI_span_days = 0L)
  } else {
    second_date <- eligible_second[1]
    list(n_distinct_MRI_dates = length(dates), first_MRI_date = dates[1],
         second_MRI_date = second_date, last_MRI_date = dates[length(dates)],
         first_to_second_days = as.integer(second_date - dates[1]),
         total_MRI_span_days = as.integer(dates[length(dates)] - dates[1]))
  }
}, by = RID]

endpoint_cols <- intersect(c("MMSE", "ADAS13", "CDRSB", "FAQ", "MOCA",
                             "mPACCdigit", "mPACCtrailsB", "DX"), names(pgs))
if (!length(endpoint_cols)) stop("No ADNI future clinical endpoint columns found.", call. = FALSE)
pgs[, valid_future_endpoint := FALSE]
for (v in endpoint_cols) pgs[, valid_future_endpoint := valid_future_endpoint | present(get(v))]

future <- merge(
  pgs[genetic_complete == TRUE & !is.na(clinical_date) & valid_future_endpoint == TRUE,
      c("RID", "clinical_date", endpoint_cols), with = FALSE],
  adni_span[!is.na(second_MRI_date)], by = "RID", allow.cartesian = TRUE
)
future[, days_after_second_MRI := as.integer(clinical_date - second_MRI_date)]
future_summary <- future[days_after_second_MRI >= 180L,
                         .(n_future_assessments = .N,
                           first_future_date = min(clinical_date),
                           last_future_date = max(clinical_date),
                           future_followup_days = max(days_after_second_MRI)), by = RID]

adni_people <- merge(adni_span, future_summary, by = "RID", all.x = TRUE)
adni_people[, acquisition_eligible := n_distinct_MRI_dates >= 2L &
              is.finite(first_to_second_days) & first_to_second_days >= 270L &
              is.finite(n_future_assessments) & n_future_assessments >= 1L]
write_csv(adni_people, "01_ADNI_longitudinal_eligibility_flow_by_person.csv")

adni_eligible <- adni_people[acquisition_eligible == TRUE]
adni_targets <- merge(fs_visit, adni_eligible[, .(RID, first_MRI_date, second_MRI_date,
                                                  n_distinct_MRI_dates, last_future_date)],
                      by = "RID")
adni_targets <- adni_targets[mri_date == first_MRI_date | mri_date == second_MRI_date]
adni_targets[, window_position := fifelse(mri_date == first_MRI_date, "MRI1", "MRI2")]
adni_targets <- adni_targets[, .(
  RID, PTID, PHASE, VISCODE, VISCODE2, window_position, mri_date,
  processed_FreeSurfer_IMAGEUID = IMAGEUID, FIELD_STRENGTH, OVERALLQC, FSVER,
  n_distinct_MRI_dates, last_future_date,
  repository_action = "Locate and download the original T1 for this PTID/phase/visit/date; do not assume the processed IMAGEUID is the raw T1 ID"
)]
setorder(adni_targets, RID, mri_date)
write_csv(adni_targets, "02_ADNI_original_T1_visit_lookup_targets.csv")
write_batches(unique(adni_targets$PTID), "02B_ADNI_subject_inventory_search", 80L)

adni_all_visits <- merge(fs_visit, adni_eligible[, .(RID)], by = "RID")
adni_all_visits <- adni_all_visits[, .(RID, PTID, PHASE, VISCODE, VISCODE2, mri_date,
                                       processed_FreeSurfer_IMAGEUID = IMAGEUID,
                                       FIELD_STRENGTH, OVERALLQC, FSVER)]
write_csv(adni_all_visits, "02C_ADNI_all_longitudinal_visits_sensitivity_lookup.csv")

# Deterministic pilot: broad phase/scanner representation without outcome use.
pilot_people <- unique(adni_targets[, .(RID, PTID, PHASE, FIELD_STRENGTH)])
setorder(pilot_people, PHASE, FIELD_STRENGTH, RID)
pilot_people[, stratum_index := seq_len(.N), by = .(PHASE, FIELD_STRENGTH)]
pilot_people[, stratum_n := .N, by = .(PHASE, FIELD_STRENGTH)]
pilot_people[, pilot_rank := (stratum_index - 0.5) / stratum_n]
setorder(pilot_people, pilot_rank, RID)
adni_pilot_ids <- head(unique(pilot_people$RID), 50L)
adni_pilot <- adni_targets[RID %in% adni_pilot_ids]
write_csv(adni_pilot, "03_ADNI_T1_pilot50_visit_targets.csv")

# -----------------------------------------------------------------------------
# 3. AIBL baseline T1 and subsequent clinical follow-up.
# -----------------------------------------------------------------------------
aibl_root <- study_path("aibl_data")
aibl_manifest_file <- need(study_path("aibl_mri_manifest"))
aibl_dx_file <- need(file.path(aibl_root, "aibl_pdxconv_01-Jun-2018.csv"))
aibl_mmse_file <- need(file.path(aibl_root, "aibl_mmse_01-Jun-2018.csv"))
aibl_cdr_file <- need(file.path(aibl_root, "aibl_cdr_01-Jun-2018.csv"))
aibl_apoe_file <- need(file.path(aibl_root, "aibl_apoeres_01-Jun-2018.csv"))

aibl_img <- fread(aibl_manifest_file, showProgress = FALSE)
aibl_img[, RID_key := norm_id(`Subject ID`)]
aibl_img[, image_id_num := num(`Image ID`)]
aibl_img[, is_3T := grepl("3\\.0|Field Strength=3", as.character(`Imaging Protocol`))]
aibl_img[, is_MPRAGE := grepl("MPRAGE", paste(Description, `Imaging Protocol`), ignore.case = TRUE)]
aibl_img <- aibl_img[present(RID_key) & present(`Image ID`) & is_MPRAGE == TRUE]
setorder(aibl_img, RID_key, -is_3T, image_id_num)
aibl_img <- aibl_img[, .SD[1], by = RID_key]

dx <- fread(aibl_dx_file, showProgress = FALSE); dx[, RID_key := norm_id(RID)]
mmse <- fread(aibl_mmse_file, showProgress = FALSE); mmse[, RID_key := norm_id(RID)]
cdr <- fread(aibl_cdr_file, showProgress = FALSE); cdr[, RID_key := norm_id(RID)]
apoe <- fread(aibl_apoe_file, showProgress = FALSE); apoe[, RID_key := norm_id(RID)]

is_nonbaseline <- function(x) !tolower(trimws(as.character(x))) %in% c("", "bl", "baseline", "sc", "screening")
dx_future <- dx[is_nonbaseline(VISCODE) &
                  (present(DXCURREN) | present(DXNORM) | present(DXMCI) | present(DXAD)),
                .(has_future_diagnosis = TRUE, n_future_diagnosis = .N), by = RID_key]
mmse_future <- mmse[is_nonbaseline(VISCODE) & is.finite(num(MMSCORE)),
                    .(has_future_MMSE = TRUE, n_future_MMSE = .N), by = RID_key]
cdr_future <- cdr[is_nonbaseline(VISCODE) & is.finite(num(CDGLOBAL)),
                  .(has_future_CDR = TRUE, n_future_CDR = .N), by = RID_key]

aibl_follow <- Reduce(function(x, y) merge(x, y, by = "RID_key", all = TRUE),
                      list(dx_future, mmse_future, cdr_future))
for (v in c("has_future_diagnosis", "has_future_MMSE", "has_future_CDR")) {
  if (!v %in% names(aibl_follow)) aibl_follow[, (v) := FALSE]
  aibl_follow[is.na(get(v)), (v) := FALSE]
}
for (v in c("n_future_diagnosis", "n_future_MMSE", "n_future_CDR")) {
  if (!v %in% names(aibl_follow)) aibl_follow[, (v) := 0L]
  aibl_follow[is.na(get(v)), (v) := 0L]
}
aibl_follow[, has_any_future_clinical := has_future_diagnosis | has_future_MMSE | has_future_CDR]

apoe_base <- apoe[, .SD[1], by = RID_key]
apoe_base[, APOE4_carrier := num(APGEN1) == 4 | num(APGEN2) == 4]
aibl_targets <- merge(aibl_img, aibl_follow[has_any_future_clinical == TRUE], by = "RID_key")
aibl_targets <- merge(aibl_targets, apoe_base[, .(RID_key, APGEN1, APGEN2, APOE4_carrier)],
                      by = "RID_key", all.x = TRUE)
aibl_targets <- aibl_targets[, .(
  RID = RID_key, Project, Sex, `Research Group`, Visit, `Study Date`, Age,
  Image_ID = `Image ID`, Description, Imaging_Protocol = `Imaging Protocol`,
  has_future_diagnosis, n_future_diagnosis, has_future_MMSE, n_future_MMSE,
  has_future_CDR, n_future_CDR, APGEN1, APGEN2, APOE4_carrier
)]
setorder(aibl_targets, RID)
write_csv(aibl_targets, "04_AIBL_locked_baseline_T1_validation_download_manifest.csv")
write_batches(aibl_targets$Image_ID, "04B_AIBL_T1_image_ids", 200L)

# Deterministic AIBL pilot, selected without future values or outcomes.
aibl_pilot_people <- copy(aibl_targets)
setorder(aibl_pilot_people, Sex, Age, RID)
aibl_pilot_people[, stratum_index := seq_len(.N), by = Sex]
aibl_pilot_people[, stratum_n := .N, by = Sex]
aibl_pilot_people[, pilot_rank := (stratum_index - 0.5) / stratum_n]
setorder(aibl_pilot_people, pilot_rank, RID)
aibl_pilot <- head(aibl_pilot_people, 50L)
write_csv(aibl_pilot, "05_AIBL_T1_pilot50_download_manifest.csv")
write_batches(aibl_pilot$Image_ID, "05B_AIBL_pilot50_image_ids", 100L)

# -----------------------------------------------------------------------------
# 4. HABS-HD: document the current derived replication and optional raw pilot.
# -----------------------------------------------------------------------------
habs_manifest_file <- need(file.path(ad_root, "LBA_MRI", "download_manifests",
                                     "HABS_HD_T1_download_manifest.csv"))
habs_pilot_file <- need(file.path(ad_root, "LBA_MRI", "download_manifests",
                                  "HABS_HD_T1_pilot100.csv"))
habs_manifest <- fread(habs_manifest_file, showProgress = FALSE)
habs_pilot <- fread(habs_pilot_file, showProgress = FALSE)
write_csv(habs_pilot, "06_HABSHD_optional_voxel_replication_pilot100.csv")
write_batches(habs_pilot$Image_ID, "06B_HABSHD_pilot100_image_ids", 100L)

# -----------------------------------------------------------------------------
# 5. Acquisition sample counts.
# -----------------------------------------------------------------------------
adni_n <- uniqueN(adni_targets$RID)
adni_n3 <- adni_people[acquisition_eligible == TRUE & n_distinct_MRI_dates >= 3L, .N]
aibl_n <- uniqueN(aibl_targets$RID)
flow <- data.table(
  cohort = c("ADNI", "ADNI", "ADNI", "AIBL", "AIBL", "HABS-HD", "HABS-HD"),
  stage = c(
    "genetically eligible people", "longitudinal MRI plus future outcome eligible people",
    "eligible people with >=3 MRI dates", "baseline T1 manifest people",
    "baseline T1 plus future clinical follow-up", "existing baseline raw-image manifest",
    "optional voxel pilot"
  ),
  n = c(length(eligible_rids), adni_n, adni_n3, uniqueN(aibl_img$RID_key), aibl_n,
        uniqueN(norm_id(habs_manifest$Med_ID)), uniqueN(norm_id(habs_pilot$Med_ID)))
)
write_csv(flow, "07_acquisition_sample_flow.csv")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
log_msg("Step 13B completed: ADNI n=", adni_n, "; AIBL n=", aibl_n, ".")
