# DICOM series identity and visit alignment
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

project_dir <- study_path("project")
adni_dir <- study_path("cohorts", "LBA_MRI/ADNI")
aibl_dir <- study_path("cohorts", "LBA_MRI/AIBL")
out_dir <- file.path(project_dir, "results", "13C_ADNI_AIBL_ZIP_DICOM_identity_visit_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  message(sprintf("[%s] [INFO] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}
write_csv <- function(x, name) {
  x <- as.data.table(x)
  fwrite(x, file.path(out_dir, name), bom = TRUE)
  log_msg("Wrote ", name, " (", nrow(x), " rows)")
}
need <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  path
}
norm_id <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("^0+([0-9]+)$", "\\1", y)
  y[y %in% c("", "NA", "NaN", ".", "-4")] <- NA_character_
  y
}
date_mdy <- function(x) suppressWarnings(as.IDate(as.character(x), format = "%m/%d/%Y"))
date_iso <- function(x) suppressWarnings(as.IDate(as.character(x)))

log_msg("Step 13C started: ZIP-level integrity, DICOM-series identity, visit pairing and blinded QC preparation.")
log_msg("No ZIP is extracted and no association or outcome model is fitted.")

# -----------------------------------------------------------------------------
# 1. Frozen expected image manifests.
# -----------------------------------------------------------------------------
adni_primary_file <- need(file.path(
  project_dir, "results", "13B1_ADNI_original_T1_ImageID_matching",
  "06_download_ready_complete_pairs.csv"
))
adni_backup_file <- need(file.path(
  project_dir, "results", "13B1_ADNI_original_T1_ImageID_matching",
  "06B_QC_backup_T1_candidates_for_ambiguous_pairs.csv"
))
aibl_file <- need(file.path(
  project_dir, "results", "13B_longitudinal_MRI_acquisition_cohort_freeze",
  "04_AIBL_locked_baseline_T1_validation_download_manifest.csv"
))

adni_primary <- fread(adni_primary_file, showProgress = FALSE)
adni_backup <- fread(adni_backup_file, showProgress = FALSE)
aibl <- fread(aibl_file, showProgress = FALSE)

if (nrow(adni_primary) != 1278L || uniqueN(adni_primary$RID) != 639L) {
  stop("ADNI primary manifest is not the frozen 1,278 images / 639 people.", call. = FALSE)
}
if (uniqueN(adni_backup$Image_ID) != 89L) {
  stop("ADNI backup manifest does not contain 89 unique Image IDs.", call. = FALSE)
}
if (nrow(aibl) != 338L || uniqueN(aibl$Image_ID) != 338L) {
  stop("AIBL manifest is not the frozen 338-image set.", call. = FALSE)
}

expected_primary <- rbindlist(list(
  adni_primary[, .(
    cohort = "ADNI", role = "PRIMARY", image_id = as.character(Image_ID),
    participant_id = as.character(PTID), RID = as.character(RID),
    window_position = as.character(window_position),
    frozen_target_date = date_iso(target_date),
    expected_date = date_mdy(Study_Date), expected_phase = as.character(target_phase),
    expected_visit = as.character(target_visit), expected_description = as.character(Description),
    expected_field_strength = suppressWarnings(as.numeric(field_strength)),
    expected_slice_thickness = suppressWarnings(as.numeric(slice_thickness)),
    expected_sequence_family = as.character(sequence_family)
  )],
  aibl[, .(
    cohort = "AIBL", role = "PRIMARY", image_id = as.character(Image_ID),
    participant_id = norm_id(RID), RID = norm_id(RID),
    window_position = "BASELINE", frozen_target_date = date_mdy(`Study Date`),
    expected_date = date_mdy(`Study Date`),
    expected_phase = "AIBL", expected_visit = as.character(Visit),
    expected_description = as.character(Description),
    expected_field_strength = suppressWarnings(as.numeric(sub(
      ".*Field Strength=([^;]+).*", "\\1", as.character(Imaging_Protocol)
    ))),
    expected_slice_thickness = suppressWarnings(as.numeric(sub(
      ".*Slice Thickness=([^;]+).*", "\\1", as.character(Imaging_Protocol)
    ))),
    expected_sequence_family = "MPRAGE_STANDARD"
  )]
), fill = TRUE)

expected_backup <- unique(adni_backup[, .(
  cohort = "ADNI", role = "QC_BACKUP", image_id = as.character(Image_ID),
  participant_id = as.character(PTID), RID = as.character(RID),
  window_position = as.character(window_position),
  frozen_target_date = date_iso(target_date),
  expected_date = date_mdy(Study_Date), expected_phase = as.character(target_phase),
  expected_visit = as.character(target_visit), expected_description = as.character(Description),
  expected_field_strength = suppressWarnings(as.numeric(field_strength)),
  expected_slice_thickness = suppressWarnings(as.numeric(slice_thickness)),
  expected_sequence_family = as.character(sequence_family)
)], by = "image_id")

expected <- rbindlist(list(expected_primary, expected_backup), fill = TRUE)
if (anyDuplicated(expected[, .(cohort, image_id)])) {
  stop("An Image ID occurs in more than one frozen role.", call. = FALSE)
}
write_csv(expected, "00_frozen_expected_image_manifest.csv")

# -----------------------------------------------------------------------------
# 2. Read ZIP central directories directly. No extraction is performed.
# -----------------------------------------------------------------------------
zip_files <- rbindlist(list(
  data.table(cohort = "ADNI", zip_path = list.files(adni_dir, "\\.zip$", full.names = TRUE, ignore.case = TRUE)),
  data.table(cohort = "AIBL", zip_path = list.files(aibl_dir, "\\.zip$", full.names = TRUE, ignore.case = TRUE))
))
if (!nrow(zip_files)) stop("No ZIP files found.", call. = FALSE)
zip_files[, zip_name := basename(zip_path)]
zip_files[, is_metadata := grepl("metadata", zip_name, ignore.case = TRUE)]
zip_files[, archive_role := fifelse(
  cohort == "ADNI" & grepl("backup", zip_name, ignore.case = TRUE), "QC_BACKUP",
  fifelse(cohort == "ADNI" & grepl("956358|ADNI_T0", zip_name, ignore.case = TRUE), "PRIMARY",
  fifelse(cohort == "AIBL" & !is_metadata, "PRIMARY", "METADATA"))
)]

zip_list_safe <- function(path) {
  tryCatch({
    z <- as.data.table(utils::unzip(path, list = TRUE))
    if ("Name" %in% names(z)) setnames(z, "Name", "entry_name")
    if ("Length" %in% names(z)) setnames(z, "Length", "entry_bytes")
    if ("Date" %in% names(z)) setnames(z, "Date", "entry_date")
    # Normalize names explicitly because unzip() column names can vary by R build.
    if (!"entry_name" %in% names(z)) setnames(z, names(z)[1], "entry_name")
    if (!"entry_bytes" %in% names(z)) setnames(z, names(z)[2], "entry_bytes")
    z[, list_error := NA_character_]
    z
  }, error = function(e) {
    data.table(entry_name = NA_character_, entry_bytes = NA_real_,
               entry_date = as.POSIXct(NA), list_error = conditionMessage(e))
  })
}

entry_list <- rbindlist(lapply(seq_len(nrow(zip_files)), function(i) {
  z <- zip_list_safe(zip_files$zip_path[i])
  z[, `:=`(
    cohort = zip_files$cohort[i], zip_path = zip_files$zip_path[i],
    zip_name = zip_files$zip_name[i], is_metadata = zip_files$is_metadata[i],
    archive_role = zip_files$archive_role[i]
  )]
  z
}), fill = TRUE, use.names = TRUE)

parse_path <- function(x) {
  x <- gsub("\\\\", "/", as.character(x))
  data.table(
    image_id = ifelse(grepl("/I[0-9]+/", x), sub(".*/I([0-9]+)/.*", "\\1", x), NA_character_),
    path_cohort = ifelse(grepl("/", x), sub("/.*", "", x), NA_character_),
    path_participant = ifelse(grepl("^[^/]+/[^/]+/", x), sub("^[^/]+/([^/]+)/.*", "\\1", x), NA_character_),
    path_description = ifelse(grepl("^[^/]+/[^/]+/[^/]+/", x), sub("^[^/]+/[^/]+/([^/]+)/.*", "\\1", x), NA_character_),
    path_datetime = ifelse(grepl("^[^/]+/[^/]+/[^/]+/[^/]+/", x), sub("^[^/]+/[^/]+/[^/]+/([^/]+)/.*", "\\1", x), NA_character_),
    extension = tolower(sub(".*\\.", "", basename(x)))
  )
}
parsed <- parse_path(entry_list$entry_name)
entry_list <- cbind(entry_list, parsed)
entry_list[, path_date := date_iso(substr(path_datetime, 1, 10))]
entry_list[, is_dicom := extension %in% c("dcm", "ima")]
entry_list[, is_xml := extension == "xml"]

archive_audit <- entry_list[, .(
  central_directory_readable = all(is.na(list_error)),
  list_error = paste(unique(na.omit(list_error)), collapse = " | "),
  n_entries = .N,
  n_dicom_entries = sum(is_dicom, na.rm = TRUE),
  n_xml_entries = sum(is_xml, na.rm = TRUE),
  n_unique_image_ids = uniqueN(image_id[!is.na(image_id)]),
  uncompressed_GB = sum(as.numeric(entry_bytes), na.rm = TRUE) / 1024^3
), by = .(cohort, zip_name, zip_path, is_metadata, archive_role)]
archive_audit[, zip_bytes := file.info(zip_path)$size]
archive_audit[, zip_GB := zip_bytes / 1024^3]
archive_audit[, md5 := unname(tools::md5sum(zip_path))]
write_csv(archive_audit, "01_ZIP_integrity_and_MD5_audit.csv")

if (any(!archive_audit$central_directory_readable)) {
  stop("At least one ZIP central directory is unreadable. See 01_ZIP_integrity_and_MD5_audit.csv", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 3. Series-level inventory and duplicate-volume handling.
# -----------------------------------------------------------------------------
dicom_entries <- entry_list[is_metadata == FALSE & is_dicom == TRUE & !is.na(image_id)]
xml_entries <- entry_list[is_metadata == TRUE & is_xml == TRUE & !is.na(image_id)]

series_copy <- dicom_entries[, .(
  n_dicom = .N,
  total_dicom_bytes = sum(as.numeric(entry_bytes), na.rm = TRUE),
  min_dicom_bytes = min(as.numeric(entry_bytes), na.rm = TRUE),
  max_dicom_bytes = max(as.numeric(entry_bytes), na.rm = TRUE),
  n_path_participants = uniqueN(path_participant),
  n_path_dates = uniqueN(path_date[!is.na(path_date)]),
  path_participant = paste(unique(path_participant), collapse = " | "),
  path_date = paste(unique(as.character(path_date[!is.na(path_date)])), collapse = " | "),
  path_description = paste(unique(path_description), collapse = " | ")
), by = .(cohort, archive_role, image_id, zip_name, zip_path)]

series <- series_copy[, .(
  n_archive_copies = .N,
  archive_names = paste(zip_name, collapse = " | "),
  archive_roles = paste(unique(archive_role), collapse = " | "),
  n_dicom_min_copy = min(n_dicom),
  n_dicom_max_copy = max(n_dicom),
  duplicate_copy_slice_count_concordant = min(n_dicom) == max(n_dicom),
  total_dicom_bytes_min_copy = min(total_dicom_bytes),
  total_dicom_bytes_max_copy = max(total_dicom_bytes),
  duplicate_copy_byte_count_concordant = min(total_dicom_bytes) == max(total_dicom_bytes),
  path_participant = paste(unique(path_participant), collapse = " | "),
  path_date = paste(unique(path_date), collapse = " | "),
  path_description = paste(unique(path_description), collapse = " | ")
), by = .(cohort, image_id)]

metadata <- xml_entries[, .(
  n_metadata_copies = .N,
  metadata_archives = paste(unique(zip_name), collapse = " | ")
), by = .(cohort, image_id)]

series <- merge(series, metadata, by = c("cohort", "image_id"), all.x = TRUE)
series[is.na(n_metadata_copies), n_metadata_copies := 0L]
series[, has_metadata := n_metadata_copies > 0L]
series[, structural_T1_name := grepl(
  "MPRAGE|MP[-_ ]?RAGE|IR[-_ ]?SPGR|IR[-_ ]?FSPGR|FSPGR|SPGR", path_description,
  ignore.case = TRUE
)]
series[, excluded_sequence_name := grepl(
  "CALIB|B1|LOCALIZER|SCOUT|FIELD.?MAP|PHANTOM|FLAIR|DTI|BOLD|REST|ASL",
  path_description, ignore.case = TRUE
)]
series[, enhanced_multiframe_candidate :=
         n_dicom_max_copy == 1L & total_dicom_bytes_max_copy >= 5 * 1024^2]
series[, plausible_slice_count :=
         (n_dicom_max_copy >= 80L & n_dicom_max_copy <= 400L) |
         enhanced_multiframe_candidate]
write_csv(series_copy, "02A_series_inventory_by_archive_copy.csv")
write_csv(series, "02B_unique_DICOM_series_inventory.csv")

# -----------------------------------------------------------------------------
# 4. Exact expected-versus-downloaded identity and visit audit.
# -----------------------------------------------------------------------------
audit <- merge(expected, series, by = c("cohort", "image_id"), all.x = TRUE)
audit[, downloaded := !is.na(n_dicom_max_copy)]
audit[, participant_match := downloaded &
        vapply(seq_len(.N), function(i) {
          expected_id <- norm_id(participant_id[i])
          actual_ids <- norm_id(trimws(unlist(strsplit(path_participant[i], "\\|"))))
          expected_id %in% actual_ids
        }, logical(1))]
audit[, date_match := downloaded &
        vapply(seq_len(.N), function(i) {
          as.character(expected_date[i]) %in%
            trimws(unlist(strsplit(path_date[i], "\\|")))
        }, logical(1))]
audit[, role_match := downloaded & (
  (role == "PRIMARY" & grepl("PRIMARY", archive_roles)) |
    (role == "QC_BACKUP" & grepl("QC_BACKUP", archive_roles))
)]
audit[, series_pass := downloaded & participant_match & date_match & role_match &
        has_metadata & structural_T1_name & !excluded_sequence_name & plausible_slice_count]
audit[, failure_reason := paste(
  fifelse(!downloaded, "missing_image", ""),
  fifelse(downloaded & !participant_match, "participant_mismatch", ""),
  fifelse(downloaded & !date_match, "date_mismatch", ""),
  fifelse(downloaded & !role_match, "role_mismatch", ""),
  fifelse(downloaded & !has_metadata, "metadata_missing", ""),
  fifelse(downloaded & !structural_T1_name, "not_structural_T1_name", ""),
  fifelse(downloaded & excluded_sequence_name, "excluded_sequence", ""),
  fifelse(downloaded & !plausible_slice_count, "implausible_slice_count", ""),
  sep = ";"
)]
audit[, failure_reason := gsub("(^;+|;+$)", "", gsub(";+", ";", failure_reason))]
write_csv(audit, "03_expected_downloaded_identity_visit_series_audit.csv")
write_csv(audit[series_pass == FALSE], "03B_failed_series_audit.csv")

extra <- series[!paste(cohort, image_id) %in% paste(expected$cohort, expected$image_id)]
write_csv(extra, "03C_unexpected_downloaded_series.csv")

duplicate_audit <- series[n_archive_copies > 1L, .(
  cohort, image_id, n_archive_copies, archive_names,
  n_dicom_min_copy, n_dicom_max_copy, duplicate_copy_slice_count_concordant,
  total_dicom_bytes_min_copy, total_dicom_bytes_max_copy,
  duplicate_copy_byte_count_concordant
)]
write_csv(duplicate_audit, "04_duplicate_series_across_split_ZIPs.csv")

# -----------------------------------------------------------------------------
# 5. Participant-level longitudinal pairing and AIBL baseline identity.
# -----------------------------------------------------------------------------
adni_pair <- audit[cohort == "ADNI" & role == "PRIMARY", .(
  n_expected = .N,
  n_downloaded = sum(downloaded),
  n_pass = sum(series_pass),
  has_MRI1 = any(window_position == "MRI1" & series_pass),
  has_MRI2 = any(window_position == "MRI2" & series_pass),
  MRI1_date = paste(expected_date[window_position == "MRI1"], collapse = " | "),
  MRI2_date = paste(expected_date[window_position == "MRI2"], collapse = " | "),
  MRI1_image_id = paste(image_id[window_position == "MRI1"], collapse = " | "),
  MRI2_image_id = paste(image_id[window_position == "MRI2"], collapse = " | ")
), by = .(RID, participant_id)]
adni_pair[, complete_longitudinal_pair := n_expected == 2L & n_pass == 2L & has_MRI1 & has_MRI2]
write_csv(adni_pair, "05_ADNI_longitudinal_pair_QC.csv")

aibl_qc <- audit[cohort == "AIBL" & role == "PRIMARY", .(
  RID, participant_id, image_id, expected_date, downloaded, series_pass,
  n_dicom = n_dicom_max_copy, path_description, failure_reason
)]
write_csv(aibl_qc, "06_AIBL_baseline_T1_QC.csv")

# -----------------------------------------------------------------------------
# 6. Outcome-blinded primary/backup adjudication manifest.
# -----------------------------------------------------------------------------
ambiguous_rids <- unique(expected_backup$RID)
blind_pool <- audit[cohort == "ADNI" & RID %in% ambiguous_rids &
                      role %in% c("PRIMARY", "QC_BACKUP")]
setorder(blind_pool, RID, window_position, role, image_id)
blind_pool[, blind_series_code := sprintf("BQC_%04d", .GRP),
           by = .(RID, window_position, image_id)]
blind_manifest <- blind_pool[, .(
  blind_series_code, RID, participant_id, window_position, image_id,
  source_role = role, expected_date, downloaded, series_pass,
  n_dicom = n_dicom_max_copy, path_description,
  qc_identity = NA_character_, qc_orientation = NA_character_,
  qc_brain_coverage = NA_character_, qc_motion = NA_character_,
  qc_conversion = NA_character_, qc_final = NA_character_,
  blinded_selection_reason = NA_character_
)]
write_csv(blind_manifest, "07_outcome_blinded_primary_backup_image_QC_manifest.csv")

blind_key <- blind_pool[, .(
  blind_series_code, RID, participant_id, window_position, image_id,
  source_role = role
)]
write_csv(blind_key, "07B_blind_code_key_DO_NOT_USE_DURING_IMAGE_QC.csv")

# -----------------------------------------------------------------------------
# 7. Archive identity and series checks.
# -----------------------------------------------------------------------------
flow <- rbindlist(list(
  audit[, .(
    expected_images = .N, downloaded_images = sum(downloaded),
    passed_ZIP_series = sum(series_pass), failed_ZIP_series = sum(!series_pass)
  ), by = .(cohort, role)],
  data.table(
    cohort = "ADNI", role = "LONGITUDINAL_PRIMARY_PAIR",
    expected_images = nrow(adni_pair),
    downloaded_images = sum(adni_pair$complete_longitudinal_pair),
    passed_ZIP_series = sum(adni_pair$complete_longitudinal_pair),
    failed_ZIP_series = sum(!adni_pair$complete_longitudinal_pair)
  )
), fill = TRUE)
write_csv(flow, "08_sample_and_image_flow.csv")

zip_gate <- all(archive_audit$central_directory_readable) &&
  nrow(extra) == 0L && all(audit$series_pass) &&
  all(adni_pair$complete_longitudinal_pair) && all(aibl_qc$series_pass)

gate <- data.table(
  gate = c(
    "ZIP_CENTRAL_DIRECTORY_AND_MD5_FROZEN",
    "EXPECTED_IMAGE_IDENTITY_AND_METADATA_COMPLETE",
    "ADNI_639_LONGITUDINAL_PAIRS_COMPLETE",
    "AIBL_338_BASELINE_IMAGES_COMPLETE",
    "BLINDED_PRIMARY_BACKUP_QC_FRAME_READY",
    "PIXEL_LEVEL_PREPROCESSING_UNLOCKED"
  ),
  pass = c(
    all(archive_audit$central_directory_readable),
    nrow(extra) == 0L && all(audit$series_pass),
    nrow(adni_pair) == 639L && all(adni_pair$complete_longitudinal_pair),
    nrow(aibl_qc) == 338L && all(aibl_qc$series_pass),
    nrow(blind_manifest) > 0L && all(blind_manifest$downloaded),
    zip_gate
  ),
  interpretation = c(
    "Every ZIP central directory was read and a reproducibility checksum was generated; CRC/pixel decoding is tested during 13D extraction/conversion.",
    "All expected IDs, subjects, dates, metadata, structural sequence names and slice-count ranges pass.",
    "Both frozen MRI visits are present for every ADNI participant.",
    "Every locked AIBL baseline T1 is present.",
    "Primary and alternative scans are prepared without clinical outcomes for image-level adjudication.",
    "PASS permits automatic selective extraction, DICOM-to-NIfTI conversion and pixel-level QC in 13D."
  )
)
write_csv(gate, "09_13C_readiness_gates.csv")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
log_msg("Step 13C completed: ", ifelse(zip_gate, "PASS", "REQUIRES_REVIEW"))
message("Gate: ", file.path(out_dir, "09_13C_readiness_gates.csv"))
