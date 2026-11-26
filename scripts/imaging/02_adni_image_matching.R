# ADNI original-image identifiers and repeat pairs
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

project_dir <- study_path("project")
export_dir <- study_path("cohorts", "LBA_MRI/IDA_exports")
out_dir <- file.path(project_dir, "results", "13B1_ADNI_original_T1_ImageID_matching")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] [INFO] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
write_csv <- function(x, n) {
  fwrite(as.data.table(x), file.path(out_dir, n), bom = TRUE)
  log_msg("Wrote ", n, " (", nrow(x), " rows)")
}
present <- function(x) !is.na(x) & trimws(as.character(x)) != "" &
  !trimws(as.character(x)) %in% c("NA", "NaN", ".", "-4")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
date_mdy <- function(x) suppressWarnings(as.IDate(as.character(x), format = "%m/%d/%Y"))
date_iso <- function(x) suppressWarnings(as.IDate(as.character(x)))
norm_phase <- function(x) {
  y <- toupper(gsub("[^A-Z0-9]", "", as.character(x)))
  fifelse(y == "ADNI1", "ADNI1",
  fifelse(y %in% c("ADNIGO", "GO"), "ADNIGO",
  fifelse(y == "ADNI2", "ADNI2",
  fifelse(y == "ADNI3", "ADNI3",
  fifelse(y == "ADNI4", "ADNI4", y)))))
}
extract_protocol <- function(x, key) {
  z <- sub(paste0(".*", key, "=([^;]+).*"), "\\1", as.character(x), ignore.case = TRUE)
  z[!grepl(paste0(key, "="), as.character(x), ignore.case = TRUE)] <- NA_character_
  trimws(z)
}
sequence_family <- function(desc) {
  x <- toupper(as.character(desc))
  accelerated <- grepl("ACCEL|GRAPPA|SENSE", x)
  mprage <- grepl("MP.?RAGE|MPRAGE", x)
  irspgr <- grepl("IR.?F?SPGR|FSPGR|SPGR", x)
  fifelse(mprage & accelerated, "MPRAGE_ACCELERATED",
  fifelse(mprage, "MPRAGE_STANDARD",
  fifelse(irspgr & accelerated, "IRSPGR_ACCELERATED",
  fifelse(irspgr, "IRSPGR_STANDARD", "OTHER_T1"))))
}
target_visit_class <- function(phase, viscode, viscode2) {
  p <- norm_phase(phase)
  v <- tolower(ifelse(present(viscode2), as.character(viscode2), as.character(viscode)))
  out <- rep(NA_character_, length(v))
  out[v %in% c("sc", "scmri")] <- "SCREENING"
  out[v %in% c("bl", "init")] <- "BASELINE_OR_INITIAL"
  m <- suppressWarnings(as.integer(sub("^m", "", v)))
  out[grepl("^m[0-9]+$", v)] <- paste0("MONTH_", m[grepl("^m[0-9]+$", v)])
  out[p == "ADNI2" & v == "scmri"] <- "ADNI2_SCREENING_MRI"
  out[p == "ADNI3" & v == "init"] <- "ADNI3_INITIAL"
  out
}
ida_visit_class <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("SCREENING", y) & !grepl("ADNI2|ADNIGO", y)] <- "SCREENING"
  out[grepl("BASELINE", y)] <- "BASELINE_OR_INITIAL"
  out[grepl("ADNI2 SCREENING MRI", y)] <- "ADNI2_SCREENING_MRI"
  out[grepl("ADNI3 INITIAL", y)] <- "ADNI3_INITIAL"
  month <- suppressWarnings(as.integer(sub(".*MONTH ([0-9]+).*", "\\1", y)))
  has_month <- grepl("MONTH [0-9]+", y)
  out[has_month] <- paste0("MONTH_", month[has_month])
  year <- suppressWarnings(as.integer(sub(".*YEAR ([0-9]+).*", "\\1", y)))
  has_year <- grepl("YEAR [0-9]+", y)
  out[has_year] <- paste0("MONTH_", 12L * year[has_year])
  out
}
write_batches <- function(ids, prefix, size = 200L) {
  ids <- unique(as.character(ids[present(ids)]))
  if (!length(ids)) return(0L)
  groups <- split(ids, ceiling(seq_along(ids) / size))
  for (i in seq_along(groups)) {
    writeLines(paste(groups[[i]], collapse = ","),
               file.path(out_dir, sprintf("%s_batch_%02d.txt", prefix, i)), useBytes = TRUE)
  }
  length(groups)
}

log_msg("Step 13B1 started: match IDA original T1 Image IDs to the 1,320 frozen ADNI target visits.")

target_file <- file.path(project_dir, "results", "13B_longitudinal_MRI_acquisition_cohort_freeze",
                         "02_ADNI_original_T1_visit_lookup_targets.csv")
if (!file.exists(target_file)) stop("Missing target file: ", target_file, call. = FALSE)
target <- fread(target_file, showProgress = FALSE)
if (nrow(target) != 1320L || uniqueN(target$PTID) != 660L) {
  stop("Frozen ADNI target is not 1,320 visits/660 people.", call. = FALSE)
}
target[, target_id := .I]
target[, target_date := date_iso(mri_date)]
target[, target_phase := norm_phase(PHASE)]
target[, target_visit := target_visit_class(PHASE, VISCODE, VISCODE2)]
target[, target_field := num(gsub("[^0-9.]", "", FIELD_STRENGTH))]

# Private source-availability exclusions apply only to raw-image processing.
exclusion_file <- study_path("adni_unavailable_ptids")
source_unavailable_ptid <- character()
if (file.exists(exclusion_file)) {
  exclusions <- fread(exclusion_file, colClasses = "character", showProgress = FALSE)
  if (!"PTID" %in% names(exclusions)) stop("Source-exclusion CSV requires a PTID column.", call. = FALSE)
  source_unavailable_ptid <- unique(trimws(exclusions$PTID[present(exclusions$PTID)]))
  if (any(!source_unavailable_ptid %in% target$PTID))
    stop("Source-exclusion CSV contains participants outside the acquisition target.", call. = FALSE)
} else {
  message("No private source-exclusion CSV supplied; retaining all acquisition targets.")
}
source_unavailable <- target[PTID %in% source_unavailable_ptid]
source_unavailable[, exclusion_reason :=
  "No image returned by unrestricted IDA participant search; raw source cannot be verified or downloaded"]
write_csv(source_unavailable, "00B_source_unavailable_prespecified_exclusions.csv")
analysis_target <- target[!PTID %in% source_unavailable_ptid]
expected_people <- uniqueN(analysis_target$RID)
expected_visits <- nrow(analysis_target)
if (file.exists(exclusion_file) && (expected_people != 657L || expected_visits != 1314L)) {
  stop("Unexpected post-exclusion target size: ", expected_people, " people / ", expected_visits, " visits.", call. = FALSE)
}

files <- list.files(export_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No IDA CSV exports in: ", export_dir, call. = FALSE)
exports <- rbindlist(lapply(files, function(f) {
  z <- fread(f, showProgress = FALSE, fill = TRUE)
  required <- c("Subject ID", "Phase", "Visit", "Study Date", "Description", "Imaging Protocol", "Image ID")
  if (!all(required %in% names(z))) return(NULL)
  z[, source_export := basename(f)]
  z
}), fill = TRUE, use.names = TRUE)
if (!nrow(exports)) stop("No export contains the required full IDA fields.", call. = FALSE)
exports <- unique(exports, by = "Image ID")

exports[, candidate_date := date_mdy(`Study Date`)]
exports[, candidate_phase := norm_phase(Phase)]
exports[, candidate_visit := ida_visit_class(Visit)]
exports[, field_strength := num(extract_protocol(`Imaging Protocol`, "Field Strength"))]
exports[, slice_thickness := num(extract_protocol(`Imaging Protocol`, "Slice Thickness"))]
exports[, acquisition_plane := toupper(extract_protocol(`Imaging Protocol`, "Acquisition Plane"))]
exports[, acquisition_type := toupper(extract_protocol(`Imaging Protocol`, "Acquisition Type"))]
exports[, weighting := toupper(extract_protocol(`Imaging Protocol`, "Weighting"))]
exports[, manufacturer := extract_protocol(`Imaging Protocol`, "Manufacturer")]
exports[, model := extract_protocol(`Imaging Protocol`, "Mfg Model")]
exports[, sequence_family := sequence_family(Description)]

desc_up <- toupper(as.character(exports$Description))
exports[, excluded_nonstructural := grepl(
  "CALIB|B1|LOCALIZER|SCOUT|FIELD.?MAP|PHANTOM|COIL|SURVEY|ASL|DTI|FLAIR|T2|PD|FMRI|BOLD|REST",
  desc_up
)]
exports[, structural_name := grepl("MP.?RAGE|MPRAGE|IR.?F?SPGR|FSPGR|SPGR", desc_up)]
exports[, protocol_3d := acquisition_type == "3D"]
exports[, protocol_sagittal := acquisition_plane == "SAGITTAL"]
exports[, protocol_t1 := weighting == "T1"]
exports[, repeat_named := grepl("REPEAT", desc_up)]

coverage <- data.table(
  metric = c("full_export_files", "unique_export_rows", "target_people", "covered_target_people", "missing_target_people"),
  n = c(length(unique(exports$source_export)), nrow(exports), uniqueN(target$PTID),
        uniqueN(exports[`Subject ID` %in% analysis_target$PTID, `Subject ID`]),
        uniqueN(analysis_target[!PTID %in% exports$`Subject ID`, PTID]))
)
write_csv(coverage, "00_export_coverage.csv")

# Candidate target joins are restricted to the same participant and a 30-day
# window. Exact date, correct phase/visit and structural protocol dominate rank.
cand <- merge(
  analysis_target[, .(target_id, RID, PTID, target_phase, VISCODE, VISCODE2, target_visit,
             window_position, target_date, target_field, processed_FreeSurfer_IMAGEUID)],
  exports[, .(`Subject ID`, Phase, candidate_phase, Visit, candidate_visit,
              `Study Date`, candidate_date, Description, `Imaging Protocol`, `Image ID`,
              field_strength, slice_thickness, acquisition_plane, acquisition_type,
              weighting, manufacturer, model, sequence_family, excluded_nonstructural,
              structural_name, protocol_3d, protocol_sagittal, protocol_t1,
              repeat_named, source_export)],
  by.x = "PTID", by.y = "Subject ID", allow.cartesian = TRUE
)
cand[, date_gap_days := abs(as.integer(candidate_date - target_date))]
cand <- cand[is.finite(date_gap_days) & date_gap_days <= 30L &
               structural_name == TRUE & excluded_nonstructural == FALSE]
cand[, phase_match := candidate_phase == target_phase]
cand[, visit_match := !is.na(candidate_visit) & !is.na(target_visit) & candidate_visit == target_visit]
cand[, field_match := !is.na(field_strength) & !is.na(target_field) & abs(field_strength - target_field) < 0.2]
cand[, individual_score :=
       10000 * (!protocol_3d) + 10000 * (!protocol_t1) + 5000 * (!protocol_sagittal) +
       1000 * (!phase_match) + 200 * (!visit_match) + 100 * (!field_match) +
       20 * repeat_named + 10 * date_gap_days +
       fifelse(is.finite(slice_thickness) & slice_thickness <= 1.5, 0, 25)]

candidate_counts <- cand[, .(
  n_candidates = .N,
  n_exact_date = sum(date_gap_days == 0L),
  n_exact_date_phase_visit = sum(date_gap_days == 0L & phase_match & visit_match),
  best_individual_score = min(individual_score)
), by = target_id]
target_audit <- merge(analysis_target, candidate_counts, by = "target_id", all.x = TRUE)
for (v in c("n_candidates", "n_exact_date", "n_exact_date_phase_visit")) target_audit[is.na(get(v)), (v) := 0L]
write_csv(target_audit, "01_target_candidate_counts.csv")

# Jointly choose MRI1/MRI2 per person. Prefer the same sequence family and field
# strength across time after the individual target-match criteria are satisfied.
c1 <- cand[window_position == "MRI1"]
c2 <- cand[window_position == "MRI2"]
pairs <- merge(c1, c2, by = c("RID", "PTID"), suffixes = c("_MRI1", "_MRI2"), allow.cartesian = TRUE)
pairs <- pairs[`Image ID_MRI1` != `Image ID_MRI2`]
pairs[, pair_sequence_match := sequence_family_MRI1 == sequence_family_MRI2]
pairs[, pair_field_match := is.finite(field_strength_MRI1) & is.finite(field_strength_MRI2) &
        abs(field_strength_MRI1 - field_strength_MRI2) < 0.2]
pairs[, pair_manufacturer_match := !is.na(manufacturer_MRI1) & !is.na(manufacturer_MRI2) &
        toupper(manufacturer_MRI1) == toupper(manufacturer_MRI2)]
pairs[, pair_score := individual_score_MRI1 + individual_score_MRI2 +
        40 * (!pair_sequence_match) + 20 * (!pair_field_match) + 10 * (!pair_manufacturer_match)]
setorder(pairs, RID, pair_score, date_gap_days_MRI1, date_gap_days_MRI2,
         `Image ID_MRI1`, `Image ID_MRI2`, na.last = TRUE)
pairs[, tied_best_pairs := sum(pair_score == min(pair_score)), by = RID]
best <- pairs[, .SD[1], by = RID]

selected <- rbindlist(list(
  best[, .(
    RID, PTID, window_position = "MRI1", target_date = target_date_MRI1,
    target_phase = target_phase_MRI1, target_visit = target_visit_MRI1,
    Image_ID = `Image ID_MRI1`, Study_Date = `Study Date_MRI1`, Phase = Phase_MRI1,
    Visit = Visit_MRI1, Description = Description_MRI1,
    Imaging_Protocol = `Imaging Protocol_MRI1`, field_strength = field_strength_MRI1,
    slice_thickness = slice_thickness_MRI1, sequence_family = sequence_family_MRI1,
    date_gap_days = date_gap_days_MRI1, phase_match = phase_match_MRI1,
    visit_match = visit_match_MRI1, field_match = field_match_MRI1,
    pair_sequence_match, pair_field_match, pair_manufacturer_match,
    pair_score, tied_best_pairs, source_export = source_export_MRI1
  )],
  best[, .(
    RID, PTID, window_position = "MRI2", target_date = target_date_MRI2,
    target_phase = target_phase_MRI2, target_visit = target_visit_MRI2,
    Image_ID = `Image ID_MRI2`, Study_Date = `Study Date_MRI2`, Phase = Phase_MRI2,
    Visit = Visit_MRI2, Description = Description_MRI2,
    Imaging_Protocol = `Imaging Protocol_MRI2`, field_strength = field_strength_MRI2,
    slice_thickness = slice_thickness_MRI2, sequence_family = sequence_family_MRI2,
    date_gap_days = date_gap_days_MRI2, phase_match = phase_match_MRI2,
    visit_match = visit_match_MRI2, field_match = field_match_MRI2,
    pair_sequence_match, pair_field_match, pair_manufacturer_match,
    pair_score, tied_best_pairs, source_export = source_export_MRI2
  )]
), use.names = TRUE)
setorder(selected, RID, window_position)
selected[, match_grade := fifelse(
  date_gap_days == 0L & phase_match & visit_match & field_match &
    tied_best_pairs == 1L, "A_LOCKED",
  fifelse(date_gap_days <= 30L & phase_match & field_match &
            tied_best_pairs <= 2L, "B_ACCEPTABLE",
          "C_REVIEW_REQUIRED"))]
write_csv(selected, "02_selected_original_T1_ImageIDs.csv")

selected_people <- unique(selected$RID)
unresolved <- analysis_target[!RID %in% selected_people]
unresolved[, exclusion_reason :=
  "No verifiable paired structural T1 within 30 days of both frozen target visits in current IDA export"]
write_csv(unresolved, "03_unresolved_target_visits.csv")

review <- selected[match_grade == "C_REVIEW_REQUIRED" | tied_best_pairs > 2L]
write_csv(review, "04_selected_matches_requiring_review.csv")

# Preserve every close structural candidate for transparent adjudication.
candidate_export <- cand[, .(
  target_id, RID, PTID, window_position, target_date, target_phase, target_visit,
  Image_ID = `Image ID`, Study_Date = `Study Date`, Phase, Visit, Description,
  Imaging_Protocol = `Imaging Protocol`, field_strength, slice_thickness,
  sequence_family, date_gap_days, phase_match, visit_match, field_match,
  protocol_3d, protocol_sagittal, protocol_t1, repeat_named, individual_score,
  source_export
)]
setorder(candidate_export, RID, window_position, individual_score, Image_ID)
write_csv(candidate_export, "05_all_close_structural_candidates.csv")

# Selected pairs meet the participant and acquisition-window criteria.
# C grade denotes metadata ambiguity. Blinded image QC resolves scan choice.
downloadable <- copy(selected)
complete_people <- downloadable[, .N, by = RID][N == 2L, RID]
downloadable <- downloadable[RID %in% complete_people]
write_csv(downloadable, "06_download_ready_complete_pairs.csv")
review_people <- unique(review$RID)
backup_candidates <- candidate_export[RID %in% review_people & !Image_ID %in% downloadable$Image_ID]
write_csv(backup_candidates, "06B_QC_backup_T1_candidates_for_ambiguous_pairs.csv")
old_batches <- list.files(out_dir, pattern = "^07B?_ADNI_.*ImageIDs_batch_[0-9]+\\.txt$",
                          full.names = TRUE)
if (length(old_batches)) file.remove(old_batches)
n_batches <- write_batches(downloadable$Image_ID, "07_ADNI_original_T1_download_ImageIDs", 200L)
n_backup_batches <- write_batches(backup_candidates$Image_ID, "07B_ADNI_QC_backup_T1_ImageIDs", 200L)

summary <- data.table(
  metric = c(
    "frozen_target_people", "frozen_target_visits", "export_covered_people",
    "people_with_selected_pairs", "A_locked_visits", "B_acceptable_visits",
    "C_review_visits", "unresolved_people", "download_ready_complete_people",
    "download_ready_images", "download_batches", "scanner_or_field_change_people",
    "QC_backup_images", "QC_backup_batches"
  ),
  n = c(
    expected_people, expected_visits, uniqueN(exports[`Subject ID` %in% analysis_target$PTID, `Subject ID`]),
    uniqueN(selected$RID), selected[match_grade == "A_LOCKED", .N],
    selected[match_grade == "B_ACCEPTABLE", .N], selected[match_grade == "C_REVIEW_REQUIRED", .N],
    uniqueN(unresolved$RID), length(complete_people), nrow(downloadable), n_batches,
    uniqueN(best[pair_field_match == FALSE | pair_manufacturer_match == FALSE, RID]),
    uniqueN(backup_candidates$Image_ID), n_backup_batches
  )
)
write_csv(summary, "08_matching_summary.csv")

gate_pass <- length(complete_people) >= 600L && nrow(downloadable) == 2L * length(complete_people)
gate <- data.table(
  gate = "ADNI_SOURCE_AVAILABLE_PAIRED_T1_IMAGE_IDS_FROZEN",
  pass = gate_pass,
  status = if (gate_pass) "READY_TO_DOWNLOAD_07_PRIMARY_AND_07B_BACKUP_BATCHES_FOR_BLINDED_IMAGE_QC" else "REVIEW_03_AND_04_BEFORE_FULL_DOWNLOAD",
  note = if (gate_pass) {
    paste0("All ", length(complete_people), " source-available participants have two primary T1 Image IDs; ambiguous metadata alternatives are placed in 07B for blinded image QC. Portal-unavailable target visits are documented separately.")
  } else {
    paste0("Download batches contain A/B complete pairs. Unresolved or C-grade matches require image-level adjudication. Expected sample size: ", expected_visits, "-image source-available set.")
  }
)
write_csv(gate, "09_13B1_gate.csv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
log_msg("Step 13B1 completed: ", gate$status)
message("Review: ", file.path(out_dir, "08_matching_summary.csv"))
message("Gate: ", file.path(out_dir, "09_13B1_gate.csv"))
