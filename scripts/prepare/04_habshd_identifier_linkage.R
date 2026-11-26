# Consistent identifiers and nearest-date pathology linkage
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

pick_root <- function(ascii, native) if (file.exists(ascii)) ascii else native
project_dir <- pick_root(study_path("project"), study_path("project"))
habs_root <- pick_root(study_path("habs_release"), study_path("cohorts", "0MRI"))
out_dir <- file.path(project_dir, "results", "13F_HABS_linkage_repair")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  path
}
num <- function(x) {
  z <- suppressWarnings(as.numeric(as.character(x)))
  z[z <= -9000] <- NA_real_
  z
}
normalize_id <- function(x) {
  z <- trimws(as.character(x))
  n <- suppressWarnings(as.numeric(z))
  z[is.finite(n)] <- format(n[is.finite(n)], scientific = FALSE, trim = TRUE, digits = 15)
  z <- sub("\\.0+$", "", z)
  z[z %in% c("", "NA", "N/A", ".")] <- NA_character_
  z
}
clean_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  make.unique(x, sep = "_")
}
month_num <- function(x) {
  y <- suppressWarnings(as.integer(x))
  bad <- !is.finite(y)
  if (any(bad)) {
    y[bad] <- match(tolower(substr(trimws(as.character(x[bad])), 1, 3)),
                    tolower(month.abb))
  }
  y
}
make_date <- function(year, month) {
  y <- suppressWarnings(as.integer(year))
  m <- month_num(month)
  out <- as.IDate(rep(NA_character_, length(y)))
  ok <- is.finite(y) & is.finite(m) & m >= 1L & m <= 12L
  out[ok] <- as.IDate(sprintf("%04d-%02d-15", y[ok], m[ok]))
  out
}
read_inclusions <- function(path, pattern = "inclusion|included") {
  sheets <- excel_sheets(path)
  keep <- sheets[grepl(pattern, sheets, ignore.case = TRUE) &
                   !grepl("excl|pending", sheets, ignore.case = TRUE)]
  if (!length(keep)) stop("No inclusion sheets found in ", basename(path), call. = FALSE)
  rbindlist(lapply(keep, function(sheet) {
    x <- suppressWarnings(as.data.table(read_excel(
      path, sheet = sheet, col_types = "text", na = c("", "NA", "N/A", "."),
      .name_repair = "minimal"
    )))
    setnames(x, clean_names(names(x)))
    x[, source_sheet := sheet]
    x
  }), fill = TRUE)
}
find_col <- function(d, patterns) {
  for (pattern in patterns) {
    hit <- grep(pattern, names(d), value = TRUE, ignore.case = TRUE)
    if (length(hit)) return(hit[1])
  }
  stop("Column not found: ", paste(patterns, collapse = " | "), call. = FALSE)
}

ct_path <- need(file.path(habs_root, "HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx"))
hip_path <- need(file.path(habs_root, "HABS_HD_RP_8_Hippocampal_Volume_HippoDeep_Data_BraskieLab_20260501.xlsx"))
amy_path <- need(file.path(habs_root, "HABS_HD_RP_8_FBB_PET_data_BraskieLab_20260501.xlsx"))
tau_path <- need(file.path(habs_root, "HABS_HD_RP_8_PI2620_PET_Data_BraskieLab_20260501.xlsx"))
freeze_path <- file.path(project_dir, "results", "12B_ADNI_HABSHD_structural_MRI_feature_freeze", "01_HABSHD_frozen_cortical_features.rds")

ct_source <- read_inclusions(ct_path)
ct_id <- find_col(ct_source, c("^Med_ID$", "^MED_ID$"))
ct_tp <- find_col(ct_source, c("^Time_point_MRI$", "^Time_point$"))
ct_month <- find_col(ct_source, c("^MRI_month$"))
ct_year <- find_col(ct_source, c("^MRI_year$"))
ct_dates <- unique(ct_source[, .(
  participant_id = normalize_id(get(ct_id)),
  mri_timepoint = as.character(get(ct_tp)),
  mri_date = make_date(get(ct_year), get(ct_month))
)])
ct_dates <- ct_dates[!is.na(participant_id)]

ct_features <- as.data.table(readRDS(need(freeze_path)))
ct_features[, participant_id := normalize_id(participant_id)]
ct_features[, mri_timepoint := as.character(mri_timepoint)]
ct_features <- merge(ct_features, ct_dates,
                      by = c("participant_id", "mri_timepoint"),
                      all.x = TRUE, sort = FALSE)

hip <- read_inclusions(hip_path)
hip_id <- find_col(hip, c("^Med_ID$"))
hip_month <- find_col(hip, c("^MRI_month$"))
hip_year <- find_col(hip, c("^MRI_year$"))
hip_l <- find_col(hip, c("^L_hippocampus$"))
hip_r <- find_col(hip, c("^R_hippocampus$"))
hip_icv <- find_col(hip, c("^ICV_HD$"))
hip <- hip[, .(
  participant_id = normalize_id(get(hip_id)),
  endpoint_date = make_date(get(hip_year), get(hip_month)),
  hippocampus_icv_ratio = (num(get(hip_l)) + num(get(hip_r))) / num(get(hip_icv))
)]
hip <- hip[!is.na(participant_id) & is.finite(hippocampus_icv_ratio)]

brain_path <- need(file.path(habs_root, "HABS_HD_RP_8_Whole_Brain_Volume_Data_BraskieLab_20260501.xlsx"))
brain <- read_inclusions(brain_path)
brain_id <- find_col(brain, c("^Med_ID$"))
brain_month <- find_col(brain, c("^MRI_month$"))
brain_year <- find_col(brain, c("^MRI_year$"))
brain_vol <- find_col(brain, c("^BrainSegNotVent$"))
brain_icv <- find_col(brain, c("^ICV_HD$"))
brain <- brain[, .(
  participant_id = normalize_id(get(brain_id)),
  endpoint_date = make_date(get(brain_year), get(brain_month)),
  brain_icv_ratio = num(get(brain_vol)) / num(get(brain_icv))
)]
brain <- brain[!is.na(participant_id) & is.finite(brain_icv_ratio)]

amy <- read_inclusions(amy_path)
amy_id <- find_col(amy, c("^Med_ID$"))
amy_month <- find_col(amy, c("FBB_PET_month"))
amy_year <- find_col(amy, c("FBB_PET_year"))
amy_value <- find_col(amy, c("Centiloid_Units_GlobalAmyloid", "Global_FBB_SUVR"))
amy <- amy[, .(
  participant_id = normalize_id(get(amy_id)),
  endpoint_date = make_date(get(amy_year), get(amy_month)),
  amyloid_value = num(get(amy_value))
)]
amy <- amy[!is.na(participant_id) & is.finite(amyloid_value)]

tau <- read_inclusions(tau_path)
tau_id <- find_col(tau, c("^Med_ID$"))
tau_month <- find_col(tau, c("PI2620_PET_month"))
tau_year <- find_col(tau, c("PI2620_PET_year"))
tau_value <- find_col(tau, c("Median_Medial_Temporal_PI2620_SUVR", "Medial_Temporal_PI2620_SUVR"))
tau <- tau[, .(
  participant_id = normalize_id(get(tau_id)),
  endpoint_date = make_date(get(tau_year), get(tau_month)),
  tau_value = num(get(tau_value))
)]
tau <- tau[!is.na(participant_id) & is.finite(tau_value)]

audit_one <- function(endpoint, x, value_col) {
  ct_ids <- unique(ct_features[!is.na(mri_date), participant_id])
  endpoint_ids <- unique(x[is.finite(get(value_col)), participant_id])
  overlap_ids <- intersect(ct_ids, endpoint_ids)
  z <- merge(
    ct_features[!is.na(mri_date), .(participant_id, mri_date)],
    x[is.finite(get(value_col)), .(participant_id, endpoint_date, value = get(value_col))],
    by = "participant_id", allow.cartesian = TRUE
  )
  z[, gap_days := abs(as.integer(mri_date - endpoint_date))]
  near <- z[is.finite(gap_days), .(nearest_gap_days = min(gap_days)), by = participant_id]
  data.table(
    endpoint = endpoint,
    raw_rows = nrow(x),
    raw_value_rows = x[is.finite(get(value_col)), .N],
    raw_people = uniqueN(x[is.finite(get(value_col)), participant_id]),
    cortical_rows_with_valid_date = ct_features[!is.na(mri_date), .N],
    cortical_people_with_valid_date = uniqueN(ct_features[!is.na(mri_date), participant_id]),
    id_overlap_people = length(overlap_ids),
    date_valid_pairs = z[!is.na(mri_date) & !is.na(endpoint_date), .N],
    within_548_people = uniqueN(z[is.finite(gap_days) & gap_days <= 548, participant_id]),
    within_1095_people = uniqueN(z[is.finite(gap_days) & gap_days <= 1095, participant_id]),
    median_nearest_gap_days = if (nrow(near)) median(near$nearest_gap_days) else NA_real_,
    max_nearest_gap_days = if (nrow(near)) max(near$nearest_gap_days) else NA_real_
  )
}

audits <- rbindlist(list(
  audit_one("amyloid", amy, "amyloid_value"),
  audit_one("hippocampus", hip, "hippocampus_icv_ratio"),
  audit_one("brain", brain, "brain_icv_ratio"),
  audit_one("tau", tau, "tau_value")
), fill = TRUE)

nearest_link <- function(index, candidate, value_col, prefix, max_days = 548L) {
  a <- copy(index); b <- copy(candidate)
  a[, row_id__ := .I]
  setnames(b, "endpoint_date", "candidate_date__")
  z <- merge(a, b, by = "participant_id", allow.cartesian = TRUE, all.x = TRUE)
  z[, gap_days__ := abs(as.integer(mri_date - candidate_date__))]
  setorder(z, row_id__, gap_days__)
  z <- z[, head(.SD, 1L), by = row_id__]
  keep_value <- value_col
  z[!is.finite(gap_days__) | gap_days__ > max_days, (keep_value) := NA_real_]
  setnames(z, c("candidate_date__", keep_value), c(paste0(prefix, "date"), paste0(prefix, value_col)))
  z[, c("row_id__", "gap_days__") := NULL]
  z
}

linked <- nearest_link(ct_features, hip, "hippocampus_icv_ratio", "hip_", 548L)
linked <- nearest_link(linked, brain, "brain_icv_ratio", "brain_", 548L)
linked <- nearest_link(linked, amy, "amyloid_value", "amyloid_", 548L)
linked <- nearest_link(linked, tau, "tau_value", "tau_", 548L)

fwrite(audits, file.path(out_dir, "01_repaired_linkage_summary.csv"), bom = TRUE)
saveRDS(
  list(
    cortical_features = ct_features,
    amyloid = amy,
    tau = tau,
    hippocampus = hip,
    brain = brain,
    linked = linked
  ),
  file.path(out_dir, "02_repaired_HABSHD_linked_frames.rds"), compress = "xz"
)
gate <- data.table(
  gate = "HABSHD_LINKAGE_REPAIR_WITH_SHARED_ID_DATE_RULES",
  id_rule = "normalize_id used for cortical and every endpoint table",
  date_rule = "15th of reported MRI/PET month and year",
  amyloid_overlap_people = audits[endpoint == "amyloid", id_overlap_people],
  hippocampus_overlap_people = audits[endpoint == "hippocampus", id_overlap_people],
  tau_overlap_people = audits[endpoint == "tau", id_overlap_people],
  within_548_amyloid_people = audits[endpoint == "amyloid", within_548_people],
  within_548_hippocampus_people = audits[endpoint == "hippocampus", within_548_people],
  status = if (all(audits$id_overlap_people > 0 & audits$date_valid_pairs > 0)) "PASS" else "REQUIRES_REVIEW"
)
fwrite(gate, file.path(out_dir, "03_repaired_linkage_gate.csv"), bom = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
message("HABS linkage repair completed: ", out_dir)
