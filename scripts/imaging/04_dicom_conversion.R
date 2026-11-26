# DICOM-to-NIfTI conversion and geometry
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

project_dir <- study_path("project")
dcm2niix_exe <- study_path("dcm2niix_executable")
stage_root <- study_path("mri_stage", "13D")
processed_root <- study_path("mri_processed", "13D")
out_dir <- file.path(project_dir, "results", "13D_ADNI_AIBL_DICOM_to_NIfTI_geometry_QC")
dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CHECKPOINT_EVERY <- 10L
COMPUTE_OUTPUT_MD5 <- TRUE

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
sanitize <- function(x) {
  y <- gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  y <- gsub("_+", "_", y)
  gsub("^_|_$", "", y)
}
norm_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
safe_clear_stage <- function(path) {
  p <- norm_path(path)
  root <- paste0(norm_path(stage_root), "/")
  if (!startsWith(paste0(p, "/"), root)) stop("Unsafe stage deletion target: ", p, call. = FALSE)
  if (dir.exists(p)) unlink(p, recursive = TRUE, force = TRUE)
}

log_msg("Step 13D started: selective extraction, DICOM decoding, NIfTI conversion and geometry QC.")
log_msg("Original ZIP files are read-only and are never modified.")

# -----------------------------------------------------------------------------
# 1. Verify dcm2niix. Version query may return status 3 on Windows despite
#    printing a valid version, so success is based on executable output.
# -----------------------------------------------------------------------------
need(dcm2niix_exe)
version_output <- suppressWarnings(system2(
  dcm2niix_exe, "-v", stdout = TRUE, stderr = TRUE
))
if (!length(version_output) || !any(grepl("dcm2niiX version|v1\\.0", version_output, ignore.case = TRUE))) {
  stop("dcm2niix executable did not return a recognizable version.", call. = FALSE)
}
writeLines(version_output, file.path(out_dir, "00_dcm2niix_version.txt"), useBytes = TRUE)
log_msg(paste(version_output, collapse = " | "))

# -----------------------------------------------------------------------------
# 2. Load the frozen 13C manifest and archive-copy inventory.
# -----------------------------------------------------------------------------
qc13c_dir <- file.path(project_dir, "results", "13C_ADNI_AIBL_ZIP_DICOM_identity_visit_QC")
gate13c <- fread(need(file.path(qc13c_dir, "09_13C_readiness_gates.csv")), showProgress = FALSE)
if (!isTRUE(gate13c[gate == "PIXEL_LEVEL_PREPROCESSING_UNLOCKED", as.logical(pass)])) {
  stop("13C pixel-level preprocessing gate is not PASS.", call. = FALSE)
}
expected <- fread(need(file.path(qc13c_dir, "00_frozen_expected_image_manifest.csv")), showProgress = FALSE)
copies <- fread(need(file.path(qc13c_dir, "02A_series_inventory_by_archive_copy.csv")), showProgress = FALSE)
blind <- fread(need(file.path(qc13c_dir, "07_outcome_blinded_primary_backup_image_QC_manifest.csv")), showProgress = FALSE)

expected[, image_id := as.character(image_id)]
copies[, image_id := as.character(image_id)]
blind[, image_id := as.character(image_id)]
copies[, zip_path := norm_path(zip_path)]

# Choose one canonical ZIP copy for duplicated split archives. The copies were
# already shown to have concordant Image IDs in 13C; prefer the largest complete
# copy and then a stable lexical archive name.
candidate <- merge(
  expected[, .(cohort, role, image_id)],
  copies[, .(cohort, archive_role, image_id, zip_name, zip_path, n_dicom,
             total_dicom_bytes, path_participant, path_date, path_description)],
  by = c("cohort", "image_id"), allow.cartesian = TRUE
)
candidate <- candidate[archive_role == role]
setorder(candidate, cohort, role, image_id, -n_dicom, -total_dicom_bytes, zip_name)
canonical <- candidate[, .SD[1], by = .(cohort, role, image_id)]
if (nrow(canonical) != nrow(expected)) {
  missing_key <- expected[!paste(cohort, role, image_id) %in%
                            paste(canonical$cohort, canonical$role, canonical$image_id)]
  write_csv(missing_key, "00B_images_without_canonical_ZIP.csv")
  stop("Not every expected image has a canonical ZIP copy.", call. = FALSE)
}

work <- merge(expected, canonical, by = c("cohort", "role", "image_id"), all.x = TRUE)
blind_key <- unique(blind[, .(cohort = "ADNI", image_id, blind_series_code)])
work <- merge(work, blind_key, by = c("cohort", "image_id"), all.x = TRUE)

work[, output_stem := fifelse(
  !is.na(blind_series_code) & nzchar(blind_series_code),
  paste0("ADNI_", sanitize(blind_series_code)),
  fifelse(
    cohort == "ADNI",
    paste0("ADNI_RID", sanitize(RID), "_", sanitize(window_position), "_I", image_id),
    paste0("AIBL_RID", sanitize(RID), "_BASELINE_I", image_id)
  )
)]
work[, output_class := fifelse(
  !is.na(blind_series_code) & nzchar(blind_series_code), "blinded_adjudication",
  fifelse(cohort == "ADNI", "primary", "validation")
)]
work[, output_dir := file.path(processed_root, cohort, output_class, output_stem)]
write_csv(work, "01_frozen_conversion_work_manifest.csv")

# -----------------------------------------------------------------------------
# 3. Build an entry index from ZIP central directories. No extraction yet.
# -----------------------------------------------------------------------------
zip_paths <- unique(work$zip_path)
log_msg("Indexing ", length(zip_paths), " canonical ZIP archives.")
entry_index <- rbindlist(lapply(seq_along(zip_paths), function(i) {
  zp <- zip_paths[i]
  z <- as.data.table(utils::unzip(zp, list = TRUE))
  if ("Name" %in% names(z)) setnames(z, "Name", "entry_name")
  if ("Length" %in% names(z)) setnames(z, "Length", "entry_bytes")
  if (!"entry_name" %in% names(z)) setnames(z, names(z)[1], "entry_name")
  z[, image_id := ifelse(
    grepl("/I[0-9]+/", gsub("\\\\", "/", entry_name)),
    sub(".*/I([0-9]+)/.*", "\\1", gsub("\\\\", "/", entry_name)),
    NA_character_
  )]
  z[, extension := tolower(sub(".*\\.", "", basename(entry_name)))]
  z <- z[!is.na(image_id) & extension %in% c("dcm", "ima")]
  z[, zip_path := zp]
  z[, .(zip_path, image_id, entry_name, entry_bytes)]
}), fill = TRUE)
setkey(entry_index, zip_path, image_id)

entry_counts <- entry_index[, .(indexed_dicom_entries = .N), by = .(zip_path, image_id)]
work <- merge(work, entry_counts, by = c("zip_path", "image_id"), all.x = TRUE)
work[, indexed_dicom_entries := fifelse(is.na(indexed_dicom_entries), 0L, indexed_dicom_entries)]
if (any(work$indexed_dicom_entries == 0L)) {
  write_csv(work[indexed_dicom_entries == 0L], "01B_images_without_indexed_DICOM.csv")
  stop("At least one image has no indexed DICOM entry.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 4. Base-R NIfTI-1 header reader. This avoids requiring an additional R package.
# -----------------------------------------------------------------------------
read_raw_num <- function(raw_header, offset, n, what, size, endian) {
  first <- offset + 1L
  last <- offset + n * size
  if (last > length(raw_header)) stop("Header offset exceeds available bytes.")
  con <- rawConnection(raw_header[first:last], "rb")
  on.exit(close(con))
  readBin(con, what = what, n = n, size = size, endian = endian)
}

read_nifti_header <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rb") else file(path, "rb")
  on.exit(close(con))
  h <- readBin(con, what = "raw", n = 352L)
  if (length(h) < 348L) stop("NIfTI header contains fewer than 348 bytes.")
  hdr_little <- read_raw_num(h, 0L, 1L, integer(), 4L, "little")
  hdr_big <- read_raw_num(h, 0L, 1L, integer(), 4L, "big")
  endian <- if (hdr_little %in% c(348L, 540L)) "little" else if (hdr_big %in% c(348L, 540L)) "big" else NA_character_
  if (is.na(endian)) stop("Unrecognized NIfTI header size.")
  sizeof_hdr <- if (endian == "little") hdr_little else hdr_big
  if (sizeof_hdr != 348L) stop("NIfTI-2 is not expected from this dcm2niix pipeline.")
  dim8 <- read_raw_num(h, 40L, 8L, integer(), 2L, endian)
  ndim <- max(0L, min(7L, dim8[1]))
  dims <- if (ndim > 0L) dim8[seq_len(ndim) + 1L] else integer()
  datatype <- read_raw_num(h, 70L, 1L, integer(), 2L, endian)
  bitpix <- read_raw_num(h, 72L, 1L, integer(), 2L, endian)
  pixdim8 <- read_raw_num(h, 76L, 8L, numeric(), 4L, endian)
  vox_offset <- read_raw_num(h, 108L, 1L, numeric(), 4L, endian)
  qform_code <- read_raw_num(h, 252L, 1L, integer(), 2L, endian)
  sform_code <- read_raw_num(h, 254L, 1L, integer(), 2L, endian)
  srow_x <- read_raw_num(h, 280L, 4L, numeric(), 4L, endian)
  srow_y <- read_raw_num(h, 296L, 4L, numeric(), 4L, endian)
  srow_z <- read_raw_num(h, 312L, 4L, numeric(), 4L, endian)
  spatial_dims <- c(dims, rep(NA_integer_, 3L))[1:3]
  spatial_vox <- abs(c(pixdim8[-1], rep(NA_real_, 3L))[1:3])
  fov <- spatial_dims * spatial_vox
  srow <- rbind(srow_x, srow_y, srow_z)
  affine_det <- if (sform_code > 0L && all(is.finite(srow[, 1:3]))) {
    det(srow[, 1:3, drop = FALSE])
  } else NA_real_
  list(
    sizeof_hdr = sizeof_hdr, endian = endian, ndim = ndim,
    dim1 = spatial_dims[1], dim2 = spatial_dims[2], dim3 = spatial_dims[3],
    datatype = datatype, bitpix = bitpix,
    voxel1 = spatial_vox[1], voxel2 = spatial_vox[2], voxel3 = spatial_vox[3],
    fov1 = fov[1], fov2 = fov[2], fov3 = fov[3], vox_offset = vox_offset,
    qform_code = qform_code, sform_code = sform_code, affine_det = affine_det
  )
}

read_json_scalar <- function(path, key, numeric_value = FALSE) {
  if (!file.exists(path)) return(if (numeric_value) NA_real_ else NA_character_)
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  pattern <- if (numeric_value) {
    paste0('"', key, '"\\s*:\\s*(-?[0-9]+(?:\\.[0-9]+)?)')
  } else {
    paste0('"', key, '"\\s*:\\s*"([^"]*)"')
  }
  m <- regexec(pattern, txt, perl = TRUE)
  hit <- regmatches(txt, m)[[1]]
  if (length(hit) < 2L) return(if (numeric_value) NA_real_ else NA_character_)
  if (numeric_value) suppressWarnings(as.numeric(hit[2])) else hit[2]
}

# -----------------------------------------------------------------------------
# 5. Resumable one-series-at-a-time conversion.
# -----------------------------------------------------------------------------
checkpoint_file <- file.path(out_dir, "02_conversion_geometry_QC_checkpoint.csv")
previous <- if (file.exists(checkpoint_file)) fread(checkpoint_file, showProgress = FALSE) else data.table()
if (nrow(previous)) previous[, image_id := as.character(image_id)]

result_rows <- vector("list", nrow(work))

for (i in seq_len(nrow(work))) {
  w <- work[i]
  out_series_dir <- w$output_dir
  dir.create(out_series_dir, recursive = TRUE, showWarnings = FALSE)
  complete_flag <- file.path(out_series_dir, "13D_conversion_complete.flag")

  existing_nifti <- list.files(out_series_dir, pattern = "\\.nii(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  existing_json <- list.files(out_series_dir, pattern = "\\.json$", full.names = TRUE, ignore.case = TRUE)
  resume_ok <- file.exists(complete_flag) && length(existing_nifti) == 1L && length(existing_json) >= 1L

  stage_dir <- file.path(stage_root, paste0("I", w$image_id))
  safe_clear_stage(stage_dir)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  conversion_log <- file.path(out_series_dir, "dcm2niix.log")
  extracted_count <- NA_integer_
  conversion_status <- NA_integer_
  conversion_message <- NA_character_

  if (!resume_ok) {
    # Remove incomplete outputs from an earlier interrupted attempt, but retain
    # the log until it is overwritten by this attempt.
    stale <- list.files(out_series_dir, full.names = TRUE, all.files = FALSE)
    stale <- stale[!basename(stale) %in% c("dcm2niix.log")]
    if (length(stale)) unlink(stale, recursive = TRUE, force = TRUE)

    entries <- entry_index[.(w$zip_path, w$image_id), entry_name]
    entries <- unique(entries[!is.na(entries)])
    extracted <- tryCatch(
      utils::unzip(w$zip_path, files = entries, exdir = stage_dir, junkpaths = TRUE),
      error = function(e) structure(character(), extraction_error = conditionMessage(e))
    )
    extracted_count <- length(extracted)
    extraction_error <- attr(extracted, "extraction_error")

    if (!is.null(extraction_error) || extracted_count == 0L) {
      conversion_status <- 9001L
      conversion_message <- if (!is.null(extraction_error)) extraction_error else "No DICOM extracted"
      writeLines(conversion_message, conversion_log, useBytes = TRUE)
    } else {
      args <- c(
        "-b", "y", "-ba", "y", "-z", "y", "-p", "y", "-x", "n",
        "-f", w$output_stem, "-o", norm_path(out_series_dir), norm_path(stage_dir)
      )
      output <- tryCatch(
        suppressWarnings(system2(dcm2niix_exe, args = shQuote(args), stdout = TRUE, stderr = TRUE)),
        error = function(e) structure(conditionMessage(e), status = 9002L)
      )
      conversion_status <- attr(output, "status")
      if (is.null(conversion_status)) conversion_status <- 0L
      conversion_message <- paste(output, collapse = "\n")
      writeLines(as.character(output), conversion_log, useBytes = TRUE)
    }
  } else {
    conversion_status <- 0L
    conversion_message <- "RESUMED_FROM_COMPLETE_OUTPUT"
    extracted_count <- w$indexed_dicom_entries
  }

  nifti_files <- list.files(out_series_dir, pattern = "\\.nii(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  json_files <- list.files(out_series_dir, pattern = "\\.json$", full.names = TRUE, ignore.case = TRUE)
  one_nifti <- length(nifti_files) == 1L
  one_or_more_json <- length(json_files) >= 1L
  nifti_path <- if (one_nifti) nifti_files[1] else NA_character_
  json_path <- if (one_or_more_json) json_files[1] else NA_character_

  header <- tryCatch(
    if (one_nifti) read_nifti_header(nifti_path) else NULL,
    error = function(e) structure(NULL, header_error = conditionMessage(e))
  )
  header_error <- attr(header, "header_error")
  header_ok <- !is.null(header)

  if (header_ok) {
    dimension_ok <- header$ndim >= 3L &&
      all(c(header$dim1, header$dim2, header$dim3) >= 64L) &&
      all(c(header$dim1, header$dim2, header$dim3) <= 512L)
    voxel_ok <- all(c(header$voxel1, header$voxel2, header$voxel3) >= 0.4) &&
      all(c(header$voxel1, header$voxel2, header$voxel3) <= 3.0)
    fov_ok <- all(c(header$fov1, header$fov2, header$fov3) >= 120) &&
      all(c(header$fov1, header$fov2, header$fov3) <= 350)
    spatial_transform_ok <- header$qform_code > 0L || header$sform_code > 0L
    affine_ok <- is.na(header$affine_det) || abs(header$affine_det) > 1e-6
  } else {
    dimension_ok <- voxel_ok <- fov_ok <- spatial_transform_ok <- affine_ok <- FALSE
  }

  nifti_bytes <- if (one_nifti) file.info(nifti_path)$size else NA_real_
  output_size_ok <- is.finite(nifti_bytes) && nifti_bytes >= 500 * 1024
  conversion_pass <- conversion_status == 0L && one_nifti && one_or_more_json &&
    header_ok && dimension_ok && voxel_ok && fov_ok && spatial_transform_ok &&
    affine_ok && output_size_ok

  json_series <- if (one_or_more_json) read_json_scalar(json_path, "SeriesDescription") else NA_character_
  json_manufacturer <- if (one_or_more_json) read_json_scalar(json_path, "Manufacturer") else NA_character_
  json_model <- if (one_or_more_json) read_json_scalar(json_path, "ManufacturersModelName") else NA_character_
  json_field <- if (one_or_more_json) read_json_scalar(json_path, "MagneticFieldStrength", TRUE) else NA_real_
  json_slice <- if (one_or_more_json) read_json_scalar(json_path, "SliceThickness", TRUE) else NA_real_

  if (conversion_pass && !file.exists(complete_flag)) {
    writeLines(c(
      paste0("cohort=", w$cohort), paste0("image_id=", w$image_id),
      paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("dcm2niix=", paste(version_output, collapse = " | "))
    ), complete_flag, useBytes = TRUE)
  }

  nifti_md5 <- if (COMPUTE_OUTPUT_MD5 && one_nifti) unname(tools::md5sum(nifti_path)) else NA_character_
  json_md5 <- if (COMPUTE_OUTPUT_MD5 && one_or_more_json) unname(tools::md5sum(json_path)) else NA_character_

  result_rows[[i]] <- data.table(
    cohort = w$cohort, role = w$role, image_id = w$image_id,
    RID = w$RID, participant_id = w$participant_id,
    window_position = w$window_position, blind_series_code = w$blind_series_code,
    output_stem = w$output_stem, zip_name = w$zip_name, zip_path = w$zip_path,
    indexed_dicom_entries = w$indexed_dicom_entries,
    extracted_dicom_entries = extracted_count,
    resumed = resume_ok, conversion_status = conversion_status,
    conversion_message = conversion_message,
    n_nifti_outputs = length(nifti_files), n_json_outputs = length(json_files),
    nifti_path = if (one_nifti) norm_path(nifti_path) else NA_character_,
    json_path = if (one_or_more_json) norm_path(json_path) else NA_character_,
    nifti_bytes = nifti_bytes, nifti_md5 = nifti_md5, json_md5 = json_md5,
    header_ok = header_ok, header_error = if (is.null(header_error)) NA_character_ else header_error,
    ndim = if (header_ok) header$ndim else NA_integer_,
    dim1 = if (header_ok) header$dim1 else NA_integer_,
    dim2 = if (header_ok) header$dim2 else NA_integer_,
    dim3 = if (header_ok) header$dim3 else NA_integer_,
    voxel1 = if (header_ok) header$voxel1 else NA_real_,
    voxel2 = if (header_ok) header$voxel2 else NA_real_,
    voxel3 = if (header_ok) header$voxel3 else NA_real_,
    fov1 = if (header_ok) header$fov1 else NA_real_,
    fov2 = if (header_ok) header$fov2 else NA_real_,
    fov3 = if (header_ok) header$fov3 else NA_real_,
    datatype = if (header_ok) header$datatype else NA_integer_,
    bitpix = if (header_ok) header$bitpix else NA_integer_,
    qform_code = if (header_ok) header$qform_code else NA_integer_,
    sform_code = if (header_ok) header$sform_code else NA_integer_,
    affine_det = if (header_ok) header$affine_det else NA_real_,
    json_series_description = json_series, json_manufacturer = json_manufacturer,
    json_model = json_model, json_field_strength = json_field,
    json_slice_thickness = json_slice,
    dimension_ok = dimension_ok, voxel_ok = voxel_ok,
    geometric_FOV_coverage_ok = fov_ok,
    spatial_transform_ok = spatial_transform_ok, affine_ok = affine_ok,
    output_size_ok = output_size_ok, conversion_pass = conversion_pass
  )

  safe_clear_stage(stage_dir)

  if (i %% CHECKPOINT_EVERY == 0L || i == nrow(work)) {
    current <- rbindlist(result_rows[seq_len(i)], fill = TRUE)
    fwrite(current, checkpoint_file, bom = TRUE)
    log_msg("Processed ", i, "/", nrow(work),
            "; passed=", current[conversion_pass == TRUE, .N],
            "; failed=", current[conversion_pass != TRUE, .N])
  }
}

results <- rbindlist(result_rows, fill = TRUE)
write_csv(results, "02_conversion_geometry_QC_final.csv")
write_csv(results[conversion_pass != TRUE], "02B_failed_conversion_or_geometry_QC.csv")

# -----------------------------------------------------------------------------
# 6. Geometry manifest for outcome-blinded image inspection.
# -----------------------------------------------------------------------------
blind_qc <- results[!is.na(blind_series_code) & nzchar(blind_series_code), .(
  blind_series_code, RID, participant_id, window_position, image_id,
  conversion_pass, n_nifti_outputs, dim1, dim2, dim3,
  voxel1, voxel2, voxel3, fov1, fov2, fov3,
  geometric_FOV_coverage_ok, spatial_transform_ok, affine_ok,
  nifti_bytes, nifti_path
)]
setorder(blind_qc, RID, window_position, blind_series_code)
write_csv(blind_qc, "03_outcome_blinded_geometry_and_visual_QC_manifest.csv")

blind_key <- results[!is.na(blind_series_code) & nzchar(blind_series_code), .(
  blind_series_code, image_id, source_role = role, RID, participant_id, window_position
)]
write_csv(blind_key, "03B_blind_key_DO_NOT_OPEN_DURING_VISUAL_QC.csv")

# -----------------------------------------------------------------------------
# 7. Cohort flow and readiness gates.
# -----------------------------------------------------------------------------
flow <- results[, .(
  expected_images = .N,
  converted_images = sum(n_nifti_outputs >= 1L),
  passed_geometry_QC = sum(conversion_pass == TRUE),
  failed_geometry_QC = sum(conversion_pass != TRUE)
), by = .(cohort, role)]

adni_pairs <- results[cohort == "ADNI" & role == "PRIMARY", .(
  expected_visits = .N,
  passed_visits = sum(conversion_pass == TRUE),
  MRI1_pass = any(window_position == "MRI1" & conversion_pass == TRUE),
  MRI2_pass = any(window_position == "MRI2" & conversion_pass == TRUE)
), by = .(RID, participant_id)]
adni_pairs[, complete_converted_pair := expected_visits == 2L & passed_visits == 2L & MRI1_pass & MRI2_pass]
write_csv(adni_pairs, "04_ADNI_converted_longitudinal_pair_QC.csv")

flow <- rbind(
  flow,
  data.table(
    cohort = "ADNI", role = "COMPLETE_LONGITUDINAL_PAIR",
    expected_images = nrow(adni_pairs),
    converted_images = sum(adni_pairs$complete_converted_pair),
    passed_geometry_QC = sum(adni_pairs$complete_converted_pair),
    failed_geometry_QC = sum(!adni_pairs$complete_converted_pair)
  ), fill = TRUE
)
write_csv(flow, "05_conversion_sample_flow.csv")

all_conversion_pass <- all(results$conversion_pass)
gate <- data.table(
  gate = c(
    "DCM2NIIX_VERSION_FROZEN",
    "ZIP_SELECTIVE_EXTRACTION_COMPLETE",
    "ALL_EXPECTED_DICOM_SERIES_DECODED",
    "NIFTI_GEOMETRY_AND_FOV_PLAUSIBLE",
    "ADNI_LONGITUDINAL_NIFTI_PAIRS_COMPLETE",
    "AIBL_BASELINE_NIFTI_COMPLETE",
    "BLINDED_AMBIGUOUS_SCAN_QC_FRAME_READY",
    "VISUAL_MOTION_AND_BRAIN_COVERAGE_QC_COMPLETE"
  ),
  pass = c(
    TRUE,
    all(results$extracted_dicom_entries > 0L),
    all(results$n_nifti_outputs == 1L),
    all_conversion_pass,
    nrow(adni_pairs) == 639L && all(adni_pairs$complete_converted_pair),
    results[cohort == "AIBL" & role == "PRIMARY", .N] == 338L &&
      all(results[cohort == "AIBL" & role == "PRIMARY", conversion_pass]),
    nrow(blind_qc) > 0L && all(blind_qc$conversion_pass),
    FALSE
  ),
  interpretation = c(
    paste(version_output, collapse = " | "),
    "Every frozen Image ID was selectively extracted from its canonical ZIP; original ZIPs remain unchanged.",
    "dcm2niix produced exactly one NIfTI volume for every selected DICOM series.",
    "NIfTI headers, spatial dimensions, voxel sizes, field of view and transforms passed automated checks.",
    "Both longitudinal MRI visits converted and passed for every retained ADNI participant.",
    "All locked AIBL baseline images converted and passed automated geometry QC.",
    "Ambiguous primary/backup scans are coded without source role for blinded image review.",
    "Conversion checks cover image identity and geometry. Raw-T1 metrics are calculated in imaging/05_t1_qc.R."
  )
)
write_csv(gate, "06_13D_readiness_gates.csv")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_sessionInfo.txt"), useBytes = TRUE)
log_msg("Step 13D completed: ", ifelse(all_conversion_pass, "AUTOMATED_QC_PASS", "REQUIRES_REVIEW"))
message("Gate: ", file.path(out_dir, "06_13D_readiness_gates.csv"))
