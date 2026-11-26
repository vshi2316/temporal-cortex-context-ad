# Automated T1 flags and processing manifest
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)

# Compute protocol-stratified automatic QC metrics for conversion-passing
# primary images. Flagged images remain available for blinded adjudication.

suppressPackageStartupMessages(library(data.table))

project_dir <- study_path("project")
input_file <- file.path(
  project_dir, "results", "13D_ADNI_AIBL_DICOM_to_NIfTI_geometry_QC",
  "02_conversion_geometry_QC_final.csv"
)
out_dir <- file.path(project_dir, "results", "13F_automated_raw_T1_QC_and_processing_manifest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# A checkpoint makes a long read of 1,614 compressed T1 volumes resumable.
# Set FORCE_RECOMPUTE <- TRUE only if the metric implementation is changed.
CHECKPOINT_EVERY <- 20L
FORCE_RECOMPUTE <- FALSE
ROBUST_Z_THRESHOLD <- 4.5
MIN_STRATUM_N <- 10L

need <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  path
}
write_csv <- function(x, file_name) {
  fwrite(as.data.table(x), file.path(out_dir, file_name), bom = TRUE, na = "")
}
norm_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
safe_chr <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | !nzchar(y)] <- "UNKNOWN"
  y
}
is_true <- function(x) as.character(x) %in% c("TRUE", "True", "true", "1")
finite_quantile <- function(x, probs) {
  z <- x[is.finite(x)]
  if (!length(z)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(z, probs = probs, names = FALSE, type = 8, na.rm = TRUE))
}
safe_mean <- function(x) {
  z <- x[is.finite(x)]
  if (!length(z)) NA_real_ else mean(z)
}
safe_sd <- function(x) {
  z <- x[is.finite(x)]
  if (length(z) < 2L) NA_real_ else stats::sd(z)
}
safe_mad <- function(x) {
  z <- x[is.finite(x)]
  if (length(z) < 3L) NA_real_ else stats::mad(z, center = stats::median(z), constant = 1, na.rm = TRUE)
}

if (!requireNamespace("RNifti", quietly = TRUE)) {
  stop("Package 'RNifti' is required. Install once with install.packages('RNifti').", call. = FALSE)
}

# Otsu thresholding defines a coarse foreground mask for anomaly screening.
otsu_threshold <- function(x, bins = 128L) {
  z <- x[is.finite(x)]
  if (length(z) < 100L) return(NA_real_)
  q <- finite_quantile(z, c(0.005, 0.995))
  if (!all(is.finite(q)) || q[2] <= q[1]) return(NA_real_)
  z <- pmin(q[2], pmax(q[1], z))
  h <- hist(z, breaks = bins, plot = FALSE)
  p <- h$counts / sum(h$counts)
  omega <- cumsum(p)
  mu <- cumsum(p * h$mids)
  mu_t <- mu[length(mu)]
  denom <- omega * (1 - omega)
  sigma_b <- rep(-Inf, length(p))
  ok <- denom > 0
  sigma_b[ok] <- (mu_t * omega[ok] - mu[ok])^2 / denom[ok]
  h$mids[which.max(sigma_b)]
}

entropy_from_values <- function(x, bins = 128L) {
  z <- x[is.finite(x)]
  if (length(z) < 100L) return(NA_real_)
  q <- finite_quantile(z, c(0.005, 0.995))
  if (!all(is.finite(q)) || q[2] <= q[1]) return(NA_real_)
  h <- hist(pmin(q[2], pmax(q[1], z)), breaks = bins, plot = FALSE)$counts
  p <- h[h > 0] / sum(h)
  -sum(p * log(p)) / log(bins)
}

sample_indices <- function(n, maximum = 750000L) {
  if (n <= maximum) seq_len(n) else unique(as.integer(round(seq(1, n, length.out = maximum))))
}

read_volume_metrics <- function(nifti_path) {
  result <- data.table(
    metric_status = "ERROR", metric_error = NA_character_,
    dim1 = NA_integer_, dim2 = NA_integer_, dim3 = NA_integer_,
    finite_fraction = NA_real_, nonzero_fraction = NA_real_,
    intensity_p01 = NA_real_, intensity_p50 = NA_real_, intensity_p99 = NA_real_,
    intensity_iqr = NA_real_, intensity_entropy = NA_real_,
    foreground_fraction = NA_real_, foreground_intensity_cv = NA_real_,
    foreground_bbox_fraction = NA_real_, foreground_center_offset = NA_real_,
    slice_coverage_iqr = NA_real_, slice_dropout_fraction = NA_real_,
    low_resolution_gradient = NA_real_
  )
  tryCatch({
    img <- RNifti::readNifti(nifti_path, internal = FALSE)
    a <- as.array(img)
    d <- dim(a)
    if (length(d) < 3L || any(d[1:3] < 16L)) stop("Volume is not a usable 3D image.")
    # Some valid dcm2niix files carry a trailing singleton time dimension.
    # Collapse it explicitly so all subsequent array indexing is strictly 3D.
    if (length(d) > 3L) a <- a[, , , 1L, drop = TRUE]
    d <- dim(a)[1:3]
    n_all <- length(a)
    idx <- sample_indices(n_all)
    values <- as.numeric(a)[idx]
    finite <- is.finite(values)
    pos <- values[finite & values > 0]
    q <- finite_quantile(pos, c(0.01, 0.25, 0.50, 0.75, 0.99))
    thr <- otsu_threshold(values)

    # A modest regularly spaced grid keeps foreground/bounding-box measures
    # stable while avoiding a second full-volume copy in memory.
    gx <- unique(as.integer(round(seq(1, d[1], length.out = min(96L, d[1])))))
    gy <- unique(as.integer(round(seq(1, d[2], length.out = min(96L, d[2])))))
    gz <- unique(as.integer(round(seq(1, d[3], length.out = min(96L, d[3])))))
    g <- a[gx, gy, gz, drop = FALSE]
    gm <- if (is.finite(thr)) is.finite(g) & g > thr else matrix(FALSE, nrow = length(g), ncol = 1L)
    gm <- array(gm, dim = dim(g))
    fg_fraction <- mean(gm)

    bbox_fraction <- center_offset <- NA_real_
    if (sum(gm) >= 20L) {
      coordinates <- arrayInd(which(gm), dim(gm))
      spans <- apply(coordinates, 2L, function(z) max(z) - min(z) + 1L)
      bbox_fraction <- prod(spans / dim(gm))
      centre <- colMeans(coordinates) / dim(gm)
      center_offset <- sqrt(sum((centre - 0.5)^2))
    }

    # Compute slice foreground fractions on at most 96 evenly spaced slices.
    slice_id <- unique(as.integer(round(seq(1, d[3], length.out = min(96L, d[3])))))
    slice_cov <- vapply(slice_id, function(k) {
      z <- a[, , k]
      if (!is.finite(thr)) return(NA_real_)
      mean(is.finite(z) & z > thr)
    }, numeric(1))
    slice_cov <- slice_cov[is.finite(slice_cov)]

# A low-resolution gradient summarizes image texture within each protocol.
    g_num <- as.numeric(g)
    grad_parts <- c(
      as.numeric(abs(g[-1L, , , drop = FALSE] - g[-dim(g)[1L], , , drop = FALSE])),
      as.numeric(abs(g[, -1L, , drop = FALSE] - g[, -dim(g)[2L], , drop = FALSE])),
      as.numeric(abs(g[, , -1L, drop = FALSE] - g[, , -dim(g)[3L], drop = FALSE]))
    )
    fg_values <- values[finite & is.finite(thr) & values > thr]

    result[, `:=`(
      metric_status = "OK",
      dim1 = as.integer(d[1]), dim2 = as.integer(d[2]), dim3 = as.integer(d[3]),
      finite_fraction = mean(finite),
      nonzero_fraction = mean(finite & values != 0),
      intensity_p01 = q[1], intensity_p50 = q[3], intensity_p99 = q[5],
      intensity_iqr = q[4] - q[2],
      intensity_entropy = entropy_from_values(values),
      foreground_fraction = fg_fraction,
      foreground_intensity_cv = safe_sd(fg_values) / pmax(1e-12, safe_mean(fg_values)),
      foreground_bbox_fraction = bbox_fraction,
      foreground_center_offset = center_offset,
      slice_coverage_iqr = if (length(slice_cov) >= 2L) stats::IQR(slice_cov) else NA_real_,
      slice_dropout_fraction = if (length(slice_cov)) mean(slice_cov < 0.01) else NA_real_,
      low_resolution_gradient = safe_mean(grad_parts[is.finite(grad_parts)])
    )]
    invisible(g_num)
  }, error = function(e) {
    result[, metric_error := conditionMessage(e)]
  })
  result
}

message("13F step 1 started: automatic raw-T1 QC and processing-manifest freeze.")
qc <- fread(need(input_file), showProgress = FALSE)
required_columns <- c("cohort", "role", "image_id", "participant_id", "nifti_path", "conversion_pass")
missing_columns <- setdiff(required_columns, names(qc))
if (length(missing_columns)) stop("13D input is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)

qc[, cohort := safe_chr(cohort)]
qc[, role := safe_chr(role)]
qc[, image_id := as.character(image_id)]
qc[, participant_id := as.character(participant_id)]
qc[, nifti_path := norm_path(nifti_path)]
qc[, conversion_pass_logical := is_true(conversion_pass)]
qc[, nifti_exists := file.exists(nifti_path)]

# Frozen main sample: ADNI primary longitudinal scans and AIBL primary baseline
# scans. ADNI QC_BACKUP scans are retained in 13D but cannot enter a main model
# unless a separately frozen replacement decision is made after processing QC.
selected <- qc[
  conversion_pass_logical & nifti_exists &
    ((cohort == "ADNI" & role == "PRIMARY") | (cohort == "AIBL" & role == "PRIMARY"))
]
if (!nrow(selected)) stop("No conversion-passing primary images were available for 13F.", call. = FALSE)
if (anyDuplicated(selected[, .(cohort, image_id)])) stop("Duplicate cohort/image_id entries in the 13F primary manifest.", call. = FALSE)

if (!all(c("ADNI", "AIBL") %in% selected$cohort)) stop("13F requires both ADNI and AIBL primary images.", call. = FALSE)
if (selected[cohort == "ADNI", .N] != 1278L) warning("ADNI primary image count is not the frozen expected value of 1278.")
if (selected[cohort == "AIBL", .N] != 336L) warning("AIBL primary image count is not the expected value after 13D exclusions (336).")

setorderv(selected, c("cohort", "participant_id", "window_position", "image_id"), na.last = TRUE)
selected[, processing_image_code := sprintf("PROC_%05d", seq_len(.N))]
selected[, processing_subject_code := sprintf("SUBJ_%05d", match(paste(cohort, participant_id), unique(paste(cohort, participant_id))))]
selected[, analysis_timepoint := fifelse(
  cohort == "ADNI", paste0("ADNI_", safe_chr(window_position)), "AIBL_BASELINE"
)]

# Use sufficiently populated scanner/protocol strata for automatic comparisons.
# Sparse protocol cells fall back to cohort × field strength, then cohort only.
selected[, scanner_protocol := paste(
  safe_chr(json_manufacturer), safe_chr(json_model), safe_chr(json_field_strength),
  sep = " | "
)]
selected[, field_protocol := paste(cohort, safe_chr(json_field_strength), sep = " | ")]
selected[, n_scanner_protocol := .N, by = scanner_protocol]
selected[, n_field_protocol := .N, by = field_protocol]
selected[, qc_stratum := fifelse(
  n_scanner_protocol >= MIN_STRATUM_N, scanner_protocol,
  fifelse(n_field_protocol >= MIN_STRATUM_N, field_protocol, cohort)
)]

manifest <- selected[, .(
  processing_image_code, processing_subject_code, cohort, role, image_id, RID,
  participant_id, window_position, analysis_timepoint, nifti_path, nifti_md5,
  json_path, json_manufacturer, json_model, json_field_strength,
  json_series_description, dim1, dim2, dim3, voxel1, voxel2, voxel3,
  qc_stratum, scanner_protocol, field_protocol
)]
write_csv(manifest, "01_13F_frozen_primary_processing_manifest.csv")

checkpoint_path <- file.path(out_dir, "02_raw_T1_automated_metrics_checkpoint.csv")
old <- if (file.exists(checkpoint_path) && !FORCE_RECOMPUTE) {
  fread(checkpoint_path, showProgress = FALSE)
} else data.table()
if (nrow(old) && !"processing_image_code" %in% names(old)) old <- data.table()
if (nrow(old)) old <- old[processing_image_code %in% manifest$processing_image_code]
# Resume completed metrics and retry metric-read errors.
if (nrow(old) && "metric_status" %in% names(old)) {
  retry_codes <- old[metric_status != "OK", processing_image_code]
  old <- old[!processing_image_code %in% retry_codes]
}

unprocessed_images <- manifest[!processing_image_code %in% old$processing_image_code]
rows <- vector("list", nrow(unprocessed_images))
if (nrow(unprocessed_images)) {
  for (i in seq_len(nrow(unprocessed_images))) {
    w <- unprocessed_images[i]
    rows[[i]] <- cbind(w[, .(processing_image_code)], read_volume_metrics(w$nifti_path))
    if (i %% CHECKPOINT_EVERY == 0L || i == nrow(unprocessed_images)) {
      current <- rbindlist(c(list(old), rows[seq_len(i)]), fill = TRUE)
      fwrite(current, checkpoint_path, bom = TRUE, na = "")
      message("13F automatic raw QC: ", i, "/", nrow(unprocessed_images), " newly read; errors=", current[metric_status != "OK", .N])
    }
  }
}
metrics <- rbindlist(c(list(old), rows), fill = TRUE)
metrics <- metrics[match(manifest$processing_image_code, processing_image_code)]
if (nrow(metrics) != nrow(manifest)) stop("Automated-metric checkpoint does not match frozen manifest.", call. = FALSE)
full <- merge(manifest, metrics, by = "processing_image_code", all.x = TRUE, sort = FALSE)
setorderv(full, "processing_image_code")

metric_columns <- c(
  "intensity_iqr", "intensity_entropy", "foreground_fraction",
  "foreground_intensity_cv", "foreground_bbox_fraction",
  "foreground_center_offset", "slice_coverage_iqr",
  "slice_dropout_fraction", "low_resolution_gradient"
)
for (metric in metric_columns) {
  stat_name <- paste0(metric, "_robust_z")
  stats <- full[metric_status == "OK" & is.finite(get(metric)), .(
    stratum_n = .N,
    center = stats::median(get(metric), na.rm = TRUE),
    spread = safe_mad(get(metric))
  ), by = qc_stratum]
  full <- merge(full, stats, by = "qc_stratum", all.x = TRUE, sort = FALSE)
  full[, (stat_name) := fifelse(
    is.finite(get(metric)) & stratum_n >= MIN_STRATUM_N & is.finite(spread) & spread > 0,
    abs(get(metric) - center) / spread,
    NA_real_
  )]
  full[, c("stratum_n", "center", "spread") := NULL]
}
z_columns <- paste0(metric_columns, "_robust_z")
full[, n_robust_metric_flags := rowSums(as.matrix(.SD) > ROBUST_Z_THRESHOLD, na.rm = TRUE), .SDcols = z_columns]
full[, hard_technical_flag :=
       metric_status != "OK" |
       !is.finite(finite_fraction) | finite_fraction < 0.999 |
       !is.finite(nonzero_fraction) | nonzero_fraction < 0.01 |
       !is.finite(intensity_iqr) | intensity_iqr <= 0 |
       !is.finite(foreground_fraction) | foreground_fraction < 0.02 | foreground_fraction > 0.95]
full[, auto_QC_status := fifelse(
  hard_technical_flag | n_robust_metric_flags >= 2L,
  "AUTO_FLAGGED_BLINDED_REVIEW", "AUTO_CLEAR"
)]
full[, auto_QC_reason := vapply(seq_len(.N), function(i) {
  reasons <- character()
  if (metric_status[i] != "OK") reasons <- c(reasons, "volume_read_error")
  if (is.finite(finite_fraction[i]) && finite_fraction[i] < 0.999) reasons <- c(reasons, "nonfinite_voxels")
  if (is.finite(nonzero_fraction[i]) && nonzero_fraction[i] < 0.01) reasons <- c(reasons, "near_empty_volume")
  if (is.finite(intensity_iqr[i]) && intensity_iqr[i] <= 0) reasons <- c(reasons, "zero_intensity_range")
  if (is.finite(foreground_fraction[i]) && (foreground_fraction[i] < 0.02 || foreground_fraction[i] > 0.95)) reasons <- c(reasons, "implausible_foreground_proxy")
  z_value <- as.numeric(unlist(full[i, z_columns, with = FALSE], use.names = FALSE))
  hits <- metric_columns[is.finite(z_value) & z_value > ROBUST_Z_THRESHOLD]
  if (length(hits)) reasons <- c(reasons, paste0("robust_outlier:", hits))
  if (!length(reasons)) "none" else paste(reasons, collapse = ";")
}, character(1))]

write_csv(full, "03_automated_raw_T1_QC_metrics_and_flags.csv")
write_csv(full[auto_QC_status == "AUTO_FLAGGED_BLINDED_REVIEW"], "04_automatic_QC_flagged_images.csv")

flagged <- full[auto_QC_status == "AUTO_FLAGGED_BLINDED_REVIEW"]

failed_13d <- qc[!conversion_pass_logical | !nifti_exists, .(
  cohort, role, image_id, RID, participant_id, nifti_path, conversion_pass, nifti_exists,
  n_nifti_outputs, header_error
)]
write_csv(failed_13d, "06_excluded_before_13F_from_13D_technical_QC.csv")

gates <- data.table(
  gate = c(
    "FROZEN_PRIMARY_PROCESSING_MANIFEST_READY",
    "ADNI_PRIMARY_LONGITUDINAL_INPUT_COMPLETE",
    "AIBL_PRIMARY_BASELINE_INPUT_COMPLETE_AFTER_13D_EXCLUSIONS",
    "AUTOMATED_RAW_T1_QC_COMPLETE",
    "NO_AUTOMATIC_OUTLIERS"
  ),
  pass = c(
    nrow(manifest) == 1614L,
    manifest[cohort == "ADNI", .N] == 1278L,
    manifest[cohort == "AIBL", .N] == 336L,
    nrow(full) == nrow(manifest) && all(full$metric_status == "OK"),
    nrow(flagged) == 0L
  ),
  n = c(nrow(manifest), manifest[cohort == "ADNI", .N], manifest[cohort == "AIBL", .N],
        nrow(full), nrow(flagged)),
  interpretation = c(
    "Only conversion-passing ADNI PRIMARY and AIBL PRIMARY T1 images are admitted; ADNI QC_BACKUP images are not used in the main processing manifest.",
    "Expected 639 ADNI participants x 2 retained primary scans.",
    "Expected 338 AIBL baseline candidates minus the two 13D multiple-output technical failures.",
    "Metrics summarize protocol-stratified image anomalies.",
    "TRUE indicates zero automatic outliers. FALSE identifies images for inspection."
  )
)
write_csv(gates, "07_13F_step1_readiness_gates.csv")

runtime_stamp <- c(
  paste0("R.version=", R.version.string),
  paste0("platform=", R.version$platform),
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("data.table=", as.character(utils::packageVersion("data.table"))),
  paste0("RNifti=", as.character(utils::packageVersion("RNifti"))),
  paste0("input=", input_file),
  paste0("retained_primary_T1=", nrow(manifest)),
  paste0("automatic_flags=", nrow(flagged))
)
writeLines(runtime_stamp, file.path(out_dir, "99_runtime_stamp.txt"), useBytes = TRUE)

message("13F automatic raw-T1 QC complete: ", nrow(manifest), " primary images; ", nrow(flagged), " flagged for limited blinded review.")
