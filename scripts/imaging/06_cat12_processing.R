# Compiled CAT12 batch generation, execution and XML/ROI parsing
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
# Use ASCII staging paths for compiled CAT12.
suppressPackageStartupMessages({
  library(data.table)
  library(processx)
  library(xml2)
  library(survival)
})

# CAT12 runner for AIBL baseline images and the ADNI calibration sample.
# Modes are audit, prepare, run and parse.
# Set CAT12_R_WORKERS for process concurrency.
# CAT12_R_LIMIT = 0 selects the complete configured image set.

env_or_option <- function(env_name, option_name, default) {
  z <- Sys.getenv(env_name, unset = "")
  if (nzchar(z)) z else getOption(option_name, default)
}
MODE <- as.character(env_or_option("CAT12_R_MODE", "CAT12_R_MODE", "audit"))
COHORT <- as.character(env_or_option("CAT12_R_COHORT", "CAT12_R_COHORT", "AIBL"))
WORKERS <- as.integer(env_or_option("CAT12_R_WORKERS", "CAT12_R_WORKERS", 2L))
LIMIT <- as.integer(env_or_option("CAT12_R_LIMIT", "CAT12_R_LIMIT", 0L))
TIMEOUT_MIN <- as.numeric(env_or_option("CAT12_R_TIMEOUT_MIN", "CAT12_R_TIMEOUT_MIN", 120))
SURFACE_OUTPUT <- as.logical(env_or_option("CAT12_R_SURFACE", "CAT12_R_SURFACE", FALSE))
KEEP_VBM_MAPS <- as.logical(env_or_option("CAT12_R_KEEP_VBM_MAPS", "CAT12_R_KEEP_VBM_MAPS", FALSE))
ADNI_CALIBRATION_SUBJECTS <- as.integer(env_or_option("CAT12_R_ADNI_CALIBRATION_SUBJECTS", "CAT12_R_ADNI_CALIBRATION_SUBJECTS", 25L))

if (!MODE %in% c("audit", "prepare", "run", "parse")) {
  stop("CAT12_R_MODE must be audit, prepare, run or parse.", call. = FALSE)
}
if (!COHORT %in% c("AIBL", "ADNI_CALIBRATION", "ADNI")) {
  stop("CAT12_R_COHORT must be AIBL, ADNI_CALIBRATION or ADNI.", call. = FALSE)
}
if (!is.finite(WORKERS) || WORKERS < 1L || WORKERS > 4L) {
  stop("Use 1-4 CAT12 processes; start with 2.", call. = FALSE)
}

PROJECT_ROOT <- study_path("project")
CAT_ROOT <- study_path("cat12")
CAT_EXE <- study_path("cat12_executable")
CAT_TEMPLATE <- file.path(CAT_ROOT, "standalone", "cat_standalone_segment.m")
MCR_ROOT <- study_path("matlab_runtime")
MANIFEST_CANDIDATES <- c(
  file.path(PROJECT_ROOT, "results", "13F_automated_raw_T1_QC_and_processing_manifest",
            "01_13F_frozen_primary_processing_manifest.csv"),
  study_path("scratch", "13F_aibl_manifest.csv")
)
# Allow an explicit manifest path from the caller. This avoids failures when
# a mounted drive is visible in Windows Explorer but not resolved identically
# by the R session. The override is deliberately checked before the defaults.
MANIFEST_OVERRIDE <- Sys.getenv("CAT12_R_MANIFEST", unset = "")
if (nzchar(MANIFEST_OVERRIDE)) {
  MANIFEST <- normalizePath(MANIFEST_OVERRIDE, winslash = "/", mustWork = FALSE)
} else {
  MANIFEST <- MANIFEST_CANDIDATES[file.exists(MANIFEST_CANDIDATES)][1L]
}
RESULT_ROOT <- env_or_option(
  "CAT12_R_RESULT_ROOT", "CAT12_R_RESULT_ROOT",
  file.path(PROJECT_ROOT, "results", "13F_CAT12_R_direct_minimal")
)
STAGE_ROOT <- env_or_option(
  "CAT12_R_STAGE_ROOT", "CAT12_R_STAGE_ROOT",
  study_path("cat12_stage")
)
XML_ARCHIVE <- file.path(STAGE_ROOT, "xml")
BATCH_ROOT <- file.path(STAGE_ROOT, "batches")
RUNTIME_LOG_ROOT <- file.path(STAGE_ROOT, "logs")
RUNTIME_RECORD_ROOT <- file.path(STAGE_ROOT, "records")
SUMMARY_ROOT <- file.path(RESULT_ROOT, "tables")

LOG_ROOT <- RUNTIME_LOG_ROOT
RECORD_ROOT <- RUNTIME_RECORD_ROOT
for (p in c(RESULT_ROOT, XML_ARCHIVE, LOG_ROOT, RECORD_ROOT, BATCH_ROOT, SUMMARY_ROOT,
            STAGE_ROOT)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

need <- function(path, label = path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("Missing ", label, ": ", path, call. = FALSE)
  }
  path
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
norm_id <- function(x) {
  z <- sub("^.*?([0-9]+)$", "\\1", as.character(x))
  z <- sub("^0+", "", z); z[z == ""] <- "0"; z
}
md5 <- function(path) unname(tools::md5sum(path))
matlab_path <- function(path) gsub("\\\\", "/", path)
escape_regex <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)

set_matlab_field <- function(lines, field, value, only_if_present = FALSE) {
  pat <- paste0("^\\s*matlabbatch\\{1\\}\\.", escape_regex(field), "\\s*=")
  ix <- grep(pat, lines, perl = TRUE)
  repl <- paste0("matlabbatch{1}.", field, " = ", value, ";")
  if (length(ix)) {
    lines[ix[1L]] <- repl
  } else if (!only_if_present) {
    lines <- c(lines, repl)
  }
  lines
}

read_manifest <- function() {
  x <- fread(file = need(MANIFEST, "frozen manifest"), showProgress = FALSE)
  # Normalize ADNI PTID to the common participant_id field.
  if (!"participant_id" %in% names(x)) {
    participant_col <- intersect(c("PTID", "ptid", "subject_id", "participant"), names(x))
    if (!length(participant_col)) {
      stop("Manifest must contain participant_id or PTID.", call. = FALSE)
    }
    x[, participant_id := as.character(get(participant_col[1L]))]
  }
  x[, `:=`(cohort = as.character(cohort), role = as.character(role),
           RID = norm_id(RID), processing_image_code = as.character(processing_image_code),
           nifti_path = as.character(nifti_path))]
  x <- x[role == "PRIMARY"]
  # The frozen full manifest contains 1,614 rows. A caller may provide a
  # deliberately smaller, explicit bridge-subset manifest for ADNI; in that
  # case only the non-empty PRIMARY requirement is enforced.
  if (nzchar(MANIFEST_OVERRIDE)) {
    if (!nrow(x)) stop("Explicit manifest contains no PRIMARY rows.", call. = FALSE)
  } else if (nrow(x) != 1614L) {
    stop("Frozen PRIMARY manifest must contain 1,614 rows; found ", nrow(x), call. = FALSE)
  }
  x
}

select_adni_calibration <- function(x, n_subjects) {
  a <- x[cohort == "ADNI"]
  counts <- a[, .(n_dates = uniqueN(analysis_timepoint)), by = RID]
  eligible <- counts[n_dates >= 2L, RID]
  eligible <- sort(eligible)
  if (!length(eligible)) stop("No paired ADNI calibration subjects found.", call. = FALSE)
  n_subjects <- min(n_subjects, length(eligible))
  # Deterministic evenly spaced selection across the frozen RID range.
  ix <- unique(round(seq(1, length(eligible), length.out = n_subjects)))
  chosen <- eligible[ix]
  a[RID %in% chosen][order(RID, analysis_timepoint)]
}

select_rows <- function(x) {
  if (COHORT == "AIBL") out <- x[cohort == "AIBL"][order(RID, image_id)]
  else if (COHORT == "ADNI") out <- x[cohort == "ADNI"][order(RID, analysis_timepoint)]
  else out <- select_adni_calibration(x, ADNI_CALIBRATION_SUBJECTS)
  if (LIMIT > 0L) out <- out[seq_len(min(LIMIT, nrow(out)))]
  if (!nrow(out)) stop("No rows selected for cohort ", COHORT, call. = FALSE)
  out[, job_id := sprintf("%s_%s", COHORT, processing_image_code)]
  out
}

runtime_environment <- function(cache_dir, temp_dir) {
  e <- Sys.getenv()
  runtime_dirs <- c(file.path(MCR_ROOT, "runtime", "win64"),
                    file.path(MCR_ROOT, "bin", "win64"),
                    file.path(MCR_ROOT, "sys", "os", "win64"),
                    file.path(MCR_ROOT, "sys", "opengl", "lib", "win64"))
  runtime_dirs <- runtime_dirs[dir.exists(runtime_dirs)]
  e["PATH"] <- paste(c(runtime_dirs, e[["PATH"]]), collapse = ";")
  e["MCR_CACHE_ROOT"] <- cache_dir
  e["TMP"] <- temp_dir
  e["TEMP"] <- temp_dir
  e
}

make_batch <- function(stage_input, batch_path) {
  lines <- readLines(need(CAT_TEMPLATE, "CAT12 batch template"), encoding = "UTF-8", warn = FALSE)
  lines <- lines[!grepl("<UNDEFINED>", lines, fixed = TRUE)]
  lines <- c(paste0("matlabbatch{1}.spm.tools.cat.estwrite.data = {'",
                    matlab_path(stage_input), "'};"), lines)

# Run voxel segmentation with the fixed atlas and configured surface outputs.
  lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.output.surface", if (SURFACE_OUTPUT) "1" else "0")
  lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.extopts.admin.print", "0")
  lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.output.label.native", "1")
  lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics", "1")
  atlas_fields <- c(
    "lpba40", "cobra", "hammers", "ibsr", "aal3", "mori", "thalamus", "anatomy3",
    "julichbrain", "Schaefer2018_100Parcels_17Networks_order",
    "Schaefer2018_200Parcels_17Networks_order",
    "Schaefer2018_400Parcels_17Networks_order",
    "Schaefer2018_600Parcels_17Networks_order"
  )
  for (a in atlas_fields) {
    lines <- set_matlab_field(lines, paste0("spm.tools.cat.estwrite.output.ROImenu.atlases.", a), "0", only_if_present = TRUE)
  }
  if (!KEEP_VBM_MAPS) {
    for (tissue in c("GM", "WM", "CSF", "WMH")) {
      for (out in c("native", "warped", "mod", "dartel")) {
        lines <- set_matlab_field(lines, paste0("spm.tools.cat.estwrite.output.", tissue, ".", out), "0", only_if_present = TRUE)
      }
    }
    lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.output.bias.native", "0", only_if_present = TRUE)
    lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.output.bias.warped", "0", only_if_present = TRUE)
    lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.output.warps", "[0 0]", only_if_present = TRUE)
    lines <- set_matlab_field(lines, "spm.tools.cat.estwrite.output.jacobianwarped", "0", only_if_present = TRUE)
  }
  writeLines(lines, batch_path, useBytes = TRUE)
  batch_path
}

stage_job <- function(row) {
  job_dir <- file.path(STAGE_ROOT, row$job_id)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  input_name <- basename(row$nifti_path)
  stage_input <- file.path(job_dir, input_name)
  if (!file.exists(stage_input)) {
    ok <- tryCatch(file.link(row$nifti_path, stage_input), warning = function(w) FALSE, error = function(e) FALSE)
    if (!isTRUE(ok) || !file.exists(stage_input)) {
      ok <- file.copy(row$nifti_path, stage_input, overwrite = FALSE)
      if (!isTRUE(ok)) stop("Could not stage input: ", row$nifti_path, call. = FALSE)
    }
  }
  batch_path <- file.path(BATCH_ROOT, paste0(row$job_id, ".m"))
  make_batch(stage_input, batch_path)
  list(job_dir = job_dir, stage_input = stage_input, batch_path = batch_path)
}

xml_paths <- function(stage_input) {
  stem <- tools::file_path_sans_ext(basename(stage_input))
  if (grepl("\\.nii$", stem, ignore.case = TRUE)) stem <- tools::file_path_sans_ext(stem)
  root <- dirname(stage_input)
  list(
    report = file.path(root, "report", paste0("cat_", stem, ".xml")),
    roi = file.path(root, "label", paste0("catROI_", stem, ".xml")),
    surface_roi = file.path(root, "label", paste0("catROIs_", stem, ".xml"))
  )
}

is_complete <- function(stage_input) {
  p <- xml_paths(stage_input)
  file.exists(p$report) && if (SURFACE_OUTPUT) file.exists(p$surface_roi) else file.exists(p$roi)
}

run_one <- function(row, worker_id) {
  staged <- stage_job(row)
  paths <- xml_paths(staged$stage_input)
  rec_path <- file.path(RECORD_ROOT, paste0(row$job_id, ".rds"))
  if (is_complete(staged$stage_input)) {
    return(data.table(job_id = row$job_id, worker = worker_id, status = "ALREADY_COMPLETE",
                      started = NA_character_, finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                      elapsed_minutes = 0, exit_code = 0L, input_path = row$nifti_path,
                      stage_input = staged$stage_input, batch_path = staged$batch_path,
    report_xml = paths$report, roi_xml = paths$roi, surface_roi_xml = paths$surface_roi))
  }
  worker_root <- file.path(STAGE_ROOT, "workers", sprintf("worker_%02d", worker_id))
  cache_dir <- file.path(worker_root, "mcr_cache")
  temp_dir <- file.path(worker_root, "tmp")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  stdout <- file.path(LOG_ROOT, paste0(row$job_id, ".stdout.log"))
  stderr <- file.path(LOG_ROOT, paste0(row$job_id, ".stderr.log"))
  started <- Sys.time()
  p <- process$new(CAT_EXE, c("batch", staged$batch_path),
                   stdout = stdout, stderr = stderr,
                   wd = CAT_ROOT, env = runtime_environment(cache_dir, temp_dir),
                   windows_hide_window = TRUE, cleanup_tree = FALSE)
  finished <- FALSE
  while (p$is_alive()) {
    p$wait(timeout = min(30000L, as.integer(TIMEOUT_MIN * 60 * 1000)))
    if (as.numeric(difftime(Sys.time(), started, units = "mins")) > TIMEOUT_MIN) {
      p$kill_tree(); break
    }
  }
  exit_code <- p$get_exit_status()
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  done <- is_complete(staged$stage_input)
  status <- if (done) "SUCCESS" else if (elapsed > TIMEOUT_MIN) "TIMEOUT" else paste0("INCOMPLETE_EXIT_", ifelse(is.na(exit_code), "NA", exit_code))
  result <- data.table(
    job_id = row$job_id, worker = worker_id, status = status,
    started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_minutes = elapsed, exit_code = exit_code,
    input_path = row$nifti_path, input_md5 = md5(row$nifti_path),
    stage_input = staged$stage_input, batch_path = staged$batch_path,
    stdout_log = stdout, stderr_log = stderr, report_xml = paths$report, roi_xml = paths$roi,
    surface_roi_xml = paths$surface_roi,
    cat_exe = CAT_EXE, mcr_root = MCR_ROOT
  )
  saveRDS(result, rec_path)
  result
}

run_queue <- function(jobs) {
  # CAT12 standalone is not internally parallelized.  This is explicitly
  # bounded process-level concurrency, with a separate MCR cache/temp root per
  # worker.  Start with two benchmark images and inspect XML completion before
  # increasing to four.
  remaining <- seq_len(nrow(jobs)); active <- list(); results <- list(); worker_ids <- seq_len(WORKERS)
  next_worker <- 1L
  while (length(remaining) || length(active)) {
    while (length(remaining) && length(active) < WORKERS) {
      idx <- remaining[1L]; remaining <- remaining[-1L]
      worker <- worker_ids[next_worker]; next_worker <- if (next_worker == WORKERS) 1L else next_worker + 1L
      job <- jobs[idx]
      # Run each CAT12 worker through a clean R subprocess so its MCR cache,
      # launch and timeout are isolated; this outer R process is the scheduler.
      worker_script <- file.path(BATCH_ROOT, paste0(job$job_id, "_worker.R"))
      saveRDS(job, paste0(worker_script, ".job.rds"))
      worker_lines <- c(
        "source('scripts/imaging/06_cat12_processing.R', encoding='UTF-8')"
      )
      # The scheduler does not use this recursive path; direct launch is used
      # below to avoid re-running the global script.
      active[[as.character(idx)]] <- list(index = idx, worker = worker, started = Sys.time())
      # A direct external CAT process is launched here, while the R loop keeps
      # worker-specific state.  Stage/batch preparation is deterministic.
      staged <- stage_job(job)
      paths <- xml_paths(staged$stage_input)
      if (is_complete(staged$stage_input)) {
        results[[length(results) + 1L]] <- data.table(job_id = job$job_id, worker = worker, status = "ALREADY_COMPLETE",
                                                       started = NA_character_, finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                                                       elapsed_minutes = 0, exit_code = 0L, input_path = job$nifti_path,
                                                       stage_input = staged$stage_input, batch_path = staged$batch_path,
                                                       report_xml = paths$report, roi_xml = paths$roi,
                                                       surface_roi_xml = paths$surface_roi)
        active[[as.character(idx)]] <- NULL
        next
      }
      worker_root <- file.path(STAGE_ROOT, "workers", sprintf("worker_%02d", worker))
      cache_dir <- file.path(worker_root, "mcr_cache")
      temp_dir <- file.path(worker_root, "tmp")
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE); dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
      stdout <- file.path(LOG_ROOT, paste0(job$job_id, ".stdout.log")); stderr <- file.path(LOG_ROOT, paste0(job$job_id, ".stderr.log"))
      proc <- process$new(CAT_EXE, c("batch", staged$batch_path), stdout = stdout, stderr = stderr, wd = CAT_ROOT,
                          env = runtime_environment(cache_dir, temp_dir), windows_hide_window = TRUE, cleanup_tree = FALSE)
      active[[as.character(idx)]]$process <- proc
      active[[as.character(idx)]]$job <- job
      active[[as.character(idx)]]$staged <- staged
      active[[as.character(idx)]]$paths <- paths
      active[[as.character(idx)]]$stdout <- stdout
      active[[as.character(idx)]]$stderr <- stderr
    }
    if (!length(active)) next
    Sys.sleep(2)
    for (key in names(active)) {
      a <- active[[key]]; elapsed <- as.numeric(difftime(Sys.time(), a$started, units = "mins"))
      alive <- a$process$is_alive()
      if (alive && elapsed <= TIMEOUT_MIN) next
      if (alive && elapsed > TIMEOUT_MIN) a$process$kill_tree()
      exit_code <- a$process$get_exit_status()
      complete <- is_complete(a$staged$stage_input)
      status <- if (complete) "SUCCESS" else if (elapsed > TIMEOUT_MIN) "TIMEOUT" else paste0("INCOMPLETE_EXIT_", ifelse(is.na(exit_code), "NA", exit_code))
      result <- data.table(job_id = a$job$job_id, worker = a$worker, status = status,
                           started = format(a$started, "%Y-%m-%d %H:%M:%S %Z"), finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                           elapsed_minutes = elapsed, exit_code = exit_code, input_path = a$job$nifti_path,
                           input_md5 = md5(a$job$nifti_path), stage_input = a$staged$stage_input,
                           batch_path = a$staged$batch_path, stdout_log = a$stdout, stderr_log = a$stderr,
                           report_xml = a$paths$report, roi_xml = a$paths$roi,
                           surface_roi_xml = a$paths$surface_roi, cat_exe = CAT_EXE, mcr_root = MCR_ROOT)
      saveRDS(result, file.path(RECORD_ROOT, paste0(a$job$job_id, ".rds")))
      results[[length(results) + 1L]] <- result
      active[[key]] <- NULL
    }
  }
  rbindlist(results, fill = TRUE)
}

# -----------------------------------------------------------------------------
# XML parsers.  CAT12 report/ROI XML is the primary quantitative interface;
# MAT files are retained by CAT12 but are not required for this R analysis.
# -----------------------------------------------------------------------------
parse_num_vector <- function(x) {
# Parse CAT12 numeric vectors after converting XML separators to whitespace.
  if (!length(x) || all(is.na(x))) return(numeric())
  txt <- paste(as.character(x), collapse = " ")
  txt <- gsub("[^0-9eE+.-]+", " ", txt, perl = TRUE)
  txt <- trimws(txt)
  if (!nzchar(txt)) return(numeric())
  vals <- suppressWarnings(as.numeric(strsplit(txt, "\\s+")[[1L]]))
  vals[is.finite(vals)]
}
xml_text1 <- function(doc, xpath) {
  node <- xml_find_first(doc, xpath)
  # xml2 versions differ in how a relative XPath is evaluated when `doc` is
  # itself a child node. Retry without the leading `./` and then with a
  # descendant search so nested CAT12 <data><Vgm>/<Vwm> nodes are recovered.
  if (inherits(node, "xml_missing") && startsWith(xpath, "./")) {
    node <- xml_find_first(doc, substring(xpath, 3L))
  }
  if (inherits(node, "xml_missing") && grepl("data/", xpath, fixed = TRUE)) {
    leaf <- sub(".*/", "", xpath)
    node <- xml_find_first(doc, paste0(".//*[local-name()='", leaf, "']"))
  }
  if (inherits(node, "xml_missing")) return(NA_character_)
  xml_text(node)
}
parse_report_xml <- function(path, job_id) {
  doc <- read_xml(need(path, "CAT12 report XML"))
  cgw <- parse_num_vector(xml_text1(doc, "/S/subjectmeasures/vol_abs_CGW"))
  data.table(
    job_id = job_id, report_xml = path,
    cat12_version = xml_text1(doc, "/S/parameter/extopts/cat12_version"),
    iqr = num(xml_text1(doc, "/S/qualityratings/IQR")),
    siqr = num(xml_text1(doc, "/S/qualityratings/SIQR")),
    ncr = num(xml_text1(doc, "/S/qualitymeasures/NCR")),
    icr = num(xml_text1(doc, "/S/qualitymeasures/ICR")),
    fec = num(xml_text1(doc, "/S/qualitymeasures/FEC")),
    resolution = num(xml_text1(doc, "/S/qualitymeasures/resolution")),
    csf_ml = if (length(cgw) >= 1L) cgw[1L] else NA_real_,
    gm_ml = if (length(cgw) >= 2L) cgw[2L] else NA_real_,
    wm_ml = if (length(cgw) >= 3L) cgw[3L] else NA_real_,
    wmh_ml = if (length(cgw) >= 4L) cgw[4L] else NA_real_,
    tiv_ml = num(xml_text1(doc, "/S/subjectmeasures/vol_TIV")),
    wmh_abs_ml = num(xml_text1(doc, "/S/subjectmeasures/vol_abs_WMH")),
    report_status = "PARSED"
  )
}
parse_roi_xml <- function(path, job_id) {
  doc <- read_xml(need(path, "CAT12 ROI XML"))
  atlas_nodes <- xml_children(xml_root(doc))
  xml_vector_from_node <- function(node, tag) {
    # CAT12 writes atlas vectors as literal text inside <data><Vgm>...</Vgm>
    # and <data><Vwm>...</Vwm>. Use XML lookup first, then a serialized-node
    # fallback so this remains robust across xml2 relative-XPath behavior.
    x <- xml_text1(node, paste0("./data/", tag))
    if (is.na(x) || !nzchar(trimws(x))) {
      raw <- as.character(node)
      hit <- regmatches(
        raw,
        regexpr(
          paste0("<", tag, ">.*?</", tag, ">"),
          raw,
          perl = TRUE
        )
      )
      if (length(hit) && !is.na(hit) && nzchar(hit)) {
        x <- sub(paste0("^.*?<", tag, ">"), "", hit, perl = TRUE)
        x <- sub(paste0("</", tag, ">.*$"), "", x, perl = TRUE)
      }
    }
    parse_num_vector(x)
  }
  out <- list()
  for (node in atlas_nodes) {
    atlas <- xml_name(node)
    names_vec <- xml_text(xml_find_all(node, "./names/item"))
    ids_vec <- parse_num_vector(xml_text1(node, "./ids"))
    vgm <- xml_vector_from_node(node, "Vgm")
    vwm <- xml_vector_from_node(node, "Vwm")
    # The XML node serialization can expose the values, but the parser's
    # numeric fallback may still fail on long semicolon vectors. If so,
    # extract directly from the document text using the atlas tag boundary.
    if (!length(vgm)) {
      whole <- paste(as.character(doc), collapse = "")
      hit <- regmatches(whole, regexpr("<neuromorphometrics>.*?</neuromorphometrics>", whole, perl = TRUE))
      if (length(hit) && !is.na(hit)) {
        vv <- regmatches(hit, regexpr("<Vgm>.*?</Vgm>", hit, perl = TRUE))
        if (length(vv) && !is.na(vv)) vgm <- parse_num_vector(sub("^.*?<Vgm>|</Vgm>.*$", "", vv, perl = TRUE))
      }
    }
    if (!length(vwm)) {
      whole <- paste(as.character(doc), collapse = "")
      hit <- regmatches(whole, regexpr("<neuromorphometrics>.*?</neuromorphometrics>", whole, perl = TRUE))
      if (length(hit) && !is.na(hit)) {
        vv <- regmatches(hit, regexpr("<Vwm>.*?</Vwm>", hit, perl = TRUE))
        if (length(vv) && !is.na(vv)) vwm <- parse_num_vector(sub("^.*?<Vwm>|</Vwm>.*$", "", vv, perl = TRUE))
      }
    }
    thickness <- parse_num_vector(xml_text1(node, "./data/thickness"))
    if (!length(thickness)) thickness <- parse_num_vector(xml_text1(node, "./thickness"))
    gyrification <- parse_num_vector(xml_text1(node, "./data/gyrification"))
    n <- max(length(names_vec), length(ids_vec), length(vgm), length(vwm), length(thickness), length(gyrification))
    if (!n) next
    out[[length(out) + 1L]] <- data.table(
      job_id = job_id, roi_xml = path, atlas = atlas,
      roi_id = if (length(ids_vec) == n) ids_vec else rep(NA_real_, n),
      roi_name = if (length(names_vec) == n) names_vec else rep(NA_character_, n),
      Vgm_ml = if (length(vgm) == n) vgm else rep(NA_real_, n),
      Vwm_ml = if (length(vwm) == n) vwm else rep(NA_real_, n),
      thickness_mm = if (length(thickness) == n) thickness else rep(NA_real_, n),
      gyrification = if (length(gyrification) == n) gyrification else rep(NA_real_, n),
      measure_type = if (length(thickness) == n && any(is.finite(thickness))) "surface_or_volume" else "volume"
    )
  }
  rbindlist(out, fill = TRUE)
}

parse_existing <- function(jobs = NULL) {
  if (is.null(jobs)) {
    recs <- list.files(RECORD_ROOT, pattern = "\\.rds$", full.names = TRUE)
    if (length(recs)) jobs <- rbindlist(lapply(recs, readRDS), fill = TRUE)
  }
  if (is.null(jobs) || !nrow(jobs)) stop("No CAT12 records available for parsing.", call. = FALSE)
  jobs <- as.data.table(jobs)
  jobs <- jobs[status %in% c("SUCCESS", "ALREADY_COMPLETE")]
  report_list <- list(); roi_list <- list()
  for (i in seq_len(nrow(jobs))) {
    roi_path <- if (SURFACE_OUTPUT && file.exists(jobs$surface_roi_xml[i])) jobs$surface_roi_xml[i] else jobs$roi_xml[i]
    if (!file.exists(jobs$report_xml[i]) || !file.exists(roi_path)) next
    report_list[[length(report_list) + 1L]] <- parse_report_xml(jobs$report_xml[i], jobs$job_id[i])
    roi_list[[length(roi_list) + 1L]] <- parse_roi_xml(roi_path, jobs$job_id[i])
  }
  if (!length(report_list)) stop("No complete report/ROI XML pairs found.", call. = FALSE)
  report <- rbindlist(report_list, fill = TRUE)
  roi <- rbindlist(roi_list, fill = TRUE)
  manifest <- read_manifest()
  manifest_meta <- data.table(
    processing_image_code = manifest[["processing_image_code"]],
    manifest_cohort = manifest[["cohort"]],
    manifest_RID = manifest[["RID"]],
    manifest_analysis_timepoint = manifest[["analysis_timepoint"]],
    manifest_image_id = manifest[["image_id"]],
    manifest_participant_id = manifest[["participant_id"]],
    manifest_nifti_path = manifest[["nifti_path"]]
  )
  # Parsed and input records may share
  # metadata columns. Remove those conflicting names before the manifest
  # join, then restore one canonical field set.
  drop_meta <- intersect(
    c("cohort", "cohort.x", "cohort.y", "RID", "RID.x", "RID.y",
      "analysis_timepoint", "image_id", "participant_id", "nifti_path"),
    names(jobs)
  )
  jobs_clean <- copy(jobs)
  if (length(drop_meta)) jobs_clean[, (drop_meta) := NULL]
  jobs_clean[, processing_image_code := sub("^[^_]+_", "", job_id)]
  jobs2 <- merge(jobs_clean, manifest_meta, by = "processing_image_code",
                 all.x = TRUE, sort = FALSE)
  jobs2[, `:=`(
    cohort = as.character(manifest_cohort),
    RID = as.character(manifest_RID),
    analysis_timepoint = as.character(manifest_analysis_timepoint),
    image_id = as.character(manifest_image_id),
    participant_id = as.character(manifest_participant_id),
    nifti_path = as.character(manifest_nifti_path)
  )]
  jobs2[, c("manifest_cohort", "manifest_RID", "manifest_analysis_timepoint",
            "manifest_image_id", "manifest_participant_id", "manifest_nifti_path") := NULL]
  report <- merge(report, unique(jobs2[, .(job_id, cohort, RID, analysis_timepoint, image_id, participant_id, nifti_path)]), by = "job_id", all.x = TRUE)
  roi <- merge(roi, unique(jobs2[, .(job_id, cohort, RID, analysis_timepoint, image_id, participant_id)]), by = "job_id", all.x = TRUE)
  fwrite(report, file.path(SUMMARY_ROOT, "01_CAT12_report_metrics.csv"), bom = TRUE)
  fwrite(roi, file.path(SUMMARY_ROOT, "02_CAT12_roi_values.csv"), bom = TRUE)
  fwrite(jobs2, file.path(SUMMARY_ROOT, "03_CAT12_job_manifest.csv"), bom = TRUE)
  list(report = report, roi = roi, jobs = jobs2)
}

derive_aibl_endpoint <- function(parsed) {
  report <- as.data.table(parsed$report)
  roi <- as.data.table(parsed$roi)

  # Some resumed .rds records contain only job_id and XML paths. In that
  # case the manifest merge can leave cohort absent or NA. Recover it from
  # the deterministic job_id prefix before applying the AIBL filter.
  if (!"cohort" %in% names(report)) report[["cohort"]] <- sub("_.*$", "", report[["job_id"]])
  if (!"cohort" %in% names(roi)) roi[["cohort"]] <- sub("_.*$", "", roi[["job_id"]])
  report[["cohort"]][is.na(report[["cohort"]]) | !nzchar(report[["cohort"]])] <-
    sub("_.*$", "", report[["job_id"]][is.na(report[["cohort"]]) | !nzchar(report[["cohort"]])])
  roi[["cohort"]][is.na(roi[["cohort"]]) | !nzchar(roi[["cohort"]])] <-
    sub("_.*$", "", roi[["job_id"]][is.na(roi[["cohort"]]) | !nzchar(roi[["cohort"]])])

  report <- report[report[["cohort"]] == "AIBL"]
  roi <- roi[roi[["cohort"]] == "AIBL" & roi[["atlas"]] == "neuromorphometrics"]
  targets <- c("Right Ent entorhinal area", "Left Ent entorhinal area",
               "Right FuG fusiform gyrus", "Left FuG fusiform gyrus",
               "Right ITG inferior temporal gyrus", "Left ITG inferior temporal gyrus",
               "Right MTG middle temporal gyrus", "Left MTG middle temporal gyrus",
               "Right TMP temporal pole", "Left TMP temporal pole",
               "Right PHG parahippocampal gyrus", "Left PHG parahippocampal gyrus")
  roi_target <- roi[roi_name %in% targets]
  temporal <- roi_target[, .(cat12_temporal_gm_ml = sum(Vgm_ml, na.rm = TRUE),
                             n_target_rois = sum(is.finite(Vgm_ml))), by = job_id]
  out <- merge(report, temporal, by = "job_id", all.x = TRUE)
  # data.table evaluates RHS expressions within the same := call before
  # newly-created columns are guaranteed to exist. Create the ratio first,
  # then derive the standardized AD-vulnerable GM score in a second step.
  out[, cat12_temporal_gm_tiv_ratio := cat12_temporal_gm_ml / tiv_ml]
  out[, cat12_ad_vulnerable_gm := -as.numeric(
    scale(log(cat12_temporal_gm_tiv_ratio))
  )]
  out
}

run_aibl_event_bridge <- function(parsed) {
  aibl_path <- c(
    study_path("aibl_clinical"),
    study_path("aibl_clinical")
  )
  aibl_path <- aibl_path[file.exists(aibl_path)][1L]
  if (is.na(aibl_path)) stop("AIBL clinical event file not found.", call. = FALSE)
  clinical <- fread(file = aibl_path, showProgress = FALSE)
  clinical[, RID := norm_id(RID)]
  clinical[, `:=`(event = as.integer(Event) == 1L,
                  followup_years = pmax(num(Time), 0) / 12,
                  age = num(Age), apoe4 = num(APOE4_Carrier),
                  sex = factor(Gender), baseline_dx = factor(DX_Baseline))]
  endpoint <- derive_aibl_endpoint(parsed)
  endpoint[, RID := norm_id(RID)]
  x <- merge(endpoint, clinical[, .(RID, event, followup_years, age, apoe4, sex, baseline_dx, MMSE)],
              by = "RID", all = FALSE)
  x <- x[is.finite(cat12_ad_vulnerable_gm) & is.finite(followup_years) & followup_years > 0]
  if (nrow(x) < 100L || sum(x$event) < 15L) {
    stop("AIBL CAT12 endpoint has insufficient complete clinical events: n=", nrow(x), ", events=", sum(x$event), call. = FALSE)
  }
  x[, `:=`(age_z = as.numeric(scale(age)), apoe4_z = as.numeric(scale(apoe4)),
           mmse_z = as.numeric(scale(num(MMSE)) ))]
  x[, event := as.integer(event)]
  fit0 <- coxph(Surv(followup_years, event) ~ cat12_ad_vulnerable_gm + age_z + apoe4_z + sex,
                data = x, ties = "efron", robust = TRUE)
  fit1 <- coxph(Surv(followup_years, event) ~ cat12_ad_vulnerable_gm + age_z + apoe4_z + sex + mmse_z,
                data = x[complete.cases(x[, .(cat12_ad_vulnerable_gm, age_z, apoe4_z, sex, mmse_z)])],
                ties = "efron", robust = TRUE)
  extract <- function(fit, model_id) {
    b <- coef(fit); se <- sqrt(diag(vcov(fit))); keep <- is.finite(b) & is.finite(se)
    z <- b[keep] / se[keep]
    data.table(model_id = model_id, term = names(b)[keep], estimate = b[keep],
               hazard_ratio = exp(b[keep]), conf_low = exp(b[keep] - 1.96 * se[keep]),
               conf_high = exp(b[keep] + 1.96 * se[keep]), p_value = 2 * pnorm(-abs(z)),
               n = nobs(fit), events = fit$nevent)
  }
  event_results <- rbindlist(list(extract(fit0, "AIBL_CAT12_event_core"), extract(fit1, "AIBL_CAT12_event_plus_MMSE")), fill = TRUE)
  fwrite(endpoint, file.path(SUMMARY_ROOT, "04_AIBL_CAT12_endpoint.csv"), bom = TRUE)
  fwrite(x, file.path(SUMMARY_ROOT, "05_AIBL_CAT12_clinical_linked_frame.csv"), bom = TRUE)
  fwrite(event_results, file.path(SUMMARY_ROOT, "06_AIBL_CAT12_event_cox_results.csv"), bom = TRUE)
  data.table(n_endpoint = nrow(endpoint), n_linked = nrow(x), events = sum(x$event),
             mean_followup_years = mean(x$followup_years),
             CAT12_event_core_HR = event_results[model_id == "AIBL_CAT12_event_core" & term == "cat12_ad_vulnerable_gm", hazard_ratio],
             CAT12_event_core_p = event_results[model_id == "AIBL_CAT12_event_core" & term == "cat12_ad_vulnerable_gm", p_value],
             CAT12_event_plus_MMSE_HR = event_results[model_id == "AIBL_CAT12_event_plus_MMSE" & term == "cat12_ad_vulnerable_gm", hazard_ratio],
             CAT12_event_plus_MMSE_p = event_results[model_id == "AIBL_CAT12_event_plus_MMSE" & term == "cat12_ad_vulnerable_gm", p_value])
}

run_adni_calibration <- function(parsed) {
  roi <- parsed$roi[cohort == "ADNI" & atlas == "neuromorphometrics"]
  targets <- c("Right Ent entorhinal area", "Left Ent entorhinal area",
               "Right FuG fusiform gyrus", "Left FuG fusiform gyrus",
               "Right ITG inferior temporal gyrus", "Left ITG inferior temporal gyrus",
               "Right MTG middle temporal gyrus", "Left MTG middle temporal gyrus",
               "Right TMP temporal pole", "Left TMP temporal pole",
               "Right PHG parahippocampal gyrus", "Left PHG parahippocampal gyrus")
  roi <- roi[roi_name %in% targets]
  cat12 <- roi[, .(cat12_temporal_gm_ml = sum(Vgm_ml, na.rm = TRUE)), by = .(job_id, RID, image_id)]
  fs_path <- file.path(PROJECT_ROOT, "results", "13F_R_only_structural_analysis", "04_ADNI_R_only_longitudinal_structural_frame.rds")
  if (!file.exists(fs_path)) stop("ADNI structural frame unavailable for calibration.", call. = FALSE)
  fs <- as.data.table(readRDS(fs_path)); fs[, `:=`(RID = norm_id(RID), image_id = as.character(IMAGEUID))]
  cal <- merge(cat12, fs[, .(RID, image_id, ad_vulnerable_index)], by = c("RID", "image_id"), all = FALSE)
  cal <- cal[is.finite(cat12_temporal_gm_ml) & is.finite(ad_vulnerable_index)]
  if (nrow(cal) < 20L) stop("Fewer than 20 ADNI calibration images linked to the frozen structural frame.", call. = FALSE)
  cal[, `:=`(cat12_z = as.numeric(scale(cat12_temporal_gm_ml)), fs_z = as.numeric(scale(ad_vulnerable_index)))]
  fit <- lm(fs_z ~ cat12_z, data = cal)
  cor_pearson <- cor(cal$cat12_z, cal$fs_z)
  cor_spearman <- cor(cal$cat12_z, cal$fs_z, method = "spearman")
  out <- data.table(n_images = nrow(cal), n_people = uniqueN(cal$RID), pearson_r = cor_pearson,
                    spearman_rho = cor_spearman, slope = coef(fit)[["cat12_z"]],
                    p_value = summary(fit)$coefficients["cat12_z", "Pr(>|t|)"],
                    direction = ifelse(cor_pearson > 0, "concordant", "discordant"))
  fwrite(cal, file.path(SUMMARY_ROOT, "07_ADNI_CAT12_cross_pipeline_calibration_frame.csv"), bom = TRUE)
  fwrite(out, file.path(SUMMARY_ROOT, "08_ADNI_CAT12_cross_pipeline_calibration.csv"), bom = TRUE)
  out
}

# -----------------------------------------------------------------------------
# Main control flow.
# -----------------------------------------------------------------------------
if (.Platform$OS.type != "windows" && MODE %in% c("prepare", "run"))
  stop("This compiled CAT12 launcher uses the Windows MATLAB Runtime. Use completed CAT12 derivatives for statistical analyses on other platforms.", call. = FALSE)
manifest <- read_manifest()
selected <- select_rows(manifest)
audit <- data.table(
  item = c("project_root", "cat_executable", "cat_template", "runtime_root", "manifest", "manifest_primary_rows",
           "selected_cohort", "selected_rows", "raw_input_files_exist", "stage_root", "workers", "keep_vbm_maps"),
  value = c(PROJECT_ROOT, CAT_EXE, CAT_TEMPLATE, MCR_ROOT, MANIFEST, nrow(manifest), COHORT, nrow(selected),
            all(file.exists(selected$nifti_path)), STAGE_ROOT, WORKERS, KEEP_VBM_MAPS)
)
fwrite(audit, file.path(RESULT_ROOT, "00_CAT12_R_direct_audit.csv"), bom = TRUE)
if (!file.exists(CAT_EXE) || !file.exists(CAT_TEMPLATE) || !dir.exists(MCR_ROOT)) stop("CAT12/MCR installation audit failed.", call. = FALSE)
if (!all(file.exists(selected$nifti_path))) stop("At least one selected raw NIfTI is missing.", call. = FALSE)

if (MODE == "prepare") {
  prepared <- rbindlist(lapply(seq_len(nrow(selected)), function(i) {
    st <- stage_job(selected[i]); data.table(job_id = selected$job_id[i], stage_input = st$stage_input, batch_path = st$batch_path)
  }), fill = TRUE)
  fwrite(prepared, file.path(RESULT_ROOT, "04_prepared_jobs.csv"), bom = TRUE)
  message("CAT12 R direct prepare completed: ", RESULT_ROOT)
}
if (MODE == "run") {
  queue_result <- run_queue(selected)
  fwrite(queue_result, file.path(RESULT_ROOT, "05_CAT12_run_summary.csv"), bom = TRUE)
  message("CAT12 R direct run completed: ", RESULT_ROOT)
}
if (MODE == "parse") {
  parsed <- parse_existing()
  if (COHORT == "AIBL") {
    event_summary <- run_aibl_event_bridge(parsed)
    fwrite(event_summary, file.path(SUMMARY_ROOT, "09_AIBL_CAT12_event_summary.csv"), bom = TRUE)
  }
  if (COHORT == "ADNI_CALIBRATION") {
    calibration <- run_adni_calibration(parsed)
    fwrite(calibration, file.path(SUMMARY_ROOT, "09_ADNI_CAT12_calibration_summary.csv"), bom = TRUE)
  }
  message("CAT12 R direct parse completed: ", RESULT_ROOT)
}
if (MODE == "audit") message("CAT12 R direct audit completed: ", RESULT_ROOT)
