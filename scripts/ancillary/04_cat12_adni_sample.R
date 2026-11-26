# Prespecified CAT12 baseline sample and image manifest
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

# CAT12 input sample requires both genetic scores and APOE4.
# Eligible participants have an A/T/V/N measurement and cognitive follow-up.
# Each participant is linked to a baseline T1 image.

suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

out_root <- study_path("project", "results/13F_ADNI_CAT12_mechanistic_bridge_subset")
dir.create(file.path(out_root, "tables"), recursive = TRUE, showWarnings = FALSE)

pgs_file <- paste0(
  study_path("project", "results/10A_ADNI_gMAV_PGS_assembly/"),
  "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"
)
manifest_file <- paste0(
  study_path("project", "results/13F_automated_raw_T1_QC_and_processing_manifest/"),
  "01_13F_frozen_primary_processing_manifest.csv"
)
endpoint_dir <- study_path("project", "results/10B_ADNI_ATVN_PGS_clinical_expression")
endpoint_files <- c(
  A = file.path(endpoint_dir, "05_paired_A_baseline_data.tsv"),
  T = file.path(endpoint_dir, "05_paired_T_baseline_data.tsv"),
  V = file.path(endpoint_dir, "05_paired_V_baseline_data.tsv"),
  N = file.path(endpoint_dir, "05_paired_N_baseline_data.tsv")
)

stopifnot(file.exists(pgs_file), file.exists(manifest_file))

# -----------------------------------------------------------------------------
# 1. PGS and cognitive-follow-up eligibility
# -----------------------------------------------------------------------------
pgs <- fread(pgs_file, sep = "\t", showProgress = FALSE)
need_pgs <- c("PTID", "RID", "VISCODE", "EXAMDATE", "gMAV_noAPOE",
              "gAD_noAPOE", "MMSE")
miss <- setdiff(need_pgs, names(pgs))
if (length(miss)) stop("PGS file missing: ", paste(miss, collapse = ", "))

pgs[, PTID := as.character(PTID)]
pgs[, RID := as.integer(RID)]
pgs[, EXAMDATE := as.IDate(EXAMDATE)]
pgs[, gMAV_noAPOE := as.numeric(gMAV_noAPOE)]
pgs[, gAD_noAPOE := as.numeric(gAD_noAPOE)]

pgs[, cognitive_available := is.finite(MMSE) | is.finite(ADAS13) |
      is.finite(CDRSB) | is.finite(FAQ) | is.finite(MOCA) |
      is.finite(mPACCdigit) | is.finite(mPACCtrailsB)]

pgs_person <- pgs[, .(
  n_pgs_rows = .N,
  n_cognitive_visits = sum(cognitive_available),
  n_total_visits = uniqueN(EXAMDATE),
  gMAV_noAPOE = gMAV_noAPOE[1L],
  gAD_noAPOE = gAD_noAPOE[1L],
  age_bl = age_bl[1L],
  sex = sex[1L],
  baseline_date = suppressWarnings(min(EXAMDATE, na.rm = TRUE))
), by = .(PTID, RID)]
pgs_person[, cognitive_followup_ok :=
  n_total_visits >= 2L & n_cognitive_visits >= 2L]

# -----------------------------------------------------------------------------
# 2. A/T/V/N and APOE availability
# -----------------------------------------------------------------------------
endpoint_rows <- list()
for (nm in names(endpoint_files)) {
  f <- endpoint_files[[nm]]
  if (!file.exists(f) || file.info(f)$size == 0) next
  z <- fread(f, sep = "\t", showProgress = FALSE)
  if (!all(c("PTID", "RID") %in% names(z))) next
  z[, PTID := as.character(PTID)]
  z[, RID := as.integer(RID)]
  apoe_present <- if ("APOE4" %in% names(z)) any(is.finite(as.numeric(z$APOE4))) else FALSE
  endpoint_rows[[length(endpoint_rows) + 1L]] <- z[, .(
    PTID, RID,
    endpoint_available = TRUE,
    endpoint_type = nm,
    APOE4 = if ("APOE4" %in% names(z)) as.numeric(APOE4) else NA_real_
  )]
}
if (!length(endpoint_rows)) stop("No A/T/V/N endpoint files were found.")
ep <- rbindlist(endpoint_rows, fill = TRUE)

ep_person <- ep[, .(
  n_ATVN_endpoints = uniqueN(endpoint_type),
  ATVN_endpoint_types = paste(sort(unique(endpoint_type)), collapse = ";"),
  any_ATVN = TRUE,
  APOE4 = {
    x <- APOE4[is.finite(APOE4)]
    if (length(x)) x[1L] else NA_real_
  }
), by = .(PTID, RID)]

# -----------------------------------------------------------------------------
# 3. Matched baseline T1 from the frozen primary manifest
# -----------------------------------------------------------------------------
man <- fread(manifest_file, showProgress = FALSE)
need_man <- c("processing_image_code", "cohort", "role", "participant_id",
              "RID", "analysis_timepoint", "nifti_path", "image_id")
miss <- setdiff(need_man, names(man))
if (length(miss)) stop("Manifest missing: ", paste(miss, collapse = ", "))

man[, participant_id := as.character(participant_id)]
man[, RID := as.integer(RID)]
man_adni <- man[
  cohort == "ADNI" & role == "PRIMARY" &
    analysis_timepoint %in% c("ADNI_MRI1", "ADNI_BASELINE", "MRI1") &
    file.exists(nifti_path)
]

# Keep one baseline image per participant; prefer MRI1 and the earliest image.
setorder(man_adni, participant_id, image_id)
man_adni <- man_adni[, .SD[1L], by = participant_id]
setnames(man_adni, "participant_id", "PTID")
man_adni[, PTID := as.character(PTID)]

# -----------------------------------------------------------------------------
# 4. Eligibility gate and frozen bridge manifest
# -----------------------------------------------------------------------------
bridge <- Reduce(
  function(x, y) merge(x, y, by = c("PTID", "RID"), all = TRUE),
  list(pgs_person, ep_person)
)
bridge <- merge(bridge, man_adni, by = "PTID", all = TRUE,
                suffixes = c("", "_t1"))
if ("RID_t1" %in% names(bridge)) {
  bridge[is.na(RID) & !is.na(RID_t1), RID := RID_t1]
}
bridge[, has_gMAV := is.finite(gMAV_noAPOE)]
bridge[, has_gAD_noAPOE := is.finite(gAD_noAPOE)]
bridge[, has_APOE4 := is.finite(APOE4)]
bridge[, has_ATVN := !is.na(any_ATVN) & any_ATVN == TRUE]
bridge[, has_cognitive_followup := !is.na(cognitive_followup_ok) &
      cognitive_followup_ok == TRUE]
bridge[, has_baseline_T1 := !is.na(nifti_path) & file.exists(nifti_path)]
bridge[, eligible_bridge := has_gMAV & has_gAD_noAPOE & has_APOE4 &
      has_ATVN & has_cognitive_followup & has_baseline_T1]

eligible <- bridge[eligible_bridge == TRUE]
if (!nrow(eligible)) stop("No ADNI participants passed the mechanistic bridge gate.")

setorder(eligible, PTID)
fwrite(bridge, file.path(out_root, "tables", "01_ADNI_bridge_eligibility_audit.csv"))
fwrite(eligible, file.path(out_root, "tables", "02_ADNI_CAT12_bridge_baseline_manifest.csv"))
dir.create(study_path("cat12_config"), recursive = TRUE, showWarnings = FALSE)
fwrite(eligible, study_path("cat12_config", "02_ADNI_CAT12_bridge_baseline_manifest.csv"))

flow <- data.table(
  gate = c("all_PGS_people", "has_gMAV", "has_gAD_noAPOE", "has_APOE4",
           "has_ATVN", "has_cognitive_followup", "has_baseline_T1",
           "final_bridge_eligible"),
  n_people = c(
    uniqueN(bridge$PTID), sum(bridge$has_gMAV, na.rm = TRUE),
    sum(bridge$has_gAD_noAPOE, na.rm = TRUE), sum(bridge$has_APOE4, na.rm = TRUE),
    sum(bridge$has_ATVN, na.rm = TRUE), sum(bridge$has_cognitive_followup, na.rm = TRUE),
    sum(bridge$has_baseline_T1, na.rm = TRUE), nrow(eligible)
  )
)
fwrite(flow, file.path(out_root, "tables", "03_ADNI_bridge_flow.csv"))

writeLines(c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Eligible participants: ", nrow(eligible)),
  paste0("Eligible baseline T1 images: ", sum(file.exists(eligible$nifti_path))),
  "This is a CAT12 input-manifest step; CAT12 was not run."
), file.path(out_root, "analysis_run_log.txt"))

cat("ADNI CAT12 mechanistic bridge subset frozen: ", out_root, "\n", sep = "")
cat("Eligible participants: ", nrow(eligible), "\n", sep = "")
cat("Manifest: tables/02_ADNI_CAT12_bridge_baseline_manifest.csv\n")
