# Official MRI derivatives and prespecified anatomical indices
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

project_dir <- study_path("project")
adni_root <- study_path("cohorts", "ADNI/data")
habs_root <- study_path("habs_release")
out_dir <- file.path(project_dir, "results", "12B_ADNI_HABSHD_structural_MRI_feature_freeze")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] [INFO] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
need <- function(p) { if (!file.exists(p)) stop("Missing input: ", p, call. = FALSE); p }
clean_names <- function(x) {
  x <- trimws(as.character(x)); x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x); make.unique(x, sep = "_")
}
to_num <- function(x) suppressWarnings(as.numeric(x))
write_csv <- function(x, n) { fwrite(as.data.table(x), file.path(out_dir, n), bom = TRUE); log_msg("Wrote ", n, " (", nrow(x), " rows)") }

adni_fs_file <- need(file.path(adni_root, "MRI", "MRI_FreeSurfer.csv"))
adni_dict_file <- need(file.path(adni_root, "MRI", "MRI_FreeSurfer_Data_Dictionary.csv"))
adni_hip_file <- need(file.path(adni_root, "MRI", "MRI_Hippocampus_ASHS.csv"))
adni_wmh_file <- need(file.path(adni_root, "MRI", "UCD_WMH_03Jan2026.csv"))
habs_ct_file <- need(file.path(habs_root, "HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx"))
habs_meta_file <- need(file.path(habs_root, "HABS_HD_RP_8_CorticalThickness_ATN_MetaROI_Data_BraskieLab_20260501.xlsx"))
habs_hip_file <- need(file.path(habs_root, "HABS_HD_RP_8_Hippocampal_Volume_HippoDeep_Data_BraskieLab_20260501.xlsx"))
habs_brain_file <- need(file.path(habs_root, "HABS_HD_RP_8_Whole_Brain_Volume_Data_BraskieLab_20260501.xlsx"))
habs_gen_file <- need(file.path(habs_root, "RP_HD_7_Genomics.xlsx"))

log_msg("Step 12B started: freeze regional structural MRI features; no outcome model is fitted.")

# Freeze biological regions before looking at genetic or cognitive associations.
region_sets <- list(
  AD_vulnerable = c("entorhinal", "fusiform", "inferiortemporal", "middletemporal", "temporalpole", "parahippocampal"),
  frontal = c("caudalmiddlefrontal", "rostralmiddlefrontal", "superiorfrontal", "lateralorbitofrontal", "medialorbitofrontal", "frontalpole"),
  posterior = c("precuneus", "posteriorcingulate", "isthmuscingulate", "inferiorparietal", "superiorparietal", "lateraloccipital"),
  primary_sensorimotor = c("precentral", "postcentral", "pericalcarine", "transversetemporal")
)
region_manifest <- rbindlist(lapply(names(region_sets), function(s) data.table(index=s, region=region_sets[[s]])))
write_csv(region_manifest, "00_frozen_region_sets.csv")

# Read HABS-HD inclusion sheets and standardize headers.
read_inclusion_sheets <- function(path, pattern) {
  sheets <- excel_sheets(path)
  keep <- sheets[grepl(pattern, sheets, ignore.case = TRUE) & !grepl("excl|pending", sheets, ignore.case = TRUE)]
  if (!length(keep)) stop("No inclusion sheet found in ", basename(path), call. = FALSE)
  rbindlist(lapply(seq_along(keep), function(i) {
    # Read mixed-type release sheets as text first. Literal "NA" cells in
    # numeric columns are converted explicitly below and should not trigger
    # readxl type-guessing warnings.
    x <- suppressWarnings(as.data.table(read_excel(
      path, sheet = keep[i], col_types = "text", na = c("", "NA", "N/A", "."),
      .name_repair = "minimal"
    )))
    setnames(x, clean_names(names(x))); x[, source_sheet := keep[i]]; x
  }), fill = TRUE, use.names = TRUE)
}

habs_ct <- read_inclusion_sheets(habs_ct_file, "inclusion")
habs_meta <- read_inclusion_sheets(habs_meta_file, "inclusion")
habs_hip <- read_inclusion_sheets(habs_hip_file, "inclusion")
habs_brain <- read_inclusion_sheets(habs_brain_file, "included|inclusion")
habs_gen <- suppressWarnings(as.data.table(read_excel(
  habs_gen_file, col_types = "text", na = c("", "NA", "N/A", "."),
  .name_repair = "minimal"
)))
setnames(habs_gen, clean_names(names(habs_gen)))

# Region-level HABS-HD thickness indices.
ct_cols <- grep("^[LR]_.*_thickavg$", names(habs_ct), value = TRUE)
if (length(ct_cols) < 60) stop("HABS-HD cortical thickness columns <60; check release format.", call. = FALSE)
for (v in ct_cols) set(habs_ct, j=v, value=to_num(habs_ct[[v]]))

cols_for <- function(regions) unique(unlist(lapply(regions, function(r) grep(paste0("^[LR]_", r, "_thickavg$"), names(habs_ct), value=TRUE))))
row_mean_safe <- function(dt, cols) { if (!length(cols)) return(rep(NA_real_, nrow(dt))); x=as.matrix(dt[, cols, with=FALSE]); out=rowMeans(x,na.rm=TRUE); out[rowSums(is.finite(x))==0]=NA_real_; out }
habs_ct[, global_cortical_thickness := row_mean_safe(.SD, ct_cols)]
habs_ct[, ad_vulnerable_index := row_mean_safe(.SD, cols_for(region_sets$AD_vulnerable))]
habs_ct[, frontal_index := row_mean_safe(.SD, cols_for(region_sets$frontal))]
habs_ct[, posterior_index := row_mean_safe(.SD, cols_for(region_sets$posterior))]
habs_ct[, frontoposterior_gradient := frontal_index - posterior_index]
habs_ct[, spatial_dispersion := {x=as.matrix(.SD[, ct_cols, with=FALSE]); apply(x,1,sd,na.rm=TRUE)}]

paired_regions <- sub("^[LR]_", "", ct_cols[grepl("^L_",ct_cols)])
asym_cols <- intersect(paste0("L_",paired_regions),ct_cols)
asym_r <- paste0("R_",sub("^L_", "", asym_cols))
ok_pair <- asym_r %in% ct_cols; asym_cols <- asym_cols[ok_pair]; asym_r <- asym_r[ok_pair]
if (length(asym_cols)) {
  L=as.matrix(habs_ct[,..asym_cols]); R=as.matrix(habs_ct[,..asym_r])
  habs_ct[, mean_absolute_asymmetry := rowMeans(abs(L-R)/pmax(1e-6,(L+R)/2),na.rm=TRUE)]
} else habs_ct[, mean_absolute_asymmetry := NA_real_]

id_candidates <- c("Med_ID","MED_ID"); id_h <- id_candidates[id_candidates %in% names(habs_ct)][1]
if (is.na(id_h)) stop("Med_ID missing in HABS-HD cortical thickness.", call. = FALSE)
time_h <- c("Time_point_MRI","Time_point"); time_h <- time_h[time_h %in% names(habs_ct)][1]
habs_features <- habs_ct[, .(
  cohort="HABS-HD", participant_id=as.character(get(id_h)),
  mri_timepoint=if(!is.na(time_h))as.character(get(time_h)) else NA_character_,
  age=to_num(Age), scanner=as.character(Scanner_MRI),
  global_cortical_thickness, ad_vulnerable_index, frontal_index, posterior_index,
  frontoposterior_gradient, spatial_dispersion, mean_absolute_asymmetry
)]
saveRDS(habs_features, file.path(out_dir,"01_HABSHD_frozen_cortical_features.rds"),compress="xz")

# ADNI variable dictionary maps opaque ST codes to anatomical labels.
adni_dict <- fread(adni_dict_file, showProgress=FALSE)
write_csv(adni_dict, "02_ADNI_FreeSurfer_dictionary_snapshot.csv")
fs <- fread(adni_fs_file, showProgress=FALSE)
hip <- fread(adni_hip_file, showProgress=FALSE)
wmh <- fread(adni_wmh_file, showProgress=FALSE)

# Read the official ADNI structural derivatives.
hip[, hippocampus_total := to_num(LEFT_HIPP_VOL)+to_num(RIGHT_HIPP_VOL)]
hip[, hippocampus_icv_ratio := hippocampus_total/to_num(ICV)]
wmh[, wmh_brain_ratio := to_num(TOTAL_WMH)/to_num(TOTAL_BRAIN)]
adni_nonregional <- merge(
  hip[,.(RID,PTID,VISCODE,VISCODE2,EXAMDATE,hippocampus_total,hippocampus_icv_ratio)],
  wmh[,.(RID,PTID_wmh=PTID,VISCODE_wmh=VISCODE,VISCODE2_wmh=VISCODE2,EXAMDATE_wmh=EXAMDATE,wmh_brain_ratio,total_brain=to_num(TOTAL_BRAIN))],
  by="RID",all=TRUE,allow.cartesian=TRUE
)
saveRDS(adni_nonregional,file.path(out_dir,"03_ADNI_hippocampus_WMH_features.rds"),compress="xz")

# Availability and cross-cohort construct mapping.
mapping <- data.table(
  construct=c("global cortical thickness","AD-vulnerable cortical index","frontoposterior gradient","spatial dispersion","mean absolute asymmetry","hippocampal volume/ICV","whole brain volume/ICV","WMH burden","APOE4","genome-wide PGS"),
  ADNI=c("Requires ST-code dictionary mapping","Requires ST-code dictionary mapping","Requires ST-code dictionary mapping","Requires ST-code dictionary mapping","Requires ST-code dictionary mapping","AVAILABLE","AVAILABLE","AVAILABLE","AVAILABLE","AVAILABLE"),
  HABSHD=c("AVAILABLE","AVAILABLE","AVAILABLE","AVAILABLE","AVAILABLE","AVAILABLE","AVAILABLE","NOT AVAILABLE in current release","AVAILABLE","NOT AVAILABLE in current release"),
  role=c(rep("common structural MRI construct",7),"ADNI-only vascular sensitivity","common genetic modifier","ADNI-only genetic mechanism")
)
write_csv(mapping,"04_cross_cohort_construct_mapping.csv")

gate <- data.table(
  gate="ADNI_HABSHD_REGIONAL_STRUCTURAL_MRI_FEATURE_FREEZE",
  habshd_cortical_rows=nrow(habs_features),
  habshd_unique_people=uniqueN(habs_features$participant_id),
  habshd_thickness_regions=length(ct_cols),
  adni_freesurfer_rows=nrow(fs),
  adni_hippocampus_rows=nrow(hip),
  adni_wmh_rows=nrow(wmh),
  raw_MRI_required=FALSE,
  status=if(length(ct_cols)>=60 && nrow(fs)>1000)"READY_FOR_12C" else "REQUIRES_REPAIR"
)
write_csv(gate,"05_12B_gate.csv")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"99_sessionInfo.txt"),useBytes=TRUE)
log_msg("Step 12B completed: ",gate$status)
