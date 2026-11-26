# Four-trait gMAV common-factor GWAS
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

#!/usr/bin/env Rscript

# Step 09F: GWAS of the frozen shared MAV genetic factor. The latent factor
# was defined in Step 09E without AD; AD is not used to estimate SNP effects.

options(stringsAsFactors = FALSE, scipen = 999, width = 180, warn = 1)

PROJECT_DIR <- study_path("project")
RESULT_09D <- file.path(PROJECT_DIR, "results", "09D_GWAS_QC_LDSC")
RESULT_09E <- file.path(PROJECT_DIR, "results", "09E_GenomicSEM_model_comparison")
RESULT_DIR <- file.path(PROJECT_DIR, "results", "09F_MAV_commonfactorGWAS")
STAGE_DIR <- study_path("gwas_stage", "09F")
INPUT_DIR <- file.path(STAGE_DIR, "sumstats_input")
VROOM_TEMP_DIR <- file.path(STAGE_DIR, "vroom_temp")
# GenomicSEM 0.0.5 requires a MAF column in `ref`; the local LDSC HM3 list
# contains SNP/A1/A2 only. A project-side adapter is created below.
REF_FILE_RAW <- study_path("ldsc", "w_hm3.snplist")
REF_FILE <- file.path(STAGE_DIR, "w_hm3.genomicsem_ref.tsv")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(VROOM_TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
if (file.access(VROOM_TEMP_DIR, mode = 2) != 0) {
  stop("Cannot write to the staged vroom temporary directory: ", VROOM_TEMP_DIR, call. = FALSE)
}
# Use the configured scratch directory for vroom temporary files.
Sys.setenv(VROOM_TEMP_PATH = normalizePath(VROOM_TEMP_DIR, winslash = "/", mustWork = TRUE))
# Workers inherit TMP/TEMP/TMPDIR. Recreate tempdir() if it has been removed.
Sys.setenv(
  TMPDIR = normalizePath(VROOM_TEMP_DIR, winslash = "/", mustWork = TRUE),
  TMP = normalizePath(VROOM_TEMP_DIR, winslash = "/", mustWork = TRUE),
  TEMP = normalizePath(VROOM_TEMP_DIR, winslash = "/", mustWork = TRUE)
)
SESSION_TEMP_DIR <- tempdir(check = TRUE)
if (!dir.exists(SESSION_TEMP_DIR)) {
  dir.create(SESSION_TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(SESSION_TEMP_DIR) || file.access(SESSION_TEMP_DIR, mode = 2) != 0) {
  stop("The current R session temporary directory cannot be recreated: ", SESSION_TEMP_DIR, call. = FALSE)
}
LOG_FILE <- file.path(RESULT_DIR, "09F_MAV_commonfactorGWAS.log")
if (file.exists(LOG_FILE)) unlink(LOG_FILE, force = TRUE)

log_line <- function(level = "INFO", ...) {
  line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  level, paste0(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("GenomicSEM is not installed.", call. = FALSE)
}
if (!file.exists(REF_FILE_RAW)) stop("Missing HM3 reference file: ", REF_FILE_RAW, call. = FALSE)
M1_RDS <- file.path(RESULT_09E, "04_M1_single_gMAV.rds")
LDSC_RDS <- file.path(RESULT_09D, "05_primary_five_trait_covstruc.rds")
if (!file.exists(M1_RDS)) stop("Run Step 09E first: missing ", M1_RDS, call. = FALSE)
if (!file.exists(LDSC_RDS)) stop("Run Step 09D first: missing ", LDSC_RDS, call. = FALSE)

source_files <- c(
  CMM = file.path(RESULT_09D, "01_formatted", "CMM", "CMM_TwosampleMR.txt"),
  INSOMNIA = file.path(RESULT_09D, "01_formatted", "INSOMNIA", "INSOMNIA_TwosampleMR.txt"),
  HEARING = file.path(RESULT_09D, "01_formatted", "HEARING", "HEARING_TwosampleMR.txt"),
  GRIP = file.path(RESULT_09D, "01_formatted", "GRIP", "GRIP_TwosampleMR.txt")
)
N <- c(CMM = 367147, INSOMNIA = 462341, HEARING = 330759, GRIP = 256523)
stop_missing <- source_files[!file.exists(source_files)]
if (length(stop_missing)) stop("Missing formatted GWAS files:\n", paste(stop_missing, collapse = "\n"), call. = FALSE)

REF_MAF_SOURCE_DIR <- study_path("ldsc", "eur_w_ld_chr")
ref_needs_rebuild <- !file.exists(REF_FILE)
if (!ref_needs_rebuild) {
  ref_probe <- data.table::fread(REF_FILE, nrows = 1000L, showProgress = FALSE)
  ref_needs_rebuild <- !("MAF" %in% names(ref_probe)) ||
    (nrow(ref_probe) > 0L && all(abs(as.numeric(ref_probe$MAF) - 0.05) < 1e-12))
}
if (ref_needs_rebuild) {
  ref_raw <- data.table::fread(REF_FILE_RAW, showProgress = FALSE)
  required_ref <- c("SNP", "A1", "A2")
  if (!all(required_ref %in% names(ref_raw))) {
    stop("HM3 reference must contain SNP, A1, and A2 columns.", call. = FALSE)
  }
  ld_files <- file.path(REF_MAF_SOURCE_DIR, paste0(1:22, ".l2.ldscore.gz"))
  if (any(!file.exists(ld_files))) {
    stop("Missing LDSC files needed for real MAF values:\n", paste(ld_files[!file.exists(ld_files)], collapse = "\n"), call. = FALSE)
  }
  log_line("INFO", "Building GenomicSEM reference with real 1000G EUR MAF from LDSC files")
  ref_maf <- data.table::rbindlist(lapply(ld_files, function(x) {
    data.table::fread(x, select = c("SNP", "MAF"), showProgress = FALSE)
  }), use.names = TRUE)
  ref_maf <- unique(ref_maf, by = "SNP")
  ref_out <- ref_raw[, required_ref, with = FALSE]
  ref_out <- merge(ref_out, ref_maf, by = "SNP", all = FALSE, sort = FALSE)
  ref_out <- ref_out[is.finite(MAF) & MAF >= 0.01 & MAF <= 0.5]
  if (!nrow(ref_out)) stop("No HM3 SNPs retained after merging real MAF values.", call. = FALSE)
  data.table::fwrite(ref_out, REF_FILE, sep = "\t", quote = FALSE)
  write.csv(data.frame(reference_file = REF_FILE, maf_source = REF_MAF_SOURCE_DIR,
                       retained_snps = nrow(ref_out), maf_min = min(ref_out$MAF),
                       maf_max = max(ref_out$MAF)),
            file.path(RESULT_DIR, "01C_reference_MAF_audit.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
  rm(ref_raw, ref_maf, ref_out)
  gc(verbose = FALSE)
}

make_sumstats_input <- function(path, trait, n) {
  log_line("INFO", "Reading ", trait, ": ", path)
  dat <- data.table::fread(path, showProgress = FALSE)
  pick <- function(...) {
    candidates <- c(...)
    hit <- candidates[candidates %in% names(dat)]
    if (!length(hit)) stop("Missing required column for ", trait, ": ", paste(candidates, collapse = ","), call. = FALSE)
    dat[[hit[[1]]]]
  }
  out <- data.frame(
    SNP = as.character(pick("SNP")),
    A1 = toupper(as.character(pick("effect_allele.outcome"))),
    A2 = toupper(as.character(pick("other_allele.outcome"))),
    beta = suppressWarnings(as.numeric(pick("beta.outcome"))),
    se = suppressWarnings(as.numeric(pick("se.outcome"))),
    p = suppressWarnings(as.numeric(pick("pval.outcome"))),
    N = rep(as.numeric(n), nrow(dat)),
    stringsAsFactors = FALSE
  )
  keep <- is.finite(out$beta) & is.finite(out$se) & out$se > 0 &
    is.finite(out$p) & out$p > 0 & out$p <= 1 &
    nzchar(out$SNP) & nchar(out$A1) == 1L & nchar(out$A2) == 1L & out$A1 != out$A2
  out <- out[keep, , drop = FALSE]
  out$p[out$p == 0] <- 1e-300
  if (!nrow(out)) stop("No valid rows remain for ", trait, call. = FALSE)
  output <- file.path(INPUT_DIR, paste0(trait, ".sumstats.txt"))
  data.table::fwrite(out, output, sep = "\t", quote = FALSE, na = "NA")
  log_line("INFO", trait, " retained ", nrow(out), " rows")
  output
}

input_files <- vapply(names(source_files), function(trait) {
  make_sumstats_input(source_files[[trait]], trait, N[[trait]])
}, character(1))
write.csv(data.frame(trait = names(input_files), file = unname(input_files), N = unname(N)),
          file.path(RESULT_DIR, "01_sumstats_input_manifest.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

# Preflight the reference merge explicitly. GenomicSEM otherwise fails later
# with an opaque effect[[1]] subscript error when no SNP survives the merge.
reference_ids <- data.table::fread(REF_FILE, select = "SNP", showProgress = FALSE)$SNP
overlap_audit <- do.call(rbind, lapply(names(input_files), function(trait) {
  input_ids <- data.table::fread(input_files[[trait]], select = "SNP", showProgress = FALSE)$SNP
  data.frame(trait = trait, input_snps = length(unique(input_ids)),
             reference_snps = length(unique(reference_ids)),
             overlap_snps = length(intersect(input_ids, reference_ids)))
}))
write.csv(overlap_audit, file.path(RESULT_DIR, "01B_reference_overlap_audit.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
if (any(overlap_audit$overlap_snps < 1000L)) {
  stop("Too few SNPs overlap the GenomicSEM reference. See 01B_reference_overlap_audit.csv.",
       call. = FALSE)
}

log_line("INFO", "Running GenomicSEM::sumstats()")
harmonized <- GenomicSEM::sumstats(
  files = unname(input_files),
  ref = REF_FILE,
  trait.names = names(input_files),
  se.logit = c(TRUE, TRUE, TRUE, FALSE),
  OLS = c(FALSE, FALSE, FALSE, TRUE),
  N = unname(N),
  info.filter = 0.6,
  maf.filter = 0.01,
  keep.indel = FALSE,
  parallel = FALSE,
  cores = NULL,
  ambig = FALSE,
  direct.filter = FALSE
)
beta_cols <- grep("^beta\\.", names(harmonized), value = TRUE)
se_cols <- grep("^se\\.", names(harmonized), value = TRUE)
if (length(beta_cols) != 4L || length(se_cols) != 4L) {
  stop("Expected four beta.* and four se.* columns after GenomicSEM::sumstats().", call. = FALSE)
}
se_matrix <- as.matrix(harmonized[, se_cols, drop = FALSE])
beta_matrix <- as.matrix(harmonized[, beta_cols, drop = FALSE])
finite_keep <- apply(se_matrix, 1L, function(x) all(is.finite(x) & x > 1e-6)) &
  apply(beta_matrix, 1L, function(x) all(is.finite(x)))
commonfactor_qc <- data.frame(
  input_rows = nrow(harmonized),
  retained_rows = sum(finite_keep),
  removed_rows = sum(!finite_keep),
  minimum_se_threshold = 1e-6
)
write.csv(commonfactor_qc, file.path(RESULT_DIR, "02C_commonfactor_SNP_QC.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
if (!all(finite_keep)) {
  harmonized <- harmonized[finite_keep, , drop = FALSE]
}
if (!nrow(harmonized)) {
  stop("No SNPs remain after commonfactorGWAS numerical QC.", call. = FALSE)
}
saveRDS(harmonized, file.path(RESULT_DIR, "02_harmonized_sumstats.rds"))
log_line("INFO", "Harmonized summary statistics saved")

# commonfactorGWAS() uses the native four-trait ldsc() object in V, S, I order.
ld_dir <- study_path("ldsc", "eur_w_ld_chr/")
ldsc_sumstats <- c(
  CMM = study_path("gwas_stage", "09D/sumstats/CMM/CMM.sumstats.gz"),
  INSOMNIA = study_path("gwas_stage", "09D/sumstats/INSOMNIA/INSOMNIA.sumstats.gz"),
  HEARING = study_path("gwas_stage", "09D/sumstats/HEARING/HEARING.sumstats.gz"),
  GRIP = study_path("gwas_stage", "09D/sumstats/GRIP/GRIP.sumstats.gz")
)
if (any(!file.exists(ldsc_sumstats))) {
  stop("Missing Step 09D LDSC sumstats:\n", paste(ldsc_sumstats[!file.exists(ldsc_sumstats)], collapse = "\n"), call. = FALSE)
}
log_line("INFO", "Creating native four-trait GenomicSEM LDSC covariance structure")
ldsc_mav <- GenomicSEM::ldsc(
  traits = unname(ldsc_sumstats),
  sample.prev = rep(NA_real_, length(ldsc_sumstats)),
  population.prev = rep(NA_real_, length(ldsc_sumstats)),
  ld = ld_dir,
  wld = ld_dir,
  trait.names = names(ldsc_sumstats),
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = file.path(RESULT_DIR, "02B_MAV4_ldsc"),
  stand = FALSE,
  select = FALSE,
  chisq.max = NA
)
saveRDS(ldsc_mav, file.path(RESULT_DIR, "02B_MAV4_ldsc_covstruc.rds"))

log_line("INFO", "Running GenomicSEM::commonfactorGWAS()")
factor_gwas <- GenomicSEM::commonfactorGWAS(
  covstruc = ldsc_mav,
  SNPs = harmonized,
  estimation = "DWLS",
  cores = NULL,
  # A single SNP can produce an extremely ill-conditioned but non-zero
  # sampling covariance matrix. The prefilter above removes invalid SEs;
  # this smaller solve tolerance permits valid near-singular matrices to be
  # handled without changing the model or package internals.
  toler = 1e-40,
  SNPSE = FALSE,
  parallel = TRUE,
  GC = "standard",
  MPI = FALSE,
  TWAS = FALSE,
  smooth_check = FALSE
)
saveRDS(factor_gwas, file.path(RESULT_DIR, "03_MAV_commonfactorGWAS.rds"))
if (is.data.frame(factor_gwas)) {
  data.table::fwrite(factor_gwas, file.path(RESULT_DIR, "03_MAV_commonfactorGWAS.tsv"), sep = "\t")
}
writeLines(capture.output(sessionInfo()), file.path(RESULT_DIR, "04_sessionInfo.txt"), useBytes = TRUE)
log_line("INFO", "Step 09F completed")
