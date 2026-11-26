# Clumping, allele alignment, PLINK scoring and clinical linkage
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

#!/usr/bin/env Rscript

# Step 10A: Outcome-blind construction of gMAV and AD-noAPOE PGS in ADNI,
# followed by longitudinal clinical/PC assembly. No association model is fit here.

options(stringsAsFactors = FALSE, scipen = 999, width = 180, warn = 1)

PROJECT_DIR <- study_path("project")
RESULT_09G <- file.path(PROJECT_DIR, "results", "09G_MAV_factorGWAS_QC")
GWAS_DIR <- study_path("gwas_input")
ADNI_ROOT <- study_path("cohorts", "ADNI")
PC_FILE <- study_path("genetic_pcs")
PLINK <- study_path("plink")
OUTDIR <- file.path(PROJECT_DIR, "results", "10A_ADNI_gMAV_PGS_assembly")
STAGE <- study_path("gwas_stage", "10A")

GMAV_FILE <- file.path(RESULT_09G, "02_gMAV_factorGWAS_valid_rows.tsv")
AD_FILE <- file.path(GWAS_DIR, "GCST007511.Late-onset Alzheimer's disease.Kunkle_etal_Stage1_results.txt")
CLINICAL_FILE <- file.path(ADNI_ROOT, "data/Clinical/ADNIMERGE_03Jan2026.csv")
GENO_PREFIX <- study_path("genotype_prefix")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAGE, recursive = TRUE, showWarnings = FALSE)
for (d in c("weights", "clump", "score")) dir.create(file.path(STAGE, d), recursive = TRUE, showWarnings = FALSE)

required <- c(GMAV_FILE, AD_FILE, CLINICAL_FILE, PC_FILE, PLINK,
              paste0(GENO_PREFIX, ".bed"), paste0(GENO_PREFIX, ".bim"), paste0(GENO_PREFIX, ".fam"))
if (any(!file.exists(required))) stop("Missing required inputs:\n", paste(required[!file.exists(required)], collapse = "\n"), call. = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.", call. = FALSE)

LOG <- file.path(OUTDIR, "10A_MAV_ADNI_PGS_build_and_longitudinal_assembly.log")
if (file.exists(LOG)) unlink(LOG, force = TRUE)
log_line <- function(level = "INFO", ...) {
  z <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(..., collapse = ""))
  cat(z, "\n"); cat(z, "\n", file = LOG, append = TRUE)
}
run_plink <- function(args, label) {
  log_line("INFO", "PLINK: ", label)
  status <- system2(PLINK, args = shQuote(args), stdout = TRUE, stderr = TRUE)
  if (!identical(attr(status, "status"), NULL) && attr(status, "status") != 0L)
    stop("PLINK failed for ", label, ". See staged .log file.", call. = FALSE)
  invisible(status)
}
is_nonambiguous <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  !(pair %in% c("AT", "TA", "CG", "GC"))
}
standardize <- function(x) {
  if (sum(is.finite(x)) < 2L || is.na(stats::sd(x, na.rm = TRUE)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE))
}

log_line("INFO", "Step 10A started: fixed PGS scoring only; no outcome model is fit.")

# Final ADNI genotype BIM after the prior QC/relatedness pipeline.
bim <- data.table::fread(paste0(GENO_PREFIX, ".bim"), header = FALSE,
                         col.names = c("CHR", "SNP", "CM", "BP", "BIM_A1", "BIM_A2"))
bim[, SNP := as.character(SNP)]
bim <- bim[!duplicated(SNP) & BIM_A1 %in% c("A", "C", "G", "T") & BIM_A2 %in% c("A", "C", "G", "T")]

g <- data.table::fread(GMAV_FILE, sep = "\t", showProgress = TRUE)
need_g <- c("SNP", "A1", "A2", "est", "Pval_Estimate")
if (length(setdiff(need_g, names(g)))) stop("gMAV GWAS file misses required columns.", call. = FALSE)
g <- g[is.finite(est) & is.finite(Pval_Estimate) & Pval_Estimate >= 0 & Pval_Estimate <= 1]
g[, `:=`(SNP = as.character(SNP), A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)),
         BETA = as.numeric(est), P = as.numeric(Pval_Estimate))]
g <- g[!duplicated(SNP) & A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T") & is_nonambiguous(A1, A2)]
g <- merge(g[, .(SNP, A1, A2, BETA, P)], bim[, .(SNP, CHR, BP, BIM_A1, BIM_A2)], by = "SNP")
# Allele matching must be performed within the same SNP row, not against the
# pooled set of all alleles in the BIM file.
g <- g[(A1 == BIM_A1 | A1 == BIM_A2) & !(CHR == 6L & BP >= 25000000 & BP <= 34000000)]

ad <- data.table::fread(AD_FILE, showProgress = TRUE)
need_ad <- c("Chromosome", "Position", "MarkerName", "Effect_allele", "Non_Effect_allele", "Beta", "Pvalue")
if (length(setdiff(need_ad, names(ad)))) stop("AD GWAS file misses expected Kunkle columns.", call. = FALSE)
ad <- ad[, .(CHR = suppressWarnings(as.integer(Chromosome)), BP = suppressWarnings(as.integer(Position)),
             SNP = as.character(MarkerName), A1 = toupper(as.character(Effect_allele)),
             A2 = toupper(as.character(Non_Effect_allele)), BETA = suppressWarnings(as.numeric(Beta)),
             P = suppressWarnings(as.numeric(Pvalue)))]
ad <- ad[is.finite(CHR) & is.finite(BP) & is.finite(BETA) & is.finite(P) & P >= 0 & P <= 1 &
           !duplicated(SNP) & A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T") & is_nonambiguous(A1, A2)]
ad <- merge(ad, bim[, .(SNP, BIM_A1, BIM_A2)], by = "SNP")
ad <- ad[(A1 == BIM_A1 | A1 == BIM_A2) & !(CHR == 6L & BP >= 25000000 & BP <= 34000000)]

# The AD primary score deliberately excludes APOE (GRCh37 chr19:44--46 Mb),
# so it is separable from direct APOE4 adjustment in all clinical models.
weights <- list(
  gMAV = g,
  gMAV_noAPOE = g[!(CHR == 19L & BP >= 44000000 & BP <= 46000000)],
  gAD_noAPOE = ad[!(CHR == 19L & BP >= 44000000 & BP <= 46000000)]
)

score_audit <- list()
scores <- list()
for (nm in names(weights)) {
  w <- weights[[nm]]
  if (!nrow(w)) stop("No usable score weights for ", nm, call. = FALSE)
  pfile <- file.path(STAGE, "weights", paste0(nm, "_clump_input.tsv"))
  data.table::fwrite(w[, .(SNP, P)], pfile, sep = "\t")
  pref <- file.path(STAGE, "clump", nm)
  run_plink(c("--bfile", GENO_PREFIX, "--clump", pfile, "--clump-snp-field", "SNP", "--clump-field", "P",
              "--clump-p1", "1", "--clump-r2", "0.1", "--clump-kb", "500", "--out", pref), paste0(nm, " fixed C+T LD clumping"))
  clump_file <- paste0(pref, ".clumped")
  if (!file.exists(clump_file)) stop("Missing PLINK clump output for ", nm, call. = FALSE)
  cl <- data.table::fread(clump_file, fill = TRUE, showProgress = FALSE)
  if (!"SNP" %in% names(cl)) stop("PLINK clump output lacks SNP column for ", nm, call. = FALSE)
  keep <- unique(as.character(cl$SNP[!is.na(cl$SNP) & cl$SNP != "NONE"]))
  wscore <- w[SNP %in% keep, .(SNP, A1, BETA)]
  if (!nrow(wscore)) stop("No clumped score variants remained for ", nm, call. = FALSE)
  wfile <- file.path(STAGE, "weights", paste0(nm, "_plink_score_weights.tsv"))
  data.table::fwrite(wscore, wfile, sep = "\t")
  spref <- file.path(STAGE, "score", nm)
  run_plink(c("--bfile", GENO_PREFIX, "--score", wfile, "1", "2", "3", "header", "sum", "--out", spref), paste0(nm, " PGS"))
  pf <- paste0(spref, ".profile")
  if (!file.exists(pf)) stop("Missing PLINK score output for ", nm, call. = FALSE)
  x <- data.table::fread(pf, showProgress = FALSE)
  # PLINK 1.9 writes SCORESUM; other supported builds may write SCORE or
  # SCORE1_SUM. Select the score column provided by the installed PLINK version.
  score_col <- intersect(c("SCORESUM", "SCORE", "SCORE1_SUM", "SCORE1_AVG"), names(x))
  if (!length(score_col)) stop("PLINK score output lacks a recognized score column for ", nm,
                               ". Found: ", paste(names(x), collapse = ", "), call. = FALSE)
  score_col <- score_col[[1L]]
  scores[[nm]] <- x[, c("FID", "IID", score_col), with = FALSE]
  data.table::setnames(scores[[nm]], score_col, nm)
  score_audit[[nm]] <- data.frame(score = nm, input_weights = nrow(w), clumped_weights = nrow(wscore),
                                  scored_participants = nrow(x), stringsAsFactors = FALSE)
}
score_audit <- do.call(rbind, score_audit)
write.csv(score_audit, file.path(OUTDIR, "01_PGS_weight_and_scoring_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")

pgs <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID"), all = TRUE), scores)
pc <- data.table::fread(PC_FILE, showProgress = FALSE)
pc[, `:=`(FID = as.character(FID), IID = as.character(IID))]
pgs[, `:=`(FID = as.character(FID), IID = as.character(IID))]
pgs <- merge(pgs, pc, by = c("FID", "IID"), all.x = TRUE)

# RID is not used for PGS linkage: genotype IIDs map directly to ADNI PTIDs.
clinical <- data.table::fread(CLINICAL_FILE, showProgress = TRUE)
must_clin <- c("PTID", "VISCODE", "EXAMDATE", "DX_bl", "DX", "AGE", "PTGENDER", "PTEDUCAT", "APOE4",
               "MMSE", "ADAS13", "CDRSB", "FAQ", "MOCA", "mPACCdigit", "mPACCtrailsB", "PTETHCAT", "PTRACCAT")
absent <- setdiff(must_clin, names(clinical))
if (length(absent)) stop("Clinical file missing: ", paste(absent, collapse = ", "), call. = FALSE)
clinical[, PTID := as.character(PTID)]
clinical[, EXAMDATE := as.Date(EXAMDATE)]
baseline <- clinical[tolower(VISCODE) %in% c("bl", "sc") & !is.na(PTID)][order(PTID, EXAMDATE)]
baseline <- baseline[, .SD[1L], by = PTID]
baseline <- baseline[, .(PTID, age_bl = as.numeric(AGE), sex = as.character(PTGENDER), education = as.numeric(PTEDUCAT),
                         dx_bl = as.character(DX_bl), APOE4 = as.numeric(APOE4),
                         ethnicity = as.character(PTETHCAT), race = as.character(PTRACCAT), baseline_date = EXAMDATE)]
baseline[, primary_european_ancestry := ethnicity == "Not Hisp/Latino" & race == "White"]

long <- merge(clinical, pgs, by.x = "PTID", by.y = "IID", all = FALSE)
long <- merge(long, baseline, by = "PTID", all.x = TRUE)
long[, time_years := as.numeric(EXAMDATE - baseline_date) / 365.25]
long <- long[is.finite(time_years) & time_years >= 0]
for (nm in c("gMAV", "gMAV_noAPOE", "gAD_noAPOE")) {
  long[, paste0(nm, "_z_all") := standardize(get(nm))]
  long[, paste0(nm, "_z_eur") := {
    z <- rep(NA_real_, .N); ii <- which(primary_european_ancestry %in% TRUE); z[ii] <- standardize(get(nm)[ii]); z
  }]
}

keep <- c("PTID", "RID", "VISCODE", "EXAMDATE", "time_years", "age_bl", "sex", "education", "dx_bl", "DX", "APOE4",
          "ethnicity", "race", "primary_european_ancestry", "gMAV", "gMAV_noAPOE", "gAD_noAPOE",
          "gMAV_z_all", "gMAV_noAPOE_z_all", "gAD_noAPOE_z_all", "gMAV_z_eur", "gMAV_noAPOE_z_eur", "gAD_noAPOE_z_eur",
          paste0("PC", 1:10), "MMSE", "ADAS13", "CDRSB", "FAQ", "MOCA", "mPACCdigit", "mPACCtrailsB")
keep <- intersect(keep, names(long))
long <- long[, ..keep]
data.table::fwrite(long, file.path(OUTDIR, "02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"), sep = "\t", na = "NA")

flow <- data.frame(
  stage = c("genotype_QC_participants", "scored_all_three_PGS", "with_PC1_PC10", "with_clinical_followup", "primary_European_ancestry_followup"),
  n_unique_participants = c(data.table::uniqueN(pgs$IID),
                            data.table::uniqueN(pgs[stats::complete.cases(pgs[, .(gMAV, gMAV_noAPOE, gAD_noAPOE)]), IID]),
                            data.table::uniqueN(pgs[stats::complete.cases(pgs[, paste0("PC", 1:10), with = FALSE]), IID]),
                            data.table::uniqueN(long$PTID),
                            data.table::uniqueN(long[primary_european_ancestry %in% TRUE, PTID])),
  stringsAsFactors = FALSE
)
write.csv(flow, file.path(OUTDIR, "03_ADNI_PGS_assembly_flow.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c(
  "10A completed: outcome-blind PGS construction and longitudinal clinical assembly.",
  "Primary score set: gMAV, gMAV-noAPOE, and gAD-noAPOE; all use fixed C+T (P<=1, r2=0.1, 500 kb).",
  "MHC is excluded from every score; APOE (GRCh37 chr19:44-46 Mb) is excluded from noAPOE scores.",
  "Clumping and scoring use fixed genotype/reference parameters selected independently of ADNI outcomes.",
  "Primary genetic analyses must use the European-ancestry subset with PC1-PC10 adjustment; all-ancestry results are sensitivity only."
), file.path(OUTDIR, "04_INTERPRETATION.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "05_sessionInfo.txt"), useBytes = TRUE)
log_line("INFO", "Step 10A completed")
