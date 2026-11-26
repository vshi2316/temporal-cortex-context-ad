# gMAV and AD genetic correlation including APOE exclusion
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

#!/usr/bin/env Rscript

# Estimate gMAV/AD genetic correlation and APOE-region sensitivity.

options(stringsAsFactors = FALSE, scipen = 999, width = 180, warn = 1)

PROJECT_DIR <- study_path("project")
IN_09G <- file.path(PROJECT_DIR, "results", "09G_MAV_factorGWAS_QC",
                    "02_gMAV_factorGWAS_valid_rows.tsv")
AD_FORMATTED <- file.path(PROJECT_DIR, "results", "09D_GWAS_QC_LDSC",
                          "01_formatted", "AD", "AD_TwosampleMR.txt")
AD_RAW <- study_path("gwas_input",
                    "GCST007511.Late-onset Alzheimer's disease.Kunkle_etal_Stage1_results.txt")
LD_ROOT <- study_path("ldsc")
LD_DIR <- file.path(LD_ROOT, "eur_w_ld_chr")
HM3_FILE <- file.path(LD_ROOT, "w_hm3.noMHC.snplist")
OUTDIR <- file.path(PROJECT_DIR, "results", "09H_gMAV_AD_LDSC")
STAGE <- study_path("gwas_stage", "09H")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAGE, recursive = TRUE, showWarnings = FALSE)
for (x in c("ldsc_input", "sumstats", "h2", "rg"))
  dir.create(file.path(STAGE, x), recursive = TRUE, showWarnings = FALSE)

need <- c(IN_09G, AD_FORMATTED, AD_RAW, HM3_FILE)
if (any(!file.exists(need))) stop("Missing required input(s):\n", paste(need[!file.exists(need)], collapse = "\n"), call. = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.", call. = FALSE)
source(file.path(REPO_ROOT, "scripts/genetics/ldsc_commands.R"), local = TRUE)
if (!requireNamespace("GenomicSEM", quietly = TRUE)) stop("GenomicSEM is required.", call. = FALSE)

LOG <- file.path(OUTDIR, "09H_gMAV_AD_LDSC.log")
if (file.exists(LOG)) unlink(LOG, force = TRUE)
log_line <- function(level = "INFO", ...) {
  z <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(..., collapse = ""))
  cat(z, "\n"); cat(z, "\n", file = LOG, append = TRUE)
}

# Primary LDSC uses the smallest contributing GWAS sample size.
# The sensitivity analysis uses the median contributing sample size.
N_TRAITS <- c(CMM = 367147, INSOMNIA = 462341, HEARING = 330759, GRIP = 256523)
N_GMAV_PRIMARY <- min(N_TRAITS)
N_GMAV_SENS <- median(N_TRAITS)
N_AD <- 63926

g <- data.table::fread(IN_09G, sep = "\t", showProgress = TRUE)
req <- c("SNP", "A1", "A2", "MAF", "est", "se_c", "Pval_Estimate")
if (length(setdiff(req, names(g)))) stop("09G file lacks required columns.", call. = FALSE)
g <- g[is.finite(MAF) & MAF >= .01 & MAF <= .5 & is.finite(est) & is.finite(se_c) & se_c > 0 &
         is.finite(Pval_Estimate) & Pval_Estimate >= 0 & Pval_Estimate <= 1]
g <- g[!duplicated(SNP)]

gmr <- g[, .(SNP = as.character(SNP), effect_allele.outcome = toupper(A1),
             other_allele.outcome = toupper(A2), beta.outcome = as.numeric(est),
             se.outcome = as.numeric(se_c), pval.outcome = as.numeric(Pval_Estimate),
             samplesize.outcome = N_GMAV_PRIMARY, Phenotype = "gMAV", id.outcome = "gMAV")]
gmr_sens <- data.table::copy(gmr); gmr_sens[, samplesize.outcome := N_GMAV_SENS]
gmr_path <- file.path(STAGE, "ldsc_input", "gMAV_TwosampleMR.txt")
gmr_sens_path <- file.path(STAGE, "ldsc_input", "gMAV_Nmedian_TwosampleMR.txt")
ad_stage_path <- file.path(STAGE, "ldsc_input", "AD_TwosampleMR.txt")
if (!file.copy(AD_FORMATTED, ad_stage_path, overwrite = TRUE))
  stop("Failed to stage AD TwoSampleMR input at ASCII path.", call. = FALSE)
data.table::fwrite(gmr, gmr_path, sep = "\t", na = "NA")
data.table::fwrite(gmr_sens, gmr_sens_path, sep = "\t", na = "NA")

# Identify the APOE region in the original GRCh37 AD file and remove matching
# rsIDs from both traits for a transparent no-APOE sensitivity analysis.
ad_raw <- data.table::fread(AD_RAW, select = c("Chromosome", "Position", "MarkerName"), showProgress = TRUE)
apoe <- ad_raw[as.character(Chromosome) %in% c("19", "chr19") & Position >= 44000000 & Position <= 46000000, MarkerName]
apoe <- unique(as.character(apoe))
if (!length(apoe)) stop("No APOE-region SNPs found in AD raw file; check genome build/columns.", call. = FALSE)
ad <- data.table::fread(AD_FORMATTED, showProgress = TRUE)
if (!all(c("SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "se.outcome", "pval.outcome") %in% names(ad)))
  stop("AD formatted file is not in expected TwoSampleMR format.", call. = FALSE)
ad_no <- ad[!SNP %in% apoe]
g_no <- gmr[!SNP %in% apoe]
ad_no_path <- file.path(STAGE, "ldsc_input", "AD_noAPOE_TwosampleMR.txt")
g_no_path <- file.path(STAGE, "ldsc_input", "gMAV_noAPOE_TwosampleMR.txt")
data.table::fwrite(ad_no, ad_no_path, sep = "\t", na = "NA")
data.table::fwrite(g_no, g_no_path, sep = "\t", na = "NA")

write.csv(data.frame(
  analysis = c("primary", "N_sensitivity", "noAPOE"),
  gMAV_N = c(N_GMAV_PRIMARY, N_GMAV_SENS, N_GMAV_PRIMARY),
  AD_N = N_AD,
  apoe_region = c("retained", "retained", "chr19:44-46 Mb GRCh37 excluded"),
  gMAV_rows = c(nrow(gmr), nrow(gmr_sens), nrow(g_no)),
  AD_rows = c(nrow(ad), nrow(ad), nrow(ad_no)),
  stringsAsFactors = FALSE
), file.path(OUTDIR, "01_analysis_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")

make_sumstats <- ldsc_munge

log_line("INFO", "Formatting LDSC files")
ss_g <- make_sumstats(gmr_path, "gMAV", N_GMAV_PRIMARY, file.path(STAGE, "sumstats", "gMAV"))
ss_g_med <- make_sumstats(gmr_sens_path, "gMAV_Nmedian", N_GMAV_SENS, file.path(STAGE, "sumstats", "gMAV_Nmedian"))
ss_ad <- make_sumstats(ad_stage_path, "AD", N_AD, file.path(STAGE, "sumstats", "AD"))
ss_g_no <- make_sumstats(g_no_path, "gMAV_noAPOE", N_GMAV_PRIMARY, file.path(STAGE, "sumstats", "gMAV_noAPOE"))
ss_ad_no <- make_sumstats(ad_no_path, "AD_noAPOE", N_AD, file.path(STAGE, "sumstats", "AD_noAPOE"))

run_rg <- function(a, b, name) {
  ldsc_fit(c(a, b), name, file.path(STAGE, "rg", name), "rg")
}
run_h2 <- function(a, name) {
  ldsc_fit(a, name, file.path(STAGE, "h2", name), "h2")
}

run_h2(ss_g, "gMAV_h2"); run_h2(ss_ad, "AD_h2"); run_h2(ss_g_no, "gMAV_noAPOE_h2"); run_h2(ss_ad_no, "AD_noAPOE_h2")
run_rg(ss_g, ss_ad, "gMAV_AD_rg_primary")
run_rg(ss_g_med, ss_ad, "gMAV_AD_rg_Nmedian_sensitivity")
run_rg(ss_g_no, ss_ad_no, "gMAV_AD_rg_noAPOE")

# Native two-trait GenomicSEM covariance objects are retained for reproducibility.
ld_path <- paste0(normalizePath(LD_DIR, winslash = "/", mustWork = TRUE), "/")
run_native <- function(a, b, nameset, tag) {
  x <- GenomicSEM::ldsc(traits = c(a, b), sample.prev = c(NA_real_, NA_real_),
                        population.prev = c(NA_real_, NA_real_), ld = ld_path, wld = ld_path,
                        trait.names = nameset, sep_weights = FALSE, chr = 22, n.blocks = 200,
                        ldsc.log = file.path(OUTDIR, paste0(tag, "_GenomicSEM_ldsc")),
                        stand = FALSE, select = FALSE, chisq.max = NA)
  saveRDS(x, file.path(OUTDIR, paste0(tag, "_covstruc.rds")))
  if (!is.null(x$S)) write.csv(x$S, file.path(OUTDIR, paste0(tag, "_S.csv")), row.names = TRUE)
  if (!is.null(x$V)) write.csv(x$V, file.path(OUTDIR, paste0(tag, "_V.csv")), row.names = TRUE)
  x
}
run_native(ss_g, ss_ad, c("gMAV", "AD"), "02_primary")
run_native(ss_g_no, ss_ad_no, c("gMAV_noAPOE", "AD_noAPOE"), "03_noAPOE")

writeLines(c(
  "09H completed.",
  "Primary analysis: gMAV versus AD using conservative minimum contributing N.",
  "Sensitivity: median contributing N and exclusion of GRCh37 chr19:44-46 Mb APOE region."

), file.path(OUTDIR, "04_INTERPRETATION.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "05_sessionInfo.txt"), useBytes = TRUE)
log_line("INFO", "Step 09H completed")
