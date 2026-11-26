# Variant quality control and scoring weights
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

#!/usr/bin/env Rscript

# Filter variant-level factor-GWAS results and export scoring weights.

options(stringsAsFactors = FALSE, scipen = 999, width = 180, warn = 1)

PROJECT_DIR <- study_path("project")
INFILE <- file.path(PROJECT_DIR, "results", "09F_MAV_commonfactorGWAS",
                    "03_MAV_commonfactorGWAS.tsv")
OUTDIR <- file.path(PROJECT_DIR, "results", "09G_MAV_factorGWAS_QC")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE))
  stop("Package data.table is required.", call. = FALSE)
if (!file.exists(INFILE)) stop("Missing 09F output: ", INFILE, call. = FALSE)

log_file <- file.path(OUTDIR, "09G_MAV_factorGWAS_QC.log")
if (file.exists(log_file)) unlink(log_file, force = TRUE)
log_line <- function(level = "INFO", ...) {
  z <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               level, paste0(..., collapse = ""))
  cat(z, "\n"); cat(z, "\n", file = log_file, append = TRUE)
}

DT <- data.table::fread(file = INFILE, sep = "\t", showProgress = TRUE,
                        na.strings = c("", "NA", "NaN", "NULL"))
req <- c("SNP", "A1", "A2", "MAF", "est", "se_c", "Pval_Estimate",
         "Q", "Q_df", "Q_pval", "fail")
miss <- setdiff(req, names(DT))
if (length(miss)) stop("09F output is missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

num_cols <- intersect(c("MAF", "est", "se_c", "Z_Estimate", "Pval_Estimate",
                        "Q", "Q_df", "Q_pval"), names(DT))
for (j in num_cols) DT[, (j) := suppressWarnings(as.numeric(get(j)))]
DT[, SNP := as.character(SNP)]
DT[, A1 := toupper(as.character(A1))]
DT[, A2 := toupper(as.character(A2))]

# Variant-level quality checks.
DT[, duplicate_snp := duplicated(SNP)]
DT[, valid_alleles := A1 %in% c("A", "C", "G", "T") &
     A2 %in% c("A", "C", "G", "T") & A1 != A2]
DT[, valid_maf := is.finite(MAF) & MAF >= 0.01 & MAF <= 0.5]
DT[, valid_effect := is.finite(est) & is.finite(se_c) & se_c > 0 & is.finite(Pval_Estimate) &
     Pval_Estimate >= 0 & Pval_Estimate <= 1]
DT[, p_recomputed := 2 * pnorm(-abs(est / se_c))]
DT[, p_log10 := -log10(pmax(Pval_Estimate, .Machine$double.xmin))]
DT[, genomewide_sig := Pval_Estimate < 5e-8]
DT[, heterogeneity_sig := is.finite(Q_pval) & Q_pval < 5e-8]
DT[, fail_zero := as.character(fail) %in% c("0", "FALSE", "")]

qc <- data.table::data.table(
  metric = c("input_rows", "unique_snps", "duplicate_rows", "invalid_alleles",
             "invalid_maf", "invalid_effect", "fail_nonzero", "warnings_nonempty",
             "optimizer_not_converged_warning_rows",
             "genomewide_significant", "heterogeneity_significant"),
  value = c(nrow(DT), data.table::uniqueN(DT$SNP), sum(DT$duplicate_snp, na.rm = TRUE),
            sum(!DT$valid_alleles, na.rm = TRUE), sum(!DT$valid_maf, na.rm = TRUE),
            sum(!DT$valid_effect, na.rm = TRUE), sum(!DT$fail_zero, na.rm = TRUE),
            if ("warning" %in% names(DT)) sum(!is.na(DT$warning) & trimws(DT$warning) != "") else NA_integer_,
            if ("warning" %in% names(DT)) sum(grepl("optimizer warns|solution has NOT been found",
                                                     DT$warning, ignore.case = TRUE, fixed = FALSE), na.rm = TRUE) else NA_integer_,
            sum(DT$genomewide_sig, na.rm = TRUE), sum(DT$heterogeneity_sig, na.rm = TRUE))
)
data.table::fwrite(qc, file.path(OUTDIR, "01_gMAV_factorGWAS_QC_summary.csv"))

# Valid variants retain their Q-heterogeneity annotations.
VALID <- DT[!duplicate_snp & valid_alleles & valid_maf & valid_effect & fail_zero]
data.table::fwrite(VALID, file.path(OUTDIR, "02_gMAV_factorGWAS_valid_rows.tsv"), sep = "\t", na = "NA")

LEAD <- VALID[genomewide_sig == TRUE][order(Pval_Estimate)]
data.table::fwrite(LEAD, file.path(OUTDIR, "03_gMAV_genomewide_significant_SNPs.tsv"), sep = "\t", na = "NA")
HET <- VALID[heterogeneity_sig == TRUE][order(Q_pval)]
data.table::fwrite(HET, file.path(OUTDIR, "04_gMAV_heterogeneity_flagged_SNPs.tsv"), sep = "\t", na = "NA")

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "08_sessionInfo.txt"), useBytes = TRUE)
log_line("INFO", "Step 09G completed; retained ", nrow(VALID), " SNP rows and ", nrow(LEAD), " genome-wide significant SNPs.")
