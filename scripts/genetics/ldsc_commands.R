# Direct calls to the original LDSC command-line programs.
# Shared by GWAS preparation and the gMAV/AD correlation analysis.
ldsc_command <- function(program, args) {
  script <- study_path("ldsc_code", program)
  python <- study_path("ldsc_python")
  if (!file.exists(script) || !file.exists(python))
    stop("Configure the LDSC checkout and its Python interpreter in config.local.R.")
  status <- system2(python, args = shQuote(c(script, args)))
  if (status != 0L) stop("LDSC failed: ", program, " (exit ", status, ")")
  invisible(status)
}

ldsc_munge <- function(input, trait, N, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  x <- data.table::fread(input, showProgress = FALSE)
  needed <- c("SNP", "effect_allele.outcome", "other_allele.outcome",
              "beta.outcome", "se.outcome", "pval.outcome")
  if (!all(needed %in% names(x))) stop("Missing GWAS columns in ", input)
  d <- x[, .(SNP, A1 = effect_allele.outcome, A2 = other_allele.outcome,
             Z = beta.outcome / se.outcome, P = pval.outcome, N = as.numeric(N))]
  if ("eaf.outcome" %in% names(x) &&
      data.table::uniqueN(x$eaf.outcome, na.rm = TRUE) > 1L)
    d[, MAF := pmin(x$eaf.outcome, 1 - x$eaf.outcome)]
  staged <- file.path(outdir, paste0(trait, "_LDSC_input.tsv"))
  data.table::fwrite(d, staged, sep = "\t", na = "NA")
  prefix <- file.path(outdir, trait)
  args <- c("--sumstats", staged, "--merge-alleles",
            study_path("ldsc", "w_hm3.noMHC.snplist"), "--maf-min", "0.01",
            "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
            "--N-col", "N", "--signed-sumstats", "Z,0", "--out", prefix)
  if ("MAF" %in% names(d)) args <- c(args, "--frq", "MAF")
  ldsc_command("munge_sumstats.py", args)
  output <- paste0(prefix, ".sumstats.gz")
  if (!file.exists(output)) stop("Missing LDSC output: ", output)
  normalizePath(output, winslash = "/", mustWork = TRUE)
}

ldsc_fit <- function(sumstats, name, outdir, type = c("h2", "rg")) {
  type <- match.arg(type)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ld <- paste0(normalizePath(study_path("ldsc", "eur_w_ld_chr"),
                            winslash = "/", mustWork = TRUE), "/")
  ldsc_command("ldsc.py", c(paste0("--", type), paste(sumstats, collapse = ","),
               "--ref-ld-chr", ld, "--w-ld-chr", ld,
               "--out", file.path(outdir, name)))
}
