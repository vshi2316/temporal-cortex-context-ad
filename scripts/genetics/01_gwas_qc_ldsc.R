# GWAS harmonization and five-trait LDSC.
if (!exists("study_path", mode = "function")) source("config.R", encoding = "UTF-8")
source(file.path(REPO_ROOT, "scripts/genetics/ldsc_commands.R"), local = TRUE)
options(stringsAsFactors = FALSE, scipen = 999)
library(data.table)

PROJECT_DIR <- study_path("project")
PRIMARY_DIR <- study_path("gwas_input")
RESULT_DIR <- file.path(PROJECT_DIR, "results", "09D_GWAS_QC_LDSC")
FORMAT_DIR <- file.path(RESULT_DIR, "01_formatted")
STAGE <- study_path("gwas_stage", "09D")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
manifest <- data.frame(
  trait_id = c("CMM", "INSOMNIA", "HEARING", "GRIP", "AD"),
  filename = c(
    "GCST90274723.Cardiometabolic multimorbidity.stroke disorder, type 2 diabetes mellitus, coronary artery disorder.tsv",
    "ukb-b-3957.vcf.insomnia.gz",
    "GCST90012115_buildGRCh37.Adult hearing difficulty.tsv",
    "GCST90007526_buildGRCh37.tsv.Hand grip strength.gz",
    "GCST007511.Late-onset Alzheimer's disease.Kunkle_etal_Stage1_results.txt"),
  N = c(367147, 462341, 330759, 256523, 63926),
  ncase = c(14835, NA, NA, 48596, 21982),
  ncontrol = c(352312, NA, NA, 207927, 41944),
  ancestry = "European")
manifest$path <- file.path(PRIMARY_DIR, manifest$filename)
required <- c(manifest$path, study_path("ldsc", "w_hm3.noMHC.snplist"),
  unlist(lapply(1:22, function(chr) study_path("ldsc", "eur_w_ld_chr",
                   paste0(chr, c(".l2.ldscore.gz", ".l2.M_5_50"))))))
if (any(!file.exists(required)))
  stop("Missing inputs:\n", paste(required[!file.exists(required)], collapse = "\n"))
write.csv(manifest, file.path(RESULT_DIR, "00_frozen_primary_GWAS_manifest.csv"),
          row.names = FALSE)

# Column names retain the original downstream summary-statistic contract.
maps <- list(
  CMM = c(SNP="rs_id", effect_allele="effect_allele", other_allele="other_allele",
          beta="beta", se="standard_error", pval="p_value",
          eaf="effect_allele_frequency", chr="chromosome", pos="base_pair_location"),
  HEARING = c(SNP="variant_id", effect_allele="a1", other_allele="a2",
              beta="beta", se="standard_error", pval="p_value"),
  GRIP = c(SNP="variant_id", effect_allele="effect_allele", other_allele="other_allele",
           beta="beta", se="standard_error", pval="p_value",
           eaf="effect_allele_frequency", chr="chromosome", pos="base_pair_location"),
  AD = c(SNP="MarkerName", effect_allele="Effect_allele", other_allele="Non_Effect_allele",
         beta="Beta", se="SE", pval="Pvalue"))

clean_summary <- function(d, trait, N) {
  for (n in intersect(c("beta", "se", "pval", "eaf", "pos"), names(d)))
    set(d, j = n, value = suppressWarnings(as.numeric(d[[n]])))
  d[, `:=`(SNP = as.character(SNP), effect_allele = toupper(effect_allele),
            other_allele = toupper(other_allele))]
  d <- d[!is.na(SNP) & nzchar(SNP) & !is.na(effect_allele) & nzchar(effect_allele) &
    !is.na(other_allele) & nzchar(other_allele) & is.finite(beta) &
    is.finite(se) & se > 0 & is.finite(pval) & pval >= 0 & pval <= 1]
  if (trait %in% c("CMM", "GRIP", "INSOMNIA")) {
    if (trait != "INSOMNIA") {
      d <- d[grepl("rs", SNP) & !grepl(",", SNP) &
             !grepl("[ID0-9]", effect_allele)]
      setorder(d, pval)
    }
    d[, SNP := tolower(gsub("[[:space:]]", "", SNP))]
    d <- d[!duplicated(SNP)]
    d[, pval := pmax(pval, 1e-200)]
  } else d[pval == 0, pval := 1e-300]
  if (!nrow(d)) stop("No valid variants for ", trait)
  setnames(d, setdiff(names(d), "SNP"), paste0(setdiff(names(d), "SNP"), ".outcome"))
  d[, `:=`(samplesize.outcome = as.numeric(N), Phenotype = trait, id.outcome = trait)]
  d
}

# The frozen GWAS-VCF contains one study. ALT is its effect allele and
# FORMAT/LP is -log10(P). Decode FORMAT by field name in bounded chunks.
read_insomnia <- function(path, N) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("VCF header missing: ", path)
    if (startsWith(line, "#CHROM")) break
  }
  if (length(strsplit(line, "\t", fixed = TRUE)[[1]]) != 10L)
    stop("Expected the single-study insomnia GWAS-VCF.")
  parts <- list()
  repeat {
    lines <- readLines(con, n = 200000L, warn = FALSE)
    if (!length(lines)) break
    x <- as.data.table(tstrsplit(lines, "\t", fixed = TRUE))
    setnames(x, c("chr", "pos", "SNP", "other_allele", "effect_allele",
                  "QUAL", "FILTER", "INFO", "FORMAT", "STUDY"))
    x[, `:=`(beta = NA_real_, se = NA_real_, pval = NA_real_, eaf = NA_real_)]
    for (fmt in unique(x$FORMAT)) {
      ix <- which(x$FORMAT == fmt)
      keys <- strsplit(fmt, ":", fixed = TRUE)[[1]]
      if (!all(c("ES", "SE", "LP") %in% keys)) stop("Missing VCF effect fields.")
      fields <- tstrsplit(x$STUDY[ix], ":", fixed = TRUE)
      for (pair in list(c("beta", "ES"), c("se", "SE"),
                        c("pval", "LP"), c("eaf", "AF"))) {
        j <- match(pair[2], keys)
        if (!is.na(j)) {
          value <- suppressWarnings(as.numeric(fields[[j]]))
          if (pair[2] == "LP") value <- 10^(-value)
          set(x, i = ix, j = pair[1], value = value)
        }
      }
    }
    parts[[length(parts) + 1L]] <- x[, .(SNP, effect_allele, other_allele,
                                         beta, se, pval, eaf, chr, pos)]
  }
  clean_summary(rbindlist(parts), "INSOMNIA", N)
}

formatted <- setNames(character(nrow(manifest)), manifest$trait_id)
for (i in seq_len(nrow(manifest))) {
  trait <- manifest$trait_id[i]
  if (trait == "INSOMNIA") {
    d <- read_insomnia(manifest$path[i], manifest$N[i])
  } else {
    map <- maps[[trait]]
    d <- fread(manifest$path[i], select = unname(map), showProgress = TRUE)
    setnames(d, unname(map), names(map))
    d <- clean_summary(d, trait, manifest$N[i])
  }
  outdir <- file.path(FORMAT_DIR, trait)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  formatted[[trait]] <- file.path(outdir, paste0(trait, "_TwosampleMR.txt"))
  fwrite(d, formatted[[trait]], sep = "\t", na = "NA")
  rm(d)
  gc(verbose = FALSE)
}
write.csv(data.frame(trait = names(formatted), path = unname(formatted)),
          file.path(RESULT_DIR, "01_format_manifest.csv"), row.names = FALSE)

sumstats <- setNames(character(nrow(manifest)), manifest$trait_id)
for (i in seq_len(nrow(manifest))) {
  trait <- manifest$trait_id[i]
  sumstats[[trait]] <- ldsc_munge(formatted[[trait]], trait, manifest$N[i],
                                file.path(STAGE, "sumstats", trait))
  ldsc_fit(sumstats[[trait]], paste0(trait, "_h2"), file.path(STAGE, "h2", trait), "h2")
}
write.csv(data.frame(trait = names(sumstats), path = unname(sumstats), N = manifest$N),
          file.path(RESULT_DIR, "02_LDSC_sumstats_audit.csv"), row.names = FALSE)
# The original LDSC multi-trait call compares the first trait with each remaining trait.
ldsc_fit(unname(sumstats), "MAV_AD_primary_five_trait_rg",
         file.path(STAGE, "rg"), "rg")
# Keep Python command paths in the configured staging directory.
for (branch in c(h2 = "03_LDSC_h2", rg = "04_LDSC_rg")) {
  key <- if (branch == "03_LDSC_h2") "h2" else "rg"
  origin <- file.path(STAGE, key)
  target <- file.path(RESULT_DIR, branch)
  files <- list.files(origin, full.names = TRUE, recursive = TRUE)
  for (f in files) {
    to <- file.path(target, substring(f, nchar(origin) + 2L))
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(f, to, overwrite = TRUE)) stop("Failed to copy LDSC output: ", f)
  }
}
ld <- paste0(normalizePath(study_path("ldsc", "eur_w_ld_chr"), winslash = "/"), "/")
covstruc <- GenomicSEM::ldsc(traits = unname(sumstats),
  sample.prev = rep(NA_real_, 5), population.prev = rep(NA_real_, 5),
  ld = ld, wld = ld, trait.names = manifest$trait_id, sep_weights = FALSE,
  chr = 22, n.blocks = 200,
  ldsc.log = file.path(RESULT_DIR, "04_LDSC_rg", "GenomicSEM_primary_ldsc"),
  stand = FALSE, select = FALSE, chisq.max = NA)
saveRDS(covstruc, file.path(RESULT_DIR, "05_primary_five_trait_covstruc.rds"))
for (part in c("S", "V", "I")) {
  filename <- c(S="05A_genetic_covariance_S.csv", V="05B_sampling_covariance_V.csv",
                I="05C_LDSC_intercepts_I.csv")[[part]]
  if (!is.null(covstruc[[part]])) write.csv(covstruc[[part]], file.path(RESULT_DIR, filename))
}
writeLines(capture.output(sessionInfo()), file.path(RESULT_DIR, "07_sessionInfo.txt"))
