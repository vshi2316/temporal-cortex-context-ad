# Shared-SNP MAGMA gene-priority definitions
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

## Same-SNP MAGMA sensitivity for the gAD and gMAV priority programs.
## Uses the frozen shared SNP intersection and 1000 Genomes European reference.

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(data.table))
magma <- study_path("magma")
if (!file.exists(magma)) stop("Configure the MAGMA executable in config.local.R.")

root <- study_path("project")
prior_dir <- file.path(root, "results",
  "13H15_MTG_genetic_program_bridge_release")
prior_gwas_dir <- file.path(prior_dir, "01_magma_gene_results")
gwas_files <- c(
  gAD = file.path(prior_gwas_dir, "gAD_MAGMA_input.tsv"),
  gMAV = file.path(prior_gwas_dir, "gMAV_MAGMA_input.tsv")
)
ref_dir <- study_path("magma_reference")
ld_bfile <- file.path(ref_dir, "g1000_eur", "g1000_eur")
gene_loc <- file.path(ref_dir, "ENSGv110.coding.genes.txt")
required <- c(gwas_files, paste0(ld_bfile, c(".bed", ".bim", ".fam")),
              gene_loc, magma)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing input(s):\n", paste(missing, collapse = "\n"))

out_dir <- file.path(root, "results", paste0(
  "13H16_harmonized_gAD_gMAV_MAGMA_", "release"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_pval <- function(path, trait) {
  x <- fread(file = path, select = c("SNP", "P", "N"), showProgress = TRUE)
  if (!all(c("SNP", "P", "N") %in% names(x)))
    stop("Expected SNP, P and N in ", path)
  x[, SNP := as.character(SNP)]
  if (anyDuplicated(x$SNP)) stop("Duplicated SNP IDs in ", trait)
  setnames(x, c("P", "N"), paste0(c("P_", "N_"), trait))
  x
}

ad <- read_pval(gwas_files[["gAD"]], "gAD")
mav <- read_pval(gwas_files[["gMAV"]], "gMAV")
input_counts <- c(gAD = nrow(ad), gMAV = nrow(mav))
common <- merge(ad, mav, by = "SNP", all = FALSE, sort = TRUE)
rm(ad, mav)
gc(verbose = FALSE)
if (nrow(common) < 100000L)
  stop("Fewer than 100,000 shared SNP IDs; inspect GWAS variant IDs")

qc <- data.table(
  trait = names(input_counts),
  input_snps = as.integer(input_counts),
  shared_snps = nrow(common),
  retained_fraction = as.numeric(nrow(common) / input_counts),
  median_N_shared = c(median(common$N_gAD), median(common$N_gMAV))
)
fwrite(qc, file.path(out_dir, "01_shared_SNP_input_audit.csv"), bom = TRUE)

write_trait <- function(trait) {
  input_file <- file.path(out_dir, paste0(trait, "_shared_SNP_MAGMA_input.tsv"))
  p_col <- paste0("P_", trait)
  n_col <- paste0("N_", trait)
  x <- common[, .(SNP, P = get(p_col), N = get(n_col))]
  fwrite(x, input_file, sep = "\t", quote = FALSE)
  input_file
}
shared_files <- setNames(vapply(names(gwas_files), write_trait, character(1)),
                         names(gwas_files))
rm(common)
gc(verbose = FALSE)

run_magma <- function(trait) {
  trait_dir <- file.path(out_dir, trait)
  dir.create(trait_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(trait_dir, trait)
  # The frozen annotation used gene boundaries without upstream/downstream extension.
  status <- system2(magma, args = shQuote(c("--annotate", "--snp-loc",
    paste0(ld_bfile, ".bim"), "--gene-loc", gene_loc, "--out", prefix)))
  if (status != 0L) stop("MAGMA annotation failed for ", trait)
  status <- system2(magma, args = shQuote(c("--bfile", ld_bfile, "--pval",
    shared_files[[trait]], "use=1,2", "ncol=N", "--gene-annot",
    paste0(prefix, ".genes.annot"), "--out", prefix)))
  if (status != 0L) stop("MAGMA gene analysis failed for ", trait)
  result_file <- paste0(prefix, ".genes.out")
  if (!file.exists(result_file)) stop("MAGMA output missing for ", trait)
  result_file
}
gene_files <- setNames(vapply(names(gwas_files), run_magma, character(1)),
                       names(gwas_files))

read_genes <- function(path, trait) {
  x <- fread(file = path, skip = "GENE")
  if (!all(c("GENE", "P", "N") %in% names(x)))
    stop("Unexpected MAGMA gene output: ", path)
  x <- x[, .(gene_id = as.character(GENE), P = as.numeric(P),
            N = as.numeric(N))]
  x <- x[!duplicated(gene_id) & is.finite(P) & P > 0 & P <= 1]
  setnames(x, c("P", "N"), paste0(c("P_", "N_"), trait))
  x
}

ad_genes <- read_genes(gene_files[["gAD"]], "gAD")
mav_genes <- read_genes(gene_files[["gMAV"]], "gMAV")
gene <- merge(ad_genes, mav_genes, by = "gene_id", all = FALSE)
gene[, `:=`(
  rank_gAD = frank(P_gAD, ties.method = "first"),
  rank_gMAV = frank(P_gMAV, ties.method = "first"),
  q_gAD = p.adjust(P_gAD, method = "BH"),
  q_gMAV = p.adjust(P_gMAV, method = "BH")
)]
gene[, program := fifelse(rank_gAD <= 200L & rank_gMAV <= 200L,
  "shared_top200", fifelse(rank_gAD <= 200L, "gAD_priority_top200",
                   fifelse(rank_gMAV <= 200L, "gMAV_priority_top200",
                           "background")))]
fwrite(gene, file.path(out_dir, "02_shared_SNP_gene_program_membership.csv"),
       bom = TRUE)

writeLines(capture.output(sessionInfo()),file.path(out_dir,"99_sessionInfo.txt"))
