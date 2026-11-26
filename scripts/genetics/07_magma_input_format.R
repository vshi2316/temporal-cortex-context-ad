# Format MAGMA summary-statistic inputs.
if(!exists("study_path",mode="function")) source("config.R",encoding="UTF-8")
library(data.table)
prepare_magma_input <- function(path, trait, save_path) {
  x <- fread(file = path, showProgress = TRUE)
  needed <- c("SNP", "pval.outcome", "samplesize.outcome")
  if (!all(needed %in% names(x)))
    stop("Unexpected GWAS columns for ", trait)
  x <- x[, .(SNP = as.character(SNP),
            P = as.numeric(pval.outcome),
            N = as.numeric(samplesize.outcome))]
  x <- x[!is.na(SNP) & nzchar(SNP) & is.finite(P) & P > 0 & P <= 1 &
           is.finite(N) & N > 0]
  setorder(x, SNP, P)
  x <- x[!duplicated(SNP)]
  if (nrow(x) < 10000L)
    stop("Too few valid variants after deduplication for ", trait)
  path_out <- file.path(save_path, paste0(trait, "_MAGMA_input.tsv"))
  fwrite(x, path_out, sep = "\t", quote = FALSE)
  data.table(trait = trait, variants = nrow(x), median_N = median(x$N),
             input_file = path_out)
}
inputs <- c(gAD=study_path("gwas_stage","09H/ldsc_input/AD_noAPOE_TwosampleMR.txt"),
            gMAV=study_path("gwas_stage","09H/ldsc_input/gMAV_noAPOE_TwosampleMR.txt"))
out <- study_path("project","results/13H15_MTG_genetic_program_bridge_release/01_magma_gene_results")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(trait in names(inputs)) prepare_magma_input(inputs[[trait]],trait,out)
