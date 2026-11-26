# Source from the repository root. Private locations stay in config.local.R.
REPO_ROOT <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
if (!file.exists(file.path(REPO_ROOT,"config.R"))) stop("Set the working directory to the repository root.")
private <- file.path(REPO_ROOT,"data","private")
AD_PATHS <- list(
  project=file.path(private,"project"), cohorts=file.path(private,"cohorts"),
  gwas_stage=file.path(private,"gwas_stage"), ldsc=file.path(private,"ldsc"),
  magma_reference=file.path(private,"magma_reference"),
  ldsc_code=file.path(private,"tools","ldsc"),
  cat12=file.path(private,"tools","CAT12"),
  matlab_runtime=file.path(private,"tools","MATLAB_Runtime"),
  dcm2niix=file.path(private,"tools","dcm2niix"),
  mri_stage=file.path(private,"mri_stage"), mri_processed=file.path(private,"mri_processed"),
  cat12_stage=file.path(private,"cat12_stage"), scratch=file.path(private,"scratch")
)
local_config <- file.path(REPO_ROOT,"config.local.R")
if(file.exists(local_config)) source(local_config,local=TRUE,encoding="UTF-8")
if(is.null(AD_PATHS$ldsc_python)) AD_PATHS$ldsc_python <- file.path(private,"tools","ldsc-env",if(.Platform$OS.type == "windows") "python.exe" else "bin/python")
if(is.null(AD_PATHS[["magma"]])) AD_PATHS$magma <- file.path(private,"tools","magma",if(.Platform$OS.type == "windows") "magma.exe" else "magma")
if(is.null(AD_PATHS$adni_data)) AD_PATHS$adni_data <- file.path(AD_PATHS$cohorts,"ADNI","data")
if(is.null(AD_PATHS$habs_release)) AD_PATHS$habs_release <- file.path(AD_PATHS$cohorts,"HABS-HD")
if(is.null(AD_PATHS$gwas_input)) AD_PATHS$gwas_input <- file.path(AD_PATHS$project,"gwas","primary")
if(is.null(AD_PATHS$aibl_data)) AD_PATHS$aibl_data <- file.path(AD_PATHS$cohorts,"AIBL")
if(is.null(AD_PATHS$aibl_mri_manifest)) AD_PATHS$aibl_mri_manifest <- file.path(AD_PATHS$project,"mri","AIBL_BL_T1_MPRAGE_manifest.csv")
if(is.null(AD_PATHS$adni_unavailable_ptids)) AD_PATHS$adni_unavailable_ptids <- file.path(private,"exclusions","adni_unavailable_ptids.csv")
if(is.null(AD_PATHS$cat12_config)) AD_PATHS$cat12_config <- file.path(AD_PATHS$scratch,"cat12_config")
if(is.null(AD_PATHS$plink)) AD_PATHS$plink <- file.path(private,"tools",if(.Platform$OS.type == "windows") "plink.exe" else "plink")
if(is.null(AD_PATHS$dcm2niix_executable)) AD_PATHS$dcm2niix_executable <- file.path(AD_PATHS$dcm2niix,if(.Platform$OS.type == "windows") "dcm2niix.exe" else "dcm2niix")
if(is.null(AD_PATHS$cat12_executable)) AD_PATHS$cat12_executable <- file.path(AD_PATHS$cat12,"spm25.exe")
if(is.null(AD_PATHS$genotype_prefix)) AD_PATHS$genotype_prefix <- file.path(private,"genotypes","adni_qc_final")
if(is.null(AD_PATHS$genetic_pcs)) AD_PATHS$genetic_pcs <- file.path(private,"genotypes","ADNI_genetic_PCs.tsv")
if(is.null(AD_PATHS$aibl_clinical)) AD_PATHS$aibl_clinical <- file.path(AD_PATHS$aibl_data,"AIBL_Integrated_Data_FIXED.csv")
if(is.null(AD_PATHS$mtg_metadata)) AD_PATHS$mtg_metadata <- file.path(AD_PATHS$project,"GSE297004","GSE297004_metadata_all_neurons.rds")
if(is.null(AD_PATHS$mtg_counts)) AD_PATHS$mtg_counts <- file.path(AD_PATHS$project,"GSE297004","GSE297004_counts_all_neurons.rds")

study_path <- function(key, ...) {
  value <- AD_PATHS[[key]]
  if(is.null(value)||length(value)!=1L||!nzchar(value)) stop("Configure AD_PATHS$",key," in config.local.R")
  file.path(value,...)
}
