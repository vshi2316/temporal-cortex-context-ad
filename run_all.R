# This file defines the runner. Sourcing it does not start analysis.
source("config.R",encoding="UTF-8")
AD_GROUPS <- list(
genetics=c("scripts/genetics/01_gwas_qc_ldsc.R",
    "scripts/genetics/02_genomic_factor_models.R",
    "scripts/genetics/03_multisystem_factor_gwas.R",
    "scripts/genetics/04_factor_gwas_qc.R",
    "scripts/genetics/05_genetic_correlation.R",
    "scripts/genetics/06_adni_polygenic_scores.R"),
core=c("scripts/prepare/01_structural_feature_definitions.R",
    "scripts/prepare/02_adni_mri_clinical_linkage.R",
    "scripts/prepare/03_habshd_endpoint_frames.R",
    "scripts/prepare/04_habshd_identifier_linkage.R",
    "scripts/prepare/05_adni_longitudinal_frame.R",
    "scripts/analysis/01_longitudinal_models.R",
    "scripts/analysis/02_structural_models_diagnostics.R",
    "scripts/analysis/03_adni_event_frames.R",
    "scripts/analysis/04_adni_event_models.R",
    "scripts/analysis/05_complete_case_prognostic_benchmark.R",
    "scripts/analysis/06_covariance_aware_sensitivity.R",
    "scripts/analysis/07_diagnostic_time_models.R",
    "scripts/analysis/08_stage_and_ordered_longitudinal_models.R"),
regional=c("scripts/prepare/06_regional_scan_keys.R",
    "scripts/analysis/09_aligned_regional_genetic_map.R",
    "scripts/analysis/10_paired_regional_time_models.R",
    "scripts/analysis/11_frozen_pattern_event_models.R",
    "scripts/analysis/12_pattern_clinical_increment.R"),
clinical=c("scripts/analysis/13_first_second_origin_prediction.R",
    "scripts/analysis/14_mci_brief_full_prediction.R",
    "scripts/analysis/15_current_state_prior_loss.R",
    "scripts/analysis/16_six_region_whole_cortex_benchmark.R",
    "scripts/analysis/17_habshd_mci_prediction.R",
    "scripts/analysis/18_habshd_horizon_support.R",
    "scripts/analysis/19_mci_stage_interaction.R",
    "scripts/analysis/20_same_person_cognitive_prediction.R",
    "scripts/analysis/21_assessment_components_hippocampus.R",
    "scripts/analysis/22_pathology_conditioned_trajectories.R"),
pathology=c("scripts/ancillary/01_habshd_pathology_memory.R"),
cat12=c("scripts/ancillary/02_adni_pathology_baseline_inputs.R",
    "scripts/ancillary/03_adni_pathology_index_inputs.R",
    "scripts/ancillary/05_cat12_adni_models.R",
    "scripts/ancillary/06_aibl_incident_ad.R"),
tau=c("scripts/ancillary/07_tau_regional_distribution.R",
    "scripts/ancillary/08_tau_spatial_rotations.R"),
gene_programs=c("scripts/genetics/07_magma_input_format.R",
    "scripts/ancillary/09_shared_snp_gene_programs.R",
    "scripts/ancillary/10_mtg_donor_programs.R",
    "scripts/ancillary/11_mtg_pathology_selection.R"),
imaging=c("scripts/imaging/01_mri_acquisition_sample.R",
    "scripts/imaging/02_adni_image_matching.R",
    "scripts/imaging/03_dicom_identity.R",
    "scripts/imaging/04_dicom_conversion.R",
    "scripts/imaging/05_t1_qc.R",
    "scripts/ancillary/02_adni_pathology_baseline_inputs.R",
    "scripts/ancillary/03_adni_pathology_index_inputs.R",
    "scripts/ancillary/04_cat12_adni_sample.R",
    "scripts/imaging/06_cat12_processing.R")
)
AD_GROUPS$analysis <- unlist(AD_GROUPS[c("core","regional","clinical","pathology")],use.names=FALSE)
AD_GROUPS$all <- unlist(AD_GROUPS[c("genetics","analysis","cat12","tau","gene_programs")],use.names=FALSE)

run_study <- function(group="analysis", scripts=NULL) {
  if(is.null(scripts)) {
    if(!group %in% names(AD_GROUPS)) stop("Unknown group: ",group)
    scripts <- AD_GROUPS[[group]]
  }
  for(script in scripts) {
    file <- file.path(REPO_ROOT,script)
    if(!file.exists(file)) stop("Missing repository script: ",script)
    message("Running ",script)
    sys.source(file,envir=new.env(parent=globalenv()),keep.source=TRUE)
  }
  invisible(scripts)
}
