# Temporal cortex context in Alzheimer disease

Code for **Clinical stage and assessment shape genetic associations and prognostic value of temporal cortical thickness in Alzheimer disease**.

The analyses compare Alzheimer disease and multisystem genetic scores (gAD and gMAV) against a fixed six-region temporal cortical thickness measure. The Alzheimer’s Disease Neuroimaging Initiative (ADNI) links genetic thinning contrasts to subsequent dementia events. Matched assessment and anatomical benchmarks quantify the additional prognostic value of temporal thickness in ADNI and the Health and Aging Brain Study Health Disparities (HABS-HD). Reported supplementary analyses are included.

## Setup and execution

Use R 4.4.2 as the reference release. Recorded stage-specific environments are in `package-versions.txt`; `dependencies.csv` lists direct dependencies and their installation sources.

1. Obtain the cohort files and reference resources listed in [docs/INPUTS.md](docs/INPUTS.md).
2. Set the input locations in `config.R`. Private overrides can be placed in a local `config.local.R` file, which Git ignores.
3. Start R in this repository directory.

```r
source("setup.R", encoding = "UTF-8")
check_dependencies()
# Install missing CRAN/Bioconductor packages when needed.
# install_cran_dependencies()
# install_bioconductor_dependencies()
source("run_all.R", encoding = "UTF-8")
```

GenomicSEM 0.0.5 constructs the signed four-trait gMAV factor and estimates its SNP effects. The factor excludes Alzheimer disease from its indicators. Its SNP effects generate the gMAV participant score for comparison with gAD. The native functions are `GenomicSEM::ldsc()`, `usermodel()`, `sumstats()` and `commonfactorGWAS()`.

Install GenomicSEM with `install_genomicsem()` from `setup.R`. Install the original [LDSC](https://github.com/bulik/ldsc) scripts in their own Python environment following that repository's `environment.yml`. Obtain [MAGMA](https://ctg.cncr.nl/software/magma) 1.10 and the matching European reference files for the gene-priority branch. Configure `ldsc_code`, `ldsc_python` and `magma` in `config.local.R`. GWAS parsing uses `data.table`; LDSC and MAGMA are called directly through their command-line interfaces.

Loading the runner defines execution groups. Select the required group explicitly.

```r
run_study("genetics")
run_study("analysis")
```

`analysis` executes `core`, `regional`, `clinical`, then `pathology`. These groups can also be selected individually to reuse completed upstream results. The additional statistical branches are `cat12`, `tau` and `gene_programs`. `run_study("all")` runs genetics and all statistical branches from the documented inputs, including completed CAT12 derivatives. Raw-image processing is a separate `imaging` group. The `cat12` group reads completed derivatives and the prepared bridge manifest without checking original NIfTI files. The `imaging` group prepares that manifest before CAT12 processing. Previously generated raw-image QC flags can be supplied for the reported exclusion sensitivities. Without that file, those sensitivities are skipped and the primary models still run.

Each module writes its model estimates and numerical summaries under `AD_PATHS$project/results`. Publication figure and table production is outside this repository.

## Analysis and manuscript map

Script paths are relative to `scripts/`. The named outputs are written to each module's result directory under `AD_PATHS$project/results`.

| Paper content | Scripts | Principal outputs |
| --- | --- | --- |
| Fixed anatomy and cohort characteristics, Fig. 1, Table 1, S1–S2 | `prepare/01–06` | `02_derivative_sample_flow.csv`, `03_qc_sensitivity.csv` |
| Genetic construction, S3–S4, Fig. S2c–d | `genetics/01–06` | `06_all_model_parameters.csv`, `03_ADNI_PGS_assembly_flow.csv` |
| Genetic trajectories, Fig. 2, S5–S8, Figs. S1, S2a–b, S3 | `analysis/01–02`, `06–10`, `22` | `02_paired_34_region_estimates.csv`, `23_primary_two_df_joint_test.csv` |
| Event origins and spatial pattern, Fig. 3, S9–S10, Fig. S4a | `analysis/03–04`, `11–12`, `19` | `01_stage_MRI_event_interaction.csv`, origin-specific Cox outputs |
| Clinical and anatomical benchmarks, Fig. 4, Table 2, S12–S17, Figs. S6–S7 | `analysis/05`, `13–14`, `16–19`, `21` | `14_components_contrasts.csv`, `27_hippocampus_contrasts.csv` |
| Pathology and baseline cognition, Fig. 5a–d, S18–S19 | `ancillary/01`, `analysis/02` | `01_repaired_pathology_structure_models.csv`, `02_repaired_structure_cognition_models.csv` |
| Current state and prior loss, Fig. 5e–f, S11, Fig. S5 | `analysis/08`, `15`, `20` | `01_event_cox_tiers.csv`, `03_paired_increment_bootstrap.csv` |
| CAT12 and AIBL, S19–S21, Figs. S4b, S8 | `ancillary/02–06`, `imaging/01–06` | `05_AIBL_incident_AD_cox_coefficients.csv`, CAT12 report/ROI tables |
| Spatial tau comparison, S22, Fig. S9 | `ancillary/07–08` | `02_spatial_null_tests.csv`, `03_spin_null_distribution.csv` |
| Donor expression and selection, S23–S25, Fig. S9b–c | `genetics/07`, `ancillary/09–11` | `11_paired_donor_contrasts.csv`, `12_donor_level_models.csv` |

`definitions` contains the bilateral region mapping and fixed gene-program memberships. HABS-HD prediction models are fitted within that cohort. The paired clinical analyses retain the matched participants and folds. Primary assessment-component and hippocampal comparisons use 1,000 paired fixed-prediction resamples. The retained two-tier, whole-cortex and second-origin analyses use 400 resamples. ADNI CSF matching uses 90 days with post-sampling MRI records. The HABS-HD pathology matching window is 548 days. Regional spatial tests distinguish correlations from the difference between two correlations.

## Data access and release

Obtain individual-level [ADNI](https://adni.loni.usc.edu/), [HABS-HD](https://apps.unthsc.edu/itr/) and [AIBL](https://aibl.org.au/) data through their access procedures. Controlled-access data and local configuration are excluded from this repository. The required input files are described in `docs/INPUTS.md`.

Verification covers R parsing, runner groups, configurable paths and optional-input control flow. End-to-end numerical reproduction with controlled participant data has not been executed for this repository revision. Repository-authored code uses the MIT License. External software and cohort data retain their own terms.
