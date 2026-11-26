# Input contracts

## Cohort derivatives

`AD_PATHS$cohorts` holds the original cohort directory structure. The scripts read the filenames below and select columns explicitly. Columns are checked in the corresponding reader before models are fitted. Release substitutions can change linkage and eligibility.

| Study | Files used in the executed workflow | Linkage |
| --- | --- | --- |
| ADNI | `ADNI/data/Clinical/ADNIMERGE_03Jan2026.csv`; `MRI/MRI_FreeSurfer.csv` and its dictionary; `MRI/MRI_Hippocampus_ASHS.csv`; `MRI/UCD_WMH_03Jan2026.csv`; Berkeley amyloid and tau PET tables; CSF assay tables enumerated in assessment/CSF scripts | RID/PTID, examination date and IMAGEUID for MRI; nearest-date rules and tie resolution are explicit in code |
| HABS-HD | BraskieLab RP8 cortical thickness and ATN meta-ROI workbooks dated 20260501; HippoDeep hippocampal and whole-brain volume workbooks; `RP_HD_7_Genomics.xlsx`; clinical CSV release dated 12Aug2026 | Normalized shared identifier, visit and acquisition/assessment dates; endpoint-specific nearest-date matches |
| AIBL | `aibl_pdxconv_01-Jun-2018.csv`, `aibl_mrimeta_01-Jun-2018.csv` under `aibl_data`; `AIBL_Integrated_Data_FIXED.csv` at `aibl_clinical` | Normalized subject identifier, visit and MRI date; baseline exclusion and subsequent diagnosis are reconstructed from visit-level records |

The AIBL integrated clinical file is an author-prepared input used in the original execution. It supplies demographic, APOE and cognitive fields. Its required column aliases are listed in `ancillary/06_aibl_incident_ad.R` and the CAT12 parser. The file is required at the prepared-input boundary. The visit-level diagnosis table remains the event source.

The scoring input is a previously QC-filtered PLINK bed/bim/fam prefix configured by `genotype_prefix` plus the principal-component table at `genetic_pcs`. These are upstream controlled-access preprocessing inputs. The included scoring code performs clumping and allele-aware scoring; it starts after genotype QC and principal-component generation.

## GWAS and reference resources

Place the five original GWAS files in `AD_PATHS$gwas_input`. Their exact filenames are specified in `genetics/01_gwas_qc_ldsc.R`. The executed inputs were cardiometabolic multimorbidity GCST90274723, insomnia ukb-b-3957, adult hearing difficulty GCST90012115, grip strength GCST90007526 and Kunkle stage-1 AD GCST007511. GenomicSEM preserves the signed four-trait factor orientation. The main polygenic comparison excludes the APOE region using the coordinates in code.

`AD_PATHS$ldsc` supplies `w_hm3.noMHC.snplist`, the SNP/A1/A2 HapMap3 list `w_hm3.snplist`, and the `eur_w_ld_chr` directory. Each chromosome needs `.l2.ldscore.gz` and `.l2.M_5_50` files. The factor-GWAS reference also uses the original European allele-frequency files read by `genetics/03_multisystem_factor_gwas.R`.

`AD_PATHS$ldsc_code` points to an original bulik/ldsc checkout containing `munge_sumstats.py`, `ldsc.py` and the `ldscore` module. `AD_PATHS$ldsc_python` selects the Python interpreter for its environment. Use an ASCII staging path for `AD_PATHS$gwas_stage` with the original Python 2 LDSC environment. Munging uses signed Z statistics, the non-MHC HapMap3 merge list and MAF 0.01. Binary traits retain the original observed-scale specification. Formatting writes all five traits under `results/09D_GWAS_QC_LDSC/01_formatted`; the factor-GWAS and correlation scripts read that same location.

`AD_PATHS$magma` identifies the MAGMA 1.10 executable. `AD_PATHS$magma_reference` supplies the `g1000_eur/g1000_eur` bed/bim/fam prefix and GRCh37 `ENSGv110.coding.genes.txt`. Gene annotation uses 0-kb windows. Gene analysis uses the same shared SNP intersection for gAD and gMAV, the SNP/P columns and per-variant N. Its `.genes.out` files define the top-200 priority sets. Gene-set enrichment is outside the final donor comparison and is omitted.

Set the PLINK executable with `AD_PATHS$plink` and the dcm2niix executable with `AD_PATHS$dcm2niix_executable`. Their default filenames follow the operating system. The compiled CAT12 12.9 launcher in this repository uses the Windows MATLAB Runtime and the executable at `AD_PATHS$cat12_executable`. Statistical analysis of completed CAT12 tables does not use that launcher.

## Image processing and CAT12

The optional imaging scripts consume IDA exports and locally downloaded DICOM archives. They write image-identity, geometry and processing manifests. `imaging/06_cat12_processing.R` accepts `CAT12_R_MODE` as audit, prepare, run or parse. Configure `CAT12_R_RESULT_ROOT`, `CAT12_R_STAGE_ROOT` and `CAT12_R_MANIFEST` when staging the ADNI CAT12 bridge. The final ADNI analysis expects the report, ROI and job tables under `results/13F_ADNI_CAT12_bridge/tables`, plus `02_ADNI_CAT12_bridge_baseline_manifest.csv` in `AD_PATHS$cat12_config`. AIBL uses `results/13F_CAT12_R_direct_minimal/tables/04_AIBL_CAT12_endpoint.csv`.

The optional raw-image QC input is `results/13F_automated_raw_T1_QC_and_processing_manifest/04_automatic_QC_flagged_images.csv`, with columns `cohort` and `RID`. It is generated by `imaging/05_t1_qc.R`. Previously generated flags can be reused without image reprocessing. If it is absent, analyses 01, 02 and 06 skip the raw-QC exclusion models and retain their other models.

`AD_PATHS$adni_unavailable_ptids` points to a private CSV with a `PTID` column for the original source-unavailable exclusions. The CSV remains outside the repository. If absent, the raw-image matching script retains all 660 participants and 1,320 target visits. With the original exclusion list supplied, the script checks the frozen 657 participants and 1,314 visits. FreeSurfer-derived primary analyses do not use this exclusion list.

For completed CAT12 input, supply `02_ADNI_CAT12_bridge_baseline_manifest.csv` under `cat12_config` together with the report, ROI and job tables. The `cat12` group does not inspect raw-image existence. To regenerate the manifest from source images, the `imaging` group runs `ancillary/02–04` after image QC and copies the selected bridge manifest to `cat12_config`.

## Donor expression

The donor scripts read GSE297004 metadata and neuronal count RDS files at `AD_PATHS$mtg_metadata` and `AD_PATHS$mtg_counts`. Metadata requires barcode, individual, diagnosis and anno_lvl_2. Counts must be a named gene-by-nucleus matrix with column names matching barcode. The code checks diagnosis mapping and donor consistency before calculating donor-level contrasts. The prepared metadata and count matrices are required inputs.

`definitions/mtg_gene_programs.csv` preserves the 189 gAD and 191 gMAV gene memberships used in the final donor analysis. They derive from the top-200 shared-SNP MAGMA priorities after the previously fixed expression-feature coverage intersection. MAGMA input preparation and shared-SNP gene prioritization are retained in code.

## Input locations

All keys below belong to `AD_PATHS`. Set them in `config.R` or override them in a private `config.local.R` file. Existing data directories can be used directly without renaming files.

| Key | Input |
| --- | --- |
| `project` | Working project containing generated `results` and reference inputs |
| `cohorts` / `adni_data` | Cohort root / ADNI data release |
| `habs_release` | HABS-HD clinical and imaging-derivative release |
| `gwas_input` | The five original GWAS summary files |
| `genotype_prefix` / `genetic_pcs` | QC-filtered PLINK prefix / participant principal components |
| `aibl_data` / `aibl_clinical` | AIBL release tables / prepared integrated clinical table |
| `aibl_mri_manifest` | Original AIBL baseline T1 acquisition manifest |
| `mtg_metadata` / `mtg_counts` | GSE297004 metadata / neuronal counts |
| `cat12_config` | Completed ADNI CAT12 bridge manifest |
| `adni_unavailable_ptids` | Optional private source-unavailable exclusion CSV |
| `plink` / `dcm2niix_executable` / `cat12_executable` | Original executable locations |
| `ldsc` / `ldsc_code` / `ldsc_python` | LD references / original LDSC checkout / its Python interpreter |
| `magma` / `magma_reference` | MAGMA executable / European reference files |

