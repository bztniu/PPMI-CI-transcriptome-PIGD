# Code Input and Output Map

This file explains where each numbered script obtains its inputs. It is intended for manual checking before uploading the code to GitHub.

## Path Conventions Used by the Scripts

The original scripts were run in a local Windows workspace. The main root paths were:

- `E:/PPMI帕金森数据库专用`
  - local PPMI data root and derived expression objects
- `E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2`
  - analysis workspace containing derived tables used by several downstream scripts
- `E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21`
  - current manuscript package, with final tables in `04_tables/` and figures in `02_figures/`

For GitHub release, these absolute paths should be replaced by configurable project-root variables before public reuse.

## 00_install_pkgs.R

Purpose: installs or loads R package dependencies.

Inputs:

- CRAN/Bioconductor package repositories.

Outputs:

- Installed R packages in the user's R library.

## 01_prep_data.R

Purpose: core preprocessing script for CI PC1 score, DESeq2 differential expression, VST panel expression, NNLS immune deconvolution, selected LMM summaries, and an older ROC helper table.

Important: this script is not purely raw-data based. It mixes controlled-access PPMI-derived objects, manually curated gene tables, PPMI clinical CSVs, and local derived sample annotations.

Inputs:

| Input path in script | Type | Source category | What it is used for |
|---|---|---|---|
| `complex/CI_genes_converted.csv` | CSV | Curated derived file | 66 MitoCarta CI genes mapped to expression IDs and symbols |
| `dds_pd_control_BL_object.rds` | RDS | Derived expression object from PPMI RNA-seq | Baseline DESeq2 object containing PD and control whole-blood RNA-seq counts |
| `revision_v2/data/PD_all_clustering_methods.csv` | CSV | Derived sample annotation | PD sample IDs, PATNO mapping, prior clustering/group labels, baseline clinical summaries |
| `结果/免疫浸润/refer.txt` | TSV | Reference signature matrix | LM22 immune-cell signature matrix used for NNLS deconvolution |
| `运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv` | CSV | PPMI clinical export | Longitudinal MoCA total score for LMM helper outputs |
| `运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv` | CSV | PPMI clinical export | MDS-UPDRS III, PIGD items, tremor items for LMM helper outputs |
| `Demographics_03Feb2026.csv` | CSV | PPMI clinical export | Sex covariate for older ROC helper output |
| `Age_at_visit_17Mar2025.csv` | CSV | PPMI clinical export | Baseline age covariate for older ROC helper output |

Main outputs:

| Output path | What it contains |
|---|---|
| `figures_pc1/data/pca.rds` | PCA object for 66 CI genes |
| `figures_pc1/data/pc1_scores.csv` | SAMPLE_ID, PATNO, PC1 score, low/high CI group |
| `figures_pc1/data/pc1_loadings.csv` | PC1 loading per CI gene |
| `figures_pc1/data/pca_variance.csv` | Variance explained by PCs |
| `figures_pc1/data/ci_expr.csv` | log2 normalized CI expression matrix |
| `figures_pc1/data/de_full.csv` | genome-wide DESeq2 results for high-CI vs low-CI expression groups |
| `figures_pc1/data/panel_vst.csv` | VST expression of selected PD/CI/immune/synaptic genes |
| `figures_pc1/data/heatmap_anno.csv` | Heatmap annotation table |
| `figures_pc1/data/immune_fractions.csv` | NNLS-estimated immune-cell fractions |
| `figures_pc1/data/immune_cohend.csv` | Immune-cell Cohen's d and p/q values |
| `figures_pc1/data/lmm_emm.csv` | selected estimated marginal means from older LMM helper code |
| `figures_pc1/data/lmm_interactions.csv` | selected LMM interaction estimates from older helper code |
| `figures_pc1/data/raw_means.csv` | raw mean trajectories |
| `figures_pc1/data/roc_pred.csv` | older ROC helper output; not part of the current main manuscript story |

## 02_structure_assessment.R

Inputs:

- `figures_pc1/data/ci_expr.csv` from `01_prep_data.R`
- `figures_pc1/data/pc1_scores.csv` from `01_prep_data.R`

Outputs:

- `figures_pc1/data/structure_assessment.rds`
- structure assessment summary tables/plots if enabled in the script

## 03_pca_eigencor.R

Inputs:

- `figures_pc1/data/ci_expr.csv`
- `figures_pc1/data/pc1_scores.csv`
- `Demographics_03Feb2026.csv`
- `Age_at_visit_17Mar2025.csv`
- `运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv`
- `运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv`

Source category:

- first two are derived by `01_prep_data.R`
- remaining files are PPMI clinical exports

## 04_run_gsea_lowCI.R

Inputs:

- DESeq2 result table generated upstream, typically `figures_pc1/data/de_full.csv` or the matching final-table copy searched by the script.
- MSigDB/clusterProfiler gene-set resources loaded by R packages.

Outputs:

- Low-CI-ranked GSEA result tables used for Figure 2.

## 05_screen_pd_risk_genes.R

Inputs:

- `figures_pc1/data/de_full.csv`: PPMI DESeq2 results from `01_prep_data.R`
- `ppmi数据表/GSE99039/results/GSE99039_limma_DE_PC1.csv`: derived GENEPARK limma result table

Outputs:

- PD-risk gene summary tables used in Figure 2 and supplementary tables.

## 06_run_ssgsea.R

Inputs:

- `dds_pd_control_BL_object.rds`: derived PPMI RNA-seq DESeq2 object
- `revision_v2/data/PD_all_clustering_methods.csv`: derived PD sample annotation
- `complex/CI_genes_converted.csv`: curated CI gene mapping table
- `figures_pc1/data/three_scores.csv`: derived score table containing PC1 and alternative scores

Outputs:

- ssGSEA/equal-weight score sensitivity outputs.

## 07_kmeans_immune_check.R

Inputs:

- `figures_pc1/data/immune_fractions.csv`: NNLS output from `01_prep_data.R`
- `revision_v2/data/PD_all_clustering_methods.csv`: derived clustering/sample annotation

Outputs:

- k-means immune robustness tables used for Figure 2.

## 08_lmm_sweep.R

Inputs:

- `figures_pc1/data/pc1_scores.csv`: PC1 group and PATNO mapping
- PPMI clinical exports:
  - `运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv`
  - `REM_Sleep_Behavior_Disorder_Questionnaire_03Feb2026.csv`
  - `SCOPA-AUT_07Feb2026.csv`
  - other files passed through the local `read_csv(file)` helper inside the script

Outputs:

- LMM sweep outputs for selected clinical outcomes.

## 09_lmm_sweep_full.R

Inputs:

- `data/pc1_scores.csv`: derived PC1 score table
- Multiple PPMI clinical/cognitive/autonomic CSV exports loaded through the script's `rd()` helper

Outputs:

- `LMM_sweep_COMPREHENSIVE.csv` or equivalent comprehensive clinical LMM scan table.

## 10_lmm_sweep_forest.R

Inputs:

- `LMM_sweep_all_outcomes.csv`: output from `08_lmm_sweep.R`

Outputs:

- Forest-plot figure/table for selected LMM outcomes.

## 11_lmm_sweep_full_forest.R

Inputs:

- `LMM_sweep_COMPREHENSIVE.csv`: output from `09_lmm_sweep_full.R`

Outputs:

- Forest-plot figure/table for comprehensive LMM outcomes.

## 12_dat_subregion_sweep.R

Inputs:

- `data/pc1_scores.csv`: derived PC1 group table
- `运动症状数据/Xing_Core_Lab_-_Quant_SBR_23Feb2026.csv`: PPMI DAT-SPECT quantitative SBR export

Outputs:

- DAT-SPECT subregion LMM summary tables.
- Posterior dorsal putamen trajectory table used in Figure 3.

## 13_biomarker_sweep.R

Inputs:

- `data/pc1_scores.csv`: derived PC1 group table
- `data/biomarkers_long.csv`: locally curated long-format biomarker table derived from PPMI blood/CSF biomarker exports

Outputs:

- `biomarker_LMM_sweep.csv`: focused biomarker comparison table used by Figure 3.

## 14_biomarker_volcano.R

Inputs:

- `biomarker_LMM_sweep.csv`: output from `13_biomarker_sweep.R`

Outputs:

- Biomarker summary plot.

## 15_cell_fraction_sensitivity.R

Purpose: blood-cell-composition sensitivity analysis for the year-5 PIGD interaction (Supplementary Table S15).

Inputs:

- `figures_pc1/data/pc1_scores.csv`: PC1 scores and group labels
- `figures_pc1/data/immune_fractions.csv`: estimated immune-cell fractions
- PPMI MDS-UPDRS Part III clinical export (PIGD items)

Outputs:

- Supplementary Table S15 (PIGD year-5 interaction before and after adjustment for estimated leukocyte fractions).

## 16_pigd_lcmm_trajectory.R

Purpose: latent-class mixed modelling (LCMM) of PIGD trajectories; generates the supplementary LCMM trajectory figure.

Inputs:

- `04_tables/three_scores.csv`: PC1 scores and alternative partitions
- PPMI MDS-UPDRS Part III clinical export (PIGD items)

Outputs:

- `04_tables/PIGD_lcmm_long_input.csv`, `PIGD_lcmm_model_selection.csv`, `PIGD_lcmm_class_trajectories.csv`, `PIGD_lcmm_patient_labels.csv`, and CI-enrichment test tables
- `02_figures/FigS_PIGD_lcmm_trajectories.png/pdf`

Requires the `lcmm` package; run from the repository root so `04_tables/` and `02_figures/` resolve correctly.

## 17_make_figure1.R

Inputs:

- `figures_pc1/data/pca.rds`
- `figures_pc1/data/pc1_scores.csv`
- `figures_pc1/data/pca_variance.csv`

Outputs:

- Main Figure 1.

## 18_make_figure2.R

Inputs:

- Final/derived tables under `04_tables/`, including DESeq2, GSEA, GO, PD-risk validation, immune, k-means, and cell-adjustment tables.
- Some paths are searched by helper functions in the script.

Outputs:

- Main Figure 2 (9-panel layout, panels a-i).

## 19_make_figure3.R

Inputs:

- `04_tables/pc1_scores.csv` or equivalent `DATA_DIR/pc1_scores.csv`
- `04_tables/biomarkers_long.csv`
- `04_tables/biomarker_LMM_sweep.csv`
- `04_tables/PIGD_continuous_PC1_visit_sensitivity.csv`
- `04_tables/Table_LMM_factor.csv`
- `04_tables/LMM_sweep_COMPREHENSIVE.csv`
- `04_tables/DAT_posterior_dorsal_putamen_trajectory.csv`
- selected SAA raw table if available
- selected derived `lmm_emm.csv` and `lmm_interactions.csv`

Outputs:

- Main Figure 3.

## 20_make_classic_biomarker_figure.R

Inputs:

- `data/pc1_scores.csv`
- `data/biomarkers_long.csv`
- SAA raw file path specified in the script, if available
- `04_tables/biomarker_LMM_sweep.csv`

Outputs:

- Focused established-fluid-biomarker panel used as Figure 3g or helper figure.

## 21_make_figS1_structure.R

Inputs:

- `figures_pc1/data/structure_assessment.rds`
- `figures_pc1/data/pc1_scores.csv`

Outputs:

- Supplementary CI-axis structure figure.
