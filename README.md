# Peripheral blood mitochondrial complex I transcriptome continuum and year-5 PIGD progression in Parkinson's disease

This repository contains the R/Python analysis code used for the manuscript on a peripheral-blood mitochondrial complex I (CI) transcriptomic axis and late PIGD progression in Parkinson's disease.

## Structure

All analysis code is consolidated in a single file, **`main_analysis.R`**, organized into 31 sections in pipeline order. Each section corresponds to one originally standalone script and is independently runnable — later sections read intermediate tables produced by earlier sections. Reviewer-driven analyses added during the second revision round are kept as standalone scripts in **`revision2_final/`** (see below).

| Section | Content |
|---|---|
| 01 | Package dependencies |
| 02 | Data preparation and CI-group assignment |
| 03 | Distribution structure assessment (dip test, GMM, consensus clustering) |
| 04 | PCA on 66 CI genes and eigen-correlations |
| 05 | DESeq2 covariate-adjusted differential expression (`~ RIN + Plate + Age + Sex + PC1_group`) |
| 06 | KEGG-GSEA and GO enrichment on covariate-adjusted DE statistics |
| 07 | PD risk-gene screen (Nalls 2019) |
| 08 | ssGSEA whole-transcriptome validation |
| 09 | k-means partition immune robustness check |
| 10 | GENEPARK independent PC1 re-estimation and CI stratification (GSE99039) |
| 11 | GENEPARK risk-gene DE and immune deconvolution validation |
| 12 | Unified PIGD LMM with deduplicated longitudinal data (primary specification) |
| 13–16 | Pre-specified longitudinal outcome sweep (43 outcomes) and forest plots |
| 17 | DAT-SBR subregion sweep |
| 18–19 | Fluid biomarker sweep and volcano figure |
| 20 | Cell-composition sensitivity analysis |
| 21 | PIGD latent-class trajectory modelling (LCMM, supplementary) |
| 22 | Figure 1: CI axis construction and structure validation |
| 23 | Figure 2: DE, GSEA, GO, risk-gene cross-cohort validation, immune shifts |
| 24 | Figure 3: year-5 PIGD divergence and biomarker context |
| 25 | Classic biomarker figure |
| 26 | Supplementary Figure S1: structure assessment |
| 27 | Table 1: baseline clinical characteristics by CI group |
| 28 | Longitudinal sensitivity analyses on the primary dataset (covariate blocks M1–M7, RIN x visit, time-varying LEDD, group-by-LEDD, medication-state-restricted OFF/ON, PDSTATE balance, LEDD trajectories) |
| 29 | 43-outcome sweep rerun: deduplicated data + REML (supersedes Sections 13–14 as the verifiable sweep) |
| 30 | Baseline UPSIT availability and group comparison |
| 31 | Score and partition consistency metrics (kappas, consensus k=3 cluster sizes) |

Sections 05, 06, 10, 11, 12 and 23 are the final revision analyses: covariate-adjusted DESeq2 (RIN, sequencing plate, age, sex), adjusted GSEA/GO, GENEPARK independent validation, the deduplicated longitudinal LMM, and Figure 2 panel (d) showing the three cell-composition-robust PD risk genes (ITGA8, RIMS1, DNAH17).

## Second revision round (`revision2_final/`)

Standalone scripts from the second revision round (August–September 2026). Each is self-contained; input tables are the derived analysis tables referenced above.

| Script | Content |
|---|---|
| `null66_and_pdhc_test.R` | Reviewer items #3 and #5: gene-set null test (1,000 random 66-gene sets through the identical score + LMM pipeline; empirical p for the observed year-5 interaction) and PD vs HC CI-score comparison |
| `pdhc_score_only.R` | PD vs HC CI-score comparison, NA-safe standalone version |
| `subscore_pigd_rin_test.R` | Nuclear-59 vs mtDNA-7 subscores: which subscore drives the year-5 PIGD interaction, and is the mtDNA signal a RIN-degradation artifact |
| `gsea_excl66.R` | KEGG GSEA rerun excluding the 66 CI scoring genes (removes circularity) |
| `hc_pigd_specificity_test.R` | HC specificity test: PD-derived PC1 loadings projected onto healthy controls; baseline PIGD in high/low CI-score HC (Age + Sex covariates) |
| `hematology_validation.R` | Validates CIBERSORT deconvolution estimates against measured PPMI CBC differentials (neutrophils/lymphocytes/monocytes %) |
| `run_ssgsea.R` | ssGSEA enrichment of the 66-gene CI set on the full transcriptome (393 PD samples) |
| `make_figS_ssGSEA.R` | Supplementary figure: ssGSEA vs PC1 concordance |
| `make_figS3_null_test.py` | Supplementary Fig. S3: gene-set null-test distribution (random 66-gene sets) with the observed beta marked |
| `make_figure2B_adjusted_bubble.R` | Replacement panel 2B: top-8 low-CI-direction KEGG pathways by adjusted P (covariate-adjusted GSEA), bubble style |
| `add_supp_tables.py` | Builds/updates the supplementary-table workbook (Nature-style formatting) from the pipeline result tables |
| `verify_table1.py` | Row-by-row independent verification of Table 1 against raw PPMI source files |

## Important Note on Input Files

This release is a manuscript-code package, not a fully containerized one-click pipeline from raw PPMI downloads. Scripts use derived analysis tables created during local preprocessing, for example cleaned sample annotations, PC1 scores, DESeq2 outputs, immune deconvolution summaries, longitudinal LMM summaries, DAT-SPECT summaries, and focused biomarker tables:

- Some sections read original or near-original PPMI files.
- Some sections read intermediate tables generated by earlier sections.
- Some figure sections read curated result tables from the manuscript tables folder.

This structure is intentional because PPMI individual-level data are controlled access and cannot be redistributed. Before rerunning the pipeline, users should either reproduce the same intermediate tables from their authorized PPMI download or update paths and object names to match their own preprocessing workflow.

## Data Availability

PPMI individual-level clinical, imaging, biomarker, and RNA-seq data are controlled-access data (www.ppmi-info.org) and are not redistributed here. GENEPARK expression data are publicly available from GEO (accession GSE99039).

The scripts assume the local file organization used during analysis. Before running, update hard-coded paths such as:

- `E:/PPMI帕金森数据库专用/...` — PPMI downloads and analysis folders
- Manuscript table and figure output directories

## Environment

R 4.5.1. Package dependencies are listed and installed in Section 01 of `main_analysis.R`.
