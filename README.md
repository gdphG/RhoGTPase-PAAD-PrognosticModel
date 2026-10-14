# R Code for Building the 22-Gene Rho GTPase Prognostic Model in Pancreatic Cancer

This repository contains the complete R code for constructing and validating the 22-gene Rho GTPase-related prognostic model described in our manuscript submitted to *British Journal of Pharmacology*.

## Overview

The 22-gene signature comprises: **ARAP1, ARHGAP18, ARHGAP23, ARHGAP29, ARHGAP44, ARHGEF11, ARHGEF26, CALM1, CENPT, CTNNA1, CTTN, DEPDC1B, DVL1, GDI1, INPP5B, IQGAP1, NCKAP1, NUP37, PAK2, RHOBTB2, RHOV, SPDL1**

Three genes were identified as independent prognostic factors: **ARHGAP18, CENPT, NUP37**.

## Datasets

| Dataset | Role | Description |
|---------|------|-------------|
| TCGA-PAAD | Training (n=178) | RNA-seq + clinical data |
| GSE28735 | Validation (n=45) | Microarray |
| GSE62452 | Validation (n=63) | Microarray |
| ICGC-PACA-CA | Validation | Whole-genome sequencing |
| GSE57495 | Biopsy cohort (n=78) | Tissue biopsy |
| GSE78220 | Immunotherapy (n=89) | Anti-PD-1 cohort |
| IMvigor210 | Immunotherapy (n=348) | Anti-PD-L1 cohort |

## Code Structure

Scripts are organized in the recommended execution order:

### Model Construction (Steps 1–8)

| Script | Description |
|--------|-------------|
| `R/step1_differential_expression.R` | Differential gene screening (limma + Wilcoxon test, FDR < 0.05, \|logFC\| > 1) |
| `R/step2_univariate_cox.R` | Univariate Cox regression with BH correction |
| `R/step3_venn_diagram.R` | Venn diagram for gene intersection |
| `R/step4_heatmap.R` | Heatmap visualization (pheatmap) |
| `R/step5_lasso_cox.R` | **Core model**: LASSO-Cox + Bootstrap (200×) + Permutation test (500×) |
| `R/step6_multivariate_cox.R` | Multivariate Cox regression + risk score + forest plot |
| `R/step7_nomogram.R` | Nomogram construction + calibration curves (B=1000) |
| `R/step8_survival_ROC.R` | KM survival curves + time-dependent ROC (1/2/3/5-year) |

### Subgroup & Clinical Analyses (Steps 9–10)

| Script | Description |
|--------|-------------|
| `R/step9_subgroup_cox_forest.R` | Subgroup analysis + clinical variable forest plots + clinical ROC |
| `R/step10_clinical_ROC.R` | Clinical feature ROC comparison (risk score vs. individual features) |

### Functional Enrichment (Step 11)

| Script | Description |
|--------|-------------|
| `R/step11_GSEA.R` | GSEA analysis (KEGG + GO, High vs. Low risk groups) |

### Immune Microenvironment (Steps 12–14)

| Script | Description |
|--------|-------------|
| `R/step12_immune_checkpoints.R` | Immune checkpoint gene expression + Wilcoxon test |
| `R/step13_immune_correlation.R` | CIBERSORT cell composition vs. risk score correlation |
| `R/step14_cibersort_stacked.R` | Stacked bar plots of immune cell proportions |
| `R/Cibersort.R` | CIBERSORT algorithm (supporting function) |

### Genomic & Drug Analyses (Steps 15–16)

| Script | Description |
|--------|-------------|
| `R/step15_mutation_analysis.R` | Somatic mutation landscape (maftools, High vs. Low TMB) |
| `R/step16_drug_sensitivity.R` | Drug sensitivity (pRRophetic + GDSC, FDR correction) |

### Immunotherapy Prediction (Steps 17–18)

| Script | Description |
|--------|-------------|
| `R/step17_immunotherapy_prediction.R` | IPS scoring + TIDE algorithm prediction |

### Single-Cell Analysis (Step 18)

| Script | Description |
|--------|-------------|
| `R/step18_single_cell.R` | Seurat CCA integration + AddModuleScore + CellChat + NUP37 validation |

### External Validation (Step 19)

| Script | Description |
|--------|-------------|
| `R/step19_external_validation.R` | Validation in GSE28735, GSE62452, ICGC-PACA-CA, GSE57495, GSE78220, IMvigor210 |

## Required Data Files

Place the following files in the `data/` directory before running:

- `mRNA1.txt` — TCGA-PAAD expression matrix (genes × samples)
- `clinical.txt` — TCGA-PAAD clinical data
- `risk_scorecd.txt` — Risk scores with group assignments (output from step 5/6)
- `timeexpcd.txt` — Expression matrix with survival time
- `LM22.txt` — CIBERSORT LM22 signature matrix
- `c2.cp.kegg.v7.4.symbols.gmt` — KEGG gene sets
- `c5.go.v7.4.symbols.gmt` — GO gene sets
- External validation expression matrices (GSE28735, GSE62452, etc.)

## R Package Dependencies

```r
# Core survival analysis
install.packages(c("survival", "survminer", "timeROC", "rms"))

# Differential expression & enrichment
BiocManager::install(c("limma", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "maftools", "TCGAmutations"))

# Visualization
install.packages(c("pheatmap", "ggplot2", "ggpubr", "ggsci", "cowplot", "forestplot", "VennDiagram", "survminer", "rms"))

# Machine learning
install.packages(c("glmnet"))

# Single-cell
BiocManager::install(c("Seurat", "SingleCellExperiment"))
install.packages("CellChat")

# Drug sensitivity
install.packages("pRRophetic")

# Immune deconvolution
# CIBERSORT (Cibersort.R included in this repository)
```

## Risk Score Formula

```
RiskScore = Σ (Coef_i × Expression_i)
```

Where Coef_i are the multivariate Cox regression coefficients for each of the 22 signature genes.

**Median risk score** from TCGA-PAAD training cohort is used as the cutoff for High/Low risk group assignment across all cohorts.

## Citation

If you use this code, please cite our manuscript:

> [Manuscript reference to be added upon publication]

## License

This code is provided for academic research purposes.
