# Quality Control and Analysis Pipeline for scRNA-Seq Data

## Overview
This repository contains scripts and results for quality control (QC) and analysis of single-cell RNA sequencing (scRNA-seq) data used in our study of advanced-stage non-small cell lung cancer (NSCLC). Our analysis focuses on understanding tumor evolution under targeted therapy, identifying drug resistance mechanisms, and exploring tumor microenvironment interactions.

## Study Design
Our study follows the scRNA-seq analysis pipeline applied to NSCLC samples at different treatment stages:

- **TKI-naïve (TN)**: Pre-treatment baseline
- **Residual Disease (RD)**: Tumor regression or stability under therapy
- **Progressive Disease (PD)**: Tumor develops resistance and resumes growth

The study aims to identify transcriptional changes, immune interactions, and resistance pathways associated with tumor progression.

## Data Processing Pipeline
### 1. **Quality Control (QC)**
- Performed using `fastqc` and `fastp` to assess sequencing quality and remove low-quality reads.
- `multiqc` was used to aggregate and visualize QC metrics.
- `featureCounts` was applied for read quantification.

### 2. **Data Preprocessing**
- Low-quality cells filtered based on gene count and mitochondrial content.
- Normalization performed using log-normalization or `SCTransform`.
- Batch correction applied using `Harmony` or `Seurat` CCA.

### 3. **Clustering and Annotation**
- PCA used for dimensionality reduction.
- Clustering performed with Leiden/Louvain algorithms.
- Cell type annotation via marker genes, `SingleR`, and manual curation.

### 4. **Gene Expression and Pathway Analysis**
- Differential expression analysis (`DESeq2`, Wilcoxon test) between TN, RD, and PD states.
- Gene Set Enrichment Analysis (`GSEA`) for pathway discovery.
- Focus on resistance pathways: **WNT/β-catenin, kynurenine metabolism, plasminogen activation**.

### 5. **Cell-Cell Interaction and Resistance Mechanisms**
- Ligand-receptor interaction analysis (`CellPhoneDB`, `NicheNet`).
- Trajectory inference (`Monocle`, `RNA velocity`) to study cell-state transitions.

## Data Source
The data used in this study is based on **BioProject PRJNA622993**.

## Platform and Software
The analysis was conducted using **Case Western Reserve University High-Performance Computing (HPC) Platform** and **Cleveland Clinic HPC**. The following software versions were used:

- **Python**: 3.12.3
- **SRA Toolkit**: v3.2.0
- **Entrez Direct**: v23.5
- **FastQC**: v0.12.1
- **fastp**: v0.24.0
- **MultiQC**: v1.27.1
- **featureCounts**: v2.0.8

## Contributors
- **Isaac Zeng** (Department of Nutrition, Case Western Reserve University, yxz3103@case.edu, zengy2@ccf.org)
- **Max Ni** (Department of Biology, Case Western Reserve University, xxn49@case.edu)

## References
> Maynard, A., McCoach, C. E., Rotow, J. et al. (2020). *Therapy-Induced Evolution of Human Lung Cancer Revealed by Single-Cell RNA Sequencing*. Cell. DOI: [10.1016/j.cell.2020.07.013](https://doi.org/10.1016/j.cell.2020.07.013)


