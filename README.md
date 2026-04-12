# Neonatal CITE-seq atlas

Custom analysis and visualization scripts accompanying the study:

**"Human neonatal CITE-seq atlas identifies an immune transition at 32 weeks’ gestation from CD15⁺ myeloid-dominated to interferon-primed immunity"**

🔗 Preprint: https://www.biorxiv.org/content/10.64898/2026.04.01.715643v1

---

## Table of Contents
- [Overview](#overview)
- [Graphical Abstract](#graphical-abstract)
- [Abstract](#abstract)
- [Data and Materials Availability](#data-and-materials-availability)
- [Environment Setup](#environment-setup)
- [Usage](#usage)
- [Scripts](#scripts)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview
This repository provides **example analysis and visualization scripts** used in the neonatal CITE-seq atlas study.

The scripts illustrate key analysis steps and reproduce selected figures from the manuscript. 

---

## Graphical Abstract 
![Graphical Abstract](GraphicalAbstract.png)

---

## Abstract
The human neonatal immune system is developmentally specialized to balance the unique requirements of perinatal transition. Disruption of this finely tuned balance, as in preterm birth, may have profound consequences for immunity and overall health. However, the impact of prematurity on immune composition and functional responsiveness across gestational ages (GA) remains incompletely understood. Single-cell profiling has advanced our understanding of neonatal immunity, yet most studies were limited to unimodal readouts, narrow GA windows, or baseline function. Here, we present a comprehensive human neonatal CITE-seq atlas (82 samples from 25 neonates and 10 adults as controls) at the first days of life covering a wide GA range and integrating baseline and stimulated conditions. Most notably, we identify a GA-dependent immune transition point centering around 32 weeks of GA, which discriminates extremely and very preterm neonates (GA <32wks) from those of higher GA (≥32wks). In particular, early-life immunity in extremely and very preterm infants showed CD15+ granulocytic myeloid derived suppressor cell-like predominance, whereas more mature neonates exhibited interferon-primed transcriptional profiles. This was associated with divergent myeloid-to-lymphocyte signaling networks and qualitatively distinct NK- and T-cell bystander responses upon activation. Together, these findings show that intrauterine development imprints GA-specific immune programs. By defining a developmental transition around a GA of 32 weeks that regulates baseline and induced responses of neonatal immune cells, our atlas provides a framework for understanding the vulnerability of preterm infants and thus may pave the way for developing GA-adapted immunomodulatory strategies.

---

## Data and Materials Availability

All code and data associated with this study are present in the paper, the Supplementary Materials, or as indicated below:  

- **Interactive exploration and download will be possible via:** [CELLxGENE portal](https://cellxgene.cziscience.com/)  
- **Raw and processed CITE-seq data matrices will be available at:** [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXX)  
- **Custom scripts for analysis and visualization:** available in this repository

The repository contains representative scripts used for analysis and figure generation.
All other analyses were performed using publicly available software packages as described in the Methods section of the manuscript.

---

## Environment Setup

Analysis was performed in R. Below is a minimal setup:

```r
# CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr", "ggrepel", "forcats", "tibble",
                   "Matrix", "pheatmap", "svglite", "ggpattern", "hdf5r"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "DESeq2", "enrichR"))

# GitHub dependency
install.packages("remotes")
remotes::install_github("satijalab/MuDataSeurat")
```

---

## Usage

1. Download dataset from CELLxGENE 
2. Store locally (e.g., `data/`)
3. Adjust file paths within scripts 
4. Run scripts 

Example:
```r
source("scripts/PCAplot.R")
```

---

## Scripts

## Scripts

This repository contains example analysis and visualization scripts used in the study.  
The scripts implement general analysis workflows and can be applied to different subsets of the dataset depending on the use case.

### `PCAplot.R`
- **Input:** Preprocessed single-cell dataset (e.g., CELLxGENE dataset “All - neonatal and adult circulating immune cells baseline and stim”)  
- **Output:** PCA plot 
- **Description:**  
  Aggregates normalized RNA expression at the sample level and performs principal component analysis (PCA) across donors and conditions. Samples are annotated by gestational age, stimulation condition, and postnatal age.  
- **Example:** Reproduces Fig. 1B when applied to the full dataset  

### `ndensPlot.R`
- **Input:** Preprocessed single-cell dataset (e.g., CELLxGENE dataset “Myeloid cells - neonatal and adult circulating immune cells baseline (Fig 2)”)  
- **Output:** UMAP and normalized density plots  
- **Description:**  
  Generates UMAP visualizations and overlays normalized density estimates to highlight differences in cellular distributions across groups.  
- **Example:** Reproduces Fig. 2A when applied to the myeloid cell subset  

### `StackedBarPlot_cellproportions.R`
- **Input:** Preprocessed single-cell dataset (e.g., CELLxGENE dataset “Myeloid cells - neonatal and adult circulating immune cells baseline (Fig 2)”)  
- **Output:** Stacked bar plot  
- **Description:**  
  Computes relative cell type proportions per group and visualizes them using stacked bar plots, with optional pattern overlays for subpopulations.  
- **Example:** Reproduces Fig. 2B when applied to the myeloid cell subset  

### `Pheatmap.R`
- **Input:** Subsetted preprocessed single-cell dataset (e.g., CELLxGENE dataset “Myeloid cells - neonatal and adult circulating immune cells baseline (Fig 2)”, subset "CD15⁺ myeloid cells”)  
- **Output:** Gene expression heatmap  
- **Description:**  
  Calculates average gene expression per group and generates a scaled heatmap of selected marker genes to highlight transcriptional differences between groups.  
- **Example:** Reproduces Fig. 2C when applied to "CD15⁺ myeloid cells”

### `DEGs_Analysis_Visualization.R`
- **Input:** Subsetted preprocessed single-cell dataset (e.g., CELLxGENE dataset “All - neonatal and adult circulating immune cells baseline and stim”)  
- **Output:** Differential expression results and visualizations  
- **Description:**  
  Performs age-stratified and interaction-based differential expression analysis, followed by enrichment analysis and visualization of differentially expressed genes.

---

## Citation

If you use this code, please cite:

Rothämel P, Mattia A, Corey MI, Puzek B, Wiesel J, Michael-Kuschel P, Klein C, Sperandio M, Henneke P, Nussbaum C, Kim-Hellmuth S.  
**Human neonatal CITE-seq atlas identifies an immune transition at 32 weeks’ gestation from CD15⁺ myeloid-dominated to interferon-primed immunity.**  
bioRxiv (2026).  
https://www.biorxiv.org/content/10.64898/2026.04.01.715643v1

---

## Contact

- Name: Paula Rothämel, MD  
- Email: Paula.Rothaemel@med.uni-muenchen.de  
- Institution: Dr. von Hauner Children's Hospital, LMU Munich 
