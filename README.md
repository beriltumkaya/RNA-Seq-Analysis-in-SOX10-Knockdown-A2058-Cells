# RNA-Seq Analysis in SOX10 Knockdown A2058 Cells

## Bilkent University Department of Molecular Biology and Genetics Senior Project
### Author: Beril Tümkaya
### Date: January 2026

---

## Project Overview

To investigate differential gene expression between experimental conditions using transcriptome-wide
RNA-Seq analysis coupled with DESeq2, enabling the systematic identification of genes whose expression levels
change significantly. In addition, quantitative PCR (qPCR) was employed to validate the expression patterns of
selected differentially expressed genes, providing an independent and targeted method to confirm the RNA-Seq
findings and enhance the reliability of the results.

The repository contains a fully automated and reproducible RNA-Seq analysis pipeline implemented in R, covering quality control, normalization, statistical testing, visualization, and functional enrichment analysis.

---

## Experimental Design

The analysis includes the following experimental conditions:

* no_treatment
* siCNTL (negative control siRNA)
* siSOX10_1
* siSOX10_2

Differential expression analyses are performed for the following comparisons:

* no_treatment vs siCNTL
* siCNTL vs siSOX10_2
* siSOX10_1 vs siSOX10_2

---

## Repository Structure

```
RNA-Seq-Analysis-in-SOX10-Knockdown-A2058-Cells/
│
├── deseq2_pipeline.R          # Main RNA-Seq analysis pipeline
├── combined_counts.csv        # Raw RNA-Seq count matrix
├── metadata.csv               # Sample metadata
├── results/                   # Analysis outputs (auto-generated)
│   ├── QC plots
│   ├── normalized counts
│   ├── PCA and clustering
│   ├── differential expression results
│   ├── heatmaps and volcano plots
│   ├── KEGG enrichment results
│   └── annotated gene tables
└── README.md                  # Project documentation
```

---

## Input Files

### Count Matrix (`combined_counts.csv`)

* Rows represent genes (ENSEMBL IDs)
* Columns represent samples
* Values correspond to raw read counts
* Counts must be non-negative integers

### Metadata File (`metadata.csv`)

The metadata file must contain at least the following columns:

| Column    | Description                                 |
| --------- | ------------------------------------------- |
| sample_id | Must match column names of the count matrix |
| treatment | Experimental condition for each sample      |

Example:

```csv
sample_id,treatment
Sample1,no_treatment
Sample2,siCNTL
Sample3,siSOX10_1
Sample4,siSOX10_2
```

---

## Software Requirements

* R (version 4.2 or higher recommended)
* Bioconductor

The pipeline automatically installs and loads all required CRAN and Bioconductor packages, including DESeq2, ggplot2, clusterProfiler, and related dependencies.

---

## Analysis Pipeline Overview

The R script implements a complete RNA-Seq analysis workflow consisting of the following steps:

1. Validation of the raw count matrix
2. Raw count distribution quality control
3. Creation of the DESeq2 dataset
4. Filtering of low-count genes
5. Library size normalization and size factor estimation
6. Variance stabilizing transformation (VST)
7. Principal component analysis and scree plot
8. Sample correlation and hierarchical clustering
9. Differential expression analysis using DESeq2
10. Mean–variance relationship assessment
11. Dispersion estimation
12. Pairwise differential expression comparisons
13. UpSet plot for overlapping differentially expressed genes
14. Log2 fold-change shrinkage and MA plots
15. KEGG pathway enrichment analysis
16. Gene annotation (ENSEMBL ID to gene symbol)
17. Summary of significant genes
18. Heatmap visualization of top differentially expressed genes
19. Volcano plot visualization
20. Expression plots of top genes across conditions

---

## How to Run the Analysis

### 1. Prepare the Working Directory

Ensure the following files are located in the same directory:

* deseq2_pipeline.R
* combined_counts.csv
* metadata.csv

Set the working directory in R:

```r
setwd("path/to/project_directory")
```

---

### 2. Run the Full Pipeline

Execute the pipeline by sourcing the script:

```r
source("deseq2_pipeline.R")
```

The pipeline runs end-to-end automatically and generates all results in the `results/` directory.

---

## Output Files

All outputs are written to the `results/` directory and include:

* Quality control plots
* Normalized and VST-transformed expression matrices
* PCA and clustering visualizations
* Differential expression result tables
* Annotated gene lists
* Heatmaps, MA plots, and volcano plots
* KEGG pathway enrichment analysis results

---
