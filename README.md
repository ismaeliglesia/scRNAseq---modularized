# 🧬 Single-Cell RNA-Seq Analysis Toolkit (Seurat + TCR)

This repository provides a modular R toolkit for end-to-end analysis of single-cell RNA-seq data using Seurat, including:

- Loading and preprocessing 10X Genomics data (Gene Expression, Antibody Capture, TCR/BCR)
- Quality control filtering
- SCTransform normalization with optional cell cycle correction
- Dimensionality reduction (PCA, UMAP)
- Clustering at multiple resolutions
- Custom heatmap generation with annotated metadata

---

## 📦 Requirements

Before using this toolkit, install the following R packages:

```R
install.packages(c("Seurat", "dplyr", "readxl", "ggplot2", "scales", "RColorBrewer", "ComplexHeatmap"))
remotes::install_github("ncborcherding/scRepertoire")
