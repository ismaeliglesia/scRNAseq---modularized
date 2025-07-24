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
```

## 📂 Files

- **modules.R**  
  Contains reusable functions for each analysis step:  
  - `load_RNA_data()`  
  - `qc_filtering()`  
  - `normalize_SCT()`  
  - `create_UMAP()`  
  - `findclusters_resolutions()`  
  - `make_heatmap()`  

- **main.R**  
  Script to run the full analysis pipeline using the functions from `modules.R`. It sets working directories, loads data, runs QC, normalization, dimensionality reduction, and clustering.
  
- **utils.R**
  Important utils and code for making the pipeline work.

## 📝 How to use

1. Clone this repository or download the scripts
2. Prepare your project directory with the following structure:

```css
   project_folder/
├── count/
│   ├── filtered_feature_bc_matrix.h5
│   ├── aggregation.csv
├── metadata.xlsx
├── vdj/
│   ├── filtered_contig_annotations.csv
├── modules.R
├── main.R
├── utils.R
```
3. Edit main.R to set the directories, project name and function flags:
- R_directory: Absolute pathway of the folder in which the data is contained.
- project_name: Name of the seurat_obj and the saved object (in RDS format).
- Function flags: Look for all the options of each function and specify if different than default.

## 💫 Output
The output is a .RDS file with the normalized and clustered single cell experiment.
