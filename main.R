# --- Project configuration ---
R_directory <- "/path/to/your/project"  # Adjust the path
setwd(R_directory)

project_name <- "scexp_2025"  # Used only for file naming

source("utils.R")
source("modules.R")

# Load required libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)
library(scRepertoire)
library(cowplot)

# --- 1. Load data ---
tryCatch({
  seurat_obj <- load_RNA_data(
    h5_path = "./count/filtered_feature_bc_matrix.h5",
    filter_singlets = TRUE,
    aggr_dir = "./count/aggregation.csv",
    metadata_file = "./metadata.xlsx",
    contig_list = "./vdj/filtered_contig_annotations.csv"
  )
}, error = function(e) {
  stop("Error loading RNA data: ", e$message)
})

# --- 2. Quality control filtering ---
tryCatch({
  seurat_obj <- qc_filtering(
    seurat_obj = seurat_obj,
    mito_percent_qc = 10,
    min_genes_qc = 500,
    max_genes_qc = 8000
  )
}, error = function(e) {
  stop("Error in QC filtering: ", e$message)
})

# --- 3. Normalization and/or cell cycle scoring/regression ---
tryCatch({
  seurat_obj <- normalize_SCT(
    seurat_obj,
    genes_to_mask = NULL,
    cellcyclescore = TRUE,
    cellcyclenorm = TRUE,
    new_assay_name = "ccSCT"
  )
}, error = function(e) {
  stop("Error in normalization: ", e$message)
})

# --- 4. Dimensionality reduction (PCA + UMAP) ---
tryCatch({
  seurat_obj <- create_UMAP(
    seurat_obj,
    automatic_npc = TRUE,
    npcs = NULL,
    min.dist = 0.01,
    k.param = 30,  # more typical value
    new_assay_name = "ccSCT"
  )
}, error = function(e) {
  stop("Error creating UMAP: ", e$message)
})

# --- 5. Clustering ---
tryCatch({
  seurat_obj <- findclusters_resolutions(
    seurat_obj,
    resolutions = c(1.8),
    n.iter = 10,
    save_plot = TRUE,
    new_assay_name = "ccSCT"
  )
}, error = function(e) {
  stop("Error in clustering: ", e$message)
})

# --- Save final object ---
saveRDS(seurat_obj, file = paste0(project_name, "_final.rds"))
