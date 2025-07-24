R_directory <- ""  ### Set the project directory
setwd(R_directory)
project_name = "" ### Set the name of the project, i.e. scexp_2025

source("utils.R")
source("modules.R")
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)
library(scRepertoire)
library(cowplot)
library(readxl)

# --- Load data ---
tryCatch({
  # Load 10X data, filter singlets, add metadata and TCR contigs if available
  seurat_obj <- load_RNA_data(
    h5_path = "./count/filtered_feature_bc_matrix.h5",
    filter_singlets = TRUE,
    aggr_dir = "./count/aggregation.csv",
    metadata_file = "./metadata.xlsx",
    contig_list = "./vdj/filtered_contig_annotations.csv"
  )
  assign(project_name, seurat_obj)
}, error = function(e) {
  stop("Error loading RNA data: ", e$message)
})

# --- Quality control filtering ---
tryCatch({
  seurat_obj <- get(project_name)
  seurat_obj <- qc_filtering(
    seurat_obj = seurat_obj,
    mito_percent_qc = 10,
    min_genes_qc = 500,
    max_genes_qc = 8000
  )
  assign(project_name, seurat_obj)
}, error = function(e) {
  stop("Error in QC filtering: ", e$message)
})

# --- Normalization and cell cycle scoring/regression ---
tryCatch({
  seurat_obj <- get(project_name)
  seurat_obj <- normalize_SCT(
    seurat_obj,
    genes_to_mask = NULL,
    cellcyclescore = TRUE,
    cellcyclenorm = TRUE,
    new_assay_name = "ccSCT"
  )
  assign(project_name, seurat_obj)
}, error = function(e) {
  stop("Error in normalization: ", e$message)
})

# --- Dimensionality reduction (PCA + UMAP) ---
tryCatch({
  seurat_obj <- get(project_name)
  seurat_obj <- create_UMAP(
    seurat_obj,
    automatic_npc = TRUE,
    npcs = NULL,
    min.dist = 0.01,
    k.param = 300,
    new_assay_name = "ccSCT"
  )
  assign(project_name, seurat_obj)
}, error = function(e) {
  stop("Error in UMAP creation: ", e$message)
})

# --- Clustering at specified resolution ---
tryCatch({
  seurat_obj <- get(project_name)
  seurat_obj <- findclusters_resolutions(
    seurat_obj,
    resolutions = c(1.8),
    n.iter = 10,
    save_plot = TRUE,
    new_assay_name = "ccSCT"
  )
  assign(project_name, seurat_obj)
}, error = function(e) {
  stop("Error in clustering: ", e$message)
})

# --- Save final Seurat object for downstream analysis ---
saveRDS(get(project_name), file = paste0(project_name, "_final_seurat_obj.rds"))
message("Pipeline finished successfully. Seurat object saved.")
