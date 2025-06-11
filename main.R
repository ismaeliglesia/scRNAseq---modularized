# main.R

R_directory <- ""  ### Choose the directory of functions
setwd(R_directory)

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

R_directory <- ""  ### Choose the directory of data
setwd(R_directory)
project_name <- "" ### Write the project name in a brief way, i.e., project1_2025

assign(project_name, load_RNA_data("./count/filtered_feature_bc_matrix.h5", 
                                   filter_singlets = TRUE, 
                                   aggr_dir = "./count/aggregation.csv", 
                                   metadata_file = NULL,
                                   contig_list = "./count/filtered_contig_annotations.csv"))

assign(project_name, qc_filtering(seurat_obj = get(project_name), mito_percent_qc = 10, min_genes_qc = 500, max_genes_qc = 8000))

assign(project_name, normalize_SCT(get(project_name), genes_to_mask = "CD7CAR", cellcyclescore = FALSE, cellcyclenorm = FALSE, new_assay_name = "SCT"))

assign(project_name, create_UMAP(get(project_name), automatic_npc = TRUE, npcs = NULL, min.dist = 0.01, k.param = 150))

assign(project_name, findclusters_resolutions(get(project_name), resolutions = c(0.2, 0.6), n.iter = 10, save_plot = TRUE))
