# modules.R
load_RNA_data <- function(h5_path, filter_singlets = TRUE, aggr_dir = NULL, metadata_file = NULL, contig_list = NULL) {
  RNA.data <- Read10X_h5(h5_path)
  seurat_obj <- CreateSeuratObject(counts = RNA.data$`Gene Expression`)
  hto_counts <- CreateAssayObject(RNA.data$`Antibody Capture`)
  seurat_obj[['ADT']] <- hto_counts 
  if ("Antibody Capture" %in% names(RNA.data)) {
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
    seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
    seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
    seurat_obj <- HTODemux(seurat_obj, assay = "ADT", positive.quantile = 0.99)
    if (filter_singlets) {
      seurat_obj <- subset(seurat_obj, subset = ADT_classification.global == "Singlet")
    }
  }
  if (!is.null(aggr_dir)) {
    cellcodes <- as.data.frame(seurat_obj@active.ident)
    cellcodes$barcodes <- rownames(cellcodes)
    cellcodes$libcodes <- as.integer(gsub(pattern = ".+-", replacement = "", cellcodes$barcodes))
    
    samples <- read.csv(aggr_dir, stringsAsFactors = FALSE)
    cellcodes$id_prov <- as.vector(samples$sample_id[cellcodes$libcodes])
    
    id_prov <- cellcodes["id_prov"]
    seurat_obj <- AddMetaData(seurat_obj, metadata = id_prov)
  }
  if (!is.null(contig_list)) {
    contig_list <- read.csv(contig_list)
    contig_list <- createHTOContigList(contig_list, seurat_obj, group.by = "id_prov")
    contig_list <- combineTCR(contig_list)
    contig_list <- lapply(contig_list, function(df) {
      df$barcode <- sub("^[0-9]+_[A-Z]+[0-9]?_", "", df$barcode)
      return(df)
    })
    
    seurat_obj <- combineExpression(
      contig_list, seurat_obj,
      cloneCall = 'CTstrict',
      proportion = TRUE,
      cloneSize = c(Rare = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1)
    )
    
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
      mutate(cloneSize = ifelse(is.na(cloneSize), "noTCR_data", cloneSize)) %>%
      mutate(cloneSize = recode(cloneSize,
                                "Hyperexpanded (0.1 < X <= 1)" = "Hyperexpanded",
                                "2" = "High", "3" = "Medium", "4" = "Rare"))
    
    seurat_obj@meta.data$cloneSize <- factor(
      seurat_obj@meta.data$cloneSize,
      levels = c("High", "Medium", "Rare", "noTCR_data")
    )
    
  }
  if (!is.null(metadata_file)) {
    metadata <- read_xlsx(metadata_file)
    metadata <- metadata %>%
      mutate(sample_id_full_prov = paste0(id_prov, "_Sample", hash))
    seurat_obj$sample_id_full_prov <- paste0(seurat_obj$id_prov, "_", seurat_obj$hash.ID)
    cell_metadata <- as.data.frame(seurat_obj@meta.data)
    cell_metadata <- cell_metadata %>%
      left_join(metadata %>% select(sample_id_full_prov, patient_code), by = "sample_id_full_prov")
    
    seurat_obj <- AddMetaData(seurat_obj, metadata = cell_metadata)
  }
  return(seurat_obj)
}

qc_filtering <- function(seurat_obj, mito_percent_qc, min_genes_qc, max_genes_qc){
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  plot_QC(seurat_obj, mito_percent_qc, min_genes_qc, max_genes_qc, "Unfiltered", "orig.ident", "patient_code")
  seurat_obj <- subset(seurat_obj, subset= percent.mt < mito_percent_qc & nFeature_RNA > min_genes_qc & nFeature_RNA < max_genes_qc)
  plot_QC(seurat_obj, mito_percent_qc, min_genes_qc, max_genes_qc, "Filtered", "orig.ident", "patient_code")
  return(seurat_obj)
}

normalize_SCT <- function(seurat_obj, genes_to_mask = NULL, new_assay_name = "SCT", npcs = 50, cellcyclescore = FALSE, cellcyclenorm = FALSE){
  seurat_obj <- mask_var_genes(seurat_obj, assay_name = 'RNA', add.mask = genes_to_mask)
  if (cellcyclescore) {
    # Cargar listas de genes S y G2M si no están
    if (!exists("cc.genes.updated.2019")) {
      cc.genes.updated.2019 <- Seurat::cc.genes.updated.2019
    }
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, assay = "RNA")
  }
  vars_to_regress <- NULL
  if (cellcyclenorm) {
    if (!("S.Score" %in% colnames(seurat_obj@meta.data) && "G2M.Score" %in% colnames(seurat_obj@meta.data))) {
      stop("Cell cycle scores not found. Please set cellcyclescore = TRUE if you want to regress them.")
    }
    vars_to_regress <- c("S.Score", "G2M.Score")
  }
  seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", conserve.memory = T, verbose = T, new.assay.name = new_assay_name, vars.to.regress = vars_to_regress)
  seurat_obj <- RunPCA(seurat_obj, verbose = F, npcs = npcs)
  return(seurat_obj)
  }

create_UMAP <- function(seurat_obj, automatic_npc = TRUE, npcs = NULL, min.dist = 0.1, k.param = 150){
  if (automatic_npc) {
    elbow_data <- data.frame(
      PC = 1:length(seurat_obj[["pca"]]@stdev),
      StandardDeviation = seurat_obj[["pca"]]@stdev
    )
    elbow_point <- find_elbow(elbow_data$StandardDeviation)
    dims_to_use <- 1:elbow_point
    message("PCs number (elbow method): ", elbow_point)
  } else {
    if (is.null(npcs)) {
      stop("Please, stablish the number of principal components (npcs) if automatic_npc = FALSE.")
    }
    dims_to_use <- 1:npcs
    message("PCs number (manual method): ", npcs)
  }
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_to_use, verbose = TRUE, min.dist = min.dist, spread = 3, reduction = 'pca', seed.use = 041224)
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_to_use, verbose = T, k.param = k.param, reduction = 'pca')
  
  return(seurat_obj)
}

findclusters_resolutions <- function(seurat_obj, resolutions = c(0.2, 0.4, 0.8, 1.0), n.iter = 20, save_plot = TRUE){
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, verbose = TRUE, resolution = res, cluster.name = paste0("clusters_", res), algorithm = 3, n.star = 50, n.iter = n.iter)
    if (save_plot){
    p <- DimPlot(seurat_obj, group.by = paste0("clusters_", res), reduction = "umap", label = TRUE, shuffle = TRUE, raster = TRUE, repel = TRUE, label.box = TRUE) + NoLegend() +
      ggtitle(paste("UMAP - Resolution", res)) +
      theme_minimal() +
      theme(plot.title = element_text(size = 16, face = "bold"))
    ggsave(filename = file.path("./", paste0("UMAP_res", res, ".png")), plot = p, width = 7, height = 6)
    }
  }
  return(seurat_obj)
}



