# modules.R
load_RNA_data <- function(h5_path, 
                          filter_singlets = TRUE, 
                          aggr_dir = NULL, 
                          metadata_file = NULL, 
                          sample_column_name = "sample_id",
                          hash_column_name = "hash",
                          contig_list = NULL) {
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
      if (!(sample_column_name %in% colnames(metadata))) {
      stop(paste("Column", sample_column_name, "does not exist in metadata file."))
      }
      if (!(hash_column_name %in% colnames(metadata))) {
      stop(paste("Column", hash_column_name, "does not exist in metadata file."))
    }
      metadata <- metadata %>%
        mutate(sample_id_full_prov = paste0(.[[sample_column_name]], "_Sample", .[[hash_column_name]]))

      seurat_obj$sample_id_full_prov <- paste0(seurat_obj$id_sample, "_", seurat_obj[[hash_column_name]])
    
      cell_metadata <- as.data.frame(seurat_obj@meta.data)
      cell_metadata <- cell_metadata %>%
        left_join(metadata %>% select(sample_id_full_prov, patient_code), by = "sample_id_full_prov")
    
      seurat_obj <- AddMetaData(seurat_obj, metadata = cell_metadata)
  }
  
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

qc_filtering <- function(seurat_obj, 
                         mito_percent_qc = 10, 
                         min_genes_qc = 500, 
                         max_genes_qc = 8000){
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  plot_QC(seurat_obj, mito_percent_qc, min_genes_qc, max_genes_qc, "Unfiltered", "orig.ident", "patient_code")
  seurat_obj <- subset(seurat_obj, subset= percent.mt < mito_percent_qc & nFeature_RNA > min_genes_qc & nFeature_RNA < max_genes_qc)
  plot_QC(seurat_obj, mito_percent_qc, min_genes_qc, max_genes_qc, "Filtered", "orig.ident", "patient_code")
  return(seurat_obj)
}

normalize_SCT <- function(seurat_obj, 
                          genes_to_mask = NULL, 
                          new_assay_name = "SCT", 
                          npcs = 50, 
                          cellcyclescore = FALSE, 
                          cellcyclenorm = FALSE){
    if (cellcyclescore) {
    if (!exists("cc.genes.updated.2019")) {
      cc.genes.updated.2019 <- Seurat::cc.genes.updated.2019
    }
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    
    seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", conserve.memory = T, verbose = T, assay = "RNA")
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, assay = "SCT", set.ident = TRUE)
  }
  vars_to_regress <- NULL
  if (cellcyclenorm) {
    if (!("S.Score" %in% colnames(seurat_obj@meta.data) && "G2M.Score" %in% colnames(seurat_obj@meta.data))) {
      stop("Cell cycle scores not found. Please set cellcyclescore = TRUE if you want to regress them.")
    }
    vars_to_regress <- c("S.Score", "G2M.Score")
    seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", conserve.memory = T, verbose = T, new.assay.name = new_assay_name, vars.to.regress = vars_to_regress, assay = "RNA") 
  }
  mask_var_genes(seurat_obj, assay_name = new_assay_name, add.mask = genes_to_mask)
  seurat_obj <- RunPCA(seurat_obj, verbose = T, assay = new_assay_name, reduction.name = paste0("pca_", new_assay_name), features = var.gene.names)
  return(seurat_obj)
}

create_UMAP <- function(seurat_obj, 
                        automatic_npc = TRUE, 
                        npcs = NULL, 
                        min.dist = 0.1, 
                        k.param = 150, 
                        new_assay_name = NULL){
  if (automatic_npc) {
    elbow_data <- data.frame(
      PC = 1:length(seurat_obj[[paste0("pca_", new_assay_name)]]@stdev),
      StandardDeviation = seurat_obj[[paste0("pca_", new_assay_name)]]@stdev
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
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_to_use, verbose = TRUE, min.dist = min.dist, spread = 3, seed.use = 041224, reduction.name = paste0("umap_", new_assay_name), reduction = paste0("pca_", new_assay_name))
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_to_use, verbose = T, k.param = k.param, reduction = paste0("pca_", new_assay_name))
  
  return(seurat_obj)
}

findclusters_resolutions <- function(seurat_obj, 
                                     resolutions = c(0.2, 0.8, 1.0, 1.6), 
                                     n.iter = 20, 
                                     save_plot = TRUE, 
                                     new_assay_name){
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, verbose = TRUE, resolution = res, cluster.name = paste0("clusters_", res), algorithm = 3, n.star = 50, n.iter = n.iter)
    if (save_plot){
    p <- DimPlot(seurat_obj, group.by = paste0("clusters_", res), reduction = paste0("umap_", new_assay_name), label = TRUE, shuffle = TRUE, raster = TRUE, repel = TRUE, label.box = TRUE) + NoLegend() +
      ggtitle(paste("UMAP - Resolution", res)) +
      theme_minimal() +
      theme(plot.title = element_text(size = 16, face = "bold"))
    ggsave(filename = file.path("./", paste0("UMAP_res", res, ".png")), plot = p, width = 7, height = 6)
    }
  }
  return(seurat_obj)
}


#### This function is waiting to be optimized and simplified
make_heatmap <- function(seurat_obj, 
                                  markers_df, 
                                  ident_col = "ident", 
                                  sample_col = "sample", 
                                  group_col = "group", 
                                  assay = "SCT", 
                                  slot = "data", 
                                  top_n = 100, 
                                  exclude_genes = NULL,
                                  highlight_genes = NULL,
                                  annotation_cols = c("sample", "group"),
                                  annotation_palettes = list(),
                                  output_file = NULL) {
  # Set active identity
  Idents(seurat_obj) <- ident_col
  group <- seurat_obj[[group_col]][,1]
  sample <- seurat_obj[[sample_col]][,1]
  ident <- seurat_obj[[ident_col]][,1]
    
  # Get top genes by expression difference
  top_genes <- markers_df %>%
    mutate(diff = abs(pct.1 - pct.2)) %>%
    mutate(gene = row.names(markers_df)) %>%
    filter(!(gene %in% exclude_genes)) %>%
    arrange(desc(diff)) %>%
    slice_head(n = top_n) %>%
    pull(gene) %>%
    unique()
  
  # Get expression matrix
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)[top_genes, ]
  expr_matrix <- t(scale(t(as.matrix(expr_matrix))))
  
  # Apply bold styling for highlight_genes
  all_genes <- rownames(expr_matrix)
  if (is.null(highlight_genes)) highlight_genes <- character(0)
  row_font_gp <- gpar(
    fontface = ifelse(all_genes %in% highlight_genes, "bold", "plain"),
    col = ifelse(all_genes %in% highlight_genes, "black", "grey60"),
    fontsize = 6
  )
  
  # Define order
  column_order <- order(group, ident, sample)
  
  # Extraer datos de anotación
  annotation_data <- seurat_obj@meta.data[, annotation_cols, drop = FALSE]
  annotation_data[] <- lapply(annotation_data, function(x) {
    if (!is.factor(x)) factor(x) else x
  })
  
  # Paletas automáticas si no se pasan
  auto_pals <- lapply(annotation_data, function(x) {
    lvls <- levels(x)
    n <- length(lvls)
    max_n <- brewer.pal.info["Dark2", "maxcolors"]
    base_pal <- brewer.pal(min(n, max_n), "Dark2")
    pal <- colorRampPalette(base_pal)(n)
    setNames(pal, lvls)
    })
  
  # Fusionar con paletas personalizadas
  full_palette <- modifyList(auto_pals, annotation_palettes)
  
  # Crear anotaciones
  top_anno <- HeatmapAnnotation(
    df = annotation_data,
    col = full_palette,
    show_legend = rep(TRUE, length(annotation_cols)),
    gap = unit(1, "mm"),
    simple_anno_size = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 6)
  )
  
  # ---- HEATMAP ----
  ht <- Heatmap(
    matrix = expr_matrix,
    name = "Expression",
    column_split = factor(group),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_order = column_order,
    column_title_rot = 0,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_row_dend = FALSE,
    cluster_column_slices = TRUE,
    column_title_gp = gpar(fontsize = 7),
    column_gap = unit(0.5, "mm"),
    row_names_gp = row_font_gp,
    top_annotation = top_anno,
    use_raster = FALSE
  )
  
  # Output
  if (!is.null(output_file)) {
    png(filename = output_file, width = 2000, height = 1400, res = 200)
    draw(ht)
    dev.off()
  } else {
    draw(ht)
  }
}




