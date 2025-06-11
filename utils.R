# utils.R

pb_save <- function(plot, name, width.numb, height.numb, filepath.name) {
  plot.filename <- paste0(filepath.name, name, ".png")
  ggsave(plot = plot, filename = plot.filename, width = width.numb, height = height.numb)
}

plot_QC <- function(obj, mito_percent_qc = mito_percent_qc, min_genes_qc = 500, max_genes_qc = 8000,
                    title = "QC", ident.1 = "orig.ident", ident.2 = "patient_code", filepath = "./") {
  
  obj <- SetIdent(obj, value = ident.1)
  v1 <- VlnPlot(obj, features = "percent.mt", pt.size = 0, log = F, raster = T) + geom_hline(yintercept = mito_percent_qc)
  v2 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0, log = T, raster = T) +
    geom_hline(yintercept = min_genes_qc) + geom_hline(yintercept = max_genes_qc)
  v3 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0, log = T, raster = T)
  
  comb.plot <- (v1 / v2 / v3) +
    plot_annotation(title = title) & theme(plot.title = element_text(size = 18, face = "bold")) & NoLegend()
  
  obj <- SetIdent(obj, value = ident.2)
  v4 <- VlnPlot(obj, features = "percent.mt", pt.size = 0, log = F, raster = T) + geom_hline(yintercept = mito_percent_qc)
  v5 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0, log = T, raster = T) +
    geom_hline(yintercept = min_genes_qc) + geom_hline(yintercept = max_genes_qc)
  v6 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0, log = T, raster = T)
  
  pt.plot <- (v4 / v5 / v6) +
    plot_annotation(title = title) & theme(plot.title = element_text(size = 18, face = "bold")) & NoLegend()
  
  pb_save(comb.plot, paste0("QCplot_", ident.1, "_", title), 3, 10, filepath)
  pb_save(pt.plot, paste0("QCplot_", ident.2, "_", title), 10, 10, filepath)
  
  return(list(combined = comb.plot, per_sample = pt.plot))
}

mask_var_genes <- function(seurat_obj, assay_name, add.mask) {
  # Variable gene names
  var.gene.names <- VariableFeatures(seurat_obj, assay = assay_name)
  
  # Get the assay data and gene index
  counts <- GetAssayData(seurat_obj, assay = assay_name)
  gene_index <- Features(seurat_obj@assays[[assay_name]])
  
  # Index TCR and sex genes to remove
  rm_TCR_index <- grep(pattern = "^TR[AB][VJC]", gene_index)
  rm_sex_index <- unlist(sapply(c("XIST", "SRY", "RPS4Y1", "RPS4Y2"), function(gene) grep(pattern = gene, gene_index)))
  rm_index <- unique(c(rm_TCR_index, rm_sex_index))
  
  # Get the gene names to be removed
  rm_matrix <- counts[rm_index, ]
  rm_genes <- as.character(rm_matrix@Dimnames[[1]])
  rm_genes <- c(rm_genes, add.mask) # add masked gene
  
  # Update the Seurat object with the new variable gene list
  var.gene.names <- setdiff(var.gene.names, rm_genes)
  VariableFeatures(seurat_obj[[assay_name]]) <- var.gene.names
  var.gene.names <<- var.gene.names
  print(head(var.gene.names, 20))
  
  return(seurat_obj)
}

find_elbow <- function(stdev) {
  n <- length(stdev)
  line_start <- c(1, stdev[1])
  line_end <- c(n, stdev[n])
  point_line_dist <- function(x, y, start, end) {
    num <- abs((end[2] - start[2]) * x - (end[1] - start[1]) * y + end[1]*start[2] - end[2]*start[1])
    den <- sqrt((end[2] - start[2])^2 + (end[1] - start[1])^2)
    return(num / den)
  }
  distances <- sapply(1:n, function(i) {
    point_line_dist(i, stdev[i], line_start, line_end)
  })
  which.max(distances)
}
