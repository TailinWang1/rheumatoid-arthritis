library(Seurat)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)

prefixes <- c("GSE_", 
              "GSE_", 
              "GSE_")
sample_names <- c("Sample_1", "Sample_2", "Sample_3")

processed_seurat_list <- list()

marker_genes <- c("CD14", "HLA-DRA", "CD74", "CD4", "IL2RA")

for (i in 1:length(prefixes)) {
  current_sample <- sample_names[i]
  cat("\nProcessing sample:", current_sample, "\n")
  
  prefix <- prefixes[i]
  counts <- ReadMtx(mtx = paste0(prefix, "matrix.mtx.gz"), 
                    features = paste0(prefix, "features.tsv.gz"), 
                    cells = paste0(prefix, "barcodes.tsv.gz"))
  obj <- CreateSeuratObject(counts = counts, project = current_sample, min.cells = 3, min.features = 200)
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:15, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:15, verbose = FALSE)
  
  p1 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5) + 
    ggtitle(paste(current_sample, "- Clusters"))
  
  p2 <- tryCatch({
    FeaturePlot(obj, features = intersect(marker_genes, rownames(obj)), ncol = 3)
  }, error = function(e) {
    ggplot() + theme_void() + ggtitle("Some genes not detected")
  })
  
  print(p1)
  print(p2) 
  
  processed_seurat_list[[current_sample]] <- obj
  rm(counts, obj, p1, p2); gc()
}

cat("\nBasic processing completed for all samples.\n")

cat("Loading Human Primary Cell Atlas reference data...\n")
ref_data <- HumanPrimaryCellAtlasData()

for (sample_name in names(processed_seurat_list)) {
  cat("\n--------------------------------------------------\n")
  cat("Annotating:", sample_name, "...\n")
  
  obj <- processed_seurat_list[[sample_name]]
  
  test_data <- GetAssayData(obj, assay = "RNA", layer = "data") 
  
  predictions <- SingleR(test = test_data, ref = ref_data, labels = ref_data$label.main)
  
  obj$SingleR_Labels <- predictions$labels
  Idents(obj) <- "SingleR_Labels"
  
  cat("\nCell annotation statistics for", sample_name, ":\n")
  print(table(obj$SingleR_Labels))
  
  p <- DimPlot(obj, reduction = "umap", group.by = "SingleR_Labels", 
               label = TRUE, label.size = 4, repel = TRUE) + 
    ggtitle(paste(sample_name, "- SingleR Annotation")) +
    theme(legend.position = "right")
  
  print(p)
  
  processed_seurat_list[[sample_name]] <- obj
}

cat("\nSingleR annotation and UMAP plotting completed for all samples.\n")

mono_genes <- c("CD14", "HLA-DRA", "CD74")
tcell_genes <- c("CD4", "IL2RA")

for (sample_name in names(processed_seurat_list)) {
  cat("\n==================================================\n")
  cat("Plotting sample:", sample_name, "...\n")
  
  obj <- processed_seurat_list[[sample_name]]
  
  Idents(obj) <- "SingleR_Labels"
  
  available_celltypes <- unique(Idents(obj))
  
  if ("Monocyte" %in% available_celltypes) {
    mono_obj <- subset(obj, idents = "Monocyte")
    
    valid_mono_genes <- intersect(mono_genes, rownames(mono_obj))
    
    if (length(valid_mono_genes) > 0) {
      p_mono <- FeaturePlot(mono_obj, features = valid_mono_genes, ncol = 3, pt.size = 0.8) +
        patchwork::plot_annotation(title = paste(sample_name, "- Monocyte Marker Genes"))
      print(p_mono)
    } else {
      cat("Warning: Specified marker genes not detected in Monocytes of", sample_name, "\n")
    }
  } else {
    cat("Warning: Monocyte cells not found in", sample_name, "\n")
  }
  
  if ("T_cells" %in% available_celltypes) {
    tcell_obj <- subset(obj, idents = "T_cells")
    
    valid_tcell_genes <- intersect(tcell_genes, rownames(tcell_obj))
    
    if (length(valid_tcell_genes) > 0) {
      p_tcell <- FeaturePlot(tcell_obj, features = valid_tcell_genes, ncol = 2, pt.size = 0.8) +
        patchwork::plot_annotation(title = paste(sample_name, "- T Cells Marker Genes"))
      print(p_tcell)
    } else {
      cat("Warning: Specified marker genes not detected in T cells of", sample_name, "\n")
    }
  } else {
    cat("Warning: T_cells not found in", sample_name, "\n")
  }
}

cat("\nFeaturePlots for specific cell populations completed for all samples.\n")