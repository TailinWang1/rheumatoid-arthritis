library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
#Please refer to the article for the GSMID used, or replace it with your own
setwd("/Users/")
sample_info <- list(
  "WT" = "GSMID",
  "Early" = "GSMID",
  "Adv1" = "GSMID",
  "Adv2" = "GSMID"
)

seurat_list <- list()

main_temp_dir <- file.path(tempdir(), "spatial_data")
dir.create(main_temp_dir, showWarnings = FALSE, recursive = TRUE)

for(sample in names(sample_info)){
  gsm_number <- sample_info[[sample]]
  
  cat("\nProcessing sample:", sample, "\n")
  
  matrix_file <- list.files(pattern = paste0(gsm_number, "_", sample, "_matrix\\.mtx\\.gz$"), full.names = TRUE)
  features_file <- list.files(pattern = paste0(gsm_number, "_", sample, "_features\\.tsv\\.gz$"), full.names = TRUE)
  barcodes_file <- list.files(pattern = paste0(gsm_number, "_", sample, "_barcodes\\.tsv\\.gz$"), full.names = TRUE)
  image_file <- list.files(pattern = paste0(gsm_number, "_", sample, "_tissue_hires_image\\.png\\.gz$"), full.names = TRUE)
  positions_file <- list.files(pattern = paste0(gsm_number, "_", sample, "_tissue_positions_list\\.csv\\.gz$"), full.names = TRUE)
  scalefactors_file <- list.files(pattern = paste0(gsm_number, "_", sample, "_scalefactors_json\\.json\\.gz$"), full.names = TRUE)
  
  cat("Files found:\n")
  cat("Matrix file:", matrix_file, "\n")
  cat("Features file:", features_file, "\n")
  cat("Barcodes file:", barcodes_file, "\n")
  cat("Image file:", image_file, "\n")
  cat("Positions file:", positions_file, "\n")
  cat("Scalefactors file:", scalefactors_file, "\n")
  
  if(any(sapply(list(matrix_file, features_file, barcodes_file, image_file, positions_file, scalefactors_file), length) == 0)) {
    cat("Error: Missing required files for sample", sample, "\n")
    next
  }
  
  cat("All required files found for sample", sample, "\n")
  
  temp_dir <- file.path(main_temp_dir, gsm_number)
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  R.utils::gunzip(matrix_file, destname = file.path(temp_dir, "matrix.mtx"), overwrite = TRUE)
  R.utils::gunzip(features_file, destname = file.path(temp_dir, "features.tsv"), overwrite = TRUE)
  R.utils::gunzip(barcodes_file, destname = file.path(temp_dir, "barcodes.tsv"), overwrite = TRUE)
  R.utils::gunzip(image_file, destname = file.path(temp_dir, "tissue_hires_image.png"), overwrite = TRUE)
  R.utils::gunzip(positions_file, destname = file.path(temp_dir, "tissue_positions_list.csv"), overwrite = TRUE)
  R.utils::gunzip(scalefactors_file, destname = file.path(temp_dir, "scalefactors_json.json"), overwrite = TRUE)
  
  cat("Creating Seurat object...\n")
  tryCatch({
    matrix <- ReadMtx(
      mtx = file.path(temp_dir, "matrix.mtx"),
      features = file.path(temp_dir, "features.tsv"),
      cells = file.path(temp_dir, "barcodes.tsv")
    )
    

    seurat_obj <- CreateSeuratObject(counts = matrix, project = sample, assay = "Spatial")
    
    image <- Read10X_Image(
      image.dir = temp_dir,
      image.name = "tissue_hires_image.png",
      filter.matrix = TRUE
    )
    
    seurat_obj[["image"]] <- image
    
    seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
    
    seurat_list[[sample]] <- seurat_obj
 
    saveRDS(seurat_obj, file = file.path(main_temp_dir, paste0(sample, "_seurat_object.rds")))
    
    cat("Finished processing sample:", sample, "\n")
  }, error = function(e) {
    cat("Error occurred while processing sample", sample, ":\n")
    cat(as.character(e), "\n")
  })
  
  cat("Temporary files for", sample, "are saved in:", temp_dir, "\n")
}

cat("\nAll processed data and Seurat objects are saved in:", main_temp_dir, "\n")
for (sample in names(seurat_list)) {
  print(sample)
  print(head(GetTissueCoordinates(seurat_list[[sample]])))
  print(table(seurat_list[[sample]]$seurat_clusters))
}
for (sample in names(seurat_list)) {
  print(SpatialDimPlot(seurat_list[[sample]], label = TRUE, label.size = 3) +
          ggtitle(paste("Spatial clusters -", sample)))
  print(SpatialFeaturePlot(seurat_list[[sample]], features = "nCount_Spatial") +
          ggtitle(paste("UMI counts -", sample)))
}

for (sample in names(seurat_list)) {
  if ("umap" %in% names(seurat_list[[sample]]@reductions)) {
    print(DimPlot(seurat_list[[sample]], reduction = "umap", label = TRUE) +
            ggtitle(paste("UMAP -", sample)))
  } else {
    cat("UMAP not found for sample:", sample, "\n")
  }
}
for (sample in names(seurat_list)) {
  print(sample)
  print(head(seurat_list[[sample]]@meta.data))
}


seurat_list1<-seurat_list
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

seurat_list <- list()
sample_info <- list(
  "WT" = "GSM5841550",
  "Early" = "GSM5841551",
  "Adv1" = "GSM5841552",
  "Adv2" = "GSM5841553"
)

seurat_list1 <- list()

main_temp_dir <- file.path(tempdir(), "spatial_data")
dir.create(main_temp_dir, showWarnings = FALSE, recursive = TRUE)
for(sample in names(sample_info)){
  gsm_number <- sample_info[[sample]]
  
  cat("\nProcessing sample:", sample, "\n")
  

  temp_dir <- file.path(main_temp_dir, gsm_number)
  
 
  image_file <- file.path(temp_dir, "tissue_hires_image.png")
  if (!file.exists(image_file)) {
    cat("Error: Image file not found:", image_file, "\n")
    next
  }
  

  cat("Creating Seurat object...\n")
  tryCatch({
   
    matrix <- ReadMtx(
      mtx = file.path(temp_dir, "matrix.mtx"),
      features = file.path(temp_dir, "features.tsv"),
      cells = file.path(temp_dir, "barcodes.tsv")
    )
  
    seurat_obj <- CreateSeuratObject(counts = matrix, project = sample, assay = "Spatial")
    

    image <- Read10X_Image(
      image.dir = temp_dir,
      filter.matrix = TRUE,
      image.name = "tissue_hires_image.png"
    )
    

    DefaultAssay(image) <- "Spatial"
    seurat_obj[["slice1"]] <- image
    

    seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
    
    seurat_list1[[sample]] <- seurat_obj
    
    cat("Finished processing sample:", sample, "\n")
  }, error = function(e) {
    cat("Error occurred while processing sample", sample, ":\n")
    cat(as.character(e), "\n")
  })
}

for (sample in names(seurat_list1)) {
  print(sample)
  print(head(GetTissueCoordinates(seurat_list1[[sample]])))
  print(table(seurat_list1[[sample]]$seurat_clusters))
}


library(Seurat)
library(ggplot2)
library(patchwork)


check_alignment <- function(seurat_obj, sample_name) {
  img <- seurat_obj@images$slice1@image
  coords <- GetTissueCoordinates(seurat_obj)
  norm_coords <- normalize_coordinates(coords)
  
  p <- ggplot() +
    annotation_custom(grid::rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc")), 
                      -Inf, Inf, -Inf, Inf) +
    geom_point(data = norm_coords, aes(x = x, y = y), color = "red", size = 0.5, alpha = 0.5) +
    theme_void() +
    coord_fixed() +
    ggtitle(paste("Original Image -", sample_name))
  
  print(p)
}


for (sample in names(seurat_list1)) {
  cat("\n\n=== Checking original image for sample:", sample, "===\n")
  check_alignment(seurat_list1[[sample]], sample)
  check_coords(seurat_list1[[sample]], sample)
}


library(Seurat)
library(ggplot2)
library(grid)

rotate_points <- function(coords, img_height) {
  x_new <- coords[, 2]
  y_new <- img_height - coords[, 1]
  return(cbind(x_new, y_new))
}

check_alignment_combined <- function(seurat_obj, sample_name) {
  img <- seurat_obj@images$slice1@image
  coords <- GetTissueCoordinates(seurat_obj)
  
  img_height <- nrow(img)
  img_width <- ncol(img)
  
  rotated_coords <- rotate_points(coords, img_height)
  
  plot_data <- data.frame(x = rotated_coords[,1], y = rotated_coords[,2])
  

  image_scale <- 0.90  
  
  p <- ggplot() +
    annotation_custom(rasterGrob(img, width=unit(image_scale,"npc"), height=unit(image_scale,"npc")), 
                      -Inf, Inf, -Inf, Inf) +
    geom_point(data = plot_data, aes(x = x, y = y), color = "red", size = 0.5, alpha = 0.5) +
    theme_void() +
    coord_fixed() +
    scale_x_continuous(limits = c(0, img_width)) +
    scale_y_continuous(limits = c(0, img_height)) +
    ggtitle(paste(sample_name, "- Rotated Points"))
  

  p <- p + theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))
  
  print(p)
  
  cat("Sample:", sample_name, "\n")
  cat("Image dimensions:", img_width, "x", img_height, "\n")
  cat("Rotated coordinate range X:", min(rotated_coords[,1]), "-", max(rotated_coords[,1]), "\n")
  cat("Rotated coordinate range Y:", min(rotated_coords[,2]), "-", max(rotated_coords[,2]), "\n\n")
}

for (sample_name in names(seurat_list1)) {
  check_alignment_combined(seurat_list1[[sample_name]], sample_name)
}


library(Seurat)
library(ggplot2)
library(grid)
library(patchwork)
library(viridis)

rotate_points <- function(coords, img_height) {
  x_new <- coords[, 2]
  y_new <- img_height - coords[, 1]
  return(cbind(x_new, y_new))
}

create_spatial_plot <- function(seurat_obj, feature, sample_name) {
  img <- seurat_obj@images$slice1@image
  coords <- GetTissueCoordinates(seurat_obj)
  
  img_height <- nrow(img)
  img_width <- ncol(img)
  
  rotated_coords <- rotate_points(coords, img_height)
  
  plot_data <- data.frame(x = rotated_coords[,1], y = rotated_coords[,2])
  
  if (feature %in% c("nCount_Spatial", "nFeature_Spatial")) {
    plot_data[[feature]] <- FetchData(seurat_obj, vars = feature)[[feature]]
    color_scale <- scale_color_viridis_c()
  } else if (feature == "seurat_clusters") {
    plot_data[[feature]] <- Idents(seurat_obj)
    color_scale <- scale_color_discrete()
  }
  

  image_scale <- 0.92
  
  p <- ggplot() +
    annotation_custom(rasterGrob(img, width=unit(image_scale,"npc"), height=unit(image_scale,"npc")), 
                      -Inf, Inf, -Inf, Inf) +
    geom_point(data = plot_data, aes(x = x, y = y, color = .data[[feature]]), size = 1.5, alpha = 0.7) +
    color_scale +
    theme_void() +
    coord_fixed() +
    scale_x_continuous(limits = c(0, img_width)) +
    scale_y_continuous(limits = c(0, img_height)) +
    labs(color = feature) +
    ggtitle(paste(sample_name, "-", feature))
  

  p <- p + theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))
  
  return(p)
}

WT_sample <- seurat_list1$WT
Early_sample <- seurat_list1$Early
Adv1_sample <- seurat_list1$Adv1
Adv2_sample <- seurat_list1$Adv2

p1 <- create_spatial_plot(WT_sample, "nCount_Spatial", "WT")
p2 <- create_spatial_plot(WT_sample, "nFeature_Spatial", "WT")
p3 <- create_spatial_plot(WT_sample, "seurat_clusters", "WT")
p4 <- create_spatial_plot(Early_sample, "nCount_Spatial", "Early")
p5 <- create_spatial_plot(Early_sample, "nFeature_Spatial", "Early")
p6 <- create_spatial_plot(Early_sample, "seurat_clusters", "Early")
p7 <- create_spatial_plot(Adv1_sample, "nCount_Spatial", "Adv1")
p8 <- create_spatial_plot(Adv1_sample, "nFeature_Spatial", "Adv1")
p9 <- create_spatial_plot(Adv1_sample, "seurat_clusters", "Adv1")
p10 <- create_spatial_plot(Adv2_sample, "nCount_Spatial", "Adv2")
p11 <- create_spatial_plot(Adv2_sample, "nFeature_Spatial", "Adv2")
p12 <- create_spatial_plot(Adv2_sample, "seurat_clusters", "Adv2")

Assays(WT_sample)

library(Seurat)
library(dplyr)

genes_of_interest <- c('Cd74', 'H2-Ab1', 'Mif', 'Cd44',"Cxcr4")

calculate_mean_expression <- function(seurat_object, genes) {
  avg_expression <- FetchData(seurat_object, vars = genes) %>%
    summarise_all(mean)
  return(avg_expression)
}

WT_mean_expression <- calculate_mean_expression(WT_sample, genes_of_interest)
Early_mean_expression <- calculate_mean_expression(Early_sample, genes_of_interest)
Adv1_mean_expression <- calculate_mean_expression(Adv1_sample, genes_of_interest)
Adv2_mean_expression <- calculate_mean_expression(Adv2_sample, genes_of_interest)
print("WT Sample Mean Expression:")
print(WT_mean_expression)

print("Early Sample Mean Expression:")
print(Early_mean_expression)

print("Adv1 Sample Mean Expression:")
print(Adv1_mean_expression)

print("Adv2 Sample Mean Expression:")
print(Adv2_mean_expression)

perform_t_test <- function(sample1, sample2, gene) {
  expr1 <- FetchData(sample1, vars = gene)
  expr2 <- FetchData(sample2, vars = gene)
  t_test_result <- t.test(expr1, expr2)
  return(t_test_result)
}

t_test_results <- data.frame(
  Gene = character(),
  Comparison = character(),
  Statistic = numeric(),
  P.Value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in genes_of_interest) {
  comparisons <- list(
    WT_vs_Early = perform_t_test(WT_sample, Early_sample, gene),
    WT_vs_Adv1 = perform_t_test(WT_sample, Adv1_sample, gene),
    WT_vs_Adv2 = perform_t_test(WT_sample, Adv2_sample, gene),
    Early_vs_Adv1 = perform_t_test(Early_sample, Adv1_sample, gene),
    Early_vs_Adv2 = perform_t_test(Early_sample, Adv2_sample, gene),
    Adv1_vs_Adv2 = perform_t_test(Adv1_sample, Adv2_sample, gene)
  )
  
  for (comparison in names(comparisons)) {
    t_test_results <- rbind(t_test_results, data.frame(
      Gene = gene,
      Comparison = comparison,
      Statistic = comparisons[[comparison]]$statistic,
      P.Value = comparisons[[comparison]]$p.value
    ))
  }
}

print("T-test results:")
print(t_test_results)


library(monocle3)
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

perform_trajectory_analysis <- function(seurat_obj, sample_name) {

  expression_matrix <- GetAssayData(seurat_obj, assay = "SCT", slot = "counts")
  cell_metadata <- seurat_obj@meta.data
  gene_metadata <- data.frame(
    gene_short_name = rownames(expression_matrix),
    row.names = rownames(expression_matrix)
  )
  
 
  cds <- new_cell_data_set(
    expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
  

  reducedDims(cds)$PCA <- seurat_obj@reductions$pca@cell.embeddings
  reducedDims(cds)$UMAP <- seurat_obj@reductions$umap@cell.embeddings
  
 
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds)
  

  p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = TRUE,
                   label_leaves = FALSE, label_branch_points = FALSE) +
    ggtitle(paste("Cell trajectory -", sample_name))
  
  
  genes_of_interest <- c("H2-Ab1", "Cd74")
  expression_plots <- list()
  

  custom_color <- function(n) {
    viridis(n, option = "D", direction = -1, begin = 0.1, end = 0.9)
  }
  
  for (gene in genes_of_interest) {
    if (gene %in% rownames(cds)) {
      p <- plot_cells(cds, 
                      genes = gene, 
                      label_cell_groups = FALSE,
                      label_leaves = FALSE, 
                      label_branch_points = FALSE,
                      show_trajectory_graph = TRUE,
                      cell_size = 1) +
        scale_color_gradientn(colours = custom_color(100)) +
        ggtitle(paste(gene, "expression on trajectory -", sample_name)) +
        theme(plot.title = element_text(size = 14, face = "bold"),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10))
      expression_plots[[gene]] <- p
    } else {
      p <- ggplot() + 
        ggtitle(paste(gene, "not found in dataset")) +
        theme_void()
      expression_plots[[gene]] <- p
    }
  }
  
  return(c(list(trajectory = p1), expression_plots))
}

