library(hdf5r)
library(Seurat) 
library(dplyr)
library(multtest)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(future)
library(glmGamPoi)
#Please refer to the article for the GSMID used, or replace it with your own
#HC
data_seurat<-readRDS("NCpbmc3.rds")
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
#RA
data_sample <- Read10X_h5("GSM4819747_RA_filtered_feature_bc_matrix.h5")
data_seurat <- CreateSeuratObject(data_sample, project = "data_sample")
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")

VlnPlot(data_seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

plot1 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

CombinePlots(plots = list(plot1, plot2))

pbmc <- subset(data_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


ncol(GetAssayData(pbmc, assay = "RNA", slot = "counts"))

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(pbmc))

plot1 <- VariableFeaturePlot(pbmc)


plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1
plot2


all.genes <- rownames(pbmc)


pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)


JackStrawPlot(pbmc, dim = 1:20)


ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:17)


pbmc <- FindClusters(pbmc, resolution = 1.2)

head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:17)
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:17)
DimPlot(pbmc, reduction = "tsne",label = TRUE)
VlnPlot(pbmc, features = c("ROMO1"), pt.size = 0.001)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

for (i in 0:20) {
cat("Cluster", i, "markers:\n")
markers <- FindMarkers(pbmc, ident.1 = i, min.pct = 0.1, logfc.threshold = 0.1)
print(head(markers))
}

all.genes <- rownames(pbmc)


head(markers)
#RA
#monocytes，PMID: 32989127
FeaturePlot(pbmc,features = c("S100A8"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD14"),reduction = "tsne")
FeaturePlot(pbmc,features = c("S100A9"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD14"),reduction = "tsne")#cd14+
FeaturePlot(pbmc,features = c("FCGR3A"),reduction = "tsne")#CD16+
FeaturePlot(pbmc, features = c("HLA-DRA"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD74"), reduction = "tsne")
FeaturePlot(pbmc, features = c("MIF"), reduction = "tsne")
#macrophages，PMID: 32989127
FeaturePlot(RA_pbmc1,features = c("AIF1"),reduction = "tsne")
FeaturePlot(RA_pbmc1,features = c("PSAP"),reduction = "tsne")
FeaturePlot(pbmc,features = c("IFITM3"),reduction = "tsne")
FeaturePlot(pbmc,features = c("LST1"),reduction = "tsne")
FeaturePlot(pbmc,features = c("SERPINA1"),reduction = "tsne")
FeaturePlot(RA_pbmc1,features = c("CD163"),reduction = "tsne")#M2
#NK cells，PMID: 32989127
FeaturePlot(pbmc,features = c("NCR3"),reduction = "tsne")
FeaturePlot(pbmc,features = c("NKG7"),reduction = "tsne")
FeaturePlot(pbmc,features = c("KLRD1"),reduction = "tsne")
FeaturePlot(pbmc,features = c("GNLY"),reduction = "tsne")
#B cells，PMID: 38159454.
FeaturePlot(pbmc,features = c("CD79A"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD74"),reduction = "tsne")
FeaturePlot(pbmc,features = c("PAX5"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD79B"),reduction = "tsne")
#T cells，PMID: 38159454.
FeaturePlot(pbmc,features = c("CD8B"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CCR7"),reduction = "tsne")
FeaturePlot(RA_pbmc1,features = c("CD8A"),reduction = "tsne")
FeaturePlot(pbmc,features = c("IL7R"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD4"),reduction = "tsne")
FeaturePlot(pbmc,features = c("IL2RA"),reduction = "tsne")
#Dendritic cells，PMID: 38159454
FeaturePlot(pbmc,features = c("CD"),reduction = "tsne")

RA_pbmc1 <- RenameIdents(pbmc, 
                         `0` = "T cells CD8+", 
                         `1` = "T cells CD4+", 
                         `2` = "T cells CD4+",
                         `3` = "T cells CD8+", 
                         `4` = "NK cells", 
                         `5` = "Monocytes CD14+HLA-DR+CD74+",
                         `6` = "T cells CD8+", 
                         `7` = "B cells", 
                         `8` = "B cells",
                         `9` = "Monocytes CD16+HLA-DR+CD74+", 
                         `10` = "T cells CD4+IL2RA+", 
                         `11` = "T cells CD8+",
                         `12` = "T cells CD8+") 


VlnPlot(RA_pbmc1, features = c("ROMO1"), slot = "counts", log = TRUE)
DimPlot(RA_pbmc1, reduction = "tsne",label = TRUE)

FeaturePlot(RA_pbmc1, features = c("ROMO1"), reduction = "tsne",label = TRUE)
FeaturePlot(pbmc, features = c("ROMO1"), reduction = "tsne",label = TRUE)
DimPlot(pbmc, reduction = "tsne",label = TRUE)


#HC
#T cells，PMID: 38159454.
FeaturePlot(pbmc,features = c("CD8B"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CCR7"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD8A"),reduction = "tsne")
FeaturePlot(pbmc,features = c("IL7R"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD4"),reduction = "tsne")#
FeaturePlot(pbmc,features = c("IL2RA"),reduction = "tsne")
#monocytes，PMID: 32989127
FeaturePlot(pbmc,features = c("S100A8"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD14"),reduction = "tsne")
FeaturePlot(pbmc,features = c("S100A9"),reduction = "tsne")
FeaturePlot(pbmc,features = c("FCGR3A"),reduction = "tsne")#
FeaturePlot(pbmc, features = c("HLA-DRA"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD74"), reduction = "tsne")
#B cells，PMID: 38159454.
FeaturePlot(pbmc,features = c("CD79A"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD74"),reduction = "tsne")
FeaturePlot(pbmc,features = c("PAX5"),reduction = "tsne")
FeaturePlot(pbmc,features = c("CD79B"),reduction = "tsne")
#NK cells，PMID: 32989127
FeaturePlot(pbmc,features = c("GZMB"),reduction = "tsne")
FeaturePlot(pbmc,features = c("NKG7"),reduction = "tsne")
FeaturePlot(pbmc,features = c("FCGR3A"),reduction = "tsne")
FeaturePlot(pbmc,features = c("GZMH"),reduction = "tsne")
#Macrophage
FeaturePlot(pbmc,features = c("AIF1"),reduction = "tsne")
FeaturePlot(pbmc,features = c("PSAP"),reduction = "tsne")
FeaturePlot(pbmc,features = c("IFITM3"),reduction = "tsne")
FeaturePlot(pbmc,features = c("LST1"),reduction = "tsne")
FeaturePlot(pbmc,features = c("SERPINA1"),reduction = "tsne")

HC_pbmc1 <- RenameIdents(pbmc, 
                         `0` = "T cells CD4+", 
                         `1` = "Monocytes CD14+HLA-DR+CD74+", 
                         `2` = "T cells CD8+",
                         `3` = "Monocytes CD14+HLA-DR+CD74+", 
                         `4` = "Monocytes CD14+HLA-DR+CD74+", 
                         `5` = "Monocytes CD16+HLA-DR+CD74+", 
                         `6` = "T cells CD4+", 
                         `7` = "T cells CD8+", 
                         `8` = "B cells", 
                         `9` = "NK cells",
                         `10` = "Monocytes CD14+HLA-DR+CD74+", 
                         `11` = "T cells CD8+", 
                         `12` = "NK cells", 
                         `13` = "T cells CD4+IL2RA+")

VlnPlot(HC_pbmc1, features = c("ROMO1"), slot = "counts", log = TRUE)
FeaturePlot(HC_pbmc1, features = c("ROMO1"), reduction = "tsne",label = TRUE)
DimPlot(HC_pbmc1, reduction = "tsne",label = TRUE)
FeaturePlot(HC_pbmc1, features = c("CD14"), reduction = "tsne",label = TRUE)

romo1_expression <- FetchData(RA_pbmc1, vars = "ROMO1")
print(Idents(RA_pbmc1))

RA_pbmc1@meta.data$celltype <- Idents(RA_pbmc1)
HC_pbm1c@meta.data$celltype <- Idents(HC_pbmc1)

DefaultAssay(RA) <- "RNA"
DefaultAssay(HC) <- "RNA"

RA_pbmc1$celltype <- Idents(RA_pbmc1)

HC_pbmc1$celltype <- Idents(HC_pbmc1)

RA_data <- FetchData(RA_pbmc1, vars = c("ROMO1", "celltype"))
RA_data$sample <- "RA"

HC_data <- FetchData(HC_pbmc1, vars = c("ROMO1", "celltype"))
HC_data$sample <- "HC"


combined_data <- rbind(RA_data, HC_data)

library(dplyr)
library(ggpubr)

cell_types <- unique(combined_data$celltype)

results <- list()


for (cell_type in cell_types) {
  
  cell_data <- combined_data %>% filter(celltype == cell_type)

  test_result <- wilcox.test(ROMO1 ~ sample, data = cell_data)
  

  results[[cell_type]] <- list(
    p_value = test_result$p.value,
    median_RA = median(cell_data$ROMO1[cell_data$sample == "RA"]),
    median_HC = median(cell_data$ROMO1[cell_data$sample == "HC"])
  )
  

  p <- ggboxplot(cell_data, x = "sample", y = "ROMO1", 
                 color = "sample", palette = "jco",
                 add = "jitter") +
    stat_compare_means(method = "wilcox.test") +
    ggtitle(paste("ROMO1 expression in", cell_type))
  

  ggsave(filename = paste0("ROMO1_", cell_type, ".png"), plot = p, width = 6, height = 4)
}


for (cell_type in names(results)) {
  cat("Cell type:", cell_type, "\n")
  cat("P-value:", results[[cell_type]]$p_value, "\n")
  cat("Median ROMO1 in RA:", results[[cell_type]]$median_RA, "\n")
  cat("Median ROMO1 in HC:", results[[cell_type]]$median_HC, "\n")
  cat("\n")
}


pdf("ROMO1_expression_boxplots.pdf", width = 5, height = 4)

for (cell_type in cell_types) {
 
  cell_data <- combined_data %>% filter(celltype == cell_type)
  

  test_result <- wilcox.test(ROMO1 ~ sample, data = cell_data)
  

  results[[cell_type]] <- list(
    p_value = test_result$p.value,
    median_RA = median(cell_data$ROMO1[cell_data$sample == "RA"]),
    median_HC = median(cell_data$ROMO1[cell_data$sample == "HC"])
  )
  

  p <- ggboxplot(cell_data, x = "sample", y = "ROMO1", 
                 color = "sample", palette = "jco",
                 add = "jitter") +
    stat_compare_means(method = "wilcox.test") +
    ggtitle(paste("ROMO1 expression in", cell_type))
  
  
  print(p)
}

dev.off()

for (cell_type in names(results)) {
  cat("Cell type:", cell_type, "\n")
  cat("P-value:", results[[cell_type]]$p_value, "\n")
  cat("Median ROMO1 in RA:", results[[cell_type]]$median_RA, "\n")
  cat("Median ROMO1 in HC:", results[[cell_type]]$median_HC, "\n")
  cat("\n")
}


library(CellChat)
library(Seurat)

DefaultAssay(RA_pbmc1) <- "RNA"
RA_pbmc1 <- NormalizeData(RA_pbmc1, assay = "RNA", normalization.method = "LogNormalize")
RA_pbmc1 <- JoinLayers(RA_pbmc1)
cellchat <- createCellChat(object = RA_pbmc1,
                           meta = RA_pbmc1@meta.data,
                           group.by = "celltype")

CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat@DB <- CellChatDB.use

# This step is necessary even if using the whole database
cellchat <- subsetData(cellchat) 
# do parallel 
future::plan("multisession", workers = 1)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)
head(df.net)
write.csv(df.net, "RA_cell-cell_communications.all.csv")

#access the the inferred communications at the level of signaling pathways
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")

#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
levels(cellchat@idents)
df.net2 <- subsetCommunication(cellchat, sources.use = c("Epi"), targets.use = c("Fibroblast" ,"T")) 

#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))


cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
str(cellchat@net)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellchat@net$weight

cell_types_of_interest <- c("Monocytes CD14+HLA-DR+CD74+", "T cells CD4+IL2RA+")


par(mfrow = c(1,3), xpd=TRUE)


for (cell_type in cell_types_of_interest) {
  i <- which(rownames(mat) == cell_type)
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = cell_type,
                   vertex.label.cex = 0.7,  
                   edge.width.max = 8,      
                   vertex.size.max = 8)     
}
#Chord diagram
netVisual_chord_gene(cellchat, sources.use = c(4,7), 
                     targets.use = c(4,7), 
                     small.gap =3,
                     slot.name = "netP", legend.pos.x = 0.5)
netVisual_chord_gene(cellchat, sources.use = c(3,7), 
                     targets.use = c(3,7), 
                     small.gap = 5,  # 
                     slot.name = "netP", legend.pos.x = 0.5)
netVisual_chord_gene(cellchat, sources.use = c(3,7), 
                     targets.use = c(3,7), 
                     small.gap = 5,
                     slot.name = "netP", 
                     legend.pos.x = 0.5,
                     title.name = "Gene interactions")

pdf("chord_plot_large.pdf", width = 15, height = 15)


netVisual_chord_gene(cellchat, sources.use = c(4,7), 
                     targets.use = c(4,7), 
                     small.gap = 5,
                     slot.name = "netP", 
                     legend.pos.x = 0.5,
                     title.name = "Gene interactions")


dev.off()
cellchat@netP$pathways
pathways.show <- cellchat@netP$pathways
levels(cellchat@idents)   
vertex.receiver = c(4,7) 
#Chord diagram
netVisual_aggregate(cellchat, signaling = pathways.show,  
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
pdf("chord_diagram.pdf", width = 16, height = 12)
netVisual_chord_gene(cellchat, sources.use = c(1,4,7), targets.use = c(1,4,7), 
                     small.gap = 3, slot.name = "netP", legend.pos.x = 1,legend.pos.y = 30)


plotGeneExpression(cellchat, signaling = "MIF")


dev.off()
#Heatmap
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,7),
                 targets.use = c(1:5), 
                 remove.isolate = FALSE)


#Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


cell_types <- c("T cells CD8+", 
                "T cells CD4+", 
                "NK cells", 
                "Monocytes CD14+HLA-DR+CD74+", 
                "B cells", 
                "T cells CD4+IL2RA+")


mat <- cellchat@net$weight


mat <- mat[cell_types, cell_types]


groupSize <- as.numeric(table(cellchat@idents)[cell_types])


par(mfrow = c(2,3), xpd=TRUE)


for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}


pathways.show <- cellchat@netP$pathways


ht1 <- netAnalysis_contribution(cellchat, signaling = pathways.show)


print(str(ht1))


if (is(ht1, "ggplot")) {
  print(ht1)
} else {
  
  library(pheatmap)
  pheatmap(ht1, main = "Pathway Importance")
}


ht2 <- netAnalysis_contribution(cellchat, signaling = pathways.show)

print(str(ht2))


if (is(ht2, "ggplot")) {
  print(ht2)
} else {
  pheatmap(ht2, main = "Pathway Importance Between Cell Pairs")
}


pathways.show <- cellchat@netP$pathways

ht1 <- netAnalysis_contribution(cellchat, signaling = pathways.show)


print(ht1)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


cellchat <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


gg1 <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
print(gg1)


gg2 <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", signaling.role = "receiver")
print(gg2)


ht_strength <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = pathways.show, title = "Outgoing signaling pathway strength")
print(ht_strength)


ht_strength_incoming <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = pathways.show, title = "Incoming signaling pathway strength")
print(ht_strength_incoming)


pathways.show <- cellchat@netP$pathways


ht1 <- netAnalysis_contribution(cellchat, signaling = pathways.show)


print(ht1)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

print(names(cellchat))
print(names(cellchat@net))
print(names(cellchat@netP))
pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
gg_bubble <- netVisual_bubble(cellchat, pairLR.use = pairLR.use, sources.use = 1:4, targets.use = 5:8, remove.isolate = FALSE)
print(gg_bubble)

gg_heatmap <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
print(gg_heatmap)

ht_strength <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = pathways.show, title = "Outgoing signaling pathway strength")
print(ht_strength)

ht_strength_incoming <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = pathways.show, title = "Incoming signaling pathway strength")
print(ht_strength_incoming)
