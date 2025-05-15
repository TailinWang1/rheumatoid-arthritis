library(readxl)
library(dplyr)
library(limma)
library(GEOquery)
library(impute)
library(WGCNA)
library(readr)
library(preprocessCore)
#Please refer to the article for the GSMID used, or replace it with your own
RA_gset <- getGEO('GSEID', destdir=".",AnnotGPL = T,getGPL = T)
RA_exp<-exprs(RA_gset[[1]])
RA_GPL<-fData(RA_gset[[1]])
RA_gpl<- RA_GPL[, c(1, 3)]
RA_exp<-as.data.frame(RA_exp)
RA_exp$ID<-rownames(RA_exp)
RA_exp_symbol<-merge(RA_exp,RA_gpl,by="ID")
RA_exp_symbol<-na.omit(RA_exp_symbol)
table(duplicated(RA_exp_symbol$`Gene symbol`))
RA_datExpr02<-avereps(RA_exp_symbol[,-c(1,ncol(RA_exp_symbol))],ID=RA_exp_symbol$`Gene symbol`)
datExpr0<- t(RA_datExpr02)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  
  abline(h = 250000, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 250000, minSize = 0)
keepSamples = (clust==1)
datExpr = datExpr1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)
datExpr = as.data.frame(datExpr)

datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("1RA03Threshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft$powerEstimate

net = blockwiseModules(datExpr, power = 1,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)
pdf("2RA03module.pdf",width = 12, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}

RAsamples <- read_excel('RA1.xlsx')

row_names <- RAsamples[[1]]

RAsamples <- RAsamples %>% select(-1)

RAsamples_df <- as.data.frame(RAsamples)

row.names(RAsamples_df) <- row_names

print(RAsamples_df)

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
dim(MEsWW)   
dim(RAsamples_df) 
modTraitCor = cor(MEsWW, RAsamples_df, use = "p")
colnames(MEsWW)
modlues=MEsWW
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3RA03Module-trait.pdf",width = 6, height = 8)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(RAsamples_df), yLabels = names(MEsWW), cex.lab = 0.8,  yColorWidth=0.02, 
               xColorWidth = 0.04,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()
print(blue.expr)
blueGenes <- colnames(blue.expr)
brownGenes<- colnames(brown.expr)
RAcombinedGenes3 <- c(blueGenes,brownGenes)
greyGenes<-colnames(grey.expr)
RAgreyGenes02<- greyGenes
#Please refer to the article for the GSMID used, or replace it with your own
RA_gset <- getGEO('GSEID', destdir=".",AnnotGPL = T,getGPL = T)
RA_exp<-exprs(RA_gset[[1]])
RA_GPL<-fData(RA_gset[[1]])
str(RA_GPL)
RA_gpl<- RA_GPL[, c(1, 3)]
RA_exp<-as.data.frame(RA_exp)
RA_exp$ID<-rownames(RA_exp)
RA_exp_symbol<-merge(RA_exp,RA_gpl,by="ID")
RA_exp_symbol<-na.omit(RA_exp_symbol)
table(duplicated(RA_exp_symbol$`Gene symbol`))
RA_datExpr02<-avereps(RA_exp_symbol[,-c(1,ncol(RA_exp_symbol))],ID=RA_exp_symbol$`Gene symbol`)
datExpr0<- t(RA_datExpr02)
datExpr1 <- datExpr0
m.vars = apply(datExpr0, 2, var)
var_threshold = quantile(m.vars, 0.40)
expro.upper = datExpr0[, which(m.vars > var_threshold)]
datExpr1 <- data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
 
  abline(h = 250000, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 250000, minSize = 0)
keepSamples = (clust==1)
datExpr = datExpr1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)
datExpr = as.data.frame(datExpr)

datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("1RA03MUSThreshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft$powerEstimate

net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)
pdf("2RA03MUSmodule.pdf",width = 12, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}

RAsamples <- read_excel('RA2.xlsx')

row_names <- RAsamples[[1]]

RAsamples <- RAsamples %>% select(-1)

RAsamples_df <- as.data.frame(RAsamples)

row.names(RAsamples_df) <- row_names

print(RAsamples_df)

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
dim(MEsWW)   
dim(RAsamples_df) 
modTraitCor = cor(MEsWW, RAsamples_df, use = "p")
colnames(MEsWW)
modlues=MEsWW
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3RA03MUSModule-trait-60.pdf",width = 6, height = 8)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(RAsamples_df), yLabels = names(MEsWW), cex.lab = 0.8,  yColorWidth=0.02, 
               xColorWidth = 0.04,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1.0, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()

