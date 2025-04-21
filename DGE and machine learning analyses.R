library(readxl)
library(dplyr)
library(limma)
library(GEOquery)
library(impute)
library(WGCNA)
library(readr)
library(preprocessCore)
RA_top_table<- read.delim("GSE56649.top.table.tsv", header = TRUE, sep = "\t")
ra_genes_base <- RA_top_table[which(RA_top_table$adj.P.Val< 0.05),]
ra_genes<- ra_genes_base$Gene.symbol
ra_genes_DEG <- unique(ra_genes)
ra_genes_DEG_cleaned <- ra_genes_DEG[ra_genes_DEG != ""]
ROS_genes <- c("MIR675", "HDAC6", "IL18BP", "PPIF", "TRAP1", "DHRS2", "NET1", "KLF2", "RACK1", "STK25", "FBLN5", "ARL6IP5", "CAMKK2", "PPARGC1A", "PRDX3", "RIPK3", "ZNF277", "PDCD10", "TREX1", "PARK7", "CHUK", "AGAP3", "LRRK2", "CPEB2", "PPARGC1B", "NCOA7", "ATF2", "SRXN1", "ROMO1", "SIRPA", "PARP1", "CRYGD", "CYP1B1", "VKORC1L1", "DAPK1", "DHFR", "NQO1", "ECT2", "EDN1", "EDNRA", "EGFR", "EIF2S1", "AIF1", "EPAS1", "STX2", "AKT1", "ERN1", "ETV5", "EZH2", "ABCD1", "FABP1", "FANCC", "FANCD2", "ALDH3B1", "FER", "PLA2R1", "SIRT2", "SETX", "FOXO1", "FOXO3", "KDM6B", "ATP13A2", "SIRT1", "FOS", "SLC7A11", "FXN", "ALOX5", "ABL1", "FUT8", "FYN", "G6PD", "MPV17L", "PRDX5", "SIN3A", "ANKRD2", "GCH1", "GJB2", "STAU2", "FOXP1", "HTRA2", "H19", "GPR37", "GPX1", "GPX5", "GPX7", "GSR", "PYCR2", "SLC25A24", "ANXA1", "HDAC2", "HGF", "HIF1A", "HMOX1", "NUDT2", "HSF1", "HSPA1A", "HSPA1B", "APOA4", "IL6", "AQP1", "JUN", "ERCC6L2", "SUMO4", "RHOB", "TMIGD1", "ARNT", "BMAL1", "MIRLET7B", "MIR103A1", "MIR107", "MIR132", "MIR133A1", "MIR135A1", "MIR17", "MIR21", "MIR34A", "MIR92A1", "MAPT", "MDM2", "MAP3K5", "MET", "MGAT3", "MGST1", "MMP2", "MMP3", "MMP9", "MPO", "MPV17", "ABCC1", "MSRA", "MT3", "MYB", "NAGLU", "ATF4", "ATM", "NFE2L2", "NOS3", "ATP2A2", "DDR2", "NR4A2", "GPX8", "PJVK", "PRDX1", "PNPLA8", "PRKN", "PAWR", "PAX2", "GLRX2", "PCNA", "CHCHD2", "ZNF580", "NME8", "WNT16", "RWDD1", "OSER1", "PDGFRA", "PDK2", "PENK", "PEX10", "PEX12", "PEX13", "PEX14", "PKD2", "ATP7A", "RBM11", "PPIA", "TMEM161A", "ADPRS", "OXR1", "ANKZF1", "BRF2", "SMPD3", "PRKAA1", "PRKAA2", "AXL", "PRKCD", "SELENOS", "PRKD1", "MAPK1", "MAPK3", "MAPK7", "MAPK8", "MAPK9", "MAPK13", "SELENON", "CBX8", "DHFRP1", "TBC1D24", "MEAK7", "PTPRK", "PEX2", "PXN", "PEX5", "PYCR1", "RAD52", "PLEKHA1", "RELA", "GAS5", "RPS3", "MAP2K4", "PINK1", "SLC1A1", "NCF1", "BMP7", "SNCA", "BNIP3", "SOD1", "SOD2", "SOD3", "SRC", "STAT6", "STAU1", "STX4", "BTK", "PRDX2", "TNFAIP3", "TOP2B", "TP53", "TPM1", "TRPM2", "TXN", "UCP1", "VRK2", "PCGF2", "MAPKAP1", "ZFAND1", "PRR5L", "PYROXD1", "ERMP1", "ZC3H12A", "PDGFD", "HM13", "SESN2", "SLC4A11", "MAP1LC3A", "CAT", "MPV17L2", "AIFM2", "PRKRA", "BECN1", "PNPT1", "RIPK1", "IL18RAP", "KAT2B", "SPHK1", "SQSTM1", "TRPA1", "SLC25A14", "AIFM1", "GPR37L1", "LONP1", "TP53INP1", "CD36", "KEAP1", "CDK1", "CCS")
gene_names <- readLines("GOCC_MITOCHONDRION.v2024.1.Hs.grp")
head(gene_names, 10)
gene_names <- gene_names[-c(1,2)]
Mitochondrial_genes <- gene_names
head(Mitochondrial_genes)
A<-Reduce(intersect, list(ROS_genes,Mitochondrial_genes,ra_genes_DEG_cleaned))
A
#RA-LASSO
library(glmnet)
library(survival)
library(tidyverse)
RA_gset <- getGEO('GSE56649', destdir=".",AnnotGPL = T,getGPL = T)
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
selected_genes_data <- datExpr0[, A]
ch <- c(rep(1,13),rep(0,9)) %>% as.data.frame() %>%
  set_names('ch')
set.seed(1)
x <- as.matrix(selected_genes_data)
y <- model.matrix(~ch-1,ch) %>% data.matrix() %>% as.factor()
dim(x) 
length(y) 
cvla <- glmnet(x,y,family = 'binomial')
cv.fit <- cv.glmnet(x,y,family='binomial')
par(mar=c(7, 4, 4, 2) + 0.1)
pdf("RA_GSE56649_LASSO_cvla_plot.pdf")
plot(cvla, xvar = 'lambda', label = TRUE)
dev.off()
pdf("RA_GSE56649_LASSO_cv_fit_plot.pdf")
plot(cv.fit)
dev.off()
cv.fit$lambda.min
cv.fit$lambda.1se
total=cv.fit$lambda.min+cv.fit$lambda.1se
coef=coef(cvla,s=cv.fit$lambda.1se)
index=which(coef!=0)
actcoef=coef[index]
LassoGene=row.names(coef)[index]
RA_genecoef=cbind(Gene=LassoGene,coef=actcoef)
RA_genecoef
library(tidyverse)
library(glmnet)
library(VennDiagram)
library(e1071)
library(caret)
library(randomForest)
set.seed(21)
source('msvmRFE.R')
expression_data <- datExpr0[, A, drop = FALSE] 
train <- read.csv("RA_GSE56649_23genes_data.csv",row.names = 1,
                  as.is = F)
input <- train
svmRFE(input, k = 5, halve.above = 100) 
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) 
top.features = WriteFeatures(results, input, save=F) 
head(top.features)
write.csv(top.features,"RA_GSE56649_feature_svm.csv")
featsweep = lapply(1:23, FeatSweep.wrap, results, input)
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
pdf("RA_GSE56649_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()
pdf("RA_GSE56649_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info)
dev.off()
which.min(errors)
top.features[1:which.min(errors), "FeatureName"]
top_genes <- top.features[1:which.min(errors), "FeatureName"]
write.csv(top_genes, file = "RA_GSE56649_SVM_top_genes.csv", row.names = FALSE)
a1<-Reduce(intersect, list(RA_genecoef,top_genes))
a1

library(randomForest)
library(caret)
set.seed(123)
expression_data <- datExpr0[, a1, drop = FALSE] 
ch <- c(rep(1,13), rep(0,9))
ch <- as.factor(ch) 
rf_model <- randomForest(x = expression_data, y = ch, importance = TRUE, ntree = 500)

importance_df <- importance(rf_model)
importance_df <- as.data.frame(importance_df)
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
top_features_rf <- head(importance_df, 10)
importance_df
top_features_rf
pdf("RA_GSE56649_RF_feature_importance.pdf", width = 10, height = 8)
varImpPlot(rf_model, n.var = min(10, nrow(importance_df)), 
           main = "Top 10 Important Features")
dev.off()
