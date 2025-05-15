library(readxl)
library(dplyr)
library(limma)
library(GEOquery)
library(impute)
library(WGCNA)
library(readr)
library(preprocessCore)
#Please refer to the article for the GSMID used, or replace it with your own
RA_top_table<- read.delim("GSEID.top.table.tsv", header = TRUE, sep = "\t")
ra_genes_base <- RA_top_table[which(RA_top_table$adj.P.Val< 0.05),]
ra_genes<- ra_genes_base$Gene.symbol
ra_genes_DEG <- unique(ra_genes)
ra_genes_DEG_cleaned <- ra_genes_DEG[ra_genes_DEG != ""]
ROS_genes <-  read.csv("ROS.csv") 
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
selected_genes_data <- datExpr0[, A]
ch <- c(rep(1,NA),rep(0,NA)) %>% as.data.frame() %>%
  set_names('ch')
set.seed(1)
x <- as.matrix(selected_genes_data)
y <- model.matrix(~ch-1,ch) %>% data.matrix() %>% as.factor()
dim(x) 
length(y) 
cvla <- glmnet(x,y,family = 'binomial')
cv.fit <- cv.glmnet(x,y,family='binomial')
par(mar=c(7, 4, 4, 2) + 0.1)
pdf("RA_LASSO_cvla_plot.pdf")
plot(cvla, xvar = 'lambda', label = TRUE)
dev.off()
pdf("RA_LASSO_cv_fit_plot.pdf")
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
train <- read.csv("RA_23genes_data.csv",row.names = 1,
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
write.csv(top.features,"RA_feature_svm.csv")
featsweep = lapply(1:23, FeatSweep.wrap, results, input)
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
pdf("RA_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()
pdf("RA_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info)
dev.off()
which.min(errors)
top.features[1:which.min(errors), "FeatureName"]
top_genes <- top.features[1:which.min(errors), "FeatureName"]
write.csv(top_genes, file = "RA_SVM_top_genes.csv", row.names = FALSE)
a1<-Reduce(intersect, list(RA_genecoef,top_genes))
a1

library(randomForest)
library(caret)
set.seed(123)
expression_data <- datExpr0[, a1, drop = FALSE] 
ch <- c(rep(1,NA), rep(0,NA))
ch <- as.factor(ch) 
rf_model <- randomForest(x = expression_data, y = ch, importance = TRUE, ntree = 500)

importance_df <- importance(rf_model)
importance_df <- as.data.frame(importance_df)
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
top_features_rf <- head(importance_df, 10)
importance_df
top_features_rf
pdf("RA_RF_feature_importance.pdf", width = 10, height = 8)
varImpPlot(rf_model, n.var = min(10, nrow(importance_df)), 
           main = "Top 10 Important Features")
dev.off()
