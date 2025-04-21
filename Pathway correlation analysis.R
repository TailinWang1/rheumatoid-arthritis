library(tidyverse)
library(ggplot2)
library(reshape2)
library(cowplot)
library(limma)
library(Biobase)
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

RA_groups <- c(rep("RA", 12), rep("HC", 3))#mus

RA_groups <- c(rep("HC", 5), rep("RA", 9))#human

genes_of_interest <- c("MIF", "CD74", "CXCR4", "CD44","HLA-DRA","HLA-DRB4")
genes_of_interest <- c("Mif", "Cd74", "Cxcr4", "Cd44","H2-Ab1","H2-Aa")
RA_genes <- RA_datExpr02[rownames(RA_datExpr02) %in% genes_of_interest, ]

perform_differential_expression <- function(expression_data, groups) {
  design <- model.matrix(~factor(groups))
  fit <- lmFit(expression_data, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef=2, number=Inf)
  results$Gene <- rownames(results)
  return(results)
}
)

RA_de <- perform_differential_expression(RA_genes, RA_groups)
RA_de$significance <- ifelse(RA_de$adj.P.Val < 0.05, "*", "")
RA_de$significance[RA_de$adj.P.Val < 0.01] <- "**"
RA_de$significance[RA_de$adj.P.Val < 0.001] <- "***"

RA_genes_long <- RA_genes %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Group = rep(RA_groups, each = length(genes_of_interest)))


p1 <- ggplot(RA_genes_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = RA_de, aes(x = "RA", y = max(RA_genes_long$Expression), label = significance), 
            vjust = -0.5, color = "black", size = 5, inherit.aes = FALSE) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Expression of MIF Pathway Genes in HC and RA",
       x = "Group", y = "Expression")

print(p1)

calculate_pathway_score <- function(gene_data) {
  colMeans(gene_data)
}

RA_score <- calculate_pathway_score(RA_genes)

pathway_scores <- data.frame(
  Score = RA_score,
  Disease = factor(RA_groups),
  Dataset = factor(rep("RA", length(RA_score)))
)


p2 <- ggplot(pathway_scores, aes(x = Disease, y = Score, fill = Disease)) +
  geom_violin(alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 2) +
  scale_fill_manual(values = c("HC" = "#2f5c85", "RA" = "#c39797")) +
  stat_compare_means(
    aes(group = Disease),
    label.x = 1.5,
    label.y = max(pathway_scores$Score) * 1.1,
    method = "t.test",
    label = "p.signif"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "MIF Pathway Score",
    x = "",
    y = "Pathway Score"
  )
p2

p3 <- ggplot(pathway_scores, aes(x = Score, fill = Disease)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("HC" = "#2f5c85", "RA" = "#c39797")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Distribution of Pathway Scores",
    x = "Score",
    y = "Density"
  )
p3

correlation_matrix <- cor(t(RA_genes))
p_matrix <- matrix(NA, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
rownames(p_matrix) <- rownames(correlation_matrix)
colnames(p_matrix) <- colnames(correlation_matrix)


for(i in 1:nrow(correlation_matrix)) {
  for(j in 1:nrow(correlation_matrix)) {
    test <- cor.test(RA_genes[i,], RA_genes[j,], method = "pearson")
    p_matrix[i,j] <- test$p.value
  }
}


get_sig_stars <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value <= 0.001) return("***")
  if (p_value <= 0.01) return("**")
  if (p_value <= 0.05) return("*")
  return("")
}


melted_cormat <- reshape2::melt(correlation_matrix)
melted_p <- reshape2::melt(p_matrix)
melted_cormat$stars <- sapply(melted_p$value, get_sig_stars)


p4 <- ggplot(melted_cormat, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2f5c85",
    high = "#c39797",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)
  ) +
  geom_text(aes(label = sprintf("%.2f%s", value, stars)), 
            size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  ) +
  labs(title = "MIF Pathway Genes Correlations", fill = "Correlation")

library(patchwork)
combined_plot <- (p2 + p3) / p4 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "MIF Pathway Analysis in RA vs HC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

print(combined_plot)


