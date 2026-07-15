library(GEOquery)
library(pROC)
library(ggplot2)

# Download and process GSE datasets (Note the [[1]] added to extract the ExpressionSet from the list)
gse56649 <- getGEO('GSE56649', destdir=".", AnnotGPL = TRUE, getGPL = TRUE)[[1]]
gse15573 <- getGEO('GSE15573', destdir=".", AnnotGPL = TRUE, getGPL = TRUE)[[1]]
gse205962 <- getGEO('GSE205962', destdir=".", AnnotGPL = TRUE, getGPL = TRUE)[[1]]

# Define sample groups
groups_56649 <- factor(c(rep("RA", 13), rep("HC", 9)))
groups_15573 <- factor(c("RA", "RA", "Control", "RA", "Control", "RA", "Control", "Control", "RA", "Control", "Control", "RA", "Control", "RA", "RA", "Control", "RA", "Control", "RA", "RA", "Control", "RA", "Control", "RA", "RA", "Control", "Control", "RA", "RA", "Control", "RA", "Control", "RA"))
groups_205962 <- factor(c(rep("RA", 16), rep("HC", 4)))

# Function to extract gene expression data
extract_gene_expr <- function(gse, groups, gene_symbol) {
  expr <- exprs(gse)
  feature_data <- fData(gse)
  
  if ("Gene symbol" %in% colnames(feature_data)) {
    col_name <- "Gene symbol"
  } else if ("Gene Symbol" %in% colnames(feature_data)) {
    col_name <- "Gene Symbol"
  } else {
    stop("Cannot find a column for gene symbols")
  }
  
  gene_probe <- grep(gene_symbol, feature_data[[col_name]], ignore.case = TRUE)
  if (length(gene_probe) == 0) {
    stop(paste(gene_symbol, "probe not found in the dataset"))
  }
  
  gene_expr <- expr[gene_probe[1], ]
  batch <- rep(deparse(substitute(gse)), length(gene_expr))
  list(expr = gene_expr, group = groups, batch = batch)
}

# Extract expression data for ROMO1, CRP, and ADA
data56649_romo1 <- extract_gene_expr(gse56649, groups_56649, "ROMO1")
data15573_romo1 <- extract_gene_expr(gse15573, groups_15573, "ROMO1")
data205962_romo1 <- extract_gene_expr(gse205962, groups_205962, "ROMO1")

data56649_crp <- extract_gene_expr(gse56649, groups_56649, "CRP")
data15573_crp <- extract_gene_expr(gse15573, groups_15573, "CRP")
data205962_crp <- extract_gene_expr(gse205962, groups_205962, "CRP")

data56649_ada <- extract_gene_expr(gse56649, groups_56649, "ADA")
data15573_ada <- extract_gene_expr(gse15573, groups_15573, "ADA")
data205962_ada <- extract_gene_expr(gse205962, groups_205962, "ADA")

# Function to calculate ROC curve and AUC
calc_roc <- function(expr, group) {
  roc_curve <- roc(group, expr)
  auc_value <- auc(roc_curve)
  list(roc = roc_curve, auc = auc_value)
}

# Calculate ROC for ROMO1
roc_train_romo1 <- calc_roc(data56649_romo1$expr, data56649_romo1$group)
roc_valid1_romo1 <- calc_roc(data15573_romo1$expr, data15573_romo1$group)
roc_valid2_romo1 <- calc_roc(data205962_romo1$expr, data205962_romo1$group)

# Calculate ROC for CRP
roc_train_crp <- calc_roc(data56649_crp$expr, data56649_crp$group)
roc_valid1_crp <- calc_roc(data15573_crp$expr, data15573_crp$group)
roc_valid2_crp <- calc_roc(data205962_crp$expr, data205962_crp$group)

# Calculate ROC for ADA
roc_train_ada <- calc_roc(data56649_ada$expr, data56649_ada$group)
roc_valid1_ada <- calc_roc(data15573_ada$expr, data15573_ada$group)
roc_valid2_ada <- calc_roc(data205962_ada$expr, data205962_ada$group)

# Helper function to plot individual gene ROC curves
plot_single_gene_roc <- function(roc_train, roc_valid1, roc_valid2, gene_name) {
  
  # Prepare data with tissue annotations
  plot_data <- rbind(
    data.frame(specificity = roc_train$roc$specificities, sensitivity = roc_train$roc$sensitivities, 
               group = "Training: GSE56649 (PB CD4 T cells)"),
    data.frame(specificity = roc_valid1$roc$specificities, sensitivity = roc_valid1$roc$sensitivities, 
               group = "Validation 1: GSE15573 (PBMCs)"),
    data.frame(specificity = roc_valid2$roc$specificities, sensitivity = roc_valid2$roc$sensitivities, 
               group = "Validation 2: GSE205962 (Whole blood)")
  )
  
  # Define Morandi (low saturation) color palette
  morandi_colors <- c(
    "Training: GSE56649 (PB CD4 T cells)" = "#A3B18A",  # Sage Green
    "Validation 1: GSE15573 (PBMCs)"      = "#C5AFA4",  # Dusty Rose
    "Validation 2: GSE205962 (Whole blood)" = "#9EA9B1" # Dusty Blue
  )
  
  # Format AUC text
  auc_text <- sprintf(
    "Training AUC: %.3f\nValidation 1 AUC: %.3f\nValidation 2 AUC: %.3f",
    roc_train$auc, roc_valid1$auc, roc_valid2$auc
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = 1 - specificity, y = sensitivity, color = group)) +
    geom_line(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
    scale_color_manual(values = morandi_colors) +
    theme_minimal() +
    labs(
      title = paste("ROC Curves for", gene_name, "in RA Diagnosis"),
      x = "1 - Specificity",
      y = "Sensitivity",
      color = NULL # 移除图例标题让横向排版更紧凑
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal", 
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    ) +
    
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    annotate("text", x = 0.60, y = 0.15, label = auc_text, hjust = 0, vjust = 0, size = 4, color = "#4A4A4A")
  
  return(p)
}

# Generate and print the three separate plots
p_romo1 <- plot_single_gene_roc(roc_train_romo1, roc_valid1_romo1, roc_valid2_romo1, "ROMO1")
p_crp   <- plot_single_gene_roc(roc_train_crp, roc_valid1_crp, roc_valid2_crp, "CRP")
p_ada   <- plot_single_gene_roc(roc_train_ada, roc_valid1_ada, roc_valid2_ada, "ADA")

print(p_romo1)
print(p_crp)
print(p_ada)