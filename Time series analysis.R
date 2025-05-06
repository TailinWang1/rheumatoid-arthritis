library(readxl)
library(dplyr)
library(limma)
library(GEOquery)
library(impute)
library(readr)
library(preprocessCore)
library(ggplot2)
library(dplyr)
library(ggrepel)
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

R_data <- RA_datExpr02["基因", ]
R_data <- RA_datExpr02["基因", ]
group_0_3_days <- c(1,2,3)
group_1_2_weeks <- c(4,5,6)
group_3_4_weeks <- c(7,8,9)
group_remission_over_3_weeks <- c(10,11,12)
group_control <- c(13,14,15)

calculate_stats <- function(data) {
  c(mean = mean(data), 
    sd = sd(data), 
    min = min(data), 
    max = max(data))
}

stats_0_3_days <- calculate_stats(R_data[group_0_3_days])
stats_1_2_weeks <- calculate_stats(R_data[group_1_2_weeks])
stats_3_4_weeks <- calculate_stats(R_data[group_3_4_weeks])
stats_remission_over_3_weeks <- calculate_stats(R_data[group_remission_over_3_weeks])
stats_control <- calculate_stats(R_data[group_control])


results <- data.frame(
  Group = c("0-3 days", "1-2 weeks", "3-4 weeks", "Remission over 3 weeks", "Control"),
  rbind(stats_0_3_days, stats_1_2_weeks, stats_3_4_weeks, stats_remission_over_3_weeks, stats_control)
)

print("Statistical summary for R expression:")
print(results)




perform_t_test <- function(group1, group2, group1_name, group2_name) {
  t_result <- t.test(R_data[group1], R_data[group2])
  return(data.frame(
    Comparison = paste(group1_name, "vs", group2_name),
    t_statistic = t_result$statistic,
    p_value = t_result$p.value,
    mean_difference = diff(t_result$estimate)
  ))
}

t_test_results_control <- rbind(
  perform_t_test(group_0_3_days, group_control, "0-3 days", "Control"),
  perform_t_test(group_1_2_weeks, group_control, "1-2 weeks", "Control"),
  perform_t_test(group_3_4_weeks, group_control, "3-4 weeks", "Control"),
  perform_t_test(group_remission_over_3_weeks, group_control, "Remission over 3 weeks", "Control")
)

t_test_results_disease <- rbind(
  perform_t_test(group_0_3_days, group_1_2_weeks, "0-3 days", "1-2 weeks"),
  perform_t_test(group_0_3_days, group_3_4_weeks, "0-3 days", "3-4 weeks"),
  perform_t_test(group_0_3_days, group_remission_over_3_weeks, "0-3 days", "Remission over 3 weeks"),
  perform_t_test(group_1_2_weeks, group_3_4_weeks, "1-2 weeks", "3-4 weeks"),
  perform_t_test(group_1_2_weeks, group_remission_over_3_weeks, "1-2 weeks", "Remission over 3 weeks"),
  perform_t_test(group_3_4_weeks, group_remission_over_3_weeks, "3-4 weeks", "Remission over 3 weeks")
)

all_t_test_results <- rbind(t_test_results_control, t_test_results_disease)

print("T-test results for all comparisons:")
print(all_t_test_results)

library(ggplot2)
library(ggrepel)


prepare_data <- function(gene_data, gene_name) {
  data.frame(
    Time = factor(rep(c("Control", "0-3 days", "1-2 weeks", "3-4 weeks", "Remission >3 weeks"), each = 3),
                  levels = c("Control", "0-3 days", "1-2 weeks", "3-4 weeks", "Remission >3 weeks")),
    Expression = c(gene_data[group_control],
                   gene_data[group_0_3_days], 
                   gene_data[group_1_2_weeks], 
                   gene_data[group_3_4_weeks], 
                   gene_data[group_remission_over_3_weeks]),
    Gene = gene_name
  )
}

plot_data <- rbind(
  prepare_data(R_data["R",], "R"),
  prepare_data(R_data["H",], "H")
)

mean_data <- aggregate(Expression ~ Time + Gene, plot_data, mean)

d
muted_morandi_colors <- c("#D4A5A5", "#9CAF88", "#A5C0C5", "#B0A5BA", "#E2B49A")
control_level_R <- mean(R_data["R", group_control])
control_level_H <- mean(R_data["H", group_control])

p <- ggplot(plot_data, aes(x = Time, y = Expression, color = Time, shape = Gene, group = Gene)) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_hline(yintercept = control_level_R, linetype = "dashed", color = "#8B0000", size = 1) +
  geom_hline(yintercept = control_level_H, linetype = "dashed", color = "#8B0000", size = 1) +
  geom_text_repel(data = mean_data, aes(label = sprintf("%.2f", Expression)), 
                  size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = muted_morandi_colors) +
  scale_shape_manual(values = c(16, 17)) +
  theme_minimal(base_size = 14) +
  labs(title = "Gene Expression Profile in CIA Mouse Model",
       subtitle = "Temporal changes in R and H expression",
       x = "Disease Stage",
       y = "Expression Level",
       shape = "Gene",
       color = "Time Point") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  ) +
  coord_cartesian(ylim = c(min(plot_data$Expression) * 0.9, max(plot_data$Expression) * 1.1))

print(p)

