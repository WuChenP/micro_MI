# 加载必要的包
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(reshape2)

# 设置结果保存路径
results_path <- "D:/PythonProject/micro_MI/metagenome/Graphic/Diversity-Richness_Correlation_Analysis"
if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}
setwd(results_path)

# 读取数据
bacteria_data <- read_excel("D:/PythonProject/micro_MI/metagenome/Absolute_abundance_analysis/filtered_data_1percent/diversity_analysis_OTU/Bacteria_alpha_diversity.xlsx")
virus_data <- read.csv("D:/PythonProject/micro_MI/metagenome/data_figures/filtered_data_1percent_change/diversity_analysis_OTU/alpha_diversity_COMPLETE_virus.csv")
group_data <- read_excel("D:/PythonProject/micro_MI/metagenome/Absolute_abundance_analysis/filtered_data_1percent/sample_metadata.xlsx")

# 数据清洗和预处理
bacteria_clean <- bacteria_data %>%
  select(SampleID, Chao1, Shannon, Simpson) %>%
  rename(Bacterial_Richness = Chao1,
         Bacterial_Diversity = Shannon)

virus_clean <- virus_data %>%
  select(SampleID, Chao1, Shannon, Simpson) %>%
  rename(Viral_Richness = Chao1,
         Viral_Diversity = Shannon)

# 合并数据
merged_data <- group_data %>%
  inner_join(bacteria_clean, by = "SampleID") %>%
  inner_join(virus_clean, by = "SampleID")

# 定义指标顺序（保持一致的顺序）
all_columns <- c("Bacterial_Richness", "Bacterial_Diversity", 
                 "Viral_Richness", "Viral_Diversity")

# 分割数据
con_data <- merged_data %>% filter(Group == "Control")
mi_data <- merged_data %>% filter(Group == "MI")

# 计算完整相关性矩阵的函数（但只保留下三角）
calculate_lower_triangle_correlations <- function(data, columns) {
  n <- length(columns)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)
  pval_matrix <- matrix(NA, nrow = n, ncol = n)
  
  # 计算所有配对的相关性
  for (i in 1:n) {
    for (j in 1:n) {
      if (i > j) {  # 只保留下三角部分
        x <- data[[columns[i]]]
        y <- data[[columns[j]]]
        
        # 移除NA值
        complete_cases <- complete.cases(x, y)
        x_clean <- x[complete_cases]
        y_clean <- y[complete_cases]
        
        if (length(x_clean) > 3) {
          cor_test <- cor.test(x_clean, y_clean, method = "spearman", exact = FALSE)
          cor_matrix[i, j] <- cor_test$estimate
          pval_matrix[i, j] <- cor_test$p.value
        }
      }
    }
  }
  
  rownames(cor_matrix) <- columns
  colnames(cor_matrix) <- columns
  rownames(pval_matrix) <- columns
  colnames(pval_matrix) <- columns
  
  return(list(correlation = cor_matrix, pvalue = pval_matrix))
}

# 计算两组的相关性（下三角）
con_corr <- calculate_lower_triangle_correlations(con_data, all_columns)
mi_corr <- calculate_lower_triangle_correlations(mi_data, all_columns)

# 创建下三角热图函数
create_lower_triangle_heatmap <- function(corr_results, title) {
  cor_matrix <- corr_results$correlation
  pval_matrix <- corr_results$pvalue
  
  # 将矩阵转换为长格式
  cor_melted <- melt(cor_matrix, na.rm = TRUE)
  pval_melted <- melt(pval_matrix, na.rm = TRUE)
  
  # 合并相关性和p值
  plot_data <- cbind(cor_melted, pvalue = pval_melted$value)
  colnames(plot_data) <- c("RowVar", "ColVar", "Correlation", "pvalue")
  
  # 添加显著性星号
  plot_data$significance <- cut(plot_data$pvalue, 
                                breaks = c(0, 0.001, 0.01, 0.05, 1),
                                labels = c("***", "**", "*", ""),
                                include.lowest = TRUE)
  
  # 创建更友好的标签
  variable_labels <- c(
    "Bacterial_Richness" = "Bacterial\nRichness",
    "Bacterial_Diversity" = "Bacterial\nDiversity", 
    "Viral_Richness" = "Viral\nRichness",
    "Viral_Diversity" = "Viral\nDiversity"
  )
  
  plot_data$RowLabel <- variable_labels[plot_data$RowVar]
  plot_data$ColLabel <- variable_labels[plot_data$ColVar]
  
  # 设置因子水平以控制顺序
  row_levels <- c("Bacterial\nRichness", "Bacterial\nDiversity", 
                  "Viral\nRichness", "Viral\nDiversity")
  col_levels <- c("Bacterial\nRichness", "Bacterial\nDiversity", 
                  "Viral\nRichness", "Viral\nDiversity")
  
  plot_data$RowLabel <- factor(plot_data$RowLabel, levels = row_levels)
  plot_data$ColLabel <- factor(plot_data$ColLabel, levels = col_levels)
  
  # 创建热图
  p <- ggplot(plot_data, aes(x = ColLabel, y = RowLabel, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%.2f\n%s", Correlation, significance)), 
              size = 5, color = "black", fontface = "bold") +
    scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Spearman\nCorrelation") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    coord_fixed() +
    ggtitle(title) +
    # 隐藏上三角部分（通过设置不显示文本）
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE)
  
  return(p)
}

# 创建并排下三角热图
p1 <- create_lower_triangle_heatmap(con_corr, "Control Group")
p2 <- create_lower_triangle_heatmap(mi_corr, "MI Group")

# 组合图形
combined_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "right")

# 添加总标题
combined_plot <- annotate_figure(combined_plot, 
                                 top = text_grob("Lower Triangle Correlation Matrix: Bacteria vs Virus Diversity", 
                                                 face = "bold", size = 18))

# 保存ggplot图形
ggsave("Lower_Triangle_Correlations.png", combined_plot, width = 14, height = 7, dpi = 300, bg = "white")
ggsave("Lower_Triangle_Correlations.pdf", combined_plot, width = 14, height = 7, bg = "white")

# 使用corrplot包创建专业的下三角矩阵
png("Lower_Triangle_Correlations_corrplot.png", width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(0, 0, 3, 0))

# 创建友好的标签
friendly_labels <- c("Bac_Rich", "Bac_Div", "Vir_Rich", "Vir_Div")

# Control组下三角相关性图
corrplot(con_corr$correlation, 
         method = "color", 
         type = "lower",
         order = "original", 
         addCoef.col = "black",
         tl.col = "black", 
         tl.srt = 0,
         tl.cex = 1.2,
         tl.pos = "ld",  # 标签在左边和底部
         cl.cex = 0.8,
         col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(200),
         mar = c(0, 0, 2, 0),
         number.cex = 1.0,
         diag = FALSE,  # 不显示对角线
         title = "Control Group")

# MI组下三角相关性图
corrplot(mi_corr$correlation, 
         method = "color", 
         type = "lower",
         order = "original", 
         addCoef.col = "black",
         tl.col = "black", 
         tl.srt = 0,
         tl.cex = 1.2,
         tl.pos = "ld",  # 标签在左边和底部
         cl.cex = 0.8,
         col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(200),
         mar = c(0, 0, 2, 0),
         number.cex = 1.0,
         diag = FALSE,  # 不显示对角线
         title = "MI Group")

# 添加总标题
mtext("Lower Triangle Correlation Matrix: Bacteria vs Virus Diversity", 
      side = 3, outer = TRUE, line = -1.5, cex = 1.5, font = 2)

dev.off()

# 保存结果到CSV文件
write.csv(round(con_corr$correlation, 3), "Control_Group_LowerTriangle_Correlations.csv")
write.csv(round(con_corr$pvalue, 4), "Control_Group_LowerTriangle_Pvalues.csv")
write.csv(round(mi_corr$correlation, 3), "MI_Group_LowerTriangle_Correlations.csv")
write.csv(round(mi_corr$pvalue, 4), "MI_Group_LowerTriangle_Pvalues.csv")

cat("Analysis completed! All results saved to:", results_path, "\n")
cat("Generated files:\n")
cat("- Lower_Triangle_Correlations.png/pdf (ggplot版本)\n")
cat("- Lower_Triangle_Correlations_corrplot.png (corrplot版本)\n")
cat("- Control_Group_LowerTriangle_Correlations.csv\n")
cat("- Control_Group_LowerTriangle_Pvalues.csv\n")
cat("- MI_Group_LowerTriangle_Correlations.csv\n")
cat("- MI_Group_LowerTriangle_Pvalues.csv\n")

# 显示下三角矩阵的内容
cat("\nControl Group Lower Triangle Correlations:\n")
print(round(con_corr$correlation, 3))

cat("\nMI Group Lower Triangle Correlations:\n")
print(round(mi_corr$correlation, 3))