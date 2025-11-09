# =============================================================================
# ANCOM-BC2 火山图 - 基于W值百分位数的可靠性分级（图例显示数量）
# =============================================================================

library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggrepel)

# ----------------------------
# 0. 定义输出目录
# ----------------------------
output_dir <- "E:/Python/MI_Analysis/metagenome/Graphic/Covariate_Virus_Differential_Analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1. 准备数据
# ----------------------------
result_file <- "E:/Python/MI_Analysis/metagenome/Graphic/Covariate_Virus_Differential_Analysis/virus_ANCOMBC2_results_with_covariates.xlsx"
data <- read.xlsx(result_file)

# 创建火山图数据
volcano_data <- data.frame(
  taxon = data$taxon,
  log2FC = data$lfc_GroupMI,
  W_stat = abs(data$W_GroupMI),
  q_value = data$q_GroupMI,
  significant = data$diff_GroupMI
)

volcano_data <- na.omit(volcano_data)

# ----------------------------
# 2. 数据处理：基于W值百分位数分级
# ----------------------------

# 计算W值的百分位数阈值
w_threshold_70 <- quantile(volcano_data$W_stat, 0.70)
w_threshold_90 <- quantile(volcano_data$W_stat, 0.90)

# 统计各个类别的数量
high_rel_count <- sum(volcano_data$significant & volcano_data$W_stat > w_threshold_90)
medium_rel_count <- sum(volcano_data$significant & volcano_data$W_stat > w_threshold_70 & volcano_data$W_stat <= w_threshold_90)
low_rel_count <- sum(volcano_data$significant & volcano_data$W_stat <= w_threshold_70)
not_sig_count <- sum(!volcano_data$significant)
total_sig_count <- sum(volcano_data$significant)

# 创建带数量的标签
high_rel_label <- paste0("High reliability (n=", high_rel_count, ")")
medium_rel_label <- paste0("Medium reliability (n=", medium_rel_count, ")")
low_rel_label <- paste0("Low reliability (n=", low_rel_count, ")")
not_sig_label <- "Not significant"

# 基于百分位数定义可靠性，并在标签中包含数量
volcano_data$reliability <- not_sig_label
volcano_data$reliability[volcano_data$significant] <- ifelse(
  volcano_data$W_stat[volcano_data$significant] > w_threshold_90,
  high_rel_label,
  ifelse(
    volcano_data$W_stat[volcano_data$significant] > w_threshold_70,
    medium_rel_label, 
    low_rel_label
  )
)

# 设置颜色
reliability_colors <- c(
  "#E31A1C",  # 高可靠性 - 红色
  "#FF7F00",  # 中等可靠性 - 橙色
  "#1F78B4",  # 低可靠性 - 蓝色
  "grey70"    # 不显著 - 灰色
)
names(reliability_colors) <- c(high_rel_label, medium_rel_label, low_rel_label, not_sig_label)

# ----------------------------
# 3. 绘制火山图：图例显示数量
# ----------------------------

percentile_volcano <- ggplot(volcano_data, aes(x = log2FC, y = W_stat)) +
  
  # 所有点按可靠性分组
  geom_point(aes(color = reliability), alpha = 0.7, size = 2) +
  
  # 添加百分位数参考线
  geom_hline(yintercept = w_threshold_70, linetype = "dashed", color = "#FF7F00", alpha = 0.6) +
  geom_hline(yintercept = w_threshold_90, linetype = "dashed", color = "#E31A1C", alpha = 0.6) +
  
  # 添加参考线标注
  annotate("text", x = min(volcano_data$log2FC), y = w_threshold_70, 
           label = paste("W[0.7] =", round(w_threshold_70, 2)), 
           hjust = 0, vjust = -0.5, size = 3, color = "#FF7F00") +
  annotate("text", x = min(volcano_data$log2FC), y = w_threshold_90, 
           label = paste("W[0.9] =", round(w_threshold_90, 2)), 
           hjust = 0, vjust = -0.5, size = 3, color = "#E31A1C") +
  
  # 图例和标签
  labs(
    title = "ANCOM-BC2 Differential Abundance Analysis",
    subtitle = paste0("MI vs CON | Total significant: ", total_sig_count, " taxa"),
    x = "CLR mean difference",
    y = "W statistic",
    color = "Statistical Reliability"
  ) +
  
  # 颜色设置
  scale_color_manual(values = reliability_colors) +
  
  # 主题设置
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  
  # 参考线
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  
  # Y轴从0开始
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 保存PDF
ggsave(paste0(output_dir, "virus_volcano_with_counts.pdf"), 
       percentile_volcano, width = 10, height = 7)

cat("带数量显示的火山图已保存: virus_volcano_with_counts.pdf\n")

# ----------------------------
# 4. 控制台输出统计摘要
# ----------------------------

cat("\n============================================================\n")
cat("差异物种统计摘要\n")
cat("============================================================\n")
cat("总物种数:", nrow(volcano_data), "\n")
cat("显著差异物种总数:", total_sig_count, "\n")
cat("  高可靠性 (W > 90%):", high_rel_count, "\n")
cat("  中等可靠性 (W > 70%):", medium_rel_count, "\n")
cat("  低可靠性:", low_rel_count, "\n")
cat("不显著物种:", not_sig_count, "\n")
cat("W值阈值 - 第70百分位数:", round(w_threshold_70, 3), "\n")
cat("W值阈值 - 第90百分位数:", round(w_threshold_90, 3), "\n")
cat("============================================================\n")

cat("\n火山图绘制完成！图例中已显示各类别数量。\n")