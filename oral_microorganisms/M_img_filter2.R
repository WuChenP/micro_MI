# 加载必要的包
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ape)
library(ggrepel)

# 设置随机种子以保证结果可重复
set.seed(123)

## 0. 设置输出路径 ------------------------------------------------

# 设置主输出目录
main_output_dir <- "D:/机器学习与生物信息学/北京朝阳医院项目/R/M/img_filter/species"

# 创建输出目录（如果不存在）
if(!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
  cat("创建输出目录:", main_output_dir, "\n")
}

## 1. 数据预处理和筛选 ------------------------------------------------

# 读取数据
data <- read.csv("D:/机器学习与生物信息学/北京朝阳医院项目/R/M/img_filter/R_data_rename.csv")

# 查看数据结构
cat("数据维度:", dim(data), "\n")
cat("前5行前15列:\n")
print(data[1:5, 1:15])

# 筛选物种水平为"种"的行（Species列不为空）
species_data <- data[!is.na(data$Species) & data$Species != "", ]

# 查看筛选后的物种
cat("\n筛选出", nrow(species_data), "个物种水平的微生物\n")

# 彻底解决重复行名问题 - 使用name列作为唯一标识符
species_data$unique_id <- species_data$name  # name列是唯一的

# 提取样本数据（AMI和CON开头的列）
sample_cols <- grep("^AMI|^CON", names(species_data), value = TRUE)

# 创建样本分组信息
groups <- ifelse(grepl("^AMI", sample_cols), "AMI", "CON")

# 提取物种丰度矩阵（转置，行为样本，列为物种）
# 使用unique_id作为列名，确保唯一性
abundance_matrix <- t(species_data[, sample_cols])
colnames(abundance_matrix) <- species_data$unique_id
rownames(abundance_matrix) <- sample_cols

# 创建样本信息数据框
sample_info <- data.frame(
  Sample = sample_cols,
  Group = groups,
  stringsAsFactors = FALSE
)

cat("\n样本信息:\n")
print(table(sample_info$Group))

## 2. α多样性分析 ------------------------------------------------

# 计算α多样性指数
shannon_diversity <- diversity(abundance_matrix, index = "shannon")
simpson_diversity <- diversity(abundance_matrix, index = "simpson")
observed_species <- rowSums(abundance_matrix > 0)

# 创建α多样性数据框
alpha_diversity <- data.frame(
  Sample = sample_cols,
  Group = groups,
  Shannon = shannon_diversity,
  Simpson = simpson_diversity,
  Observed_Species = observed_species,
  stringsAsFactors = FALSE
)

# 绘制α多样性箱线图
p1 <- ggplot(alpha_diversity, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  labs(title = "Shannon多样性指数", x = "组别", y = "Shannon指数") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("AMI" = "#E74C3C", "CON" = "#3498DB")) +
  theme(legend.position = "none")

p2 <- ggplot(alpha_diversity, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  labs(title = "Simpson多样性指数", x = "组别", y = "Simpson指数") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("AMI" = "#E74C3C", "CON" = "#3498DB")) +
  theme(legend.position = "none")

p3 <- ggplot(alpha_diversity, aes(x = Group, y = Observed_Species, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  labs(title = "观测物种数", x = "组别", y = "物种数量") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("AMI" = "#E74C3C", "CON" = "#3498DB")) +
  theme(legend.position = "none")

# 组合图形
alpha_plots <- ggarrange(p1, p2, p3, ncol = 3, labels = "AUTO")
alpha_plots <- annotate_figure(alpha_plots, 
                             top = text_grob("α多样性分析", face = "bold", size = 16))
print(alpha_plots)

# 统计检验
shannon_test <- wilcox.test(Shannon ~ Group, data = alpha_diversity)
simpson_test <- wilcox.test(Simpson ~ Group, data = alpha_diversity)
observed_test <- wilcox.test(Observed_Species ~ Group, data = alpha_diversity)

cat("\nα多样性统计检验结果:\n")
cat("Shannon指数p值:", round(shannon_test$p.value, 4), "\n")
cat("Simpson指数p值:", round(simpson_test$p.value, 4), "\n")
cat("观测物种数p值:", round(observed_test$p.value, 4), "\n")

## 3. β多样性分析 ------------------------------------------------

# 计算Bray-Curtis距离
bray_curtis_dist <- vegdist(abundance_matrix, method = "bray")

# PCoA分析
pcoa_result <- cmdscale(bray_curtis_dist, k = 2, eig = TRUE)

# 计算方差解释比例
variance_explained <- round(pcoa_result$eig[1:2] / sum(pcoa_result$eig[pcoa_result$eig > 0]) * 100, 2)

# 创建PCoA绘图数据
pcoa_df <- data.frame(
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  Group = groups,
  Sample = sample_cols
)

# 计算各组中心点
centroids <- pcoa_df %>%
  group_by(Group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# 绘制PCoA图
pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, size = 0.8) +
  geom_point(data = centroids, size = 5, shape = 16) +
  labs(
    title = "β多样性分析 - PCoA图",
    subtitle = paste("Bray-Curtis距离 (AMI n=", sum(groups == "AMI"), 
                    ", CON n=", sum(groups == "CON"), ")", sep = ""),
    x = paste("PC1 (", variance_explained[1], "%)", sep = ""),
    y = paste("PC2 (", variance_explained[2], "%)", sep = ""),
    color = "组别"
  ) +
  scale_color_manual(values = c("AMI" = "#E74C3C", "CON" = "#3498DB")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(pcoa_plot)

# PERMANOVA检验
permanova_result <- adonis2(bray_curtis_dist ~ groups)
cat("\nPERMANOVA检验结果:\n")
print(permanova_result)

## 4. 差异物种分析（简化版，避免行名问题）--------------------------------

# 计算各组平均丰度
group_means <- data.frame(
  unique_id = species_data$unique_id,
  Species = species_data$Species,
  AMI_mean = rowMeans(species_data[, grep("^AMI", names(species_data))]),
  CON_mean = rowMeans(species_data[, grep("^CON", names(species_data))]),
  stringsAsFactors = FALSE
)

# 计算丰度差异
group_means$log2FC <- log2((group_means$AMI_mean + 1) / (group_means$CON_mean + 1))
group_means$diff <- group_means$AMI_mean - group_means$CON_mean
group_means$abs_log2FC <- abs(group_means$log2FC)

# 为每个物种计算p值（Wilcoxon秩和检验）
cat("\n正在计算每个物种的p值...\n")
p_values <- numeric(nrow(species_data))

for (i in 1:nrow(species_data)) {
  ami_values <- as.numeric(species_data[i, grep("^AMI", names(species_data))])
  con_values <- as.numeric(species_data[i, grep("^CON", names(species_data))])
  
  # 只有当两组都有非零值时才进行检验
  if (sum(ami_values) > 0 && sum(con_values) > 0) {
    test_result <- wilcox.test(ami_values, con_values, exact = FALSE)
    p_values[i] <- test_result$p.value
  } else {
    p_values[i] <- NA
  }
}

# 添加p值到结果表
group_means$p_value <- p_values

# 多重检验校正（FDR校正，计算q值）
group_means$q_value <- p.adjust(group_means$p_value, method = "fdr")

# 筛选差异物种
significant_species <- group_means[
  group_means$abs_log2FC > 1 & 
  (group_means$AMI_mean > 10 | group_means$CON_mean > 10) &
  group_means$p_value < 0.05, 
]

# 按差异大小排序
significant_species <- significant_species[order(significant_species$abs_log2FC, decreasing = TRUE), ]

cat("\n显著差异物种 (p < 0.05):\n")
print(significant_species[, c("Species", "log2FC", "AMI_mean", "CON_mean", "p_value", "q_value")])

## 5. 简化可视化（避免行名问题）---------------------------------------

# 绘制差异物种条形图（替代热图）
if(nrow(significant_species) > 0) {
  # 选择前15个差异最大的物种
  top_species <- head(significant_species, 15)
  
  # 准备数据
  bar_data <- top_species %>%
    select(Species, AMI_mean, CON_mean) %>%
    pivot_longer(cols = c(AMI_mean, CON_mean), 
                 names_to = "Group", 
                 values_to = "Abundance") %>%
    mutate(Group = ifelse(Group == "AMI_mean", "AMI", "CON"))
  
  # 添加p值和q值信息到物种标签
  bar_data <- merge(bar_data, top_species[, c("Species", "p_value", "q_value")], by = "Species")
  bar_data$Species_label <- paste0(bar_data$Species, "\n(p=", 
                                  sprintf("%.3f", bar_data$p_value), ")")
  
  bar_plot <- ggplot(bar_data, aes(x = reorder(Species_label, -Abundance), y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("AMI" = "#E74C3C", "CON" = "#3498DB")) +
    labs(
      title = "差异物种丰度比较",
      subtitle = "Top 15 显著差异物种",
      x = "物种",
      y = "平均丰度",
      fill = "组别"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    coord_flip()
  
  print(bar_plot)
  
} else {
  cat("没有找到符合筛选条件的显著差异物种\n")
}

# 绘制火山图（简化版）
if(nrow(group_means) > 0) {
  volcano_data <- group_means
  volcano_data$Significance <- ifelse(
    volcano_data$p_value < 0.05 & volcano_data$abs_log2FC > 1,
    "Significant", "Not Significant"
  )
  
  # 只标记最显著的物种
  volcano_data$Label <- ifelse(
    volcano_data$Significance == "Significant" & 
    volcano_data$p_value < 0.01 & 
    volcano_data$abs_log2FC > 2,
    as.character(volcano_data$Species), 
    ""
  )
  
  volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = -log10(p_value), 
                                          color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant" = "#E74C3C", 
                                 "Not Significant" = "#95A5A6")) +
    labs(
      title = "差异物种火山图",
      x = "Log2 Fold Change (AMI/CON)",
      y = "-log10(p-value)",
      color = "显著性"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right") +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 10, box.padding = 0.5)
  
  print(volcano_plot)
}

## 6. 保存结果 ------------------------------------------------

# 保存重要结果到指定路径
write.csv(alpha_diversity, file.path(main_output_dir, "alpha_diversity_results.csv"), row.names = FALSE)
write.csv(pcoa_df, file.path(main_output_dir, "pcoa_results.csv"), row.names = FALSE)
write.csv(group_means, file.path(main_output_dir, "species_abundance_means.csv"), row.names = FALSE)
write.csv(significant_species, file.path(main_output_dir, "significant_species.csv"), row.names = FALSE)

# 保存统计检验结果
sink(file.path(main_output_dir, "statistical_tests.txt"))
cat("α多样性统计检验结果:\n")
cat("Shannon指数p值:", round(shannon_test$p.value, 4), "\n")
cat("Simpson指数p值:", round(simpson_test$p.value, 4), "\n")
cat("观测物种数p值:", round(observed_test$p.value, 4), "\n\n")
cat("PERMANOVA检验结果:\n")
print(permanova_result)
cat("\n差异物种分析:\n")
cat("使用筛选标准: |log2FC| > 1, 丰度 > 10, p < 0.05\n")
cat("发现显著差异物种数量:", nrow(significant_species), "\n")
sink()

# 保存图形
ggsave(file.path(main_output_dir, "alpha_diversity_plots.png"), alpha_plots, width = 12, height = 5, dpi = 300)
ggsave(file.path(main_output_dir, "beta_diversity_pcoa.png"), pcoa_plot, width = 8, height = 6, dpi = 300)

if(nrow(significant_species) > 0) {
  ggsave(file.path(main_output_dir, "differential_species_barplot.png"), bar_plot, width = 10, height = 8, dpi = 300)
  ggsave(file.path(main_output_dir, "volcano_plot.png"), volcano_plot, width = 8, height = 6, dpi = 300)
}

# 生成分析报告
cat("\n=== 分析完成 ===\n")
cat("结果已保存到目录:", main_output_dir, "\n")
cat("- alpha_diversity_results.csv: α多样性结果\n")
cat("- pcoa_results.csv: PCoA坐标\n") 
cat("- species_abundance_means.csv: 物种平均丰度（包含p值和q值）\n")
cat("- significant_species.csv: 显著差异物种\n")
cat("- statistical_tests.txt: 统计检验结果\n")
cat("- 各种图形文件\n")

cat("\n主要发现:\n")
cat("- 筛选出", nrow(species_data), "个物种水平的微生物\n")
cat("- AMI组样本数:", sum(groups == "AMI"), "\n")
cat("- CON组样本数:", sum(groups == "CON"), "\n")
cat("- 发现", nrow(significant_species), "个显著差异物种 (p < 0.05)\n")

if(nrow(significant_species) > 0) {
  cat("\nTop 5 显著差异物种:\n")
  print(head(significant_species[, c("Species", "log2FC", "AMI_mean", "CON_mean", "p_value", "q_value")], 5))
}

cat("\n分析完成时间:", Sys.time(), "\n")