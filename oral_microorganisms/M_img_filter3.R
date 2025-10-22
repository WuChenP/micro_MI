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
main_output_dir <- "D:/机器学习与生物信息学/北京朝阳医院项目/R/M/img_filter/family"

# 创建输出目录（如果不存在）
if(!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
  cat("创建输出目录:", main_output_dir, "\n")
}

## 1. 数据预处理和筛选 ------------------------------------------------

# 读取数据
data <- read.csv("D:/机器学习与生物信息学/北京朝阳医院项目/R/M/img_filter/R_data_rename_family.csv")

# 查看数据结构
cat("数据维度:", dim(data), "\n")
cat("前5行前15列:\n")
print(data[1:5, 1:15])

# 筛选科水平的微生物（Species列为空或NA）
family_data <- data[is.na(data$Species) | data$Species == "", ]

# 查看筛选后的科
cat("\n筛选出", nrow(family_data), "个科水平的微生物\n")

# 彻底解决重复行名问题 - 使用name列作为唯一标识符
family_data$unique_id <- family_data$name  # name列是唯一的

# 提取样本数据（AMI和CON开头的列）
sample_cols <- grep("^AMI|^CON", names(family_data), value = TRUE)

# 创建样本分组信息
groups <- ifelse(grepl("^AMI", sample_cols), "AMI", "CON")

# 提取科丰度矩阵（转置，行为样本，列为科）
# 使用unique_id作为列名，确保唯一性
abundance_matrix <- t(family_data[, sample_cols])
colnames(abundance_matrix) <- family_data$unique_id
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
observed_families <- rowSums(abundance_matrix > 0)

# 创建α多样性数据框
alpha_diversity <- data.frame(
  Sample = sample_cols,
  Group = groups,
  Shannon = shannon_diversity,
  Simpson = simpson_diversity,
  Observed_Families = observed_families,
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

p3 <- ggplot(alpha_diversity, aes(x = Group, y = Observed_Families, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  labs(title = "观测科数", x = "组别", y = "科数量") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("AMI" = "#E74C3C", "CON" = "#3498DB")) +
  theme(legend.position = "none")

# 组合图形
alpha_plots <- ggarrange(p1, p2, p3, ncol = 3, labels = "AUTO")
alpha_plots <- annotate_figure(alpha_plots, 
                             top = text_grob("α多样性分析 - 科水平", face = "bold", size = 16))
print(alpha_plots)

# 统计检验
shannon_test <- wilcox.test(Shannon ~ Group, data = alpha_diversity)
simpson_test <- wilcox.test(Simpson ~ Group, data = alpha_diversity)
observed_test <- wilcox.test(Observed_Families ~ Group, data = alpha_diversity)

cat("\nα多样性统计检验结果:\n")
cat("Shannon指数p值:", round(shannon_test$p.value, 4), "\n")
cat("Simpson指数p值:", round(simpson_test$p.value, 4), "\n")
cat("观测科数p值:", round(observed_test$p.value, 4), "\n")

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
    title = "β多样性分析 - PCoA图 (科水平)",
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

## 4. 差异科分析 ---------------------------------

# 计算各组平均丰度
group_means <- data.frame(
  unique_id = family_data$unique_id,
  Family = family_data$Family,
  AMI_mean = rowMeans(family_data[, grep("^AMI", names(family_data))]),
  CON_mean = rowMeans(family_data[, grep("^CON", names(family_data))]),
  stringsAsFactors = FALSE
)

# 计算丰度差异
group_means$log2FC <- log2((group_means$AMI_mean + 1) / (group_means$CON_mean + 1))
group_means$diff <- group_means$AMI_mean - group_means$CON_mean
group_means$abs_log2FC <- abs(group_means$log2FC)

# 为每个科计算p值（Wilcoxon秩和检验）
cat("\n正在计算每个科的p值...\n")
p_values <- numeric(nrow(family_data))

for (i in 1:nrow(family_data)) {
  ami_values <- as.numeric(family_data[i, grep("^AMI", names(family_data))])
  con_values <- as.numeric(family_data[i, grep("^CON", names(family_data))])
  
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

# 筛选差异科 - 同时使用p值和q值
significant_families_p <- group_means[
  group_means$abs_log2FC > 1 & 
  (group_means$AMI_mean > 10 | group_means$CON_mean > 10) &
  group_means$p_value < 0.05,  # 使用p值筛选
]

significant_families_q <- group_means[
  group_means$abs_log2FC > 1 & 
  (group_means$AMI_mean > 10 | group_means$CON_mean > 10) &
  group_means$q_value < 0.05,  # 使用q值筛选
]

# 合并两种筛选结果
significant_families <- rbind(significant_families_p, significant_families_q)
significant_families <- significant_families[!duplicated(significant_families$unique_id), ]

# 添加显著性类型
significant_families$significance_type <- ifelse(
  significant_families$unique_id %in% significant_families_p$unique_id & 
  significant_families$unique_id %in% significant_families_q$unique_id,
  "p<0.05 & q<0.05",
  ifelse(
    significant_families$unique_id %in% significant_families_p$unique_id,
    "p<0.05 only",
    "q<0.05 only"
  )
)

# 按差异大小排序
significant_families <- significant_families[order(significant_families$abs_log2FC, decreasing = TRUE), ]

cat("\n使用p值筛选的显著差异科 (p < 0.05):\n")
print(significant_families_p[, c("Family", "log2FC", "AMI_mean", "CON_mean", "p_value", "q_value")])

cat("\n使用q值筛选的显著差异科 (q < 0.05):\n")
print(significant_families_q[, c("Family", "log2FC", "AMI_mean", "CON_mean", "p_value", "q_value")])

cat("\n合并后的显著差异科:\n")
print(significant_families[, c("Family", "log2FC", "AMI_mean", "CON_mean", "p_value", "q_value", "significance_type")])

## 5. 可视化 ----------------------------------------

# 绘制差异科条形图
if(nrow(significant_families) > 0) {
  # 选择前15个差异最大的科
  top_families <- head(significant_families, 15)
  
  # 准备数据
  bar_data <- top_families %>%
    select(Family, AMI_mean, CON_mean, significance_type) %>%
    pivot_longer(cols = c(AMI_mean, CON_mean), 
                 names_to = "Group", 
                 values_to = "Abundance") %>%
    mutate(Group = ifelse(Group == "AMI_mean", "AMI", "CON"))
  
  # 添加p值和q值信息到科标签
  bar_data <- merge(bar_data, top_families[, c("Family", "p_value", "q_value", "significance_type")], by = c("Family", "significance_type"))
  bar_data$Family_label <- paste0(bar_data$Family, 
                                 "\n(p=", sprintf("%.3f", bar_data$p_value), 
                                 ", q=", sprintf("%.3f", bar_data$q_value), ")")
  
  # 根据显著性类型添加颜色
  bar_plot <- ggplot(bar_data, aes(x = reorder(Family_label, -Abundance), y = Abundance, fill = Group, alpha = significance_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("AMI" = "#E74C3C", "CON" = "#3498DB")) +
    scale_alpha_manual(values = c("p<0.05 & q<0.05" = 1, "p<0.05 only" = 0.7, "q<0.05 only" = 0.7),
                      name = "显著性类型") +
    labs(
      title = "差异科丰度比较",
      subtitle = "Top 15 显著差异科 (使用p值和q值筛选)",
      x = "科",
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
  cat("没有找到符合筛选条件的显著差异科\n")
}

# 绘制火山图 - 同时显示p值和q值
if(nrow(group_means) > 0) {
  volcano_data <- group_means
  
  # 定义显著性类型
  volcano_data$Significance <- ifelse(
    volcano_data$p_value < 0.05 & volcano_data$q_value < 0.05 & volcano_data$abs_log2FC > 1,
    "p<0.05 & q<0.05",
    ifelse(
      volcano_data$p_value < 0.05 & volcano_data$abs_log2FC > 1,
      "p<0.05 only",
      ifelse(
        volcano_data$q_value < 0.05 & volcano_data$abs_log2FC > 1,
        "q<0.05 only",
        "Not Significant"
      )
    )
  )
  
  # 只标记最显著的科
  volcano_data$Label <- ifelse(
    volcano_data$Significance %in% c("p<0.05 & q<0.05", "p<0.05 only", "q<0.05 only") & 
    volcano_data$abs_log2FC > 1.5,
    as.character(volcano_data$Family), 
    ""
  )
  
  # 创建双火山图 - p值版本
  volcano_plot_p <- ggplot(volcano_data, aes(x = log2FC, y = -log10(p_value), 
                                          color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("p<0.05 & q<0.05" = "#E74C3C", 
                                 "p<0.05 only" = "#F39C12",
                                 "q<0.05 only" = "#3498DB",
                                 "Not Significant" = "#95A5A6")) +
    labs(
      title = "差异科火山图 (p值)",
      x = "Log2 Fold Change (AMI/CON)",
      y = "-log10(p-value)",
      color = "显著性"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right") +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 10, box.padding = 0.5)
  
  # 创建双火山图 - q值版本
  volcano_plot_q <- ggplot(volcano_data, aes(x = log2FC, y = -log10(q_value), 
                                          color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("p<0.05 & q<0.05" = "#E74C3C", 
                                 "p<0.05 only" = "#F39C12",
                                 "q<0.05 only" = "#3498DB",
                                 "Not Significant" = "#95A5A6")) +
    labs(
      title = "差异科火山图 (q值)",
      x = "Log2 Fold Change (AMI/CON)",
      y = "-log10(q-value)",
      color = "显著性"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right") +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 10, box.padding = 0.5)
  
  # 组合火山图
  volcano_plots <- ggarrange(volcano_plot_p, volcano_plot_q, ncol = 2, common.legend = TRUE, legend = "bottom")
  volcano_plots <- annotate_figure(volcano_plots, 
                                 top = text_grob("差异科火山图 - 同时使用p值和q值", face = "bold", size = 16))
  print(volcano_plots)
}

## 6. 热图分析 ----------------------------------------

# 选择丰度最高的20个科进行热图分析
if(nrow(family_data) > 0) {
  # 计算总丰度
  total_abundance <- rowSums(family_data[, sample_cols])
  top_families_idx <- order(total_abundance, decreasing = TRUE)[1:min(20, nrow(family_data))]
  
  # 准备热图数据
  heatmap_data <- family_data[top_families_idx, c("Family", sample_cols)]
  
  # 标准化为相对丰度
  heatmap_rel <- as.data.frame(t(apply(heatmap_data[, -1], 1, function(x) x/sum(x) * 100)))
  heatmap_rel$Family <- heatmap_data$Family
  
  # 转换为长格式
  heatmap_long <- heatmap_rel %>%
    pivot_longer(cols = -Family, names_to = "Sample", values_to = "Abundance") %>%
    mutate(Group = ifelse(grepl("^AMI", Sample), "AMI", "CON"))
  
  # 绘制热图
  heatmap_plot <- ggplot(heatmap_long, aes(x = Sample, y = reorder(Family, Abundance), fill = Abundance)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "#E74C3C", 
                        name = "相对丰度(%)") +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    labs(
      title = "Top 20 科相对丰度热图",
      x = "样本",
      y = "科"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      panel.spacing = unit(0.1, "lines"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  print(heatmap_plot)
}

## 7. 保存结果 ------------------------------------------------

# 保存重要结果到指定路径
write.csv(alpha_diversity, file.path(main_output_dir, "alpha_diversity_results.csv"), row.names = FALSE)
write.csv(pcoa_df, file.path(main_output_dir, "pcoa_results.csv"), row.names = FALSE)
write.csv(group_means, file.path(main_output_dir, "family_abundance_means.csv"), row.names = FALSE)
write.csv(significant_families, file.path(main_output_dir, "significant_families.csv"), row.names = FALSE)
write.csv(significant_families_p, file.path(main_output_dir, "significant_families_p.csv"), row.names = FALSE)
write.csv(significant_families_q, file.path(main_output_dir, "significant_families_q.csv"), row.names = FALSE)

# 保存统计检验结果
sink(file.path(main_output_dir, "statistical_tests.txt"))
cat("α多样性统计检验结果:\n")
cat("Shannon指数p值:", round(shannon_test$p.value, 4), "\n")
cat("Simpson指数p值:", round(simpson_test$p.value, 4), "\n")
cat("观测科数p值:", round(observed_test$p.value, 4), "\n\n")
cat("PERMANOVA检验结果:\n")
print(permanova_result)
cat("\n差异科分析:\n")
cat("使用p值筛选标准: |log2FC| > 1, 丰度 > 10, p < 0.05\n")
cat("发现显著差异科数量 (p值):", nrow(significant_families_p), "\n")
cat("使用q值筛选标准: |log2FC| > 1, 丰度 > 10, q < 0.05\n")
cat("发现显著差异科数量 (q值):", nrow(significant_families_q), "\n")
cat("合并后显著差异科数量:", nrow(significant_families), "\n")
sink()

# 保存图形
ggsave(file.path(main_output_dir, "alpha_diversity_plots.png"), alpha_plots, width = 12, height = 5, dpi = 300)
ggsave(file.path(main_output_dir, "beta_diversity_pcoa.png"), pcoa_plot, width = 8, height = 6, dpi = 300)

if(nrow(significant_families) > 0) {
  ggsave(file.path(main_output_dir, "differential_families_barplot.png"), bar_plot, width = 10, height = 8, dpi = 300)
}

if(nrow(group_means) > 0) {
  ggsave(file.path(main_output_dir, "volcano_plots.png"), volcano_plots, width = 16, height = 8, dpi = 300)
  ggsave(file.path(main_output_dir, "heatmap_top_families.png"), heatmap_plot, width = 12, height = 8, dpi = 300)
}

# 生成分析报告
cat("\n=== 分析完成 ===\n")
cat("结果已保存到目录:", main_output_dir, "\n")
cat("- alpha_diversity_results.csv: α多样性结果\n")
cat("- pcoa_results.csv: PCoA坐标\n") 
cat("- family_abundance_means.csv: 科平均丰度（包含p值和q值）\n")
cat("- significant_families.csv: 合并后的显著差异科\n")
cat("- significant_families_p.csv: 使用p值筛选的显著差异科\n")
cat("- significant_families_q.csv: 使用q值筛选的显著差异科\n")
cat("- statistical_tests.txt: 统计检验结果\n")
cat("- 各种图形文件\n")

cat("\n主要发现:\n")
cat("- 筛选出", nrow(family_data), "个科水平的微生物\n")
cat("- AMI组样本数:", sum(groups == "AMI"), "\n")
cat("- CON组样本数:", sum(groups == "CON"), "\n")
cat("- 使用p值发现", nrow(significant_families_p), "个显著差异科 (p < 0.05)\n")
cat("- 使用q值发现", nrow(significant_families_q), "个显著差异科 (q < 0.05)\n")
cat("- 合并后发现", nrow(significant_families), "个显著差异科\n")

if(nrow(significant_families) > 0) {
  cat("\nTop 5 显著差异科:\n")
  print(head(significant_families[, c("Family", "log2FC", "AMI_mean", "CON_mean", "p_value", "q_value", "significance_type")], 5))
}

cat("\n分析完成时间:", Sys.time(), "\n")