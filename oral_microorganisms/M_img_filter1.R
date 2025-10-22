# 加载必要包
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan")
}

library(tidyverse)
library(vegan)

# 设置输出路径
output_path <- "D:/机器学习与生物信息学/北京朝阳医院项目/R/M/"
# 如果输出路径不存在，创建它
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# 读取数据
data <- read.csv("D:/机器学习与生物信息学/北京朝阳医院项目/R/M/R_data_rename.csv", header = TRUE)

# 提取样本列
sample_cols <- grep("^(AMI|CON)", colnames(data), value = TRUE)
ami_samples <- grep("^AMI", sample_cols, value = TRUE)
con_samples <- grep("^CON", sample_cols, value = TRUE)

# 创建丰度矩阵
abundance_matrix <- data[, sample_cols]
rownames(abundance_matrix) <- data$name

# 计算各种指标的函数
calculate_metrics <- function(abundance_df) {
  metrics <- data.frame(
    ASV = rownames(abundance_df),
    # 基本丰度指标
    Total_Abundance = rowSums(abundance_df),
    Mean_Abundance = rowMeans(abundance_df),
    Prevalence = apply(abundance_df > 0, 1, mean) * 100,
    
    # 组特异性丰度
    Mean_AMI = rowMeans(abundance_df[, ami_samples]),
    Mean_CON = rowMeans(abundance_df[, con_samples]),    
    # 变异指标
    CV = apply(abundance_df, 1, function(x) sd(x)/mean(x)),
    Max_Abundance = apply(abundance_df, 1, max),
    
    # 存在/缺失模式
    Zero_Count = apply(abundance_df == 0, 1, sum),
    Presence_Count = apply(abundance_df > 0, 1, sum),
    
    stringsAsFactors = FALSE
  )
  return(metrics)
}

# 计算指标
basic_metrics <- calculate_metrics(abundance_matrix)

# 添加分类信息
taxonomy_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy_info <- data[, c("name", taxonomy_cols)]
colnames(taxonomy_info)[1] <- "ASV"

# 合并指标和分类信息
full_metrics <- basic_metrics %>%
  left_join(taxonomy_info, by = "ASV") %>%
  mutate(
    # 计算丰度比率
    AMI_CON_Ratio = ifelse(Mean_CON > 0, Mean_AMI/Mean_CON, NA),
    Log2_FoldChange = log2(ifelse(Mean_CON > 0, (Mean_AMI + 1)/(Mean_CON + 1), NA)),
    
    # 丰度排名
    Abundance_Rank = rank(-Total_Abundance),
    Prevalence_Rank = rank(-Prevalence),
    
    # 分类分组
    Dominance_Level = case_when(
      Mean_Abundance > 100 ~ "High",
      Mean_Abundance > 10 ~ "Medium", 
      Mean_Abundance > 1 ~ "Low",
      TRUE ~ "Rare"
    )
  )

# 使用非参数检验计算组间差异
calculate_differential_abundance <- function(abundance_df) {
  results <- data.frame(ASV = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  for(i in 1:nrow(abundance_df)) {
    ami_values <- as.numeric(abundance_df[i, ami_samples])
    con_values <- as.numeric(abundance_df[i, con_samples])
    
    # 只对在至少一组中有读数的ASV进行检验
    if(sum(ami_values) > 0 | sum(con_values) > 0) {
      test_result <- wilcox.test(ami_values, con_values, exact = FALSE)
      results <- rbind(results, data.frame(
        ASV = rownames(abundance_df)[i],
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  # FDR校正
  results$fdr <- p.adjust(results$p_value, method = "fdr")
  return(results)
}

# 执行差异分析
diff_results <- calculate_differential_abundance(abundance_matrix)

# 创建完整的指标表格
complete_metrics_table <- full_metrics %>%
  left_join(diff_results, by = "ASV") %>%
  mutate(
    Significant = ifelse(fdr < 0.05, "Yes", "No"),
    Change_Direction = ifelse(Log2_FoldChange > 0, "Up in AMI", 
                            ifelse(Log2_FoldChange < 0, "Up in CON", "Equal"))
  ) %>%
  select(
    # 标识信息
    ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species,
    
    # 丰度指标
    Total_Abundance, Mean_Abundance, 
    Mean_AMI, Mean_CON, AMI_CON_Ratio, Log2_FoldChange,
    Max_Abundance,
    
    # 分布指标
    Prevalence, Presence_Count, Zero_Count, CV,
    
    # 统计指标
    p_value, fdr, Significant, Change_Direction,
    
    # 分类指标
    Dominance_Level, Abundance_Rank, Prevalence_Rank
  ) %>%
  arrange(desc(Total_Abundance))

# 保存完整表格到指定路径
complete_table_path <- paste0(output_path, "Comprehensive_Microbial_Metrics_Table.csv")
write.csv(complete_metrics_table, complete_table_path, row.names = FALSE)

# 显示前20个ASV的详细指标
cat("=== 微生物多样性指标综合表格 (前20个ASV) ===\n")
top_20_metrics <- complete_metrics_table %>%
  head(20) %>%
  select(ASV, Phylum, Genus, Species, 
         Total_Abundance, Mean_Abundance, Prevalence,
         Mean_AMI, Mean_CON, AMI_CON_Ratio, Log2_FoldChange,
         fdr, Significant, Change_Direction)

print(top_20_metrics)

# 计算总体统计
summary_stats <- complete_metrics_table %>%
  summarise(
    Total_ASVs = n(),
    ASVs_with_Species = sum(!is.na(Species) & Species != ""),
    Mean_Prevalence = mean(Prevalence),
    Median_Total_Abundance = median(Total_Abundance),
    High_Abundance_ASVs = sum(Dominance_Level == "High"),
    Rare_ASVs = sum(Dominance_Level == "Rare"),
    Significant_ASVs = sum(Significant == "Yes", na.rm = TRUE)
  )

cat("\n=== 总体统计摘要 ===\n")
print(summary_stats)

# 按门水平统计
phylum_stats <- complete_metrics_table %>%
  group_by(Phylum) %>%
  summarise(
    ASV_Count = n(),
    Total_Abundance = sum(Total_Abundance),
    Mean_Prevalence = mean(Prevalence),
    Median_Abundance = median(Mean_Abundance),
    Significant_Count = sum(Significant == "Yes", na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(Total_Abundance))

cat("\n=== 按门水平统计 ===\n")
print(phylum_stats)

# 筛选显著差异的微生物
significant_microbes <- complete_metrics_table %>%
  filter(Significant == "Yes") %>%
  arrange(fdr)

cat("\n显著差异的微生物数量:", nrow(significant_microbes), "\n")

# 按效应大小排序
top_changers <- complete_metrics_table %>%
  filter(Significant == "Yes") %>%
  arrange(desc(abs(Log2_FoldChange))) %>%
  head(20)

cat("\nTop 20 变化最大的显著差异微生物:\n")
print(top_changers %>% select(ASV, Phylum, Genus, Log2_FoldChange, fdr, Change_Direction))

# 保存显著差异微生物表格
significant_table_path <- paste0(output_path, "Significant_Differential_Microbes.csv")
write.csv(significant_microbes, significant_table_path, row.names = FALSE)

# 保存门水平统计表格
phylum_table_path <- paste0(output_path, "Phylum_Level_Statistics.csv")
write.csv(phylum_stats, phylum_table_path, row.names = FALSE)

# 保存分析摘要到文本文件
sink(paste0(output_path, "Analysis_Summary.txt"))
cat("=== 微生物多样性分析摘要 ===\n\n")
cat("分析日期:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("数据文件: R_data_rename.csv\n")
cat("输出路径:", output_path, "\n\n")

cat("总体统计:\n")
print(summary_stats)
cat("\n")

cat("按门水平统计:\n")
print(phylum_stats)
cat("\n")

cat("显著差异的微生物数量:", nrow(significant_microbes), "\n")
cat("Top 20 变化最大的显著差异微生物:\n")
print(top_changers %>% select(ASV, Phylum, Genus, Log2_FoldChange, fdr, Change_Direction))
sink()

cat("\n分析完成！所有文件已保存到指定路径:\n")
cat("完整指标表格:", complete_table_path, "\n")
cat("显著差异微生物:", significant_table_path, "\n")
cat("门水平统计:", phylum_table_path, "\n")
cat("分析摘要:", paste0(output_path, "Analysis_Summary.txt"), "\n")
cat("表格包含", nrow(complete_metrics_table), "个ASV的详细指标\n")