# =============================================================================
# 细菌-病毒Family水平相关性网络分析（Cytoscape兼容修复版本）
# 要求：|Spearman rho| > 0.50 且 q < 0.01
# 结果保存到：E:\Python\MI_Analysis\metagenome\Graphic\Bacteria_and_virus_network_diagram
# =============================================================================

# 加载必要的包
library(dplyr)
library(tidyr)
library(reshape2)
library(readxl)
library(writexl)  # 新增：用于生成xlsx文件

# 设置文件路径
bacteria_file <- "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/family_level_data/bacteria_family.xlsx"
virus_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/virus_family_no_HF.xlsx"
output_dir <- "E:/Python/MI_Analysis/metagenome/Graphic/Bacteria_and_virus_network_diagram"

# 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

cat("=== 细菌-病毒Family水平相关性网络分析开始 ===\n")
cat("结果将保存到:", output_dir, "\n")

# =============================================================================
# 颜色生成函数 - 为每个family分配对比明显的颜色
# =============================================================================

generate_distinct_colors <- function(n, type = "bacteria") {
  # 预定义两组对比明显的颜色
  if (type == "bacteria") {
    # 细菌使用蓝色系，但每个family不同色调
    base_colors <- c(
      "#1f77b4", "#2ca02c", "#17becf", "#7f7f7f", "#bcbd22",
      "#8c564b", "#e377c2", "#ff7f0e", "#9467bd", "#d62728",
      "#1a9850", "#66bd63", "#a6d96a", "#d9ef8b", "#fdae61",
      "#f46d43", "#d73027", "#4575b4", "#74add1", "#abd9e9",
      "#fee090", "#e0f3f8", "#ffffbf", "#f46d43", "#d73027"
    )
  } else {
    # 病毒使用红色/暖色系，但每个family不同色调
    base_colors <- c(
      "#d62728", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2",
      "#dbdb8d", "#9edae5", "#ffbb78", "#98df8a", "#ff9896",
      "#e6550d", "#fd8d3c", "#fdae6b", "#fdd0a2", "#31a354",
      "#74c476", "#a1d99b", "#c7e9c0", "#756bb1", "#9e9ac8",
      "#bcbddc", "#dadaeb", "#3182bd", "#6baed6", "#9ecae1"
    )
  }
  
  # 如果需要的颜色多于预定义的颜色，使用颜色生成
  if (n > length(base_colors)) {
    additional_colors <- grDevices::rainbow(n - length(base_colors))
    all_colors <- c(base_colors, additional_colors)
  } else {
    all_colors <- base_colors[1:n]
  }
  
  return(all_colors)
}

# =============================================================================
# 1. 数据读取和预处理
# =============================================================================

cat("1. 读取数据...\n")

# 读取细菌数据（Excel格式）
bacteria_data <- read_excel(bacteria_file)
# 设置第一列为行名
bacteria_df <- as.data.frame(bacteria_data[, -1])
rownames(bacteria_df) <- bacteria_data[[1]]
bacteria_data <- bacteria_df

cat("   细菌数据维度: ", nrow(bacteria_data), "个Family ×", ncol(bacteria_data), "个样本\n")

# 读取病毒数据（Excel格式）
virus_data <- read_excel(virus_file)
# 设置第一列为行名
virus_df <- as.data.frame(virus_data[, -1])
rownames(virus_df) <- virus_data[[1]]
virus_data <- virus_df

cat("   病毒数据维度: ", nrow(virus_data), "个Family ×", ncol(virus_data), "个样本\n")

# 转置数据，样本为行
bacteria_t <- as.data.frame(t(bacteria_data))
virus_t <- as.data.frame(t(virus_data))

# 查看样本交集
common_samples <- intersect(rownames(bacteria_t), rownames(virus_t))
cat("   共同样本数: ", length(common_samples), "\n")

# 按共同样本筛选数据
bacteria_common <- bacteria_t[common_samples, ]
virus_common <- virus_t[common_samples, ]

cat("   处理后细菌Family数: ", ncol(bacteria_common), "\n")
cat("   处理后病毒Family数: ", ncol(virus_common), "\n")
cat("   最终样本数: ", nrow(bacteria_common), "\n")

# 显示前几个Family名称
cat("\n   细菌Family示例:\n")
print(head(colnames(bacteria_common), 5))
cat("\n   病毒Family示例:\n")
print(head(colnames(virus_common), 5))

# =============================================================================
# 2. 计算细菌-病毒Spearman相关性
# =============================================================================

cat("\n2. 计算细菌-病毒Spearman相关性...\n")

# 初始化结果存储
cor_results <- data.frame()
total_pairs <- ncol(bacteria_common) * ncol(virus_common)
cat("   需要计算的相关性对数: ", total_pairs, "\n")

# 显示计算进度
start_time <- Sys.time()
cat("   开始时间: ", format(start_time), "\n")

# 分批计算相关性
completed <- 0

for (i in 1:ncol(bacteria_common)) {
  bacteria_name <- colnames(bacteria_common)[i]
  
  for (j in 1:ncol(virus_common)) {
    virus_name <- colnames(virus_common)[j]
    
    # 计算Spearman相关性
    cor_test <- cor.test(bacteria_common[, i], virus_common[, j], 
                         method = "spearman", exact = FALSE)
    
    # 存储结果
    cor_results <- rbind(cor_results, data.frame(
      Bacteria = bacteria_name,
      Virus = virus_name,
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value
    ))
    
    completed <- completed + 1
    
    # 进度显示（每500对显示一次）
    if (completed %% 500 == 0) {
      current_time <- Sys.time()
      elapsed <- as.numeric(difftime(current_time, start_time, units = "mins"))
      cat("   已完成: ", completed, "/", total_pairs, 
          " (", round(completed/total_pairs*100, 1), "%)",
          " 耗时: ", round(elapsed, 1), "分钟\n")
    }
  }
}

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
cat("   相关性计算完成! 总耗时: ", round(total_time, 1), "分钟\n")

# =============================================================================
# 3. FDR校正和相关性过滤
# =============================================================================

cat("\n3. FDR校正和相关性过滤...\n")

# FDR校正
cor_results$Q_value <- p.adjust(cor_results$P_value, method = "fdr")

# 应用过滤条件：|r| > 0.5, q < 0.01
filtered_edges <- cor_results %>%
  filter(abs(Correlation) > 0.5 & Q_value < 0.01) %>%
  mutate(
    direction = ifelse(Correlation > 0, "positive", "negative"),
    weight = abs(Correlation),
    interaction = "correlates_with"
  )

cat("   总相关性对数: ", nrow(cor_results), "\n")
cat("   过滤后相关性边数: ", nrow(filtered_edges), "\n")
cat("   正相关边数: ", sum(filtered_edges$direction == "positive"), "\n")
cat("   负相关边数: ", sum(filtered_edges$direction == "negative"), "\n")

# 显示结果摘要
if (nrow(filtered_edges) > 0) {
  cat("\n   相关性强度分布:\n")
  print(summary(abs(filtered_edges$Correlation)))
  
  cat("\n   最强的5个正相关:\n")
  print(head(filtered_edges[filtered_edges$direction == "positive", ] %>% 
               arrange(desc(Correlation)) %>% 
               select(Bacteria, Virus, Correlation, Q_value), 5))
  
  cat("\n   最强的5个负相关:\n")  
  print(head(filtered_edges[filtered_edges$direction == "negative", ] %>% 
               arrange(Correlation) %>% 
               select(Bacteria, Virus, Correlation, Q_value), 5))
}

# =============================================================================
# 4. 生成Cytoscape输入文件（XLSX格式，每个family不同颜色）
# =============================================================================

if (nrow(filtered_edges) > 0) {
  cat("\n4. 生成Cytoscape输入文件（XLSX格式，每个family不同颜色）...\n")
  
  # 准备边文件 - 使用标准列名
  edges_file <- filtered_edges %>%
    select(Bacteria, Virus, Correlation, weight, direction, interaction) %>%
    rename(source = Bacteria, target = Virus, correlation = Correlation)
  
  # 准备节点文件
  all_bacteria <- unique(filtered_edges$Bacteria)
  all_virus <- unique(filtered_edges$Virus)
  
  # 为每个family生成独特的颜色
  bacteria_colors <- generate_distinct_colors(length(all_bacteria), "bacteria")
  virus_colors <- generate_distinct_colors(length(all_virus), "virus")
  
  # 创建颜色映射表（用于图例）
  color_mapping <- data.frame(
    family = c(all_bacteria, all_virus),
    color = c(bacteria_colors, virus_colors),
    type = c(rep("Bacteria", length(all_bacteria)), 
             rep("Virus", length(all_virus))),
    stringsAsFactors = FALSE
  )
  
  # 创建节点文件
  nodes_file <- data.frame(
    name = c(all_bacteria, all_virus),
    node_type = c(rep("bacteria_family", length(all_bacteria)), 
                  rep("virus_family", length(all_virus))),
    group = c(rep("Bacteria", length(all_bacteria)), 
              rep("Virus", length(all_virus))),
    color = c(bacteria_colors, virus_colors),
    shape = c(rep("ELLIPSE", length(all_bacteria)), 
              rep("TRIANGLE", length(all_virus))),
    size = 50,
    stringsAsFactors = FALSE
  )
  
  # 保存XLSX文件
  edges_file_path <- file.path(output_dir, "bacteria_virus_family_network_edges.xlsx")
  nodes_file_path <- file.path(output_dir, "bacteria_virus_family_network_nodes.xlsx")
  color_mapping_path <- file.path(output_dir, "family_color_mapping.xlsx")
  
  write_xlsx(edges_file, edges_file_path)
  write_xlsx(nodes_file, nodes_file_path)
  write_xlsx(color_mapping, color_mapping_path)
  
  # 保存完整的相关性结果
  full_results_path <- file.path(output_dir, "full_family_correlation_results.xlsx")
  filtered_details_path <- file.path(output_dir, "filtered_family_correlations.xlsx")
  
  write_xlsx(cor_results, full_results_path)
  write_xlsx(filtered_edges, filtered_details_path)
  
  cat("   XLSX文件已生成!\n")
  cat("   节点文件: ", nodes_file_path, "\n")
  cat("   边文件: ", edges_file_path, "\n")
  cat("   颜色映射文件: ", color_mapping_path, "\n")
  
  # 显示颜色分配信息
  cat("\n   颜色分配情况:\n")
  cat("   - 细菌Family数量: ", length(all_bacteria), "\n")
  cat("   - 病毒Family数量: ", length(all_virus), "\n")
  cat("   - 使用的细菌颜色示例: ", paste(head(bacteria_colors), collapse = ", "), "\n")
  cat("   - 使用的病毒颜色示例: ", paste(head(virus_colors), collapse = ", "), "\n")
  
  # 验证节点匹配
  edge_nodes <- unique(c(edges_file$source, edges_file$target))
  node_file_nodes <- nodes_file$name
  
  cat("\n   节点匹配验证:\n")
  cat("   - 边文件中节点数: ", length(edge_nodes), "\n")
  cat("   - 节点文件中节点数: ", length(node_file_nodes), "\n")
  cat("   - 完全匹配: ", length(intersect(edge_nodes, node_file_nodes)), "\n")
  cat("   - 只在边文件中: ", length(setdiff(edge_nodes, node_file_nodes)), "\n")
  cat("   - 只在节点文件中: ", length(setdiff(node_file_nodes, edge_nodes)), "\n")
  
} else {
  cat("\n4. 结果: 没有找到满足条件(|r|>0.5, q<0.01)的相关性关系\n")
  cat("   建议尝试以下方法:\n")
  cat("   - 降低相关系数阈值至 |r|>0.4\n")
  cat("   - 放宽q值阈值至 q<0.05\n")
  cat("   - 检查数据质量\n")
}

cat("\n=== 分析完成 ===\n")
cat("所有结果已保存到:", output_dir, "\n")
cat("生成的文件:\n")
cat("- bacteria_virus_family_network_nodes.xlsx (节点文件)\n")
cat("- bacteria_virus_family_network_edges.xlsx (边文件)\n")
cat("- family_color_mapping.xlsx (颜色映射表，用于图例)\n")
cat("- full_family_correlation_results.xlsx (完整相关性结果)\n")
cat("- filtered_family_correlations.xlsx (过滤后的相关性)\n")