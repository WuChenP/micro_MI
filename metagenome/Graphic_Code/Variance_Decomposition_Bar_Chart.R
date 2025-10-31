# 加载必要的包
library(vegan)
library(tidyverse)
library(ggplot2)
library(readxl)
library(openxlsx)
library(RColorBrewer)

# 设置结果保存路径
votu_save_path <- "E:/Python/MI_Analysis/metagenome/Graphic/Variance_Decomposition_Bar_Chart/vOTU"
family_save_path <- "E:/Python/MI_Analysis/metagenome/Graphic/Variance_Decomposition_Bar_Chart/Family"

# 创建保存目录
dir.create(votu_save_path, recursive = TRUE, showWarnings = FALSE)
dir.create(family_save_path, recursive = TRUE, showWarnings = FALSE)

# 读取数据
## 读取vOTU数据
votu_data <- read.csv("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/心梗组_病毒_filtered_1percent.csv", 
                      header = TRUE, row.names = 1, check.names = FALSE)

## 读取科水平数据
family_data <- read.xlsx("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/virus_family_no_HF.xlsx", 
                         rowNames = TRUE, check.names = FALSE)

## 读取宿主因素数据
meta_data <- read_excel("E:/Python/MI_Analysis/origin_data/样本协变量数据.xlsx")

# 数据预处理
## 转换数据格式
votu_matrix <- as.data.frame(t(votu_data))  # 转置，使样本为行
family_matrix <- as.data.frame(t(family_data))

## 整理宿主因素数据
required_meta <- meta_data %>%
  select(
    SampleID = 分析名称,
    Group = 分组,
    Gender = `1男2女`,
    Age = 年龄,
    BMI = BMI,
    Smoking = `吸烟\n（0/1）`,
    Drinking = `饮酒\n（0/1）`
  ) %>%
  mutate(
    MI = ifelse(Group == "MI", 1, 0),  # 创建MI变量
    Gender = factor(Gender, levels = c(1, 2), labels = c("Male", "Female")),
    Smoking = factor(Smoking, levels = c(0, 1), labels = c("No", "Yes")),
    Drinking = factor(Drinking, levels = c(0, 1), labels = c("No", "Yes")),
    Age = as.numeric(Age),
    BMI = as.numeric(BMI)
  ) %>%
  filter(complete.cases(.)) %>%  # 移除有缺失值的行
  column_to_rownames("SampleID")

# 数据预处理函数
preprocess_data <- function(data_matrix) {
  # 移除全为零的行
  data_clean <- data_matrix[rowSums(data_matrix, na.rm = TRUE) > 0, ]
  # 移除全为零的列
  data_clean <- data_clean[, colSums(data_clean, na.rm = TRUE) > 0]
  # Hellinger转换
  data_hellinger <- decostand(data_clean, "hellinger")
  return(data_hellinger)
}

# PERMANOVA分析函数
run_permanova_analysis <- function(community_data, metadata, level_name) {
  
  # 预处理数据
  comm_processed <- preprocess_data(community_data)
  
  # 计算Bray-Curtis距离
  dist_matrix <- vegdist(comm_processed, method = "bray")
  
  # 确保样本顺序一致
  common_samples <- intersect(rownames(comm_processed), rownames(metadata))
  if (length(common_samples) == 0) {
    stop(paste("No common samples found between", level_name, "data and metadata"))
  }
  
  dist_matrix <- as.dist(as.matrix(dist_matrix)[common_samples, common_samples])
  metadata_common <- metadata[common_samples, ]
  
  # 存储结果
  results <- data.frame()
  
  # 定义要分析的变量
  variables <- c("MI", "Age", "Gender", "BMI", "Smoking", "Drinking")
  
  for (var in variables) {
    tryCatch({
      # 无协变量模型
      formula_simple <- as.formula(paste("dist_matrix ~", var))
      permanova_simple <- adonis2(formula_simple, data = metadata_common, permutations = 999)
      
      # 有协变量模型
      other_vars <- setdiff(variables, var)
      covariates <- paste(other_vars[1:3], collapse = " + ")  # 取前3个作为协变量
      formula_covar <- as.formula(paste("dist_matrix ~", var, "+", covariates))
      permanova_covar <- adonis2(formula_covar, data = metadata_common, permutations = 999)
      
      # 提取结果
      simple_result <- data.frame(
        Level = level_name,
        Variable = var,
        Model = "No covariates",
        R2 = permanova_simple$R2[1],
        Pvalue = permanova_simple$`Pr(>F)`[1]
      )
      
      covar_result <- data.frame(
        Level = level_name,
        Variable = var,
        Model = "With covariates",
        R2 = permanova_covar$R2[1],
        Pvalue = permanova_covar$`Pr(>F)`[1]
      )
      
      results <- rbind(results, simple_result, covar_result)
    }, error = function(e) {
      message(paste("Error analyzing variable", var, "in", level_name, ":", e$message))
    })
  }
  
  return(results)
}

# 执行分析
cat("开始vOTU水平分析...\n")
votu_results <- run_permanova_analysis(votu_matrix, required_meta, "vOTU")

cat("开始科水平分析...\n")
family_results <- run_permanova_analysis(family_matrix, required_meta, "Family")

# 合并结果
all_results <- rbind(votu_results, family_results)

# 创建顶刊风格的可视化函数（解决字体错误）
create_journal_style_plot <- function(results_data, level_name) {
  
  # 定义宿主因素颜色方案
  factor_colors <- c(
    "MI" = "#2E86AB",      # 深蓝色
    "Age" = "#A23B72",     # 紫红色
    "Gender" = "#F18F01",  # 橙色
    "BMI" = "#C73E1D",     # 红色
    "Smoking" = "#3A7D34", # 绿色
    "Drinking" = "#6D4C3D" # 棕色
  )
  
  # 定义模型形状
  shape_values <- c("No covariates" = 21, "With covariates" = 22)
  
  # 处理数据
  plot_data <- results_data %>%
    filter(Level == level_name) %>%
    mutate(
      Significance = case_when(
        Pvalue < 0.001 ~ "***",
        Pvalue < 0.01 ~ "**",
        Pvalue < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Variable = factor(Variable, levels = c("MI", "Age", "Gender", "BMI", "Smoking", "Drinking")),
      R2_percent = R2 * 100,
      label_x = R2_percent + 0.5
    )
  
  # 计算横坐标范围，确保包含所有刻度
  max_r2 <- max(plot_data$R2_percent, na.rm = TRUE)
  x_breaks <- seq(0, ceiling(max_r2), by = 1)  # 每1%一个刻度
  
  # 创建图形
  p <- ggplot(plot_data, aes(x = R2_percent, y = Variable, 
                             shape = Model, fill = Variable)) +
    geom_point(size = 4, stroke = 0.8, position = position_dodge(width = 0.7)) +
    geom_text(aes(label = Significance, x = label_x), 
              position = position_dodge(width = 0.7),
              vjust = 0.7, size = 4, fontface = "bold") +
    scale_fill_manual(values = factor_colors, guide = "none") +
    scale_shape_manual(values = shape_values,
                       labels = c("No covariates" = "No covariates", 
                                  "With covariates" = "With covariates")) +
    scale_x_continuous(
      breaks = x_breaks,  # 设置详细的刻度
      limits = c(0, max(x_breaks)),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "Explained variance (%)",
      y = NULL,
      shape = "Model"
    ) +
    theme_classic(base_size = 12) +
    theme(
      # 使用默认字体避免错误
      text = element_text(family = "", color = "black"),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.key = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor.x = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# 生成并保存图形和结果
cat("生成可视化图形和保存结果...\n")

# vOTU水平结果保存
votu_plot <- create_journal_style_plot(all_results, "vOTU")
# 保存PNG文件（不使用特殊字体）
ggsave(file.path(votu_save_path, "PERMANOVA_vOTU_level.png"), 
       votu_plot, width = 8, height = 5, dpi = 300, bg = "white")
# 保存PDF文件（使用cairo_pdf避免字体问题）
ggsave(file.path(votu_save_path, "PERMANOVA_vOTU_level.pdf"), 
       votu_plot, width = 8, height = 5, bg = "white", device = cairo_pdf)

# 科水平结果保存
family_plot <- create_journal_style_plot(all_results, "Family")
ggsave(file.path(family_save_path, "PERMANOVA_Family_level.png"), 
       family_plot, width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(family_save_path, "PERMANOVA_Family_level.pdf"), 
       family_plot, width = 8, height = 5, bg = "white", device = cairo_pdf)

# 保存详细结果
write.csv(votu_results, file.path(votu_save_path, "PERMANOVA_vOTU_results.csv"), row.names = FALSE)
write.csv(family_results, file.path(family_save_path, "PERMANOVA_Family_results.csv"), row.names = FALSE)
write.csv(all_results, file.path(dirname(votu_save_path), "PERMANOVA_All_Results.csv"), row.names = FALSE)

# 显示图形
print(votu_plot)
print(family_plot)

# 显示主要结果摘要
cat("\n=== 分析结果摘要 ===\n")
cat("vOTU水平结果:\n")
print(votu_results)
cat("\n科水平结果:\n")
print(family_results)

# 样本信息统计
cat("\n=== 样本信息 ===\n")
cat("vOTU数据样本数:", nrow(votu_matrix), "\n")
cat("科水平数据样本数:", nrow(family_matrix), "\n")
cat("宿主因素数据样本数:", nrow(required_meta), "\n")

# 检查共同样本
common_votu <- intersect(rownames(votu_matrix), rownames(required_meta))
common_family <- intersect(rownames(family_matrix), rownames(required_meta))
cat("vOTU与元数据共同样本数:", length(common_votu), "\n")
cat("科水平与元数据共同样本数:", length(common_family), "\n")

cat("\n分析完成！所有结果已保存到指定目录。\n")