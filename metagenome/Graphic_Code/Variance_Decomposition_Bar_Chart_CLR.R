# 加载必要的包
library(vegan)
library(tidyverse)
library(ggplot2)
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(matrixStats)  # 用于几何均值计算

# 设置结果保存路径
votu_save_path <- "E:/Python/MI_Analysis/metagenome/Graphic/Variance_Decomposition_Bar_Chart/vOTU"
family_save_path <- "E:/Python/MI_Analysis/metagenome/Graphic/Variance_Decomposition_Bar_Chart/Family"

# 创建保存目录
dir.create(votu_save_path, recursive = TRUE, showWarnings = FALSE)
dir.create(family_save_path, recursive = TRUE, showWarnings = FALSE)

# 读取数据
votu_data <- read.csv("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/心梗组_病毒_filtered_1percent.csv", 
                      header = TRUE, row.names = 1, check.names = FALSE)
family_data <- read.xlsx("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/virus_family_no_HF.xlsx", 
                         rowNames = TRUE, check.names = FALSE)
meta_data <- read_excel("E:/Python/MI_Analysis/origin_data/样本协变量数据.xlsx")

# 数据预处理
votu_matrix <- as.data.frame(t(votu_data))  # 转置，使样本为行
family_matrix <- as.data.frame(t(family_data))

required_meta <- meta_data %>%
  select(
    SampleID = 分析名称,
    Group = 分组,
    Gender = `1男2女`,
    Age = 年龄,
    BMI = BMI
    # Smoking = `吸烟\n（0/1）`,  # 注释掉抽烟
    # Drinking = `饮酒\n（0/1）`  # 注释掉喝酒
  ) %>%
  mutate(
    MI = ifelse(Group == "MI", 1, 0),  # 创建MI变量
    Gender = factor(Gender, levels = c(1, 2), labels = c("Male", "Female")),
    # Smoking = factor(Smoking, levels = c(0, 1), labels = c("No", "Yes")),  # 注释掉抽烟
    # Drinking = factor(Drinking, levels = c(0, 1), labels = c("No", "Yes")),  # 注释掉喝酒
    Age = as.numeric(Age),
    BMI = as.numeric(BMI)
  ) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("SampleID")

# 数据验证
cat("=== Data Validation ===\n")
cat("vOTU data dimensions:", dim(votu_matrix), "\n")
cat("Family level data dimensions:", dim(family_matrix), "\n")
cat("Metadata dimensions:", dim(required_meta), "\n")

cat("\n=== Variable Distribution Check ===\n")
cat("MI distribution:\n"); print(table(required_meta$MI))
cat("Gender distribution:\n"); print(table(required_meta$Gender))
# cat("Smoking distribution:\n"); print(table(required_meta$Smoking))  # 注释掉抽烟
# cat("Drinking distribution:\n"); print(table(required_meta$Drinking))  # 注释掉喝酒

# 样本匹配检查函数
check_sample_overlap <- function(community_data, metadata, data_name) {
  comm_samples <- rownames(community_data)
  meta_samples <- rownames(metadata)
  common_samples <- intersect(comm_samples, meta_samples)
  cat(paste("\n=== ", data_name, "Sample Matching Check ===\n"))
  cat("Community data samples:", length(comm_samples), "\n")
  cat("Metadata samples:", length(meta_samples), "\n")
  cat("Common samples:", length(common_samples), "\n")
  return(common_samples)
}

# CLR 转换函数
clr_transform <- function(data_matrix, pseudocount = 1e-6) {
  # 1. 过滤掉全零行
  data_clean <- data_matrix[rowSums(data_matrix, na.rm = TRUE) > 0, ]
  # 2. 样本内相对丰度
  data_rel <- sweep(data_clean, 1, rowSums(data_clean), FUN = "/")
  # 3. 加伪计数
  data_rel <- data_rel + pseudocount
  # 4. 计算每个样本几何均值
  geom_means <- exp(rowMeans(log(data_rel)))
  # 5. CLR 转换
  data_clr <- log(data_rel) - log(geom_means)
  rownames(data_clr) <- rownames(data_rel)
  colnames(data_clr) <- colnames(data_rel)
  return(data_clr)
}

# 正确的PERMANOVA分析函数 - 使用边际检验
run_permanova_corrected <- function(community_data, metadata, level_name) {
  cat(paste("\n=== Starting", level_name, "Level Analysis ===\n"))
  common_samples <- check_sample_overlap(community_data, metadata, level_name)
  if (length(common_samples) == 0) stop(paste("No common samples found between", level_name, "data and metadata"))
  
  comm_processed <- clr_transform(community_data[common_samples, ])
  meta_common <- metadata[common_samples, ]
  
  # CLR 后用 Euclidean 距离（Aitchison 距离）
  dist_matrix <- dist(comm_processed, method = "euclidean")
  
  results <- data.frame()
  # 只保留MI、Age、Gender、BMI四个变量
  variables <- c("MI", "Age", "Gender", "BMI")  # 移除Smoking和Drinking
  
  formula_full <- as.formula(paste("dist_matrix ~", paste(variables, collapse = " + ")))
  cat("Running full model with marginal contributions...\n")
  permanova_marginal <- adonis2(formula_full, data = meta_common, permutations = 999, by = "margin")
  cat("Marginal contributions from full model:\n"); print(permanova_marginal)
  
  for (i in 1:length(variables)) {
    var <- variables[i]
    cat(paste("  Processing variable:", var, "\n"))
    tryCatch({
      formula_simple <- as.formula(paste("dist_matrix ~", var))
      permanova_simple <- adonis2(formula_simple, data = meta_common, permutations = 999)
      if (i <= nrow(permanova_marginal)) {
        marginal_R2 <- permanova_marginal$R2[i]
        marginal_p <- permanova_marginal$`Pr(>F)`[i]
        marginal_F <- permanova_marginal$F[i]
      } else {
        marginal_R2 <- NA; marginal_p <- NA; marginal_F <- NA
      }
      simple_result <- data.frame(
        Level = level_name, Variable = var, Model = "No covariates",
        R2 = permanova_simple$R2[1], Pvalue = permanova_simple$`Pr(>F)`[1], F_value = permanova_simple$F[1]
      )
      covar_result <- data.frame(
        Level = level_name, Variable = var, Model = "With covariates",
        R2 = marginal_R2, Pvalue = marginal_p, F_value = marginal_F
      )
      results <- rbind(results, simple_result, covar_result)
    }, error = function(e) {
      message(paste("Error analyzing variable", var, "in", level_name, ":", e$message))
      error_result <- data.frame(
        Level = level_name,
        Variable = var,
        Model = rep(c("No covariates", "With covariates"), each = 1),
        R2 = NA, Pvalue = NA, F_value = NA
      )
      results <- rbind(results, error_result)
    })
  }
  return(results)
}

# 执行分析
cat("\nStarting PERMANOVA analysis with marginal contributions...\n")
votu_results <- run_permanova_corrected(votu_matrix, required_meta, "vOTU")
family_results <- run_permanova_corrected(family_matrix, required_meta, "Family")
all_results <- rbind(votu_results, family_results)

# 可视化函数（更新颜色和变量）
create_journal_style_plot <- function(results_data, level_name) {
  # 只保留MI、Age、Gender、BMI的颜色定义
  factor_colors <- c("MI"="#2E86AB","Age"="#A23B72","Gender"="#F18F01","BMI"="#C73E1D")
  shape_values <- c("No covariates"=21,"With covariates"=22)
  
  plot_data <- results_data %>%
    filter(Level==level_name & !is.na(R2)) %>%
    mutate(
      Significance=case_when(Pvalue<=0.001 ~ "***", Pvalue<=0.01 ~ "**", Pvalue<=0.05 ~ "*", TRUE ~ ""),
      Variable=factor(Variable, levels=c("MI","Age","Gender","BMI")),  # 只保留四个变量
      R2_percent=R2*100
    ) %>%
    group_by(Variable) %>%
    mutate(max_r2=max(R2_percent, na.rm=TRUE), label_x=R2_percent + max_r2*0.05) %>%
    ungroup()
  
  x_max <- max(plot_data$R2_percent, na.rm=TRUE)*1.2
  p <- ggplot(plot_data, aes(x=R2_percent, y=Variable)) +
    geom_point(aes(shape=Model, fill=Variable), size=4, stroke=0.8, position=position_dodge(width=0.6)) +
    geom_text(aes(label=Significance, x=label_x, group=Model), position=position_dodge(width=0.6), vjust=0.7, size=4, fontface="bold") +
    scale_fill_manual(values=factor_colors, guide="none") +
    scale_shape_manual(values=shape_values) +
    scale_x_continuous(limits=c(0,x_max), expand=expansion(mult=c(0,0.05))) +
    labs(x="Explained variance (%)", y=NULL, shape="Model", title=paste("PERMANOVA -",level_name,"Level")) +
    theme_classic(base_size=12) +
    theme(text=element_text(color="black"), axis.text=element_text(size=11,color="black"),
          axis.title=element_text(size=12,face="bold"), axis.line=element_line(color="black",linewidth=0.5),
          axis.ticks=element_line(color="black",linewidth=0.5), legend.position="top",
          legend.title=element_text(face="bold"), panel.grid.major.x=element_line(color="grey90",linewidth=0.2),
          plot.title=element_text(hjust=0.5,face="bold"))
  return(p)
}

# 生成并保存图形和结果
cat("\n=== Generating Visualizations and Saving Results ===\n")
if (nrow(votu_results) > 0 && sum(!is.na(votu_results$R2)) > 0) {
  votu_plot <- create_journal_style_plot(all_results,"vOTU")
  ggsave(file.path(votu_save_path,"PERMANOVA_vOTU_level.png"), votu_plot, width=8,height=5,dpi=300,bg="white")
  ggsave(file.path(votu_save_path,"PERMANOVA_vOTU_level.pdf"), votu_plot, width=8,height=5,bg="white",device=cairo_pdf)
  cat("vOTU level plots saved\n")
}
if (nrow(family_results) > 0 && sum(!is.na(family_results$R2)) > 0) {
  family_plot <- create_journal_style_plot(all_results,"Family")
  ggsave(file.path(family_save_path,"PERMANOVA_Family_level.png"), family_plot, width=8,height=5,dpi=300,bg="white")
  ggsave(file.path(family_save_path,"PERMANOVA_Family_level.pdf"), family_plot, width=8,height=5,bg="white",device=cairo_pdf)
  cat("Family level plots saved\n")
}

write.csv(votu_results, file.path(votu_save_path,"PERMANOVA_vOTU_results.csv"), row.names=FALSE)
write.csv(family_results, file.path(family_save_path,"PERMANOVA_Family_results.csv"), row.names=FALSE)
write.csv(all_results, file.path(dirname(votu_save_path),"PERMANOVA_All_Results.csv"), row.names=FALSE)
cat("Result files saved\n")

cat("\n=== Analysis Results Summary ===\n")
cat("vOTU level results:\n"); print(votu_results %>% select(Level, Variable, Model, R2, Pvalue) %>% mutate(R2=round(R2,4), Pvalue=round(Pvalue,4)))
cat("\nFamily level results:\n"); print(family_results %>% select(Level, Variable, Model, R2, Pvalue) %>% mutate(R2=round(R2,4), Pvalue=round(Pvalue,4)))
cat("\nAnalysis completed! All results saved to specified directories.\n")