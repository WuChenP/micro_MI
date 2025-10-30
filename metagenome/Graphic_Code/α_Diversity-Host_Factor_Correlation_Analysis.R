# 加载必要的包
library(readxl)
library(reshape2)
library(pheatmap)
library(dplyr)

# 设置输出目录
output_dir <- "D:/PythonProject/micro_MI/metagenome/Graphic/Virus-Host_Factor_Correlation_Heatmap"

# 1. 读取数据
data <- read_excel("D:/PythonProject/micro_MI/origin_data/样本协变量数据.xlsx", sheet = "Sheet1")

# 2. 数据预处理
colnames(data) <- c("分析名称", "分组", "性别", "年龄", "BMI", "吸烟", "饮酒", 
                    "Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "ACE")

# 将分组转换为数值：CON=0, MI=1
data$分组数值 <- ifelse(data$分组 == "MI", 1, 0)

# 3. 选择用于相关性分析的变量
host_factors <- c("分组数值", "性别", "年龄", "BMI", "吸烟", "饮酒")
alpha_diversity <- c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "ACE")

# 4. 计算Spearman相关性
cor_matrix <- matrix(NA, nrow = length(alpha_diversity), ncol = length(host_factors),
                     dimnames = list(alpha_diversity, host_factors))
pval_matrix <- matrix(NA, nrow = length(alpha_diversity), ncol = length(host_factors),
                      dimnames = list(alpha_diversity, host_factors))

for (i in 1:length(alpha_diversity)) {
  for (j in 1:length(host_factors)) {
    alpha_var <- data[[alpha_diversity[i]]]
    host_var <- data[[host_factors[j]]]
    
    test_result <- cor.test(alpha_var, host_var, method = "spearman", exact = FALSE)
    cor_matrix[i, j] <- test_result$estimate
    pval_matrix[i, j] <- test_result$p.value
  }
}

# 重命名列标签为英文
colnames(cor_matrix) <- c("MI", "Gender", "Age", "BMI", "Smoking", "Drinking")

# 5. 准备热图显示内容
sig_symbols <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
sig_symbols[pval_matrix < 0.05] <- "*"
sig_symbols[pval_matrix < 0.01] <- "**"
sig_symbols[pval_matrix < 0.001] <- "***"

display_matrix <- matrix(
  paste0(sprintf("%.3f", cor_matrix), "\n", sig_symbols),
  nrow = nrow(cor_matrix),
  dimnames = dimnames(cor_matrix)
)

# 6. Lancet风格配色
lancet_style <- colorRampPalette(c("#00468B", "white", "#ED0000"))(100)

# 7. Lancet风格热图
pdf(file.path(output_dir, "AlphaDiversity_Correlation_LancetStyle.pdf"), 
    width = 9, height = 7)
pheatmap(cor_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = display_matrix,
         number_color = "black",
         fontsize_number = 8,
         main = "Spearman Correlation: Alpha-Diversity vs Host Factors",
         color = lancet_style,
         breaks = seq(-1, 1, length.out = 101),
         cellwidth = 45,
         cellheight = 35,
         fontsize = 12,
         border_color = "white",
         angle_col = 0,
         legend = TRUE)
dev.off()

# 8. 保存相关性结果
cor_results <- data.frame(
  Alpha_Diversity = rep(alpha_diversity, each = length(host_factors)),
  Host_Factor = rep(c("MI", "Gender", "Age", "BMI", "Smoking", "Drinking"), times = length(alpha_diversity)),
  Correlation = as.vector(cor_matrix),
  P_value = as.vector(pval_matrix)
)

cor_results$Significance <- ""
cor_results$Significance[cor_results$P_value < 0.05] <- "*"
cor_results$Significance[cor_results$P_value < 0.01] <- "**"
cor_results$Significance[cor_results$P_value < 0.001] <- "***"

write.csv(cor_results, file.path(output_dir, "AlphaDiversity_HostFactors_Correlation_Results.csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")

# 9. 显示统计摘要
cat("分析完成！生成Lancet风格热图：\n")
cat("文件: AlphaDiversity_Correlation_LancetStyle.pdf\n")
cat("风格: Lancet风格（深蓝-白-深红）\n")
cat("\n总样本:", nrow(data), "(CON:", sum(data$分组=="CON"), "MI:", sum(data$分组=="MI"), ")\n")