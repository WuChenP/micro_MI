# =========================================
# 🔬 Alpha多样性箱线图分析
# =========================================

setwd("D:/Rworkspace/micro_MI/oral_microorganisms")
options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# ---- 加载依赖 ----
pkgs <- c("vegan", "ggplot2", "ggpubr", "dplyr")
for(p in pkgs){
  if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ---- 读取数据 ----
file_path <- "D:/机器学习与生物信息学/北京朝阳医院项目/R/A/img_filter/R_data_filtered_species_modified.csv"
df <- read.csv(file_path, header = TRUE, check.names = FALSE)
cat("✅ 数据维度:", nrow(df), "行 ×", ncol(df), "列\n")

# ---- 数据预处理 ----
# 清理列名
colnames(df) <- gsub("[[:space:][:cntrl:]]", "", colnames(df))

# 打印列名信息
cat("前20列名:\n")
print(head(colnames(df), 20))

# ---- 提取样本数据 ----
sample_cols <- grep("^(AMI|CON)", colnames(df), value = TRUE, ignore.case = TRUE)
if(length(sample_cols) == 0){
  stop("❌ 没有匹配到 AMI/CON 样本列，请检查列名！")
}
cat("📊 样本列数:", length(sample_cols), "\n")

# 提取OTU矩阵
otu_mat <- as.matrix(df[, sample_cols, drop = FALSE])
cat("📊 OTU矩阵维度:", nrow(otu_mat), "行 ×", ncol(otu_mat), "列\n")

# ---- 转置并创建样本数据 ----
# 转置：行为样本，列为物种
otu_t <- t(otu_mat)
cat("📊 转置后矩阵维度:", nrow(otu_t), "行 ×", ncol(otu_t), "列\n")

# 创建分组信息
sample_names <- rownames(otu_t)
Group <- ifelse(grepl("^AMI", sample_names, ignore.case = TRUE), "AMI", "Control")
cat("📊 分组情况:\n")
print(table(Group))

# 转换为数据框并确保数值类型
otu_df <- as.data.frame(apply(otu_t, 2, as.numeric))
rownames(otu_df) <- sample_names

# ---- Alpha 多样性计算 ----
cat("🔬 开始计算Alpha多样性...\n")

# 初始化结果数据框
alpha_div <- data.frame(
  SampleID = rownames(otu_df),
  Group = Group,
  stringsAsFactors = FALSE
)

# 计算基本多样性指标
safe_apply <- function(data, margin, fun) {
  tryCatch({
    apply(data, margin, fun)
  }, error = function(e) {
    cat("⚠️ 计算出错:", e$message, "\n")
    rep(NA, ifelse(margin == 1, nrow(data), ncol(data)))
  })
}

# Richness (物种丰富度)
alpha_div$Richness <- safe_apply(otu_df, 1, function(x) sum(x > 0))

# Shannon多样性
alpha_div$Shannon <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  x <- x[x > 0]
  p <- x / sum(x)
  -sum(p * log(p))
})

# Simpson多样性
alpha_div$Simpson <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  x <- x[x > 0]
  p <- x / sum(x)
  1 - sum(p^2)
})

# Inverse Simpson
alpha_div$InvSimpson <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  x <- x[x > 0]
  p <- x / sum(x)
  1 / sum(p^2)
})

# Chao1和ACE计算
alpha_div$Chao1 <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  x <- x[x > 0]
  S_obs <- length(x)
  F1 <- sum(x == 1)
  F2 <- sum(x == 2)
  if(F2 == 0) F2 <- 1
  S_chao1 <- S_obs + (F1 * (F1 - 1)) / (2 * (F2 + 1))
  return(S_chao1)
})

alpha_div$ACE <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  x <- x[x > 0]
  S_obs <- length(x)
  F1 <- sum(x == 1)
  S_ace <- S_obs + F1
  return(S_ace)
})

# Coverage
alpha_div$Coverage <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  1 - sum(x == 1) / sum(x)
})

# 保存结果
write.csv(alpha_div, "Alpha_diversity_results.csv", row.names = FALSE)
cat("✅ Alpha多样性结果已保存\n")

# ---- Alpha多样性箱线图绘制 ----
cat("🔬 开始绘制Alpha多样性箱线图...\n")

# 定义颜色
colors <- c("AMI" = "#E64B35FF", "Control" = "#4DBBD5FF")

# 绘制单个指标箱线图的函数
plot_alpha_boxplot <- function(metric) {
  # 检查是否有有效数据
  if(all(is.na(alpha_div[[metric]]))) {
    cat("⚠️ 跳过", metric, "，所有值都为NA\n")
    return()
  }
  
  # 移除NA值
  plot_data <- alpha_div[!is.na(alpha_div[[metric]]), ]
  
  if(nrow(plot_data) == 0) {
    cat("⚠️ 跳过", metric, "，无有效数据\n")
    return()
  }
  
  # 绘制箱线图
  p <- ggplot(plot_data, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.x = 1.5,
                     size = 5) +
    labs(title = metric,
         x = "", 
         y = metric) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_fill_manual(values = colors)
  
  filename <- paste0(metric, "_boxplot.png")
  ggsave(filename, p, width = 5, height = 5, dpi = 300)
  cat("💾 保存箱线图:", filename, "\n")
}

# 绘制所有Alpha多样性指标的箱线图
metrics <- c("Richness", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Coverage")
for(m in metrics) {
  plot_alpha_boxplot(m)
}

cat("🎉 Alpha多样性箱线图分析完成！\n")
cat("📊 生成的箱线图文件:\n")
for(m in metrics) {
  cat("  - ", m, "_boxplot.png\n")
}