# =========================================
# 🔬 Alpha & Beta 多样性分析 - 稳健终极版
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
df <- read.csv(file_path, header = TRUE, check.names = FALSE, row.names = 1)
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

# 计算基本多样性指标（使用安全的apply）
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

# Chao1和ACE计算（简化安全版本）
alpha_div$Chao1 <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  x <- x[x > 0]
  S_obs <- length(x)
  F1 <- sum(x == 1)
  F2 <- sum(x == 2)
  if(F2 == 0) F2 <- 1  # 避免除零
  S_chao1 <- S_obs + (F1 * (F1 - 1)) / (2 * (F2 + 1))
  return(S_chao1)
})

alpha_div$ACE <- safe_apply(otu_df, 1, function(x) {
  if(sum(x) == 0) return(NA)
  x <- x[x > 0]
  S_obs <- length(x)
  # 简化的ACE估计
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

# 检查数据质量
cat("📊 Alpha多样性数据摘要:\n")
print(summary(alpha_div[, 3:ncol(alpha_div)]))

# ---- Alpha多样性绘图 ----
plot_alpha <- function(metric) {
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
  
  p <- ggplot(plot_data, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    labs(title = metric, x = "", y = metric) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    scale_fill_manual(values = c("AMI" = "#E64B35FF", "Control" = "#4DBBD5FF"))
  
  filename <- paste0(metric, "_alpha.png")
  ggsave(filename, p, width = 6, height = 5, dpi = 300)
  cat("💾 保存图:", filename, "\n")
}

# 绘制所有指标
metrics <- c("Richness", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Coverage")
for(m in metrics) {
  plot_alpha(m)
}

# ---- Beta 多样性分析 ----
cat("🔬 开始Beta多样性分析...\n")

# 筛选有效样本（至少有一个非零物种）
valid_samples <- rowSums(otu_df) > 0
otu_filtered <- otu_df[valid_samples, ]
Group_filtered <- Group[valid_samples]

cat("📊 有效样本数:", nrow(otu_filtered), "/", nrow(otu_df), "\n")

if(nrow(otu_filtered) < 3) {
  cat("⚠️ 有效样本不足，跳过Beta多样性分析\n")
} else {
  # 计算Bray-Curtis距离
  dist_bc <- vegdist(otu_filtered, method = "bray")
  
  # PCoA分析
  pcoa_res <- cmdscale(dist_bc, eig = TRUE, k = 2)
  
  # 计算解释度
  eig_pct <- round(100 * pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]), 2)
  
  # 创建绘图数据
  pcoa_df <- data.frame(
    SampleID = rownames(otu_filtered),
    PC1 = pcoa_res$points[, 1],
    PC2 = pcoa_res$points[, 2],
    Group = Group_filtered
  )
  
  # 绘制PCoA图
  p_pcoa <- ggplot(pcoa_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 2) +
    theme_bw(base_size = 14) +
    labs(
      title = "PCoA (Bray-Curtis)",
      x = paste0("PC1 (", eig_pct[1], "%)"),
      y = paste0("PC2 (", eig_pct[2], "%)")
    ) +
    scale_color_manual(values = c("AMI" = "#E64B35FF", "Control" = "#4DBBD5FF")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave("Beta_PCoA.png", p_pcoa, width = 6, height = 5, dpi = 300)
  cat("✅ Beta PCoA 图已保存\n")
  
  # PERMANOVA分析
  if(length(unique(Group_filtered)) > 1) {
    permanova <- adonis2(dist_bc ~ Group_filtered)
    cat("📊 PERMANOVA 结果:\n")
    print(permanova)
    write.csv(as.data.frame(permanova), "PERMANOVA_results.csv")
  }
}

cat("🎉 全部分析完成！\n")

# 最终总结
cat("\n📈 分析总结:\n")
cat("✅ Alpha多样性: 计算了", length(metrics), "个指标\n")
cat("✅ 样本数量: AMI =", sum(Group == "AMI"), ", Control =", sum(Group == "Control"), "\n")
cat("✅ 结果文件: Alpha_diversity_results.csv + 多样性图表\n")
if(exists("permanova")) {
  cat("✅ Beta多样性: PCoA图 + PERMANOVA分析\n")
} else {
  cat("⚠️ Beta多样性: 仅生成PCoA图\n")
}