# ======================================
# ✅ MaAsLin2：心梗 vs 正常人口腔微生物分析（修正版）
# ======================================

if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Maaslin2")
}
library(Maaslin2)

# ---- 设置工作路径 ----
setwd("D:/机器学习与生物信息学/北京朝阳医院项目/R/M")

# ---- 1. 读取数据 ----
df <- read.csv("R_data.csv", header = TRUE, row.names = 2, check.names = FALSE)

# ---- 2. 提取丰度矩阵 ----
otu_mat <- df[, grep("AMI|CON", colnames(df))]
otu_mat <- as.data.frame(otu_mat)

# 检查是否包含非数值
if (any(sapply(otu_mat, function(x) !is.numeric(x)))) {
  message("⚠️ 检测到非数值列，尝试强制转换为数值...")
  otu_mat[] <- lapply(otu_mat, function(x) as.numeric(as.character(x)))
}

# 去掉全为0的菌
otu_mat <- otu_mat[rowSums(otu_mat) > 0, ]
message("✅ 保留的特征数量: ", nrow(otu_mat))

# ---- 3. 构建 metadata ----
meta_df <- data.frame(
  sample = colnames(otu_mat),
  group = ifelse(grepl("^AMI", colnames(otu_mat)), "AMI", "CON")
)
rownames(meta_df) <- meta_df$sample

# ---- 4. 转置丰度矩阵（样本为行）----
feature_df <- as.data.frame(t(otu_mat))
feature_df <- feature_df[, colSums(feature_df) > 0, drop = FALSE]
feature_df <- as.data.frame(lapply(feature_df, as.numeric))
rownames(feature_df) <- colnames(otu_mat)

# ---- 5. 匹配样本 ----
feature_df <- feature_df[rownames(meta_df), , drop = FALSE]

# ---- 6. 检查数据完整性 ----
cat("\n=== 数据检查 ===\n")
cat("样本数:", nrow(feature_df), "\n")
cat("特征数:", ncol(feature_df), "\n")
cat("是否有NA:", any(is.na(feature_df)), "\n")

# 如果没有任何特征，直接停止
if (ncol(feature_df) == 0) {
  stop("❌ 没有有效特征列，请检查输入数据（可能丰度太低或全部为0）。")
}

# ---- 7. 运行 MaAsLin2 ----
fit_data <- Maaslin2(
  input_data = feature_df,
  input_metadata = meta_df,
  output = "Maaslin2_out",
  normalization = "NONE",        # 不做标准化
  transform = "LOG",             # 对数变换
  fixed_effects = c("group"),
  min_abundance = 0.00000001,    # 减少过滤
  min_prevalence = 0.00,         # 不过滤稀有物种
  max_significance = 0.25
)

cat("\n✅ MaAsLin2 分析完成，结果位于文件夹: Maaslin2_out\n")
