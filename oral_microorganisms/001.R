# ==================================
# ANCOM-BC2 分析
# ==================================

# ---- 设置 CRAN 镜像 ----
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# ---- 包检查与安装 ----
check_and_install <- function(pkg, github = FALSE, repo = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste(">>> Installing package:", pkg))
    if (github && !is.null(repo)) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_github(repo, upgrade = "never")
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ---- 安装并加载依赖包 ----
check_and_install("phyloseq")
check_and_install("microbiome", github = TRUE, repo = "microbiome/microbiome")
check_and_install("ANCOMBC", github = TRUE, repo = "FrederickHuangLin/ANCOMBC")
check_and_install("ggplot2")

# ---- 设置工作目录 ----
setwd("D:/机器学习与生物信息学/北京朝阳医院项目/R")  # CSV 文件路径

# ---- 1. 读取数据 ----
df <- read.csv(
  "R_data.csv",
  header = TRUE,
  row.names = 2
)  # 用 ASV ID 作为行名

# ---- 2. 构建丰度矩阵（样本 × ASV）----
otu_mat <- as.matrix(df[, grep("AMI|CON", colnames(df))])

# ---- 3. 提取分类学矩阵 ----
tax_mat <- as.matrix(
  df[, c("Kingdom", "Phylum", "Class",
         "Order", "Family", "Genus", "Species")]
)

# ---- 4. 构建 phyloseq 对象 ----
otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
tax <- tax_table(tax_mat)

sample_ids <- colnames(otu_mat)
group <- ifelse(grepl("^AMI", sample_ids), "AMI", "CON")

sample_df <- data.frame(
  group = group,
  row.names = sample_ids
)
samp <- sample_data(sample_df)

phy <- phyloseq(otu, tax, samp)

# ---- 5. 运行 ANCOM-BC2 ----
out <- ancombc2(
  data = phy,
  assay_name = "counts",
  fix_formula = "group",
  p_adj_method = "holm",
  alpha = 0.05
)

# ---- 调试：分组检查 ----
cat("\n=== 分组检查 ===\n")
print(table(sample_df$group))

# ---- 6. 提取并保存结果 ----
res <- out$res   # 新版 ANCOMBC2 输出为 data.frame

cat("\n=== ANCOMBC2 输出列 ===\n")
print(colnames(res))

# 合并 taxonomy 信息
res_df <- cbind(tax_mat[match(res$taxon, rownames(tax_mat)), ], res)

# 保存完整结果
write.csv(res_df, "ANCOMBC2_results_full.csv", row.names = FALSE)
message("✅ 完整结果已保存到 ANCOMBC2_results_full.csv")

# 保存显著差异结果 (q < 0.05 & diff = TRUE)
if ("q_groupCON" %in% colnames(res) && "diff_groupCON" %in% colnames(res)) {
  sig_df <- subset(res_df, q_groupCON < 0.05 & diff_groupCON == TRUE)
  write.csv(sig_df, "ANCOMBC2_results_sig.csv", row.names = FALSE)
  message("✅ 显著差异结果已保存到 ANCOMBC2_results_sig.csv")
} else {
  message("⚠️ 没有找到 q_groupCON 或 diff_groupCON 列，请检查输出。")
}