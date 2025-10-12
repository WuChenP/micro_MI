# =============================================================================
# 论文级 MaAsLin2 分析：病毒/古菌/细菌/真菌 family 水平
# 输出目录：E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/family_data/MaAsLin2_results
# =============================================================================

rm(list = ls())

# ----------------------------
# 1. 安装依赖及 MaAsLin2（如果没有）
# ----------------------------
if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  deps <- c("dplyr", "data.table", "ggplot2", "Rcpp", "tibble", "stringr", "reshape2", "readxl")
  for (pkg in deps) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  remotes::install_github("biobakery/Maaslin2")
}
library(Maaslin2)
library(readxl)
library(dplyr)

# ----------------------------
# 2. 文件路径
# ----------------------------
files <- list(
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/family_data/virus_family_no_HF.xlsx",
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/family_data/archaea_family.xlsx",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/family_data/bacteria_family.xlsx",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/family_data/fungi_family.xlsx"
)
metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/sample_metadata.xlsx"

# 输出目录
outdir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/MaAsLin2_family"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ----------------------------
# 3. 读取 metadata 并保存为 TSV
# ----------------------------
metadata <- read_excel(metadata_file)
colnames(metadata)[1] <- "SampleID"
metadata$SampleID <- toupper(trimws(metadata$SampleID))
metadata_tsv <- file.path(outdir, "metadata_temp.tsv")
write.table(metadata, metadata_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------
# 4. 批量分析 + 日志
# ----------------------------
all_sig_list <- list()  # 用于存储每类微生物的显著结果

for (microbe in names(files)) {
  cat("\n=========== 开始分析:", microbe, "===========\n")
  
  # 读取 Excel 文件
  feature <- as.data.frame(read_excel(files[[microbe]]))
  
  # 修改第一列列名，保留内容不变
  colnames(feature)[1] <- "Feature"
  
  # 去掉列名前后空格并统一大小写
  colnames(feature)[-1] <- toupper(trimws(colnames(feature)[-1]))
  
  # 只保留 metadata 中匹配的样本列
  common_samples <- intersect(colnames(feature)[-1], metadata$SampleID)
  if(length(common_samples) == 0) stop("❌ feature 表和 metadata 表没有匹配样本！")
  
  # 拼接第一列 Feature 与匹配的样本列
  feature_to_save <- feature[, c("Feature", common_samples)]
  
  cat("feature 行数:", nrow(feature_to_save), "列数:", ncol(feature_to_save), "\n")
  
  # 保存临时 TSV
  feature_tsv <- file.path(outdir, paste0(microbe, "_temp.tsv"))
  write.table(feature_to_save, feature_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 创建输出子目录
  outdir_microbe <- file.path(outdir, microbe)
  if (!dir.exists(outdir_microbe)) dir.create(outdir_microbe)
  
  # 运行 MaAsLin2
  Maaslin2(
    input_data = feature_tsv,
    input_metadata = metadata_tsv,
    output = outdir_microbe,
    fixed_effects = c("Group"),
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 0.05
  )
  
  # 读取显著结果并标注升降趋势
  sig_file <- file.path(outdir_microbe, "significant_results.tsv")
  if (file.exists(sig_file)) {
    sig_df <- read.delim(sig_file, check.names = FALSE)
    if (nrow(sig_df) > 0) {
      sig_df$Trend <- ifelse(sig_df$coef > 0, "Up", "Down")
      sig_df$MicrobeClass <- microbe
      all_sig_list[[microbe]] <- sig_df
      cat("✅ ", microbe, "显著 feature 数量：", nrow(sig_df), "\n")
    } else {
      cat("⚠️ ", microbe, "没有显著 feature\n")
    }
  } else {
    cat("⚠️ ", microbe, "显著结果文件不存在\n")
  }
}

# ----------------------------
# 5. 合并四类微生物显著结果并保存
# ----------------------------
if (length(all_sig_list) > 0) {
  combined_sig <- bind_rows(all_sig_list)
  combined_file <- file.path(outdir, "combined_significant_results.tsv")
  write.table(combined_sig, combined_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\n🎉 合并后的显著结果保存在：", combined_file, "\n")
} else {
  cat("\n⚠️ 所有微生物类没有显著 feature，未生成合并表\n")
}

cat("\n🎉 所有分析完成！结果保存在：", outdir, "\n")
