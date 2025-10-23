# =============================================================================
# 论文级 MaAsLin2 分析：病毒/古菌/细菌/真菌
# 自动 CSV→TSV、Excel→TSV，筛选 q<0.05 & |log2FC|>1，合并显著结果并标注升降趋势
# 修复特征名中的括号问题
# =============================================================================

rm(list = ls())

# ----------------------------
# 1. 安装依赖及 MaAsLin2（如果没有）
# ----------------------------
if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  deps <- c("dplyr", "data.table", "ggplot2", "Rcpp", "tibble", "stringr", "reshape2")
  for (pkg in deps) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  remotes::install_github("biobakery/Maaslin2")
}
library(Maaslin2)
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl)
library(dplyr)

# ----------------------------
# 2. 定义修复函数
# ----------------------------
repair_maaslin_output <- function(output_dir, original_features) {
  # 修复 significant_results.tsv
  sig_file <- file.path(output_dir, "significant_results.tsv")
  if (file.exists(sig_file)) {
    sig_df <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(sig_df) > 0) {
      # 创建特征名映射
      feature_map <- data.frame(
        cleaned = make.names(original_features, unique = TRUE),
        original = as.character(original_features),
        stringsAsFactors = FALSE
      )
      
      # 修复特征名
      sig_df <- sig_df %>%
        left_join(feature_map, by = c("feature" = "cleaned")) %>%
        mutate(feature = ifelse(!is.na(original), original, feature)) %>%
        select(-original)
      
      # 保存修复后的文件
      write.table(sig_df, sig_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("✅ 已修复", nrow(sig_df), "个特征名在", basename(sig_file), "\n")
    }
  }
  
  # 修复 all_results.tsv
  all_file <- file.path(output_dir, "all_results.tsv")
  if (file.exists(all_file)) {
    all_df <- read.delim(all_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(all_df) > 0) {
      # 创建特征名映射
      feature_map <- data.frame(
        cleaned = make.names(original_features, unique = TRUE),
        original = as.character(original_features),
        stringsAsFactors = FALSE
      )
      
      # 修复特征名
      all_df <- all_df %>%
        left_join(feature_map, by = c("feature" = "cleaned")) %>%
        mutate(feature = ifelse(!is.na(original), original, feature)) %>%
        select(-original)
      
      # 保存修复后的文件
      write.table(all_df, all_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("✅ 已修复", nrow(all_df), "个特征名在", basename(all_file), "\n")
    }
  }
}

# ----------------------------
# 3. 文件路径
# ----------------------------
files <- list(
  virus = "E:/Python/MI_Analysis/metagenome/data_figures/g/属水平.xlsx"
)
metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/g/sample_metadata.xlsx"
outdir <- "E:/Python/MI_Analysis/metagenome/data_figures/g/MaAsLin2_OTU_new/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ----------------------------
# 4. 读取 metadata 并保存为 TSV
# ----------------------------
metadata <- read_excel(metadata_file)
colnames(metadata)[1] <- "SampleID"
metadata_tsv <- file.path(outdir, "metadata_temp.tsv")
write.table(metadata, metadata_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------
# 5. 批量分析 + CSV → TSV 自动转换 + 筛选显著
# ----------------------------
all_sig_list <- list()  # 用于存储每类微生物的显著结果

for (microbe in names(files)) {
  cat("\n=========== 开始分析:", microbe, "===========\n")
  
  # 读取 Excel 文件
  feature <- read_excel(files[[microbe]])
  
  # 保存原始特征名（假设第一列是特征名）
  if (!colnames(feature)[1] %in% c("Feature", "OTU", "Taxonomy")) {
    colnames(feature)[1] <- "Feature"
  }
  original_features <- feature[[1]]
  
  # 转换为 TSV 临时文件，禁用列名检查以保持原列名
  feature_tsv <- file.path(outdir, paste0(microbe, "_temp.tsv"))
  write.table(feature, feature_tsv, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # 创建输出目录
  outdir_microbe <- file.path(outdir, microbe)
  if (!dir.exists(outdir_microbe)) dir.create(outdir_microbe)
  
  # 运行 MaAsLin2（论文建议 q ≤ 0.05）
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
  
  # ----------------------------
  # 修复 MaAsLin2 输出文件中的特征名
  # ----------------------------
  cat("🔧 修复特征名中的括号问题...\n")
  repair_maaslin_output(outdir_microbe, original_features)
  
  # ----------------------------
  # 读取显著结果并筛选 q<0.05 & |log2FC|>1
  # ----------------------------
  sig_file <- file.path(outdir_microbe, "significant_results.tsv")
  if (file.exists(sig_file)) {
    sig_df <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(sig_df) > 0) {
      # 使用 coef 近似 log2FC（如果 MaAsLin2 已经输出 log2FC，可改用 log2FC 列）
      sig_df$log2FC <- sig_df$coef
      
      # 筛选条件
      sig_df <- sig_df %>% filter(qval < 0.05 & abs(log2FC) > 1)
      
      if (nrow(sig_df) > 0) {
        # 标注升降趋势
        sig_df$Trend <- ifelse(sig_df$log2FC > 0, "Up", "Down")
        sig_df$MicrobeClass <- microbe
        all_sig_list[[microbe]] <- sig_df
        
        # 显示修复后的特征名
        cat("✅ ", microbe, "显著 feature 数量（q<0.05 & |log2FC|>1）：", nrow(sig_df), "\n")
        cat("📝 显著特征名：", paste(sig_df$feature, collapse = ", "), "\n")
      } else {
        cat("⚠️ ", microbe, "没有满足 q<0.05 & |log2FC|>1 的 feature\n")
      }
    } else {
      cat("⚠️ ", microbe, "没有显著 feature\n")
    }
  } else {
    cat("⚠️ ", microbe, "显著结果文件不存在\n")
  }
}

# ----------------------------
# 6. 合并四类微生物显著结果并保存
# ----------------------------
if (length(all_sig_list) > 0) {
  combined_sig <- bind_rows(all_sig_list)
  combined_file <- file.path(outdir, "combined_significant_results_q0.05_log2fc1.tsv")
  write.table(combined_sig, combined_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\n🎉 合并后的显著结果保存在：", combined_file, "\n")
  cat("📊 总显著特征数量：", nrow(combined_sig), "\n")
} else {
  cat("\n⚠️ 所有微生物类没有显著 feature，未生成合并表\n")
}

cat("\n🎉 所有分析完成！结果保存在：", outdir, "\n")