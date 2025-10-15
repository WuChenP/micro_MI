# ======================================================
# 将“心衰心梗对应解释.xlsx”中的 Merge 信息
# 映射到 ANCOM-BC2 病毒结果（virus_Result）中
# 按 taxon ↔ vOTUs 匹配
# ======================================================

library(openxlsx)
library(readxl)
library(dplyr)

# ----------------------------
# 文件路径
# ----------------------------
result_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"
votu_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/MaAsLin2_OTU/virus/心衰心梗对应解释.xlsx"

# ----------------------------
# 读取文件
# ----------------------------
votu <- read_excel(votu_file)
if (!all(c("vOTUs", "Merge") %in% colnames(votu))) {
  stop("❌ 心衰心梗对应解释文件缺少 'vOTUs' 或 'Merge' 列")
}

# 载入 ANCOM-BC2 结果文件
sheet_names <- getSheetNames(result_file)
if (!("virus_Result" %in% sheet_names)) stop("❌ 未找到 virus_Result sheet")

ancom_res <- read.xlsx(result_file, sheet = "virus_Result")
if (!("taxon" %in% colnames(ancom_res))) {
  colnames(ancom_res)[1] <- "taxon"
}

# ----------------------------
# 匹配与合并
# ----------------------------
merged_df <- ancom_res %>%
  left_join(votu %>% select(vOTUs, Merge), by = c("taxon" = "vOTUs")) %>%
  rename(解释信息 = Merge)

cat("✅ 匹配完成：virus_Result 中共匹配到", sum(!is.na(merged_df$解释信息)), "条解释信息\n")

# ----------------------------
# 写回 Excel（覆盖原 sheet）
# ----------------------------
wb_out <- loadWorkbook(result_file)
removeWorksheet(wb_out, "virus_Result")
addWorksheet(wb_out, "virus_Result")
writeData(wb_out, "virus_Result", merged_df)

# 保存文件
saveWorkbook(wb_out, result_file, overwrite = TRUE)
cat("🎉 已更新文件：", result_file, "\n")





# ======================================================
# 病毒相对丰度表 - 1% 流行率过滤 + CSV 输出
# 清理 ID 列中的不可见字符
# ======================================================

library(openxlsx)

# 输入文件路径（只保留病毒）
virus_file <- "E:/Python/MI_Analysis/origin_data/心梗组_病毒.xlsx"

# 输出文件夹
out_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 流行率阈值（1%）
prev_cutoff <- 0.01

# ----------------------------
# 函数：清理不可见字符 + 流行率过滤 + CSV 输出
# ----------------------------
filter_prevalence_csv <- function(file, cutoff = 0.01, out_dir) {
  # 读取 Excel
  df <- read.xlsx(file)
  
  # 清理 ID 列：去掉不可见字符并确保唯一
  df$ID <- make.unique(iconv(as.character(df$ID), from = "UTF-8", to = "UTF-8", sub = ""))
  
  # 提取物种 ID 和丰度矩阵
  taxa <- df$ID
  abund <- df[, -1]
  
  # 强制转换为数值
  abund[] <- lapply(abund, function(x) as.numeric(as.character(x)))
  
  # 计算每个物种的流行率
  prevalence <- apply(abund, 1, function(x) mean(x > 0))
  
  # 过滤低流行率物种
  keep <- prevalence >= cutoff
  df_filtered <- df[keep, ]
  
  # 构造输出 CSV 文件名
  fname <- basename(file)
  out_file <- file.path(out_dir, gsub(".xlsx", "_filtered_1percent.csv", fname))
  
  # 保存 CSV
  write.csv(df_filtered, out_file, row.names = FALSE, quote = FALSE)
  
  cat("✅ 已完成 1% 流行率过滤并生成 CSV:", fname, "=> 保留", nrow(df_filtered), "个物种\n")
}

# ----------------------------
# 执行处理（仅病毒）
# ----------------------------
filter_prevalence_csv(virus_file, cutoff = prev_cutoff, out_dir = out_dir)









# =============================================================================
# 病毒 OTU 表 ANCOM-BC2 分析（按组删除零方差 OTU，避免报错）
# =============================================================================

library(phyloseq)
library(ANCOMBC)
library(openxlsx)
library(dplyr)

# ----------------------------
# 输入文件路径（仅病毒）
# ----------------------------
virus_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/心梗组_病毒_filtered_1percent.csv"

metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
output_file   <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/ancombc2_results_OTU/virus_ANCOMBC2_results.xlsx"

if(!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

# ----------------------------
# 读取元数据
# ----------------------------
metadata <- read.xlsx(metadata_file)
rownames(metadata) <- metadata$SampleID

# ----------------------------
# 创建 Excel 工作簿
# ----------------------------
wb <- createWorkbook()

# ----------------------------
# 开始分析病毒
# ----------------------------
cat("\n开始分析: virus\n")

# 读取 OTU 表
feature_data <- read.csv(virus_file, check.names = FALSE)
rownames(feature_data) <- feature_data$ID
feature_data <- feature_data[, -1]  # 去掉 ID 列

# 对齐样本
common_samples <- intersect(colnames(feature_data), metadata$SampleID)
feature_data <- feature_data[, common_samples, drop = FALSE]
metadata_filtered <- metadata[match(common_samples, metadata$SampleID), , drop = FALSE]
rownames(metadata_filtered) <- metadata_filtered$SampleID

# ----------------------------
# 按组删除零方差 OTU
# ----------------------------
#groups <- unique(metadata_filtered$Group)
#keep_idx <- rep(TRUE, nrow(feature_data))
# for(g in groups){
#  idx <- metadata_filtered$Group == g
#  group_var <- apply(feature_data[, idx, drop = FALSE], 1, var, na.rm = TRUE)
#  keep_idx <- keep_idx & (group_var != 0)
#}
#feature_data <- feature_data[keep_idx, , drop = FALSE]
#cat("保留 OTU 数量:", nrow(feature_data), "\n")

# 构建 phyloseq 对象
otu <- otu_table(as.matrix(feature_data), taxa_are_rows = TRUE)
sam <- sample_data(metadata_filtered)
ps_obj <- phyloseq(otu, sam)

# ANCOM-BC2 分析
ancombc_res <- tryCatch({
  ancombc2(
    data = ps_obj,
    fix_formula = "Group",
    group = "Group",
    p_adj_method = "BH",
    lib_cut = 0,
    prv_cut = 0,
    struc_zero = FALSE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 6
  )
}, error = function(e) {
  cat("ANCOM-BC2 出错:", e$message, "\n")
  return(NULL)
})

# 保存结果
if (!is.null(ancombc_res)) {
  addWorksheet(wb, "virus_Result")
  writeData(wb, "virus_Result", ancombc_res$res, rowNames = FALSE)
  cat("virus 结果已写入 Excel。\n")
} else {
  cat("virus 结果为空，未写入。\n")
}

# ----------------------------
# 保存 Excel 文件
# ----------------------------
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("\n病毒 ANCOM-BC2 原始结果已保存到:", output_file, "\n")







# ==========================================================
# 从 ANCOM-BC2 结果中筛选显著病毒（q < 0.05 且 |log2FC| > 1）
# ==========================================================

library(openxlsx)
library(dplyr)

# ----------------------------
# 输入文件
# ----------------------------
input_file  <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/ancombc2_results_OTU/virus_ANCOMBC2_results.xlsx"
output_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/ancombc2_results_OTU/virus_sig_filtered.xlsx"

# ----------------------------
# 读取结果
# ----------------------------
df <- read.xlsx(input_file, sheet = 1)

# 检查列名
required_cols <- c("q_GroupMI", "lfc_GroupMI")
missing_cols <- setdiff(required_cols, colnames(df))
if(length(missing_cols) > 0){
  stop(paste("缺少列:", paste(missing_cols, collapse = ", ")))
}

# ----------------------------
# 筛选显著病毒
# ----------------------------
df_sig <- df %>%
  filter(q_GroupMI < 0.05 & abs(lfc_GroupMI) > 1)

cat("原始 OTU 数:", nrow(df), "\n")
cat("显著 OTU 数:", nrow(df_sig), "\n")

# ----------------------------
# 保存到新的 Excel 文件
# ----------------------------
wb <- createWorkbook()
addWorksheet(wb, "virus_sig")
writeData(wb, "virus_sig", df_sig, rowNames = FALSE)
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("✅ 已保存显著病毒到:", output_file, "\n")








# ======================================================
# 筛选 ANCOM-BC2 病毒结果
# - 根据富集方向和解释信息关键词分组
# - q值 < 0.05
# - |log2FC| > 1
# ======================================================

library(dplyr)
library(stringr)
library(openxlsx)

# ----------------------------
# 文件路径
# ----------------------------
result_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/ancombc2_results_OTU/virus_sig_filtered.xlsx"

# ----------------------------
# 读取 virus_Result 工作表
# ----------------------------
df <- read.xlsx(result_file, sheet = "virus_sig")

# ----------------------------
# 关键词列表
# ----------------------------
keywords_control <- c(
  "O__Unclassifed","f__Uniclassified","g__Lambdavirus","g__Hubeivirus","o__Petitvirales",
  "d__Monodnaviria","p__Phixviicota","k__Sangervirae","c__Malgrandaviricetes",
  "f__Allomimiviridae","O__Algavirales","g__Bembunaguatrovius","g__Aguilavirus",
  "f__Vilmaviridae","f__Peduoviridae","g__Alegriavirus","S__Buchavirus_coli",
  "g__Lidleunavirus","g__Muvirus","g__Brunovirus","g__Delepguintavirus","g__Lederbergvirus",
  "g__Svunavirus","f__Filixviridae","g__Glaedevirus"
)

keywords_MI <- c(
  "g__Mushuvirus","g__Spinunavirus","S__Serangoonvirus_essarone",
  "g__Serangoonvirus","g__Punavirus","S__Salacisavirus_pssm2","g__Salacisavirus"
)

# ----------------------------
# 筛选条件
# ----------------------------
filtered_df <- df %>%
  filter(
    q_GroupMI < 0.05,
    abs(lfc_GroupMI) > 1,
  ) %>%
  filter(
    (lfc_GroupMI < 0 & str_detect(taxon, str_c(keywords_control, collapse = "|"))) |
      (lfc_GroupMI > 0 & str_detect(taxon, str_c(keywords_MI, collapse = "|")))
  )

cat("✅ 共筛选出", nrow(filtered_df), "行符合条件\n")

# ----------------------------
# 输出结果
# ----------------------------
out_file <- sub("\\.xlsx$", "_筛选结果_方向关键词.xlsx", result_file)
write.xlsx(filtered_df, out_file, overwrite = TRUE)

cat("🎉 筛选结果已保存到：", out_file, "\n")








# =============================================================================
# 论文级 MaAsLin2 分析：只处理病毒
# 自动 CSV→TSV，筛选 q<0.05 & |log2FC|>1，标注升降趋势
# =============================================================================

rm(list = ls())

# ----------------------------
# 依赖包
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
# 文件路径
# ----------------------------
virus_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/心梗组_病毒_filtered_1percent.csv"
metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
outdir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/MaAsLin2_OTU/virus"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ----------------------------
# 读取 metadata 并保存为 TSV
# ----------------------------
metadata <- read_excel(metadata_file)
colnames(metadata)[1] <- "SampleID"
metadata_tsv <- file.path(outdir, "metadata_temp.tsv")
write.table(metadata, metadata_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------
# 读取病毒 CSV 并转换为 TSV 临时文件
# ----------------------------
feature <- read.csv(virus_file, check.names = FALSE)
feature_tsv <- file.path(outdir, "virus_temp.tsv")
write.table(feature, feature_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------
# 运行 MaAsLin2
# ----------------------------
Maaslin2(
  input_data = feature_tsv,
  input_metadata = metadata_tsv,
  output = outdir,
  fixed_effects = c("Group"),
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.05
)

# ----------------------------
# 读取显著结果并筛选 q<0.05 & |log2FC|>1
# ----------------------------
sig_file <- file.path(outdir, "significant_results.tsv")
if (file.exists(sig_file)) {
  sig_df <- read.delim(sig_file, check.names = FALSE)
  if (nrow(sig_df) > 0) {
    sig_df$log2FC <- sig_df$coef
    sig_df <- sig_df %>% filter(qval < 0.05 & abs(log2FC) > 1)
    if (nrow(sig_df) > 0) {
      sig_df$Trend <- ifelse(sig_df$log2FC > 0, "Up", "Down")
      sig_df$MicrobeClass <- "virus"
      out_file <- file.path(outdir, "virus_significant_q0.05_log2fc1.tsv")
      write.table(sig_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("✅ 病毒显著 feature 数量：", nrow(sig_df), "\n")
      cat("保存文件：", out_file, "\n")
    } else {
      cat("⚠️ 没有满足 q<0.05 & |log2FC|>1 的病毒 feature\n")
    }
  } else {
    cat("⚠️ 没有显著病毒 feature\n")
  }
} else {
  cat("⚠️ 显著结果文件不存在\n")
}

cat("\n🎉 病毒 MaAsLin2 分析完成，结果保存在：", outdir, "\n")








