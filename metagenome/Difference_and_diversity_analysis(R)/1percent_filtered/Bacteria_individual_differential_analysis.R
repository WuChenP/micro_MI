# =============================================================================
# 细菌 OTU 表 ANCOM-BC2 分析（结果拼接回原表并保存，不删除零方差）
# =============================================================================

library(phyloseq)
library(ANCOMBC)
library(openxlsx)
library(dplyr)

# ----------------------------
# 输入文件路径
# ----------------------------
bacteria_otu_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/心梗组_细菌_filtered_1percent.csv"
metadata_file     <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
save_dir          <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/Bacteria_individual_differential_analysis/ANCOM-BC2_OTU_results"

if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ----------------------------
# 读取元数据
# ----------------------------
metadata <- read.xlsx(metadata_file)
meta <- as.data.frame(metadata)
rownames(meta) <- meta$SampleID

# ----------------------------
# 读取 OTU 表
# ----------------------------
abundance <- read.csv(bacteria_otu_file, check.names = FALSE)
abundance$taxon <- abundance$ID  # 将 ID 列改名为 taxon
abundance <- abundance[, c("taxon", setdiff(names(abundance), c("ID", "taxon")))]

# ----------------------------
# 对齐样本顺序
# ----------------------------
common_samples <- intersect(colnames(abundance)[-1], rownames(meta))
abundance <- abundance[, c("taxon", common_samples)]
meta <- meta[common_samples, , drop = FALSE]

# 打印构建 phyloseq 前的列名和行名
cat("=== OTU 表列名（排好顺序） ===\n")
print(colnames(abundance)[-1])
cat("\n=== metadata 行名（排好顺序） ===\n")
print(rownames(meta))

# ----------------------------
# 构建 phyloseq 对象
# ----------------------------
otu_mat <- as.matrix(abundance[, -1])
otu_mat <- apply(otu_mat, 2, as.numeric)
rownames(otu_mat) <- abundance$taxon

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
SAMPLE <- sample_data(meta)
ps_obj <- phyloseq(OTU, SAMPLE)

# ----------------------------
# ANCOM-BC2 分析
# ----------------------------
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

# ----------------------------
# 拼接结果到原丰度表
# ----------------------------
if (!is.null(ancombc_res)) {
  res <- ancombc_res$res

  # 合并原始表和 ANCOM-BC2 结果
  final_res <- merge(abundance, res, by.x = "taxon", by.y = "taxon", all.x = FALSE)

  # 删除重复 taxon 列（如果存在）
  if ("taxon" %in% colnames(final_res)[(ncol(abundance)+1):ncol(final_res)]) {
    final_res$taxon <- NULL
  }

  # ----------------------------
  # 保存结果
  # ----------------------------
  save_path <- file.path(save_dir, "Bacteria_ANCOMBC2_results.xlsx")
  write.xlsx(final_res, save_path, rowNames = FALSE)

  cat("✅ 分析完成！结果已保存至：", save_path, "\n")
} else {
  cat("❌ ANCOM-BC2 分析结果为空，未生成输出文件。\n")
}
