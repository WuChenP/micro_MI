# ======================================================
# ANCOM-BC2 物种水平分析（结果拼接回原表并保存）
# ✅ 不使用 tibble 行名，严格按 SampleID 对齐
# ======================================================

library(phyloseq)
library(ANCOMBC)
library(readxl)
library(openxlsx)
library(tibble) # rownames_to_column

# ----------------------------
# 1. 路径设置
# ----------------------------
abundance_file <- "E:/Python/MI_Analysis/origin_data/物种层级.xlsx"
meta_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
save_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/Species-level_difference_analysis_results"

if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ----------------------------
# 2. 读取数据
# ----------------------------
abundance <- read_xlsx(abundance_file)
meta <- read_xlsx(meta_file, col_types = c("text", "text"))

# ----------------------------
# 3. 确保 Taxon 列存在
# ----------------------------
if(!"Taxon" %in% colnames(abundance)){
  colnames(abundance)[1] <- "Taxon"
}

# ----------------------------
# 4. 对齐 OTU 列与 metadata
# ----------------------------
common_samples <- meta$SampleID[meta$SampleID %in% colnames(abundance)[-1]]

# 按 metadata 顺序排列 OTU 表
abundance <- abundance[, c("Taxon", common_samples)]
meta <- meta[match(common_samples, meta$SampleID), , drop = FALSE]

# ----------------------------
# 4a. 打印调试信息
# ----------------------------
cat("=== OTU 表列名（排完序后） ===\n")
print(colnames(abundance)[-1])
cat("\n=== metadata SampleID 顺序 ===\n")
print(meta$SampleID)

# ----------------------------
# 5. 构建 OTU 矩阵 & phyloseq 对象
# ----------------------------
otu_mat <- as.matrix(abundance[, -1, drop = FALSE])
otu_mat <- apply(otu_mat, 2, as.numeric)
rownames(otu_mat) <- abundance$Taxon

# phyloseq 构建时使用 data.frame，并设置 row.names = SampleID
SAMPLE <- sample_data(data.frame(meta, row.names = "SampleID"))
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
ps_obj <- phyloseq(OTU, SAMPLE)

# ----------------------------
# 6. 运行 ANCOM-BC2
# ----------------------------
res <- ancombc2(
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

# ----------------------------
# 7. 拼接结果回原表
# ----------------------------
res_df <- res$res
if(!"taxon" %in% colnames(res_df)){
  res_df <- rownames_to_column(res_df, var = "taxon")
}

merged_data <- cbind(abundance, res_df[match(abundance$Taxon, res_df$taxon), -1, drop = FALSE])

# ----------------------------
# 8. 保存结果
# ----------------------------
save_path <- file.path(save_dir, "Species_level_ANCOMBC2_results.xlsx")
write.xlsx(merged_data, save_path, rowNames = FALSE)

cat("✅ 分析完成！结果已保存至：", save_path, "\n")
