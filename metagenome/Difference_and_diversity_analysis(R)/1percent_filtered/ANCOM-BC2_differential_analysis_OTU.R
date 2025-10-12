# =============================================================================
# 四类微生物 OTU 表 ANCOM-BC2 分析（按组删除零方差 OTU，避免报错）
# =============================================================================

library(phyloseq)
library(ANCOMBC)
library(openxlsx)
library(dplyr)

# ----------------------------
# 输入文件路径
# ----------------------------
otu_files <- list(
  archaea = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/心梗组_古菌_filtered_1percent.csv",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/心梗组_细菌_filtered_1percent.csv",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/心梗组_真菌_filtered_1percent.csv",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/心梗组_病毒(新)_filtered_1percent.csv"
)

metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
output_file   <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"

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
# 循环分析四类微生物
# ----------------------------
for (microbe in names(otu_files)) {
  cat("\n开始分析:", microbe, "\n")
  
  # 读取 OTU 表
  feature_data <- read.csv(otu_files[[microbe]], check.names = FALSE)
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
  groups <- unique(metadata_filtered$Group)
  keep_idx <- rep(TRUE, nrow(feature_data))
  for(g in groups){
    idx <- metadata_filtered$Group == g
    group_var <- apply(feature_data[, idx, drop = FALSE], 1, var, na.rm = TRUE)
    keep_idx <- keep_idx & (group_var != 0)
  }
  feature_data <- feature_data[keep_idx, , drop = FALSE]
  cat("保留 OTU 数量:", nrow(feature_data), "\n")
  
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
    addWorksheet(wb, paste0(microbe, "_Result"))
    writeData(wb, paste0(microbe, "_Result"), ancombc_res$res, rowNames = FALSE)
    cat(microbe, "结果已写入 Excel。\n")
  } else {
    cat(microbe, "结果为空，未写入。\n")
  }
}

# ----------------------------
# 保存 Excel 文件
# ----------------------------
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("\n四类微生物 ANCOM-BC2 原始结果已保存到:", output_file, "\n")
