# ============================
# ANCOM-BC2 family水平分析（四文件输出，保留原始 family 表全部行）
# ============================

library(phyloseq)
library(ANCOMBC)
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# ----------------------------
# 1. 文件路径
# ----------------------------
microbe_files <- list(
  archaea = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/family_data/archaea_family.xlsx",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/family_data/bacteria_family.xlsx",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/family_data/fungi_family.xlsx",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/family_data/virus_family_no_HF.xlsx"
)

meta_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/sample_metadata.xlsx"

# 输出文件夹
output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/ancombc2_results_family/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. 读取元数据
# ----------------------------
metadata <- read.xlsx(meta_file)
rownames(metadata) <- metadata$ID  # 样本 ID 对应行名

# ----------------------------
# 3. 循环分析四类微生物
# ----------------------------
for (microbe in names(microbe_files)) {
  
  cat("正在分析:", microbe, "\n")
  
  # 读取 family 水平表
  feature_data <- read.xlsx(microbe_files[[microbe]])
  rownames(feature_data) <- feature_data$Family
  feature_data_only <- feature_data[, -1]  # 删除 Family 列，只保留丰度
  
  # 确保样本顺序一致
  common_samples <- intersect(colnames(feature_data_only), rownames(metadata))
  feature_data_filtered <- feature_data_only[, common_samples, drop = FALSE]
  metadata_filtered <- metadata[common_samples, , drop = FALSE]
  
  # 构建 phyloseq 对象
  otu <- otu_table(as.matrix(feature_data_filtered), taxa_are_rows = TRUE)
  sam <- sample_data(metadata_filtered)
  ps_obj <- phyloseq(otu, sam)
  
  # ----------------------------
  # 运行 ANCOM-BC2
  # ----------------------------
  ancombc_res <- tryCatch({
    ancombc2(
      data = ps_obj,
      fix_formula = "Group",
      p_adj_method = "fdr",
      lib_cut = 0,
      group = "Group",
      struc_zero = FALSE,  # 不处理结构零
      neg_lb = TRUE,
      alpha = 0.05,
      n_cl = 6
    )
  }, error = function(e) {
    cat(paste0("ANCOM-BC2分析出错: ", microbe, " - ", e$message, "\n"))
    return(NULL)
  })
  
  # ----------------------------
  # 合并结果并保留原始 Family
  # ----------------------------
  if (!is.null(ancombc_res)) {
    res <- ancombc_res$res
    
    # 用 taxon 列创建 Family 列，确保合并正确
    res$Family <- res$taxon
    
    # 合并原始表和 ANCOM-BC2 结果，保留所有行
    final_res <- merge(feature_data, res, by = "Family", all.x = TRUE)
    
    # ----------------------------
    # 保存为 Excel
    # ----------------------------
    output_file <- paste0(output_dir, microbe, "_ANCOMBC2_results.xlsx")
    wb <- createWorkbook()
    addWorksheet(wb, "Complete_Results")
    writeData(wb, "Complete_Results", final_res)
    saveWorkbook(wb, output_file, overwrite = TRUE)
    
    cat(paste0("ANCOM-BC2完整结果已保存: ", output_file, "\n\n"))
    
  } else {
    cat(paste0("ANCOM-BC2分析结果为空: ", microbe, "\n\n"))
  }
}
