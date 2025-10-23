# ======================================================
# ANCOM-BC2 多层级分析（改进版）
# ======================================================

library(phyloseq)
library(ANCOMBC)
library(readxl)
library(openxlsx)
library(tibble)

# ----------------------------
# 1. 路径设置
# ----------------------------
# 定义四个要处理的文件路径
abundance_files <- c(
  "E:/Python/MI_Analysis/origin_data/层级reads丰度-心衰心梗/virus.class.profile.xlsx",
  "E:/Python/MI_Analysis/origin_data/层级reads丰度-心衰心梗/virus.family.profile.xlsx", 
  "E:/Python/MI_Analysis/origin_data/层级reads丰度-心衰心梗/virus.order.profile.xlsx",
  "E:/Python/MI_Analysis/origin_data/层级reads丰度-心衰心梗/virus.phylum.profile.xlsx"
)

# 定义对应的层级名称
level_names <- c("Class", "Family", "Order", "Phylum")

meta_file <- "E:/Python/MI_Analysis/origin_data/层级reads丰度-心衰心梗/sample_metadata.xlsx"
save_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/Species-level_difference_analysis_results_new"

if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ----------------------------
# 2. 读取元数据
# ----------------------------
meta <- read_xlsx(meta_file, col_types = c("text", "text"))
cat("元数据样本数量:", nrow(meta), "\n")

# 用于存储所有结果的列表
all_results <- list()

# ----------------------------
# 3. 循环处理每个文件
# ----------------------------
for(i in 1:length(abundance_files)) {
  current_file <- abundance_files[i]
  current_level <- level_names[i]
  
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("正在处理:", current_level, "层级\n")
  cat("文件:", basename(current_file), "\n")
  cat(rep("=", 50), "\n", sep = "")
  
  # ----------------------------
  # 3.1 读取丰度数据
  # ----------------------------
  tryCatch({
    abundance <- read_xlsx(current_file)
    
    # ----------------------------
    # 3.2 确保 Taxon 列存在
    # ----------------------------
    if(!"Taxon" %in% colnames(abundance)){
      colnames(abundance)[1] <- "Taxon"
    }
    
    # ----------------------------
    # 3.3 数据预处理和样本对齐
    # ----------------------------
    # 移除可能存在的空行
    abundance <- abundance[!is.na(abundance$Taxon) & abundance$Taxon != "", ]
    
    # 获取共同样本
    common_samples <- intersect(meta$SampleID, colnames(abundance)[-1])
    
    if(length(common_samples) == 0) {
      cat("❌ 警告：", current_level, "层级没有找到共同的样本！跳过此文件。\n")
      next
    }
    
    cat("共同样本数量:", length(common_samples), "\n")
    
    # 按 metadata 顺序重排
    abundance <- abundance[, c("Taxon", common_samples)]
    meta_sub <- meta[match(common_samples, meta$SampleID), ]
    
    # ----------------------------
    # 3.4 严格样本对齐验证
    # ----------------------------
    cat("\n=== 样本对齐验证 ===\n")
    cat("OTU表样本数:", length(colnames(abundance)[-1]), "\n")
    cat("metadata样本数:", nrow(meta_sub), "\n")
    cat("样本顺序一致:", identical(colnames(abundance)[-1], meta_sub$SampleID), "\n")
    
    # ----------------------------
    # 3.5 构建 phyloseq 对象
    # ----------------------------
    otu_mat <- as.matrix(abundance[, -1, drop = FALSE])
    storage.mode(otu_mat) <- "numeric"  # 更安全的数值转换
    rownames(otu_mat) <- abundance$Taxon
    
    # 移除全为零的物种（可选）
    zero_taxa <- rowSums(otu_mat) == 0
    if(any(zero_taxa)) {
      cat("移除全为零的物种数量:", sum(zero_taxa), "\n")
      otu_mat <- otu_mat[!zero_taxa, ]
      # 同时更新abundance数据框
      abundance <- abundance[!zero_taxa, ]
    }
    
    cat("有效物种数量:", nrow(otu_mat), "\n")
    
    SAMPLE <- sample_data(data.frame(meta_sub, row.names = "SampleID"))
    OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
    ps_obj <- phyloseq(OTU, SAMPLE)
    
    cat("最终phyloseq对象:\n")
    print(ps_obj)
    
    # ----------------------------
    # 3.6 运行 ANCOM-BC2
    # ----------------------------
    cat("正在运行ANCOM-BC2分析...\n")
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
    # 3.7 处理结果
    # ----------------------------
    res_df <- res$res
    
    # 确保有taxon列
    if(!"taxon" %in% colnames(res_df)){
      res_df <- rownames_to_column(res_df, var = "taxon")
    }
    
    # 安全合并
    matched_indices <- match(abundance$Taxon, res_df$taxon)
    cat("成功匹配的物种:", sum(!is.na(matched_indices)), "/", nrow(abundance), "\n")
    
    if(sum(!is.na(matched_indices)) == 0) {
      cat("❌ 警告：没有物种匹配成功，跳过此层级\n")
      next
    }
    
    merged_data <- cbind(abundance, res_df[matched_indices, -1, drop = FALSE])
    
    # 添加层级信息列
    merged_data$Level <- current_level
    
    # 存储结果
    all_results[[current_level]] <- merged_data
    
    cat("✅", current_level, "层级分析完成！\n")
    
  }, error = function(e) {
    cat("❌", current_level, "层级分析出错:", e$message, "\n")
  })
}

# ----------------------------
# 4. 保存所有结果到一个Excel文件
# ----------------------------
if(length(all_results) > 0) {
  save_path <- file.path(save_dir, "Virus_MultiLevel_ANCOMBC2_Results.xlsx")
  
  # 创建工作簿
  wb <- createWorkbook()
  
  # 为每个层级添加工作表
  for(level_name in names(all_results)) {
    addWorksheet(wb, level_name)
    writeData(wb, level_name, all_results[[level_name]], rowNames = FALSE)
  }
  
  # 保存工作簿
  saveWorkbook(wb, save_path, overwrite = TRUE)
  
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("🎉 所有分析完成！\n")
  cat("结果已保存至：", save_path, "\n")
  cat("包含的工作表：", paste(names(all_results), collapse = ", "), "\n")
  cat("总处理层级数：", length(all_results), "/", length(abundance_files), "\n")
  cat(rep("=", 50), "\n", sep = "")
  
} else {
  cat("\n❌ 没有成功生成任何结果，请检查输入文件和数据。\n")
}