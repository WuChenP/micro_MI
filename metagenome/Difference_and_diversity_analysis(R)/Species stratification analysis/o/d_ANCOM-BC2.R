# =============================================================================
# 四类微生物 OTU 表 ANCOM-BC2 分析（包含原始丰度数据）
# =============================================================================

library(phyloseq)
library(ANCOMBC)
library(openxlsx)

# ----------------------------
# 1. 文件路径
# ----------------------------
otu_files <- list(
  virus = "E:/Python/MI_Analysis/metagenome/data_figures/o/目水平.xlsx"
)

metadata_file = "E:/Python/MI_Analysis/metagenome/data_figures/o/sample_metadata.xlsx"

# 输出文件夹
output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/o/ancombc2_results_new/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. 读取元数据
# ----------------------------
metadata <- read.xlsx(metadata_file)
rownames(metadata) <- metadata$SampleID

# ----------------------------
# 循环分析四类微生物
# ----------------------------
for (microbe in names(otu_files)) {
  
  cat("\n正在分析:", microbe, "\n")
  
  # 读取 OTU 表
  feature_data <- read.xlsx(otu_files[[microbe]], check.names = FALSE)
  rownames(feature_data) <- feature_data$目
  feature_data_only <- feature_data[, -1, drop = FALSE]  # 删除 ID 列，只保留丰度
  
  cat("OTU总数:", nrow(feature_data_only), "\n")
  
  # 打印 OTU ID 列前 6 行
  cat("\n原始OTU ID前6行:\n")
  print(head(feature_data$目, 6))
  
  # ----------------------------
  # 构建 phyloseq 对象
  # ----------------------------
  common_samples <- intersect(colnames(feature_data_only), metadata$SampleID)
  feature_data_common <- feature_data_only[, common_samples, drop = FALSE]
  metadata_common <- metadata[common_samples, , drop = FALSE]
  rownames(metadata_common) <- metadata_common$SampleID
  
  cat("共同样本数量:", length(common_samples), "\n")
  
  # 检查是否有足够的样本
  if (length(common_samples) < 3) {
    cat("警告: 样本数量不足，跳过", microbe, "\n")
    next
  }
  
  # 构建phyloseq对象时保留原始OTU ID作为行名
  otu <- otu_table(as.matrix(feature_data_common), taxa_are_rows = TRUE)
  sam <- sample_data(metadata_common)
  ps_obj <- phyloseq(otu, sam)
  
  # 保存原始OTU ID和丰度数据
  original_taxon_names <- rownames(feature_data_common)
  original_abundance_data <- feature_data_common
  
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
      struc_zero = FALSE,
      neg_lb = TRUE,
      alpha = 0.05,
      n_cl = 6
    )
  }, error = function(e) {
    cat(paste0("ANCOM-BC2分析出错: ", microbe, " - ", e$message, "\n"))
    return(NULL)
  })
  
  # ----------------------------
  # 处理结果 - 仅使用名称匹配方法
  # ----------------------------
  if (!is.null(ancombc_res)) {
    res <- ancombc_res$res
    
    # 添加ANCOM-BC2结果行数统计
    cat("==========================================\n")
    cat("ANCOM-BC2 结果统计 -", microbe, "\n")
    cat("==========================================\n")
    cat("原始输入OTU数量:     ", length(original_taxon_names), "\n")
    cat("ANCOM-BC2结果OTU数量:", nrow(res), "\n")
    cat("过滤掉的OTU数量:     ", length(original_taxon_names) - nrow(res), "\n")
    cat("保留比例:            ", round(nrow(res)/length(original_taxon_names)*100, 2), "%\n")
    
    # 检查行数是否匹配
    if (nrow(res) != length(original_taxon_names)) {
      cat("警告: ANCOM-BC2内部过滤了", length(original_taxon_names) - nrow(res), "个OTU\n")
    }
    
    # 使用名称匹配方法：检查ANCOM-BC2结果中的物种名称是否在原始名称中
    matched_taxa <- intersect(original_taxon_names, res$taxon)
    cat("通过名称匹配的物种数量:", length(matched_taxa), "\n")
    
    if (length(matched_taxa) == nrow(res)) {
      # 情况1: 所有结果物种都能在原始名称中找到
      cat("情况1: 所有ANCOM-BC2结果物种都能匹配到原始名称\n")
      final_res <- res
      # 获取对应的丰度数据
      abundance_data_to_add <- original_abundance_data[final_res$taxon, , drop = FALSE]
      
    } else if (length(matched_taxa) > 0) {
      # 情况2: 只有部分匹配，保留匹配的物种
      cat("情况2: 只有部分物种匹配，保留", length(matched_taxa), "个匹配的物种\n")
      final_res <- res[res$taxon %in% matched_taxa, ]
      # 获取对应的丰度数据
      abundance_data_to_add <- original_abundance_data[final_res$taxon, , drop = FALSE]
      
    } else {
      # 情况3: 如果完全无法匹配，使用原始名称（按顺序）
      cat("情况3: 无法通过名称匹配，使用顺序映射\n")
      final_res <- res
      final_res$taxon <- original_taxon_names[1:nrow(res)]
      # 获取对应的丰度数据
      abundance_data_to_add <- original_abundance_data[final_res$taxon, , drop = FALSE]
    }
    
    # 将丰度数据插入到taxon列之后
    cat("正在合并丰度数据...\n")
    
    # 找到taxon列的位置
    taxon_col_index <- which(colnames(final_res) == "taxon")
    
    # 将结果表分成两部分：taxon列之前和之后
    before_taxon <- final_res[, 1:taxon_col_index, drop = FALSE]
    after_taxon <- final_res[, (taxon_col_index + 1):ncol(final_res), drop = FALSE]
    
    # 合并：taxon列之前 + 丰度数据 + taxon列之后
    final_res_with_abundance <- cbind(
      before_taxon,
      abundance_data_to_add,
      after_taxon
    )
    
    # 保持原始样本名称，不添加前缀
    # 丰度数据列名保持不变
    
    cat("丰度数据已成功插入，新增", ncol(abundance_data_to_add), "个丰度列\n")
    
    # 打印最终结果的前6行（只显示前几列以免输出过长）
    cat("\n最终结果前6行（显示前8列）:\n")
    print(head(final_res_with_abundance[, 1:min(8, ncol(final_res_with_abundance))], 6))
    
    # 保存为 Excel
    output_file <- paste0(output_dir, microbe, "_ANCOMBC2_results.xlsx")
    wb <- createWorkbook()
    addWorksheet(wb, "Complete_Results")
    writeData(wb, "Complete_Results", final_res_with_abundance)
    saveWorkbook(wb, output_file, overwrite = TRUE)
    
    cat("==========================================\n")
    cat(paste0("ANCOM-BC2结果已保存: ", output_file, "\n"))
    cat(paste0("结果表维度: ", nrow(final_res_with_abundance), " 行, ", ncol(final_res_with_abundance), " 列\n"))
    cat("其中包含", ncol(abundance_data_to_add), "个样本的丰度数据\n")
    cat("==========================================\n\n")
    
  } else {
    cat("==========================================\n")
    cat(paste0("ANCOM-BC2分析结果为空: ", microbe, "\n"))
    cat("==========================================\n\n")
  }
}

# 最终总结
cat("\n\n==========================================\n")
cat("所有分析完成! 总结:\n")
cat("==========================================\n")
for (microbe in names(otu_files)) {
  output_file <- paste0(output_dir, microbe, "_ANCOMBC2_results.xlsx")
  if (file.exists(output_file)) {
    wb <- loadWorkbook(output_file)
    sheet_data <- readWorkbook(wb, sheet = 1)
    # 计算丰度数据列数（统计结果列通常有固定名称，其他列为丰度数据）
    stat_columns <- c("taxon", "lfc_Group", "se_Group", "W_Group", "p_Group", "q_Group", "diff_Group")
    abundance_cols <- ncol(sheet_data) - length(intersect(colnames(sheet_data), stat_columns))
    cat(microbe, ": ", nrow(sheet_data), "行结果, ", abundance_cols, "个样本的丰度数据\n")
  } else {
    cat(microbe, ": 文件不存在\n")
  }
}
cat("==========================================\n")