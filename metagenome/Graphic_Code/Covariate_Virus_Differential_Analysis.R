# =============================================================================
# 四类微生物 OTU 表 ANCOM-BC2 分析（包含协变量校正）
# =============================================================================

library(phyloseq)
library(ANCOMBC)
library(openxlsx)
library(dplyr)

# ----------------------------
# 1. 文件路径
# ----------------------------
otu_files <- list(
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/心梗组_病毒_filtered_1percent.csv"
)

metadata_file = "E:/Python/MI_Analysis/origin_data/样本协变量数据.xlsx"

# 输出文件夹
output_dir <- "E:/Python/MI_Analysis/metagenome/Graphic/Covariate_Virus_Differential_Analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. 读取并预处理元数据
# ----------------------------
metadata <- read.xlsx(metadata_file)

# 重命名列以便于使用
colnames(metadata) <- c("分析名称", "Group", "Gender", "Age", "BMI", "Smoking", "Drinking", 
                        "Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "ACE")

# 设置行名为样本ID
rownames(metadata) <- metadata$分析名称

# 将分类变量转换为因子
metadata$Group <- factor(metadata$Group, levels = c("CON", "MI"))  # CON作为参考组
metadata$Gender <- factor(metadata$Gender, levels = c(1, 2), labels = c("Male", "Female"))  # 1=男, 2=女
metadata$Smoking <- factor(metadata$Smoking, levels = c(0, 1), labels = c("No", "Yes"))
metadata$Drinking <- factor(metadata$Drinking, levels = c(0, 1), labels = c("No", "Yes"))

# 查看元数据结构
cat("元数据变量:\n")
print(colnames(metadata))
cat("\n分组分布:\n")
print(table(metadata$Group))
cat("\n性别分布:\n") 
print(table(metadata$Gender))
cat("\n吸烟分布:\n")
print(table(metadata$Smoking))
cat("\n饮酒分布:\n")
print(table(metadata$Drinking))

# ----------------------------
# 循环分析四类微生物
# ----------------------------
for (microbe in names(otu_files)) {
  
  cat("\n正在分析:", microbe, "\n")
  cat("==========================================\n")
  
  # 读取 OTU 表
  feature_data <- read.csv(otu_files[[microbe]], check.names = FALSE)
  rownames(feature_data) <- feature_data$ID
  feature_data_only <- feature_data[, -1, drop = FALSE]  # 删除 ID 列，只保留丰度
  
  cat("OTU总数:", nrow(feature_data_only), "\n")
  
  # ----------------------------
  # 构建 phyloseq 对象
  # ----------------------------
  common_samples <- intersect(colnames(feature_data_only), metadata$分析名称)
  feature_data_common <- feature_data_only[, common_samples, drop = FALSE]
  metadata_common <- metadata[common_samples, , drop = FALSE]
  rownames(metadata_common) <- metadata_common$分析名称
  
  cat("共同样本数量:", length(common_samples), "\n")
  
  # 检查是否有足够的样本
  if (length(common_samples) < 3) {
    cat("警告: 样本数量不足，跳过", microbe, "\n")
    next
  }
  
  # 构建phyloseq对象
  otu <- otu_table(as.matrix(feature_data_common), taxa_are_rows = TRUE)
  sam <- sample_data(metadata_common)
  ps_obj <- phyloseq(otu, sam)
  
  # 保存原始OTU ID和丰度数据
  original_taxon_names <- rownames(feature_data_common)
  original_abundance_data <- feature_data_common
  
  # ----------------------------
  # 运行 ANCOM-BC2 - 包含所有协变量
  # ----------------------------
  cat("运行ANCOM-BC2，校正协变量: Group + Age + Gender + BMI + Smoking + Drinking\n")
  
  ancombc_res <- tryCatch({
    ancombc2(
      data = ps_obj,
      # 关键修改：包含所有协变量
      fix_formula = "Group + Age + Gender + BMI + Smoking + Drinking",
      p_adj_method = "fdr",
      lib_cut = 0,
      group = "Group",  # 主要分组变量
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
  # 处理结果
  # ----------------------------
  if (!is.null(ancombc_res)) {
    res <- ancombc_res$res
    
    cat("ANCOM-BC2分析成功完成!\n")
    cat("原始输入OTU数量:", length(original_taxon_names), "\n")
    cat("ANCOM-BC2结果OTU数量:", nrow(res), "\n")
    
    # 使用名称匹配方法
    matched_taxa <- intersect(original_taxon_names, res$taxon)
    cat("通过名称匹配的物种数量:", length(matched_taxa), "\n")
    
    if (length(matched_taxa) == nrow(res)) {
      final_res <- res
      abundance_data_to_add <- original_abundance_data[final_res$taxon, , drop = FALSE]
    } else if (length(matched_taxa) > 0) {
      final_res <- res[res$taxon %in% matched_taxa, ]
      abundance_data_to_add <- original_abundance_data[final_res$taxon, , drop = FALSE]
    } else {
      final_res <- res
      final_res$taxon <- original_taxon_names[1:nrow(res)]
      abundance_data_to_add <- original_abundance_data[final_res$taxon, , drop = FALSE]
    }
    
    # 将丰度数据插入到taxon列之后
    taxon_col_index <- which(colnames(final_res) == "taxon")
    before_taxon <- final_res[, 1:taxon_col_index, drop = FALSE]
    after_taxon <- final_res[, (taxon_col_index + 1):ncol(final_res), drop = FALSE]
    
    final_res_with_abundance <- cbind(
      before_taxon,
      abundance_data_to_add,
      after_taxon
    )
    
    # 打印显著差异结果
    sig_results <- final_res_with_abundance[final_res_with_abundance$diff_GroupMI == TRUE, ]
    cat("发现", nrow(sig_results), "个显著差异物种 (q < 0.05)\n")
    
    if (nrow(sig_results) > 0) {
      cat("显著差异物种前5个:\n")
      print(sig_results[1:min(5, nrow(sig_results)), 
                        c("taxon", "lfc_GroupMI", "q_GroupMI", "diff_GroupMI")])
    }
    
    # 保存为 Excel
    output_file <- paste0(output_dir, microbe, "_ANCOMBC2_results_with_covariates.xlsx")
    wb <- createWorkbook()
    addWorksheet(wb, "Complete_Results")
    writeData(wb, "Complete_Results", final_res_with_abundance)
    saveWorkbook(wb, output_file, overwrite = TRUE)
    
    cat("结果已保存:", output_file, "\n")
    cat("==========================================\n\n")
    
  } else {
    cat("ANCOM-BC2分析失败\n")
    cat("==========================================\n\n")
  }
}

# 最终总结
cat("\n\n==========================================\n")
cat("所有分析完成! 总结:\n")
cat("==========================================\n")
for (microbe in names(otu_files)) {
  output_file <- paste0(output_dir, microbe, "_ANCOMBC2_results_with_covariates.xlsx")
  if (file.exists(output_file)) {
    wb <- loadWorkbook(output_file)
    sheet_data <- readWorkbook(wb, sheet = 1)
    sig_count <- sum(sheet_data$diff_GroupMI == TRUE, na.rm = TRUE)
    cat(microbe, ": ", nrow(sheet_data), "个物种, ", sig_count, "个显著差异物种\n")
  } else {
    cat(microbe, ": 文件不存在\n")
  }
}
cat("==========================================\n")