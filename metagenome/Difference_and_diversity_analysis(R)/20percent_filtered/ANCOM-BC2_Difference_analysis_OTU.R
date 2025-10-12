library(phyloseq)
library(ANCOMBC)
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# =========================================================
# 1. 读取元数据表
# =========================================================
metadata_path <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/sample_metadata.xlsx"
metadata <- read.xlsx(metadata_path)

# 标准化列名
colnames(metadata) <- c("Sample_ID", "Group")

# 检查是否正确读取
cat("已读取样本元数据，共", nrow(metadata), "个样本。\n")

# =========================================================
# 2. 微生物丰度文件路径列表
# =========================================================
microbe_files <- list(
  archaea = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/古菌_filtered_20percent.csv",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/细菌_filtered_20percent.csv",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/真菌_filtered_20percent.csv",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/病毒_filtered_20percent.csv"
)

# 输出文件路径
output_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/ancombc2_results_OTU/four_microbes_with_abundance.xlsx"
wb <- createWorkbook()

# =========================================================
# 3. 主循环：逐一分析每类微生物
# =========================================================
for (microbe in names(microbe_files)) {
  
  cat("\n========== 开始分析:", microbe, "==========\n")
  
  # 读取丰度表
  feature_data <- read.csv(microbe_files[[microbe]], row.names = 1, check.names = FALSE)
  
  # ✅ 给第一列加 taxon 列名
  feature_data <- cbind(taxon = rownames(feature_data), feature_data)
  rownames(feature_data) <- feature_data$taxon
  
  # 筛选元数据中存在的样本
  common_samples <- intersect(colnames(feature_data)[-1], metadata$Sample_ID)
  if (length(common_samples) == 0) {
    cat("⚠️ 无匹配样本，跳过:", microbe, "\n")
    next
  }
  
  feature_data_filtered <- feature_data[, c("taxon", common_samples), drop = FALSE]
  metadata_filtered <- metadata[match(common_samples, metadata$Sample_ID), , drop = FALSE]
  rownames(metadata_filtered) <- metadata_filtered$Sample_ID
  
  # 构建 phyloseq 对象
  otu <- otu_table(as.matrix(feature_data_filtered[, -1]), taxa_are_rows = TRUE)
  sam <- sample_data(metadata_filtered)
  ps_obj <- phyloseq(otu, sam)
  
  # =========================================================
  # 运行 ANCOM-BC2
  # =========================================================
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
    cat(paste0("❌ ANCOM-BC2分析出错: ", microbe, " - ", e$message, "\n"))
    return(NULL)
  })
  
  # =========================================================
  # 结果处理
  # =========================================================
  if (!is.null(ancombc_res)) {
    res <- ancombc_res$res
    
    # 删除重复的 taxon 列（如果存在）
    if ("taxon" %in% colnames(res)) {
      res <- res[, !colnames(res) %in% "taxon", drop = FALSE]
    }
    
    # 添加显著性标记
    res$Significance <- ifelse(
      res$q_GroupMI < 0.05 & abs(res$lfc_GroupMI) > 1,
      ifelse(res$lfc_GroupMI > 1, "Up",
             ifelse(res$lfc_GroupMI < -1, "Down", "Not significant")),
      "Not significant"
    )
    
    # 添加 Prevalence 和 Mean_Abundance
    numeric_mat <- as.matrix(feature_data_filtered[, -1])
    res$Prevalence <- apply(numeric_mat, 1, function(x) sum(x > 0) / length(x))
    res$Mean_Abundance <- apply(numeric_mat, 1, mean)
    
    # 拼接原始丰度表（保留唯一的 taxon 列）
    final_res <- cbind(feature_data_filtered, res)
    
    # 写入 Excel
    addWorksheet(wb, paste0(microbe, "_Complete_Results"))
    writeData(wb, paste0(microbe, "_Complete_Results"), final_res, rowNames = FALSE)
    
    # 筛选显著结果
    sig <- final_res[final_res$Significance != "Not significant", ]
    if (nrow(sig) > 0) {
      addWorksheet(wb, paste0(microbe, "_Significant_Species"))
      writeData(wb, paste0(microbe, "_Significant_Species"), sig, rowNames = FALSE)
      
      up <- sig[sig$Significance == "Up", ]
      down <- sig[sig$Significance == "Down", ]
      if (nrow(up) > 0) {
        addWorksheet(wb, paste0(microbe, "_Up_Regulated"))
        writeData(wb, paste0(microbe, "_Up_Regulated"), up, rowNames = FALSE)
      }
      if (nrow(down) > 0) {
        addWorksheet(wb, paste0(microbe, "_Down_Regulated"))
        writeData(wb, paste0(microbe, "_Down_Regulated"), down, rowNames = FALSE)
      }
    }
  } else {
    cat(paste0("⚠️ ", microbe, " 无分析结果。\n"))
  }
}

# 确保输出目录存在（避免 “No such file or directory” 错误）
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 4. 保存 Excel
# =========================================================
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("\n✅ 四类微生物完整结果已保存至:", output_file, "\n")
