library(phyloseq)
library(ANCOMBC)
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# 四类微生物文件路径列表
microbe_files <- list(
  archaea = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/古菌_filtered_20percent.csv",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/细菌_filtered_20percent.csv",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/真菌_filtered_20percent.csv",
  virus   = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/病毒_filtered_20percent.csv"
)

# 输出文件路径
output_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/ancombc2_results/four_microbes_with_abundance.xlsx"
wb <- createWorkbook()

for (microbe in names(microbe_files)) {
  
  # 读取数据
  feature_data <- read.csv(microbe_files[[microbe]], row.names = 1, check.names = FALSE)
  
  # 确保元数据与样本匹配
  common_samples <- intersect(colnames(feature_data), metadata$Sample_ID)
  feature_data_filtered <- feature_data[, common_samples, drop = FALSE]
  metadata_filtered <- metadata[match(common_samples, metadata$Sample_ID), , drop = FALSE]
  metadata_filtered <- as.data.frame(metadata_filtered)
  rownames(metadata_filtered) <- metadata_filtered$Sample_ID
  
  # 构建 phyloseq 对象
  otu <- otu_table(as.matrix(feature_data_filtered), taxa_are_rows = TRUE)
  sam <- sample_data(metadata_filtered)
  ps_obj <- phyloseq(otu, sam)
  
  # 运行 ANCOM-BC2
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
  
  if (!is.null(ancombc_res)) {
    res <- ancombc_res$res
    
    # 添加 Significance 列
    res$Significance <- ifelse(
      (res$p_GroupMI < 0.05 | res$q_GroupMI < 0.05) & res$passed_ss_GroupMI == TRUE,
      ifelse(res$lfc_GroupMI > 0, "Up", "Down"),
      "Not significant"
    )
    
    # 添加 Prevalence 和 Mean_Abundance
    res$Prevalence <- apply(feature_data_filtered, 1, function(x) sum(x > 0) / length(x))
    res$Mean_Abundance <- apply(feature_data_filtered, 1, mean)
    
    # 拼接原始丰度表
    final_res <- cbind(feature_data_filtered, res)
    
    # 写入 Excel，每类微生物一个 Sheet
    addWorksheet(wb, paste0(microbe, "_Complete_Results"))
    writeData(wb, paste0(microbe, "_Complete_Results"), final_res, rowNames = TRUE)
    
    # 显著物种
    sig <- final_res[final_res$Significance != "Not significant", ]
    if (nrow(sig) > 0) {
      addWorksheet(wb, paste0(microbe, "_Significant_Species"))
      writeData(wb, paste0(microbe, "_Significant_Species"), sig, rowNames = TRUE)
      
      up <- sig[sig$Significance == "Up", ]
      down <- sig[sig$Significance == "Down", ]
      if (nrow(up) > 0) {
        addWorksheet(wb, paste0(microbe, "_Up_Regulated"))
        writeData(wb, paste0(microbe, "_Up_Regulated"), up, rowNames = TRUE)
      }
      if (nrow(down) > 0) {
        addWorksheet(wb, paste0(microbe, "_Down_Regulated"))
        writeData(wb, paste0(microbe, "_Down_Regulated"), down, rowNames = TRUE)
      }
    }
  } else {
    cat(paste0("ANCOM-BC2分析结果为空: ", microbe, "\n"))
  }
}

# 保存 Excel 文件
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("四类微生物完整结果已保存至:", output_file, "\n")
