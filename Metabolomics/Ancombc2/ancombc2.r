library(phyloseq)
library(ANCOMBC)
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# 文件路径
feature_file <- "E:/代谢组学/Ancombc2/abundance.tsv"
metadata_file <- "E:/代谢组学/Ancombc2/metadata.tsv"
output_file <- "E:/代谢组学/AMI_vs_CON_lev1_ANCOM-BC2.xlsx"

# 读取丰度表和元数据
feature_data <- read.delim(feature_file, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_file, check.names = FALSE)

# 检查样本匹配
common_samples <- intersect(colnames(feature_data), metadata$SampleID)
if(length(common_samples) == 0){
  stop("没有匹配的样本！请检查丰度表列名和元数据 SampleID")
}
feature_data_filtered <- feature_data[, common_samples, drop = FALSE]
metadata_filtered <- metadata[match(common_samples, metadata$SampleID), , drop = FALSE]
rownames(metadata_filtered) <- metadata_filtered$SampleID

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
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 6
  )
}, error = function(e){
  stop(paste("ANCOM-BC2分析出错:", e$message))
})

res <- ancombc_res$res


# 添加 Significance 列（只看 p/q + passed_ss，不看 diff_robust）##其中p为p值  q为fdr校正后的p值
res$Significance <- ifelse(
  res$q_groupCON < 0.05,
  ifelse(res$lfc_groupCON > 1, "up",
         ifelse(res$lfc_groupCON < -1, "down", "insignificant")),
  "insignificant"
)

# 添加 Prevalence 和 Mean_Abundance
res$Prevalence <- apply(feature_data_filtered, 1, function(x) sum(x > 0) / length(x))
res$Mean_Abundance <- apply(feature_data_filtered, 1, mean)

# 拼接原始丰度表
final_res <- cbind(feature_data_filtered, res)

# 保存到 Excel
wb <- createWorkbook()
addWorksheet(wb, "Complete_Results")
writeData(wb, "Complete_Results", final_res, rowNames = TRUE)
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("ANCOM-BC2结果已保存至:", output_file, "\n")