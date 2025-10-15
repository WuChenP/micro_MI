# =============================================================================
# 批量将 MaAsLin2 输出的 TSV 文件另存为 Excel（xlsx）
# 包含四类微生物文件夹 + 合并显著结果文件
# =============================================================================

library(openxlsx)
library(dplyr)

# ----------------------------
# 文件夹列表
# ----------------------------
folders <- c(
  "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/MaAsLin2_family/archaea",
  "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/MaAsLin2_family/bacteria",
  "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/MaAsLin2_family/fungi",
  "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/MaAsLin2_family/virus"
)

# ----------------------------
# 额外单独文件
# ----------------------------
extra_files <- c(
  "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/MaAsLin2_family/combined_significant_results.tsv"
)

# ----------------------------
# 函数：TSV 转 XLSX
# ----------------------------
convert_tsv_to_xlsx <- function(tsv_file){
  cat("读取 TSV 文件:", tsv_file, "\n")
  df <- read.delim(tsv_file, check.names = FALSE, stringsAsFactors = FALSE)
  xlsx_file <- sub("\\.tsv$", ".xlsx", tsv_file)
  write.xlsx(df, xlsx_file, rowNames = FALSE)
  cat("✅ 已保存为 Excel:", xlsx_file, "\n")
}

# ----------------------------
# 处理文件夹中的 TSV 文件
# ----------------------------
for(folder in folders){
  cat("\n处理文件夹:", folder, "\n")
  
  tsv_files <- list.files(folder, pattern = "all_results\\.tsv$|significant_results\\.tsv$", full.names = TRUE)
  
  if(length(tsv_files) == 0){
    cat("⚠️ 没有找到 TSV 文件\n")
    next
  }
  
  for(tsv in tsv_files){
    convert_tsv_to_xlsx(tsv)
  }
}

# ----------------------------
# 处理额外单独文件
# ----------------------------
for(tsv in extra_files){
  if(file.exists(tsv)){
    convert_tsv_to_xlsx(tsv)
  } else {
    cat("⚠️ 文件不存在:", tsv, "\n")
  }
}

cat("\n🎉 所有 TSV 文件已成功转换为 Excel 格式\n")
