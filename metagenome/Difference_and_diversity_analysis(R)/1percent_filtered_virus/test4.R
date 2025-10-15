library(openxlsx)

file_path <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"

# 只读取 virus_Result 表
df <- read.xlsx(file_path, sheet = "virus_Result")

# 筛选符合条件的行
row_indices <- which(
  df$`passed_ss_(Intercept)` == TRUE &
  df$`passed_ss_GroupMI` == TRUE &
  abs(df$`lfc_GroupMI`) > 1 &
  df$`q_GroupMI` < 0.1
)

filtered_rows <- df[row_indices, ]

cat("virus_Result表:", length(row_indices), "个显著差异物种\n")

# 打印筛选出来的行的"解释信息"列
if(length(row_indices) > 0 && "解释信息" %in% colnames(filtered_rows)) {
  cat("\n筛选结果的解释信息:\n")
  for(i in 1:nrow(filtered_rows)) {
    cat(i, ". ", filtered_rows$解释信息[i], "\n")
  }
} else if(length(row_indices) > 0) {
  cat("警告：未找到'解释信息'列\n")
  cat("可用的列名：\n")
  print(colnames(filtered_rows))
} else {
  cat("没有符合条件的显著差异物种\n")
}