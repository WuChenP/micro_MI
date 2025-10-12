library(openxlsx)
library(dplyr)

# ----------------------------
# 文件路径
# ----------------------------
input_file  <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1%/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"
output_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1%/ancombc2_results_OTU/四类微生物_ANCOMBC2_with_significant.xlsx"

# 获取 sheet 名称
sheet_names <- getSheetNames(input_file)

# ----------------------------
# 定义函数：标记显著性并方向
# ----------------------------
mark_significant <- function(df) {
  df %>%
    mutate(Significant = case_when(
      q_GroupMI < 0.05 & diff_GroupMI & lfc_GroupMI > 1  ~ "Up_in_MI",
      q_GroupMI < 0.05 & diff_GroupMI & lfc_GroupMI < -1 ~ "Up_in_Control",
      TRUE ~ "Not_significant"
    ))
}

# ----------------------------
# 读取、标记并保存
# ----------------------------
wb <- createWorkbook()

for (sheet in sheet_names) {
  message("处理表：", sheet)
  
  df <- read.xlsx(input_file, sheet = sheet)
  df_marked <- mark_significant(df)
  
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df_marked)
  
  message("  → 上调数量：", sum(df_marked$Significant=="Up_in_MI", na.rm=TRUE),
          "，下调数量：", sum(df_marked$Significant=="Up_in_Control", na.rm=TRUE))
}

saveWorkbook(wb, output_file, overwrite = TRUE)
message("✅ 已保存到文件：", output_file)
