# =============================================================================
# 将 ANCOM-BC2 结果追加到原丰度表后（按 taxon 匹配）
# - 保留原丰度表的 ID（改名为 taxon）
# - 仅保留匹配行，taxon 对齐后再拼接
# - 拼接后仅保留一个 taxon 列
# - 直接保存回原 Excel
# =============================================================================

library(openxlsx)
library(dplyr)

# ----------------------------
# 文件路径
# ----------------------------
result_file <- "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"

otu_files <- list(
  fungi    = "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_fungi_species_level.csv",
  archaea = "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_archaea_species_level.csv",
  bacteria = "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_bacteria_species_level.csv",
  virus    = "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_virus_species_level.csv"
)

# ----------------------------
# 读取 ANCOM-BC2 结果文件的所有 sheet
# ----------------------------
sheet_names <- getSheetNames(result_file)
wb_out <- loadWorkbook(result_file)  # 直接在原文件上操作

for (microbe in names(otu_files)) {
  sheet_name <- paste0(microbe, "_Result")
  if (!(sheet_name %in% sheet_names)) next
  
  cat("\n处理:", microbe, "\n")
  
  # 读取 ANCOM-BC2 结果
  ancom_res <- read.xlsx(result_file, sheet = sheet_name)
  if (!"taxon" %in% colnames(ancom_res)) {
    colnames(ancom_res)[1] <- "taxon"
  }
  
  # 读取原始丰度表
  abund_df <- read.csv(otu_files[[microbe]], check.names = FALSE)
  if (!"ID" %in% colnames(abund_df)) {
    stop(paste0("❌ 文件缺少 ID 列: ", otu_files[[microbe]]))
  }
  
  # 将原始 ID 改为 taxon
  abund_df <- abund_df %>% rename(taxon = ID)
  
  # 按 taxon 匹配（只保留匹配行）
  merged_df <- left_join(abund_df, ancom_res, by = "taxon")
  
  # 确保 taxon 列在第一列
  merged_df <- merged_df %>% select(taxon, everything())
  
  # 写回原 Excel（覆盖该 sheet）
  removeWorksheet(wb_out, sheet_name)
  addWorksheet(wb_out, sheet_name)
  writeData(wb_out, sheet_name, merged_df)
  
  cat("✅ 已匹配并拼接:", microbe, "（按 taxon 对齐）\n")
}

# ----------------------------
# 保存文件（直接覆盖原文件）
# ----------------------------
saveWorkbook(wb_out, result_file, overwrite = TRUE)
cat("\n🎉 已更新原结果文件:\n", result_file, "\n")
