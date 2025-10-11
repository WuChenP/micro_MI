library(openxlsx)

# 输入文件
input_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/ancombc2_results_OTU/four_microbes_with_abundance.xlsx"

# 输出目录
outdir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/ancombc2_results_OTU/split_by_microbe/"
if (!dir.exists(outdir)) dir.create(outdir)

# 读取所有 sheet 名称
sheets <- getSheetNames(input_file)

# 提取微生物类型前缀（archaea / bacteria / fungi / virus）
microbe_types <- unique(sub("_(Complete_Results|Significant_Species|Up_Regulated|Down_Regulated)$", "", sheets))

# 遍历每种微生物
for (microbe in microbe_types) {
  wb <- createWorkbook()
  
  # 找到该微生物的所有 sheet
  microbe_sheets <- sheets[grepl(paste0("^", microbe, "_"), sheets)]
  
  for (sh in microbe_sheets) {
    df <- read.xlsx(input_file, sheet = sh)
    addWorksheet(wb, sh)
    writeData(wb, sh, df)
  }
  
  # 保存一个文件
  outfile <- file.path(outdir, paste0(microbe, "_all_results.xlsx"))
  saveWorkbook(wb, outfile, overwrite = TRUE)
  message("已保存: ", outfile)
}
