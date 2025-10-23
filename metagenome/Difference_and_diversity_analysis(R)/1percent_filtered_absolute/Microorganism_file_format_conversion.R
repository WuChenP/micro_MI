library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)

# 定义一个函数，把分类路径表汇总到 Family 水平
aggregate_to_family <- function(file_path){
  df <- read_csv(file_path, show_col_types = FALSE)  # 关闭列类型提示
  
  tax_col <- colnames(df)[1]
  sample_cols <- colnames(df)[-1]
  
  df$Family <- sapply(strsplit(df[[tax_col]], "\\|"), function(x){
    fam <- grep("^f__", x, value = TRUE)
    if(length(fam) == 0) return(NA)
    gsub("^f__", "", fam)
  })
  
  df_family <- df %>%
    filter(!is.na(Family)) %>%
    group_by(Family) %>%
    summarise(across(all_of(sample_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop")  # 新写法
  
  return(df_family)
}

# 处理三个微生物类文件
archaea_family <- aggregate_to_family("E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/心梗组_古菌_filtered_1percent.csv")
bacteria_family <- aggregate_to_family("E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/心梗组_细菌_filtered_1percent.csv")
fungi_family   <- aggregate_to_family("E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/心梗组_真菌_filtered_1percent.csv")
virus_family <- aggregate_to_family("E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/心梗组_病毒_filtered_1percent.csv")

# 分别保存为 XLSX 文件
write.xlsx(archaea_family, "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/family_level_data/archaea_family.xlsx", rowNames = FALSE)
write.xlsx(bacteria_family, "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/family_level_data/bacteria_family.xlsx", rowNames = FALSE)
write.xlsx(fungi_family,   "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/family_level_data/fungi_family.xlsx",   rowNames = FALSE)
write.xlsx(virus_family, "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/family_level_data/virus_family.xlsx", rowNames = FALSE)
