library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)

# 定义一个函数，把分类路径表汇总到 Family 水平
aggregate_to_family <- function(file_path){
  df <- read_csv(file_path)
  
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
    summarise(across(all_of(sample_cols), sum, na.rm = TRUE))
  
  return(df_family)
}

# 处理三个微生物类文件
archaea_family <- aggregate_to_family("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/古菌_filtered_20percent.csv")
bacteria_family <- aggregate_to_family("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/细菌_filtered_20percent.csv")
fungi_family   <- aggregate_to_family("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/真菌_filtered_20percent.csv")

# 分别保存为 XLSX 文件
write.xlsx(archaea_family, "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/family_data/archaea_family.xlsx", rowNames = FALSE)
write.xlsx(bacteria_family, "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/family_data/bacteria_family.xlsx", rowNames = FALSE)
write.xlsx(fungi_family,   "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/family_data/fungi_family.xlsx",   rowNames = FALSE)
