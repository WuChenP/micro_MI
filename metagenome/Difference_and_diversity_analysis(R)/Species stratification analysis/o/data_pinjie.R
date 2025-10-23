# 加载必要的包
library(readxl)
library(dplyr)
library(openxlsx)

# 文件路径
file1_path <- "E:/Python/MI_Analysis/metagenome/data_figures/o/目水平.xlsx"
file2_path <- "E:/Python/MI_Analysis/metagenome/data_figures/o/ancombc2_results_new/virus_ANCOMBC2_results.xlsx"
output_path <- "E:/Python/MI_Analysis/metagenome/data_figures/o/o_ANCOM-BC2_merged_result.xlsx"

# 读取两个Excel文件
df1 <- read_excel(file1_path)
df2 <- read_excel(file2_path)

# 查找匹配的列名
find_join_column <- function(df) {
  # 查找包含vOTU的列名
  vOTU_columns <- grep("vOTU", colnames(df), value = TRUE, ignore.case = TRUE)
  if(length(vOTU_columns) > 0) {
    return(vOTU_columns[1])
  }
  
  # 如果没有找到vOTU列，返回第一列
  return(colnames(df)[1])
}

# 获取匹配列名
join_column_df1 <- find_join_column(df1)
cat("使用列进行匹配:", join_column_df1, "-> taxon\n")

# 选择需要的列
df1_selected <- df1 %>%
  select(all_of(join_column_df1))

# 使用反引号处理包含特殊字符的列名
df2_selected <- df2 %>%
  select(taxon, MI35, MI33, MI32, MI38, MI36, MI34, MI37, 
         MI27, MI28, MI29, MI31, MI17, MI21, MI20, MI24, 
         MI23, MI26, MI19, MI18, MI15, MI25, MI22, MI14, 
         MI16, MI1, MI2, MI3, MI4, MI5, MI6, MI7, MI8, 
         MI9, MI10, MI11, MI12, MI13, CON1, CON2, CON3, 
         CON4, CON5, CON6, CON7, CON8, CON9, CON10, CON11, 
         CON12, CON13, CON14, CON15, CON16, CON17, CON18, 
         CON19, CON20, CON21, CON22, CON23, CON24, CON25, 
         CON26, CON27, CON28, CON29, CON30, CON31, CON32, 
         CON33, CON34, CON35, CON36, CON37, CON38, CON39, 
         CON40, CON41, CON42, CON43, CON44, CON45, CON46, 
         CON47, `lfc_(Intercept)`, `lfc_GroupMI`, `se_(Intercept)`, 
         `se_GroupMI`, `W_(Intercept)`, `W_GroupMI`, `p_(Intercept)`, 
         `p_GroupMI`, `q_(Intercept)`, `q_GroupMI`, `diff_(Intercept)`, 
         `diff_GroupMI`, `passed_ss_(Intercept)`, `passed_ss_GroupMI`, 
         `diff_robust_(Intercept)`, `diff_robust_GroupMI`)

# 执行左连接
merged_df <- df2_selected %>%
  left_join(df1_selected, by = c("taxon" = join_column_df1))

# 保存结果
write.xlsx(merged_df, output_path)

cat("文件拼接完成！\n")
cat("结果保存路径:", output_path, "\n")
cat("原始df2行数:", nrow(df2), "\n")
cat("合并后行数:", nrow(merged_df), "\n")
cat("列数:", ncol(merged_df), "\n")











# 加载必要的包
library(readxl)
library(dplyr)
library(openxlsx)

# 文件路径
file1_path <- "E:/Python/MI_Analysis/metagenome/data_figures/o/目水平.xlsx"
file2_path <- "E:/Python/MI_Analysis/metagenome/data_figures/o/MaAsLin2_OTU_new/virus/all_results.xlsx"
output_path <- "E:/Python/MI_Analysis/metagenome/data_figures/o/o_MaAsLin2_merged_result.xlsx"

# 读取两个Excel文件
df1 <- read_excel(file1_path)
df2 <- read_excel(file2_path)

# 查找匹配的列名
find_join_column <- function(df) {
  # 查找包含vOTU的列名
  vOTU_columns <- grep("vOTU", colnames(df), value = TRUE, ignore.case = TRUE)
  if(length(vOTU_columns) > 0) {
    return(vOTU_columns[1])
  }
  
  # 如果没有找到vOTU列，返回第一列
  return(colnames(df)[1])
}

# 获取匹配列名
join_column_df1 <- find_join_column(df1)
cat("使用列进行匹配:", join_column_df1, "-> taxon\n")

# 选择需要的列
df1_selected <- df1 %>%
  select(all_of(join_column_df1), MI35, MI33, MI32, 
         MI38, MI36, MI34, MI37, MI27, MI28, MI29, MI31, 
         MI17, MI21, MI20, MI24, MI23, MI26, MI19, MI18, 
         MI15, MI25, MI22, MI14, MI16, MI1, MI2, MI3, 
         MI4, MI5, MI6, MI7, MI8, MI9, MI10, MI11, MI12, 
         MI13, CON1, CON2, CON3, CON4, CON5, CON6, CON7, 
         CON8, CON9, CON10, CON11, CON12, CON13, CON14, 
         CON15, CON16, CON17, CON18, CON19, CON20, CON21, 
         CON22, CON23, CON24, CON25, CON26, CON27, CON28, 
         CON29, CON30, CON31, CON32, CON33, CON34, CON35, 
         CON36, CON37, CON38, CON39, CON40, CON41, CON42, 
         CON43, CON44, CON45, CON46, CON47)

# 使用反引号处理包含特殊字符的列名
df2_selected <- df2 %>%
  select(feature, metadata, value, coef, stderr, N, N.not.0, pval, qval)

# 执行左连接
merged_df <- df2_selected %>%
  left_join(df1_selected, by = c("feature" = join_column_df1))

# 保存结果
write.xlsx(merged_df, output_path)

cat("文件拼接完成！\n")
cat("结果保存路径:", output_path, "\n")
cat("原始df2行数:", nrow(df2), "\n")
cat("合并后行数:", nrow(merged_df), "\n")
cat("列数:", ncol(merged_df), "\n")