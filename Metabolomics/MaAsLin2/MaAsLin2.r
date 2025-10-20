# 安装需要的包（如未安装）
# install.packages("BiocManager")
# BiocManager::install("Maaslin2")
# install.packages("readr")
# install.packages("writexl")  # 写Excel用

library(Maaslin2)
library(readr)
library(writexl)

# 1. 读入数据文件路径
feature <- read.delim("abundance.tsv", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim("metadata.tsv", header = TRUE, row.names = 1, check.names = FALSE ,stringsAsFactors=FALSE)
output_dir <- "E:/代谢组学/AMI_vs_CON_lev1_MaAsLin2"   # 输出目录

# 2. 运行 MaAsLin2
fit_data <- Maaslin2(
  input_data = feature,
  input_metadata = metadata,
  output = output_dir,
  fixed_effects = c("Group"),  # 你的分组变量名
  normalization = "TSS",       # 总和标准化
  transform = "LOG"            # log 转换
)
cat("==========样本名一致检测 ==========\n")
all(colnames(feature) %in% rownames(metadata))  # 应该返回 TRUE
# 3. 读 MaAsLin2 结果
all_results <- read_tsv(file.path(output_dir, "all_results.tsv"))

# 4. 加 type 列：根据 coef 和 qval 判断
# 4.1 新增 FC 列（基于 coef）
all_results$FC <- exp(all_results$coef)

# # 4.2 如果需要，也可以加到显著差异表里
# sig_results$FC <- exp(sig_results$coef)

all_results$type <- ifelse(
  all_results$qval < 0.05 & all_results$FC > 2, "up",
  ifelse(all_results$qval < 0.05 & all_results$FC < 0.5, "down", "insignificant")
)

# 5. 提取显著差异物（type 是 up 或 down）
sig_results <- subset(all_results, type %in% c("up", "down"))

# 6. 保存为 Excel
write_xlsx(list(
  All_Results = all_results,
  Significant_Results = sig_results
), path = file.path(output_dir, "E:/代谢组学/AMI_vs_CON_lev1_Maaslin2.xlsx"))

# 7. 打印显著差异物前几行查看
print(head(sig_results))

