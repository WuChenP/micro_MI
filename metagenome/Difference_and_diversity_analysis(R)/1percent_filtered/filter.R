# ======================================================
# 微生物相对丰度表 - 1% 流行率过滤 + CSV 输出
# 清理 ID 列中的不可见字符
# 删除病毒数据中以 HF 开头的行
# ======================================================

library(openxlsx)

# 输入文件路径
files <- c(
  archaea = "E:/Python/MI_Analysis/origin_data/心梗组_古菌.xlsx",
  bacteria = "E:/Python/MI_Analysis/origin_data/心梗组_细菌.xlsx",
  fungi = "E:/Python/MI_Analysis/origin_data/心梗组_真菌.xlsx",
  virus = "E:/Python/MI_Analysis/origin_data/心梗组_病毒(新).xlsx"
)

# 输出文件夹
out_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 流行率阈值（1%）
prev_cutoff <- 0.01

# ----------------------------
# 函数：清理不可见字符 + 流行率过滤 + CSV 输出
# ----------------------------
filter_prevalence_csv <- function(file, cutoff = 0.01, out_dir, microbe_type = NULL) {
  # 读取 Excel
  df <- read.xlsx(file)
  
  # 清理 ID 列：去掉不可见字符并确保唯一
  df$ID <- make.unique(iconv(as.character(df$ID), from = "UTF-8", to = "UTF-8", sub = ""))
  
  # 提取物种 ID 和丰度矩阵
  taxa <- df$ID
  abund <- df[, -1]
  
  # 强制转换为数值
  abund[] <- lapply(abund, function(x) as.numeric(as.character(x)))
  
  # 计算每个物种的流行率
  prevalence <- apply(abund, 1, function(x) mean(x > 0))
  
  # 过滤低流行率物种
  keep <- prevalence >= cutoff
  df_filtered <- df[keep, ]
  
  # 如果是病毒，删除 ID 以 "HF" 开头的行
  if (!is.null(microbe_type) && microbe_type == "virus") {
    df_filtered <- df_filtered[!grepl("^HF", df_filtered$ID), ]
  }
  
  # 构造输出 CSV 文件名
  fname <- basename(file)
  out_file <- file.path(out_dir, gsub(".xlsx", "_filtered_1percent.csv", fname))
  
  # 保存 CSV
  write.csv(df_filtered, out_file, row.names = FALSE, quote = FALSE)
  
  cat("✅ 已完成 1% 流行率过滤并生成 CSV:", fname, "=> 保留", nrow(df_filtered), "个物种\n")
}

# ----------------------------
# 执行处理
# ----------------------------
for (microbe_type in names(files)) {
  filter_prevalence_csv(files[[microbe_type]], cutoff = prev_cutoff, out_dir = out_dir, microbe_type = microbe_type)
}
