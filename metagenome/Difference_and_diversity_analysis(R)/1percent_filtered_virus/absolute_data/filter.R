# ======================================================
# 微生物相对丰度表 - 1% 流行率过滤 + CSV 输出
# 清理 ID 列中的不可见字符 - 修复编码问题版本
# 只对细菌进行ID过滤，其他保持原样
# ======================================================

library(readxl)
library(openxlsx)

# 输入文件路径
files <- c(
  virus = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/心梗组_病毒.xlsx"
)

# 输出文件夹
out_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 流行率阈值（1%）
prev_cutoff <- 0.01

# ----------------------------
# 函数：清理不可见字符 + 流行率过滤 + CSV 输出
# ----------------------------
filter_prevalence_csv <- function(file, cutoff = 0.01, out_dir, microbe_type = NULL) {
  
  # 使用 readxl 包读取 Excel 文件，避免编码问题
  cat("正在读取文件:", basename(file), "\n")
  df <- read_excel(file)
  df <- as.data.frame(df)  # 转换为数据框
  
  # 清理 ID 列：去掉不可见字符并确保唯一
  cat("清理 ID 列中的不可见字符...\n")
  if ("ID" %in% colnames(df)) {
    # 移除不可见字符
    df$ID <- gsub("[^[:alnum:][:punct:][:space:]]", "", df$ID)
    # 确保唯一性
    df$ID <- make.unique(as.character(df$ID))
  } else {
    # 如果第一列不是ID，假设第一列是物种ID
    colnames(df)[1] <- "ID"
    df$ID <- gsub("[^[:alnum:][:punct:][:space:]]", "", df$ID)
    df$ID <- make.unique(as.character(df$ID))
  }
  
  # 提取物种 ID 和丰度矩阵
  taxa <- df$ID
  abund <- df[, -1, drop = FALSE]
  
  # 强制转换为数值
  cat("转换数据为数值类型...\n")
  abund[] <- lapply(abund, function(x) {
    # 先转换为字符，再转换为数值
    as.numeric(as.character(x))
  })
  
  # 检查是否有非数值数据
  na_count <- sum(is.na(abund))
  if (na_count > 0) {
    cat("警告: 发现", na_count, "个NA值，将被视为0\n")
    abund[is.na(abund)] <- 0
  }
  
  # 计算每个物种的流行率
  cat("计算流行率...\n")
  prevalence <- apply(abund, 1, function(x) mean(x > 0, na.rm = TRUE))
  
  # 过滤低流行率物种
  keep <- prevalence >= cutoff
  df_filtered <- df[keep, ]
  
  # 构造输出 CSV 文件名
  fname <- basename(file)
  out_file <- file.path(out_dir, gsub(".xlsx", "_filtered_1percent.csv", fname))
  
  # 保存 CSV
  write.csv(df_filtered, out_file, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  
  cat("✅ 已完成 1% 流行率过滤并生成 CSV:", fname, "=> 保留", nrow(df_filtered), "个物种\n")
  cat("原始物种数:", nrow(df), "=> 过滤后物种数:", nrow(df_filtered), "\n")
  cat("输出文件:", out_file, "\n\n")
}

# ----------------------------
# 执行处理
# ----------------------------
cat("开始处理微生物丰度数据...\n")
cat("输出目录:", out_dir, "\n")
cat("流行率阈值:", prev_cutoff, "\n\n")

for (microbe_type in names(files)) {
  cat("=== 处理", microbe_type, "数据 ===\n")
  tryCatch({
    filter_prevalence_csv(files[[microbe_type]], cutoff = prev_cutoff, out_dir = out_dir, microbe_type = microbe_type)
  }, error = function(e) {
    cat("❌ 处理", microbe_type, "时出错:", e$message, "\n")
  })
}

cat("所有处理完成！\n")