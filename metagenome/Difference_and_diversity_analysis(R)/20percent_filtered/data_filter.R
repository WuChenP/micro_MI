# =============================================================================
# 微生物数据流行率过滤脚本
# 适用于R 4.5.1版本
# =============================================================================

# 清理环境
rm(list = ls())

# 设置工作目录（可选）
# setwd("E:/Python/MI_Analysis")

# 安装和加载必要的包
cat("正在安装和加载必要的包...\n")
if (!require("readxl", quietly = TRUE)) {
  install.packages("readxl")
  library(readxl)
} else {
  library(readxl)
}

if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
  library(dplyr)
} else {
  library(dplyr)
}

# 定义文件路径
file_paths <- c(
  virus = "E:/Python/MI_Analysis/origin_data/心梗组_病毒(新).xlsx",
  bacteria = "E:/Python/MI_Analysis/origin_data/心梗组_细菌.xlsx", 
  fungi = "E:/Python/MI_Analysis/origin_data/心梗组_真菌.xlsx",
  archaea = "E:/Python/MI_Analysis/origin_data/心梗组_古菌.xlsx"
)

# 检查文件是否存在
cat("检查文件是否存在...\n")
for (path in file_paths) {
  if (file.exists(path)) {
    cat("找到文件:", path, "\n")
  } else {
    cat("警告: 文件不存在:", path, "\n")
  }
}

# 读取Excel数据的函数
read_microbiome_data <- function(file_path) {
  cat("正在读取文件:", basename(file_path), "\n")
  
  tryCatch({
    data <- read_excel(file_path)
    
    cat("原始数据维度 - 行:", nrow(data), "列:", ncol(data), "\n")
    
    if (ncol(data) > 1) {
      row_names <- as.character(data[[1]])
      data_df <- as.data.frame(data[, -1])
      rownames(data_df) <- row_names
      colnames(data_df) <- colnames(data)[-1]
      
      cat("处理后的维度 - 物种数:", nrow(data_df), "样本数:", ncol(data_df), "\n")
      return(data_df)
    } else {
      stop("数据列数不足，无法处理")
    }
  }, error = function(e) {
    cat("读取文件时出错:", e$message, "\n")
    return(NULL)
  })
}

# 流行率过滤函数
prevalence_filter <- function(data, prevalence_threshold = 0.2) {
  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    cat("数据为空或维度异常，跳过过滤\n")
    return(NULL)
  }
  
  cat("正在进行流行率过滤...\n")
  
  data_t <- as.data.frame(t(data))
  prevalence <- apply(data_t, 2, function(x) sum(x > 0, na.rm = TRUE) / length(x[!is.na(x)]))
  
  filtered_indices <- prevalence > prevalence_threshold
  filtered_species <- names(prevalence)[filtered_indices]
  
  cat("原始物种数:", length(prevalence), "\n")
  cat("过滤后物种数:", sum(filtered_indices), "\n")
  cat("流行率阈值:", prevalence_threshold * 100, "%\n")
  
  if (sum(filtered_indices) == 0) {
    cat("警告: 没有物种通过过滤！\n")
    return(NULL)
  }
  
  filtered_data <- data[filtered_species, , drop = FALSE]
  return(filtered_data)
}

# 保存CSV数据
save_filtered_data <- function(data, data_name, output_dir = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/") {
  if (is.null(data)) {
    cat(data_name, "数据为空，跳过保存\n")
    return(FALSE)
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_file <- file.path(output_dir, paste0(data_name, "_filtered_20percent.csv"))
  
  data_to_save <- data
  data_to_save$Taxonomy <- rownames(data)
  data_to_save <- data_to_save[, c("Taxonomy", setdiff(colnames(data_to_save), "Taxonomy"))]
  
  write.csv(data_to_save, output_file, row.names = FALSE)
  cat("已保存CSV文件:", output_file, "\n")
  return(TRUE)
}

# 生成总结报告
generate_summary_report <- function(original_list, filtered_list) {
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("          流行率过滤总结报告\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  names <- c("病毒", "细菌", "真菌", "古菌")
  
  for (i in 1:4) {
    cat("\n--- ", names[i], " ---\n", sep = "")
    
    orig_data <- original_list[[i]]
    filtered_data <- filtered_list[[i]]
    
    if (!is.null(orig_data) && !is.null(filtered_data)) {
      orig_count <- nrow(orig_data)
      filtered_count <- nrow(filtered_data)
      retention_rate <- if (orig_count > 0) round(filtered_count / orig_count * 100, 1) else NA
      
      cat("原始物种数: ", orig_count, "\n", sep = "")
      cat("过滤后物种数: ", filtered_count, "\n", sep = "")
      cat("保留比例: ", retention_rate, "%\n", sep = "")
    } else {
      if (is.null(orig_data)) cat("原始数据读取失败\n")
      if (is.null(filtered_data)) cat("过滤后数据为空\n")
    }
  }
}

# 主函数
main <- function() {
  cat("开始微生物数据流行率过滤流程\n")
  cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  virus_data <- read_microbiome_data(file_paths["virus"])
  bacteria_data <- read_microbiome_data(file_paths["bacteria"])
  fungi_data <- read_microbiome_data(file_paths["fungi"])
  archaea_data <- read_microbiome_data(file_paths["archaea"])
  
  original_list <- list(virus_data, bacteria_data, fungi_data, archaea_data)
  
  cat("\n步骤2: 应用流行率过滤(阈值20%)...\n")
  virus_filtered <- prevalence_filter(virus_data, 0.2)
  bacteria_filtered <- prevalence_filter(bacteria_data, 0.2)
  fungi_filtered <- prevalence_filter(fungi_data, 0.2)
  archaea_filtered <- prevalence_filter(archaea_data, 0.2)
  
  filtered_list <- list(virus_filtered, bacteria_filtered, fungi_filtered, archaea_filtered)
  
  cat("\n步骤3: 保存过滤后的数据...\n")
  save_filtered_data(virus_filtered, "病毒")
  save_filtered_data(bacteria_filtered, "细菌")
  save_filtered_data(fungi_filtered, "真菌")
  save_filtered_data(archaea_filtered, "古菌")
  
  cat("\n步骤4: 生成总结报告...\n")
  generate_summary_report(original_list, filtered_list)
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n流程完成！\n")
}

# 执行主函数
main()
cat("脚本执行完毕！\n")
