library(openxlsx)
library(ggplot2)
library(dplyr)
library(scales)  # 可选，用于颜色渐变

# ----------------------------
# 文件路径列表
# ----------------------------
result_files <- list(
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/ancombc2_results_family/archaea_ANCOMBC2_results.xlsx",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/ancombc2_results_family/bacteria_ANCOMBC2_results.xlsx",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/ancombc2_results_family/fungi_ANCOMBC2_results.xlsx",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/ancombc2_results_family/virus_ANCOMBC2_results.xlsx"
)

# 输出 PDF 文件夹
output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/ancombc2_results_family/volcano_plots/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 循环绘图
# ----------------------------
for (microbe in names(result_files)) {
  
  cat("正在绘制火山图 (q-value) :", microbe, "\n")
  
  # 读取数据
  df <- read.xlsx(result_files[[microbe]], sheet = 1)
  
  # 检查关键列
  if(!all(c("lfc_GroupMI", "q_GroupMI", "Family", "W_GroupMI") %in% colnames(df))){
    cat(paste0("文件 ", microbe, " 缺少必要列，跳过\n"))
    next
  }
  
  # 数据处理
  volcano_data <- df %>%
    filter(!is.na(lfc_GroupMI) & !is.na(q_GroupMI) & !is.na(W_GroupMI)) %>%
    mutate(
      negLog10Q = -log10(q_GroupMI),
      Status = case_when(
        q_GroupMI < 0.05 & lfc_GroupMI > 1  ~ "MI_Enriched",
        q_GroupMI < 0.05 & lfc_GroupMI < -1 ~ "CON_Enriched",
        TRUE ~ "NotSig"
      ),
      # 使用 W 值调节颜色深浅
      color = case_when(
        Status == "MI_Enriched"  ~ scales::alpha(scales::col_numeric(
          palette = c("#FF9999", "#CC0000"),
          domain = range(W_GroupMI, na.rm = TRUE)
        )(W_GroupMI), 0.8),
        Status == "CON_Enriched" ~ scales::alpha("blue", rescale(W_GroupMI, to = c(0.4,1))),
        TRUE                     ~ "grey70"
      )
    )
  
  # 统计每类数量
  counts <- volcano_data %>% group_by(Status) %>% summarise(n = n())
  status_labels <- setNames(
    paste0(counts$Status, " (n=", counts$n, ")"),
    counts$Status
  )
  
  # 点大小
  point_size <- ifelse(nrow(volcano_data) > 500, 0.8, 1.5)
  
  # 火山图
  p <- ggplot(volcano_data, aes(x = lfc_GroupMI, y = negLog10Q)) +
    geom_point(aes(color = color, shape = Status), size = point_size) +
    scale_shape_manual(values = c("CON_Enriched"=16, "MI_Enriched"=16, "NotSig"=16),
                       labels = status_labels) +
    scale_color_identity() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.title = element_blank()
    ) +
    labs(
      title = paste("Volcano Plot -", microbe, "(q<0.05 & |log2FC|>1)"),
      x = "log2 Fold Change (MI vs CON)",
      y = "-log10(q-value)"
    )
  
  # 保存 PDF
  output_file <- paste0(output_dir, microbe, "_volcano_qvalue_Wcolor.pdf")
  ggsave(output_file, plot = p, width = 7, height = 6)
  
  cat(paste0("✅ 火山图已保存为 PDF: ", output_file, "\n\n"))
}
