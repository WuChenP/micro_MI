# ==============================================================
# 四类微生物 ANCOM-BC2 火山图（仿照 W_value 调色）
# ==============================================================

library(openxlsx)
library(dplyr)
library(ggplot2)
library(scales)  # 用于 alpha 调色

# ----------------------------
# 文件路径
# ----------------------------
input_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"
output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/volcano_plots/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

microbe_types <- c("archaea", "bacteria", "fungi", "virus")
sheet_names <- getSheetNames(input_file)

# ----------------------------
# 循环绘制
# ----------------------------
for (microbe in microbe_types) {
  
  sheet_name <- paste0(microbe, "_Result")
  if (!(sheet_name %in% sheet_names)) next
  
  df <- read.xlsx(input_file, sheet = sheet_name)
  
  # 确保所需列存在
  if (!all(c("lfc_GroupMI", "q_GroupMI", "W_GroupMI") %in% colnames(df))) next
  
  plot_df <- df %>%
    mutate(
      log2FC = lfc_GroupMI,
      negLog10Q = -log10(q_GroupMI),
      Status = case_when(
        q_GroupMI < 0.05 & log2FC > 1  ~ "MI_Enriched",
        q_GroupMI < 0.05 & log2FC < -1 ~ "CON_Enriched",
        TRUE ~ "NotSig"
      ),
      W_value = W_GroupMI
    ) %>%
    filter(!is.na(log2FC) & !is.na(negLog10Q) & is.finite(negLog10Q) & !is.na(W_value))
  
  # 统计数量
  counts <- plot_df %>% group_by(Status) %>% summarise(n = n(), .groups = "drop")
  status_labels <- setNames(
    paste0(counts$Status, " (n=", counts$n, ")"),
    counts$Status
  )
  
  # 使用 W_value 调节颜色深浅
  # 使用 W_value 调节颜色深浅
  plot_df$color <- case_when(
    plot_df$Status == "MI_Enriched"  ~ scales::alpha(scales::col_numeric(
      palette = c("#FF9999", "#990000"),
      domain = range(plot_df$W_value, na.rm=TRUE)
    )(plot_df$W_value), 
    0.8),
    plot_df$Status == "CON_Enriched" ~ alpha("blue", rescale(plot_df$W_value, to=c(0.4,1))),
    TRUE                             ~ "grey70"
  )
  
  point_size <- ifelse(nrow(plot_df) > 500, 0.8, 1.5)
  
  # 绘图
  p <- ggplot(plot_df, aes(x = log2FC, y = negLog10Q)) +
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
  pdf_file <- paste0(output_dir, microbe, "_volcano_q0.05_log2FC1_Wcolor.pdf")
  ggsave(filename = pdf_file, plot = p, width = 7, height = 6)
  message(paste("✅ 已生成火山图 PDF:", microbe))
}
