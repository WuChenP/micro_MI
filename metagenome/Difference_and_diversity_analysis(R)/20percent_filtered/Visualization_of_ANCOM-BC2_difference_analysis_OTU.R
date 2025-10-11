library(openxlsx)
library(ggplot2)
library(dplyr)
library(scales)  # 用于 rescale

file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/ancombc2_results_OTU/four_microbes_with_abundance.xlsx"
sheet_names <- getSheetNames(file)
microbe_types <- c("archaea", "bacteria", "fungi", "virus")

output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/ancombc2_results_OTU/volcano_plots/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (microbe in microbe_types) {
  
  sheet_name <- paste0(microbe, "_Complete_Results")
  if (!(sheet_name %in% sheet_names)) next
  
  df <- read.xlsx(file, sheet = sheet_name, rowNames = TRUE)
  
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
  counts <- plot_df %>% group_by(Status) %>% summarise(n = n())
  status_labels <- setNames(
    paste0(counts$Status, " (n=", counts$n, ")"),
    counts$Status
  )
  
  # 基本颜色
  base_colors <- c("MI_Enriched" = "red", "CON_Enriched" = "blue", "NotSig" = "grey70")
  
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
  
  ggsave(
    filename = paste0(output_dir, microbe, "_volcano_q0.05_log2FC1_Wcolor.pdf"),
    plot = p, width = 7, height = 6
  )
  
  message(paste("✅ 已生成火山图 PDF:", microbe))
}
