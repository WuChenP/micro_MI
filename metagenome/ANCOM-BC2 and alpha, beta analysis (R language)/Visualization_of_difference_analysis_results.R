library(openxlsx)
library(ggplot2)
library(dplyr)

file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/ancombc2_results/four_microbes_with_abundance.xlsx"
sheet_names <- getSheetNames(file)
microbe_types <- c("archaea", "bacteria", "fungi", "virus")

# ========================================
# 火山图（仅根据 log2FC，不考虑 Significance）
# ========================================
for (microbe in microbe_types) {
  
  sheet_name <- paste0(microbe, "_Complete_Results")
  if (!(sheet_name %in% sheet_names)) next
  
  df <- read.xlsx(file, sheet = sheet_name, rowNames = TRUE)
  if (!all(c("lfc_GroupMI", "p_GroupMI") %in% colnames(df))) next
  
  plot_df <- df %>%
    mutate(
      log2FC = lfc_GroupMI,
      negLog10P = -log10(p_GroupMI),
      Status = case_when(
        log2FC < -1 ~ "CON_Enriched",
        log2FC > 1  ~ "MI_Enriched",
        TRUE        ~ "NotSig"
      )
    ) %>%
    filter(!is.na(log2FC) & !is.na(negLog10P) & is.finite(negLog10P))
  
  # 统计每类数量
  counts <- plot_df %>% group_by(Status) %>% summarise(n = n())
  
  point_size <- ifelse(nrow(plot_df) > 500, 0.8, 1.5)
  
  p <- ggplot(plot_df, aes(x = log2FC, y = negLog10P, color = Status)) +
    geom_point(alpha = 0.7, size = point_size) +
    scale_color_manual(values = c("CON_Enriched" = "blue", "MI_Enriched" = "red", "NotSig" = "grey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black")) +
    labs(
      title = paste("Volcano Plot -", microbe, "(log2FC > 1 / < -1)"),
      x = "log2 Fold Change (MI vs CON)",
      y = "-log10(p-value)"
    ) +
    # 添加数量标签
    annotate("text", x = max(plot_df$log2FC), y = max(plot_df$negLog10P),
             label = paste0(counts$Status, ": ", counts$n, collapse = "\n"),
             hjust = 1, vjust = 1, size = 4)
  
  ggsave(
    filename = paste0("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/ancombc2_results/", microbe, "_volcano_log2FC_only.png"),
    plot = p, width = 7, height = 6, dpi = 300
  )
  
  message(paste("已生成火山图:", microbe))
}

# ========================================
# 火山图（考虑 Significance）
# ========================================
for (microbe in microbe_types) {
  
  sheet_name <- paste0(microbe, "_Complete_Results")
  if (!(sheet_name %in% sheet_names)) next
  
  df <- read.xlsx(file, sheet = sheet_name, rowNames = TRUE)
  if (!all(c("lfc_GroupMI", "p_GroupMI", "Significance") %in% colnames(df))) next
  
  plot_df <- df %>%
    mutate(
      log2FC = lfc_GroupMI,
      negLog10P = -log10(p_GroupMI),
      Status = case_when(
        Significance == "Up" & log2FC > 1  ~ "MI_Enriched",
        Significance == "Down" & log2FC < -1 ~ "CON_Enriched",
        TRUE ~ "NotSig"
      )
    ) %>%
    filter(!is.na(log2FC) & !is.na(negLog10P) & is.finite(negLog10P))
  
  counts <- plot_df %>% group_by(Status) %>% summarise(n = n())
  
  point_size <- ifelse(nrow(plot_df) > 500, 0.8, 1.5)
  
  p <- ggplot(plot_df, aes(x = log2FC, y = negLog10P, color = Status)) +
    geom_point(alpha = 0.7, size = point_size) +
    scale_color_manual(values = c("CON_Enriched" = "blue", "MI_Enriched" = "red", "NotSig" = "grey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black")) +
    labs(
      title = paste("Volcano Plot -", microbe, "(Significance & |log2FC|>1)"),
      x = "log2 Fold Change (MI vs CON)",
      y = "-log10(p-value)"
    ) +
    annotate("text", x = max(plot_df$log2FC), y = max(plot_df$negLog10P),
             label = paste0(counts$Status, ": ", counts$n, collapse = "\n"),
             hjust = 1, vjust = 1, size = 4)
  
  ggsave(
    filename = paste0("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/ancombc2_results/", microbe, "_volcano_log2FC_with_Significance.png"),
    plot = p, width = 7, height = 6, dpi = 300
  )
  
  message(paste("已生成火山图:", microbe))
}
