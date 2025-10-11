# =============================================================================
# MaAsLin2 火山图（仅 q<0.05 & |log2FC|>1 显著标识物）
# 图例格式与 ANCOMBC2/前面火山图一致
# =============================================================================

library(tidyverse)
library(readr)
library(ggplot2)

# ----------------------------
# 文件路径
# ----------------------------
all_results_files <- list(
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/MaAsLin2_family/archaea/all_results.tsv",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/MaAsLin2_family/bacteria/all_results.tsv",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/MaAsLin2_family/fungi/all_results.tsv",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/MaAsLin2_family/virus/all_results.tsv"
)

# 输出目录
outdir_q <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent/MaAsLin2_family/volcano_plots"
dir.create(outdir_q, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 火山图绘制函数（仅 q<0.05 & |log2FC|>1 显著点）
# ----------------------------
plot_volcano_qval <- function(df, x_col = "coef", y_col = "qval", fc_cut = 1, sig_cut = 0.05, microbe, outdir) {
  
  # 筛选显著状态
  df <- df %>%
    mutate(Status = case_when(
      .data[[y_col]] < sig_cut & .data[[x_col]] > fc_cut  ~ "MI_Enriched",
      .data[[y_col]] < sig_cut & .data[[x_col]] < -fc_cut ~ "CON_Enriched",
      TRUE ~ "NotSig"
    ),
    negLog10 = -log10(.data[[y_col]]))
  
  # 统计数量
  counts <- df %>% group_by(Status) %>% summarise(n = n())
  status_labels <- setNames(
    paste0(counts$Status, " (n=", counts$n, ")"),
    counts$Status
  )
  
  # 基础颜色
  colors <- c("MI_Enriched" = "red", "CON_Enriched" = "blue", "NotSig" = "grey70")
  
  # 火山图
  p <- ggplot(df, aes(x = .data[[x_col]], y = negLog10)) +
    geom_point(aes(color = Status), alpha = 0.8, size = 2) +
    scale_color_manual(values = colors, labels = status_labels) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(sig_cut), linetype = "dashed", color = "black") +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(microbe, " Volcano Plot (q<0.05 & |log2FC|>1)"),
      subtitle = paste0("Thresholds: |log2FC|>", fc_cut, ", q<", sig_cut),
      x = "log2FC",
      y = paste0("-log10(", y_col, ")")
    ) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    )
  
  # 保存 PDF
  ggsave(file.path(outdir, paste0(microbe, "_volcano_qval_significant.pdf")),
         p, width = 6, height = 5)
}

# ----------------------------
# 循环绘图
# ----------------------------
for (microbe in names(all_results_files)) {
  cat("\n=========== 开始绘制火山图:", microbe, "===========\n")
  
  df <- read_tsv(all_results_files[[microbe]])
  
  # 绘制 q 值火山图（仅显著点）
  plot_volcano_qval(df, x_col = "coef", y_col = "qval", fc_cut = 1, sig_cut = 0.05,
                    microbe = microbe, outdir = outdir_q)
  
  cat("✅ ", microbe, "火山图已保存 (仅 q<0.05 & |log2FC|>1)\n")
}

cat("\n🎉 所有微生物类火山图生成完成！\n")
cat("Q值火山图目录：", outdir_q, "\n")
