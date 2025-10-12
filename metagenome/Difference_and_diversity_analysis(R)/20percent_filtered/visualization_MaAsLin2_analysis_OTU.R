# ======================================================
# MaAsLin2 火山图（q<0.05 & |log2FC|>1，显著点分组，图例显示数量）
# ======================================================

library(ggplot2)
library(dplyr)
library(readr)

# ----------------------------
# 1. 路径设置
# ----------------------------
outdir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/MaAsLin2_OTU"
pdf_dir <- file.path(outdir, "volcano_grouped")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# 结果文件路径
all_results_files <- list(
  virus    = file.path(outdir, "virus/all_results.tsv"),
  archaea  = file.path(outdir, "archaea/all_results.tsv"),
  bacteria = file.path(outdir, "bacteria/all_results.tsv"),
  fungi    = file.path(outdir, "fungi/all_results.tsv")
)

# 阈值设置
fc_threshold <- 1   # |log2FC| > 1
qval_threshold <- 0.05

# ----------------------------
# 2. 火山图绘制函数
# ----------------------------
plot_volcano_grouped <- function(df, microbe, outdir, fc_cut=1, q_cut=0.05) {
  
  df <- df %>%
    mutate(Significance = case_when(
      qval <= q_cut & coef > fc_cut  ~ "MI_Enriched",
      qval <= q_cut & coef < -fc_cut ~ "CON_Enriched",
      TRUE ~ "NotSig"
    ),
    neglog10q = -log10(qval))
  
  # 统计各类别数量
  counts <- df %>% count(Significance)
  count_labels <- paste0(counts$Significance, " (n=", counts$n, ")")
  names(count_labels) <- counts$Significance
  
  # 阈值线
  thresh_y <- -log10(q_cut)
  
  # 绘图
  p <- ggplot(df, aes(x = coef, y = neglog10q, color = Significance)) +
    geom_point(size=2.5, alpha=0.8) +
    scale_color_manual(values = c("MI_Enriched"="red", "CON_Enriched"="blue", "NotSig"="grey70"),
                       labels = count_labels) +
    geom_hline(yintercept = thresh_y, linetype="dashed", color="black") +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype="dashed", color="black") +
    theme_classic(base_size=14) +
    labs(title = paste0(microbe, " Volcano Plot"),
         subtitle = "-log10(q) vs log2FC",
         x = "log2FC",
         y = "-log10(q)") +
    theme(
      plot.title = element_text(face="bold", hjust=0.5, size=16),
      axis.title = element_text(face="bold", size=14),
      axis.text = element_text(size=12, color="black"),
      legend.title = element_blank(),
      legend.text = element_text(size=12),
      legend.position = "right"
    )
  
  # 保存 PDF
  ggsave(file.path(outdir, paste0(microbe, "_volcano_grouped.pdf")),
         p, width=6, height=5)
}

# ----------------------------
# 3. 循环生成火山图
# ----------------------------
for (microbe in names(all_results_files)) {
  cat("\n=========== 开始绘制火山图:", microbe, "===========\n")
  
  file <- all_results_files[[microbe]]
  if(!file.exists(file)) {
    cat("⚠️ 文件不存在:", file, "\n")
    next
  }
  
  df <- readr::read_tsv(file)
  
  plot_volcano_grouped(df, microbe, pdf_dir, fc_cut = fc_threshold, q_cut = qval_threshold)
  
  cat("✅ ", microbe, "火山图已生成并保存至 PDF\n")
}

cat("\n🎉 所有微生物类火山图生成完成！目录：", pdf_dir, "\n")
