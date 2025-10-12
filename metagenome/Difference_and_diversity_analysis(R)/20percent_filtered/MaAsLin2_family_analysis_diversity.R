# =============================================================================
# Family 水平 α & β 多样性分析
# - α: Observed, Shannon, Simpson (小提琴图 + 箱线图，副标题加上Mean±SD)
# - β: Bray-Curtis PCoA + PERMANOVA (副标题加上R²和p值，带置信椭圆)
# - 输入: family 丰度表 + metadata
# - 输出: 每张图单独保存为 PDF
# =============================================================================

library(vegan)
library(ggplot2)
library(dplyr)
library(readxl)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(rlang)  # 用于 tidy evaluation

# ----------------------------
# 1. 文件路径
# ----------------------------
family_files <- list(
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/family_data/virus_family_no_HF.xlsx",
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/family_data/archaea_family.xlsx",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/family_data/bacteria_family.xlsx",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/family_data/fungi_family.xlsx"
)

metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/MaAsLin2/sample_metadata.xlsx"
outdir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/family_data/diversity_results_M"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ----------------------------
# 2. 读取 metadata
# ----------------------------
metadata <- read_excel(metadata_file)
colnames(metadata)[1] <- "SampleID"
metadata$Group <- factor(metadata$Group)

# ----------------------------
# 3. 定义 α 多样性计算函数
# ----------------------------
calc_alpha <- function(otu_table) {
  richness <- specnumber(t(otu_table))        # Observed
  shannon <- diversity(t(otu_table), index = "shannon")
  simpson <- diversity(t(otu_table), index = "simpson")
  data.frame(
    SampleID = rownames(t(otu_table)),
    Observed = richness,
    Shannon = shannon,
    Simpson = simpson
  )
}

# ----------------------------
# 4. 分析每类微生物
# ----------------------------
for (microbe in names(family_files)) {
  cat("\n=========== 开始分析:", microbe, "===========\n")
  
  # 读取 family 表
  df <- readxl::read_excel(family_files[[microbe]])
  df <- as.data.frame(df)
  
  # 第一列为 family 名
  rownames(df) <- df[,1]
  df <- df[,-1]
  
  # 样本对齐 metadata
  common_samples <- intersect(colnames(df), metadata$SampleID)
  df <- df[, common_samples, drop = FALSE]
  metadata_sub <- metadata %>% filter(SampleID %in% common_samples)
  
  # ---------------- α 多样性 ----------------
  alpha_df <- calc_alpha(df) %>%
    left_join(metadata_sub, by = "SampleID")
  
  for (m in c("Observed","Shannon","Simpson")) {
    stats <- alpha_df %>%
      group_by(Group) %>%
      summarise(label = paste0("Mean±SD: ",
                               round(mean(.data[[m]]),2)," ± ",
                               round(sd(.data[[m]]),2))) %>%
      pull(label) %>%
      paste(collapse = "; ")
    
    p <- ggplot(alpha_df, aes(x = Group, y = .data[[m]], fill = Group)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      geom_boxplot(width = 0.2, outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 1, alpha = 0.7) +
      theme_bw() +
      labs(title = paste0(microbe, " - ", m, " Index"),
           subtitle = stats)
    
    ggsave(file.path(outdir, paste0(microbe, "_alpha_", m, ".pdf")),
           p, width = 6, height = 5)
  }
  
  # ---------------- β 多样性 ----------------
  dist <- vegdist(t(df), method = "bray")
  pcoa <- cmdscale(dist, eig = TRUE, k = 2)
  pcoa_df <- data.frame(
    SampleID = rownames(pcoa$points),
    Axis1 = pcoa$points[,1],
    Axis2 = pcoa$points[,2]
  ) %>% left_join(metadata_sub, by = "SampleID")
  
  adonis_res <- adonis2(dist ~ Group, data = metadata_sub)
  r2 <- round(adonis_res$R2[1], 3)
  pval <- signif(adonis_res$`Pr(>F)`[1], 3)
  
  p_beta <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Group)) +
    stat_ellipse(type = "t", linetype = 2) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    labs(title = paste0(microbe, " - Bray-Curtis PCoA"),
         subtitle = paste0("PERMANOVA: R²=", r2, ", p=", pval))
  
  ggsave(file.path(outdir, paste0(microbe, "_beta_BrayCurtis_PCoA.pdf")),
         p_beta, width = 6, height = 5)
  
  cat("✅ 已保存 PDF 到: ", outdir, "\n")
}

cat("\n🎉 Family 水平 α & β 多样性分析完成！结果保存在：", outdir, "\n")
