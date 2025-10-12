# =====================================================
# 四类微生物 Family 水平 α/β 多样性分析 (小提琴+箱线，β拼图)
# =====================================================

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(patchwork) # 拼图

# ----------------------------
# 文件路径
# ----------------------------
microbe_files <- list(
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/fungi_family.xlsx",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/virus_family_no_HF.xlsx",
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/archaea_family.xlsx",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/bacteria_family.xlsx"
)

meta_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/diversity_plots_family/"

# ----------------------------
# 读取元数据
# ----------------------------
metadata <- read.xlsx(meta_file, rowNames = TRUE)
metadata$ID <- rownames(metadata)

# ----------------------------
# 函数：Excel 转 phyloseq
# ----------------------------
make_phyloseq_from_excel <- function(file, metadata){
  abundance <- read.xlsx(file, rowNames = TRUE)
  common_samples <- intersect(colnames(abundance), metadata$ID)
  abundance <- abundance[, common_samples, drop = FALSE]
  meta <- metadata[metadata$ID %in% common_samples, , drop = FALSE]
  rownames(meta) <- meta$ID
  otu <- otu_table(as.matrix(abundance), taxa_are_rows = TRUE)
  sam <- sample_data(meta)
  phyloseq(otu, sam)
}

# ----------------------------
# 循环分析四类微生物
# ----------------------------
for(microbe in names(microbe_files)){
  message("开始分析: ", microbe)
  
  ps_obj <- make_phyloseq_from_excel(microbe_files[[microbe]], metadata)
  sample_data_df <- data.frame(sample_data(ps_obj))
  
  # ----------------------------
  # α 多样性 (Shannon & Simpson)
  # ----------------------------
  alpha_df <- estimate_richness(ps_obj, measures = c("Shannon", "Simpson"))
  alpha_df$SampleID <- rownames(alpha_df)
  alpha_df <- merge(alpha_df, sample_data_df[, "Group", drop=FALSE], by.x="SampleID", by.y="row.names")
  
  alpha_means <- alpha_df %>%
    group_by(Group) %>%
    summarise(
      mean_Shannon = round(mean(Shannon), 2),
      mean_Simpson = round(mean(Simpson), 2)
    )
  
  p_alpha_shannon <- ggplot(alpha_df, aes(x=Group, y=Shannon, fill=Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha=0.6) +
    labs(title=paste(microbe, "Alpha Diversity (Shannon)"),
         subtitle=paste("Means:", paste(alpha_means$Group, alpha_means$mean_Shannon, collapse = "; "))) +
    theme_bw()
  
  p_alpha_simpson <- ggplot(alpha_df, aes(x=Group, y=Simpson, fill=Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha=0.6) +
    labs(title=paste(microbe, "Alpha Diversity (Simpson)"),
         subtitle=paste("Means:", paste(alpha_means$Group, alpha_means$mean_Simpson, collapse = "; "))) +
    theme_bw()
  
  # 保存 α-diversity 图
  ggsave(filename = file.path(output_dir, paste0(microbe, "_alpha_shannon.pdf")), p_alpha_shannon, width=8, height=6)
  ggsave(filename = file.path(output_dir, paste0(microbe, "_alpha_simpson.pdf")), p_alpha_simpson, width=8, height=6)
  
  # ----------------------------
  # β 多样性 (PCoA Bray-Curtis)
  # ----------------------------
  dist <- phyloseq::distance(ps_obj, method = "bray")
  ord <- ordinate(ps_obj, method="PCoA", distance=dist)
  ord_df <- plot_ordination(ps_obj, ord, type="samples")$data
  if(!"Group" %in% colnames(ord_df)){
    ord_df$Group <- sample_data_df[rownames(ord_df), "Group"]
  }
  
  adonis_res <- adonis2(dist ~ Group, data = sample_data_df)
  adonis_text <- paste0("PERMANOVA: R² = ", round(adonis_res$R2[1],3),
                        ", p = ", signif(adonis_res$`Pr(>F)`[1],3))
  
  p_beta <- ggplot(ord_df, aes(x=Axis.1, y=Axis.2, color=Group)) +
    geom_point(size=3, alpha=0.8) +
    stat_ellipse(level=0.95, linetype=2) +
    labs(title=paste(microbe, "Beta Diversity (PCoA-Bray)"),
         subtitle=adonis_text,
         x="PCoA1", y="PCoA2") +
    theme_bw()
  
  # Bray-Curtis 组间距离箱线图
  dist_mat <- as.matrix(dist)
  n <- nrow(dist_mat)
  group_vec <- sample_data_df$Group
  beta_box <- data.frame()
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      beta_box <- rbind(beta_box, data.frame(
        Dist = dist_mat[i,j],
        Pair = paste(group_vec[i], group_vec[j], sep="-")
      ))
    }
  }
  p_beta_box <- ggplot(beta_box, aes(x=Pair, y=Dist, fill=Pair)) +
    geom_boxplot() +
    theme_bw() +
    labs(title=paste(microbe, "Bray-Curtis distances"), x="Group Pairs", y="Distance") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  # 拼接 β-diversity 图
  p_beta_combined <- p_beta / p_beta_box + plot_layout(heights = c(2,1))
  
  # 保存 β-diversity 拼图
  ggsave(filename = file.path(output_dir, paste0(microbe, "_beta_combined.pdf")), 
         p_beta_combined, width=10, height=10)
  
  message("完成并保存: ", microbe)
}
