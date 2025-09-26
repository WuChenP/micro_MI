library(readxl)
library(dplyr)
library(ggplot2)
library(vegan)
library(ggpubr)  # 用于组合图形

# 文件路径
abundance_files <- list(
  archaea  = "E:/Python/MI_Analysis/origin_data/心梗组_古菌.xlsx",
  bacteria = "E:/Python/MI_Analysis/origin_data/心梗组_细菌.xlsx",
  fungi    = "E:/Python/MI_Analysis/origin_data/心梗组_真菌.xlsx",
  virus    = "E:/Python/MI_Analysis/origin_data/心梗组_病毒(新).xlsx"
)

metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/sample_metadata.xlsx"
plot_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/diversity_plot"

# 读取样本元数据
metadata <- read_excel(metadata_file)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$ID

# 函数：读取丰度表并处理成 样本 × 物种 数值矩阵
process_abundance <- function(file){
  df <- read_excel(file)
  df <- as.data.frame(df)
  
  # 设置行名为物种ID
  rownames(df) <- df$ID
  df$ID <- NULL
  
  # 转置，得到样本 × 物种矩阵
  mat <- t(df)
  
  # 转数值
  mat <- apply(mat, 2, as.numeric)
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat)
  
  # 修复行名和列名
  rownames(mat) <- colnames(df)
  colnames(mat) <- rownames(df)
  
  return(mat)
}

# 保存结果用列表
alpha_list <- list()
pcoa_list <- list()

for(microbe in names(abundance_files)){
  cat("Processing:", microbe, "\n")
  
  abundance <- process_abundance(abundance_files[[microbe]])
  
  # α 多样性
  shannon <- diversity(abundance, index = "shannon")
  simpson <- diversity(abundance, index = "simpson")
  alpha_df <- data.frame(
    Sample = rownames(abundance),
    Shannon = shannon,
    Simpson = simpson,
    Group = metadata[rownames(abundance), "Group"]
  )
  alpha_list[[microbe]] <- alpha_df
  
  # 绘制 α 多样性小提琴图（干净背景）
  p1 <- ggplot(alpha_df, aes(x=Group, y=Shannon, fill=Group)) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1, fill="white", alpha=0.7) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size=12),
      axis.title = element_text(size=14)
    ) +
    ggtitle(paste(microbe, "Shannon Diversity"))
  
  ggsave(filename = paste0(plot_dir, "/", microbe, "_alpha_shannon.png"), plot = p1, width=6, height=4)
  
  # β 多样性（Bray-Curtis）
  bray <- vegdist(abundance, method = "bray")
  pcoa <- cmdscale(bray, k=2, eig=TRUE)
  pcoa_df <- data.frame(
    Sample = rownames(abundance),
    PCoA1 = pcoa$points[,1],
    PCoA2 = pcoa$points[,2],
    Group = metadata[rownames(abundance), "Group"]
  )
  
  # PERMANOVA
  perm <- adonis2(bray ~ Group, data=metadata[rownames(abundance), , drop=FALSE], permutations = 999)
  r2 <- round(perm$R2[1], 3)
  pval <- perm$`Pr(>F)`[1]
  df_text <- paste0("PERMANOVA: R²=", r2, ", P=", pval)
  
  # 绘制置信椭圆 + 散点图
  ellipse_plot <- ggplot(pcoa_df, aes(x=PCoA1, y=PCoA2, color=Group)) +
    geom_point(size=3) +
    stat_ellipse(type="norm", level=0.95, linetype=2) +
    theme_classic() +
    theme(
      axis.text = element_text(size=12),
      axis.title = element_text(size=14),
      legend.title = element_text(size=12)
    ) +
    ggtitle(paste(microbe, "PCoA (Bray-Curtis)\n", df_text))
  
  # 计算每个样本到组中心的 Bray-Curtis 距离
  group_centers <- aggregate(cbind(PCoA1, PCoA2) ~ Group, data=pcoa_df, FUN=mean)
  dist_to_center <- sapply(1:nrow(pcoa_df), function(i){
    grp <- pcoa_df$Group[i]
    sqrt((pcoa_df$PCoA1[i]-group_centers$PCoA1[group_centers$Group==grp])^2 +
           (pcoa_df$PCoA2[i]-group_centers$PCoA2[group_centers$Group==grp])^2)
  })
  box_df <- data.frame(Group=pcoa_df$Group, Distance=dist_to_center)
  
  box_plot <- ggplot(box_df, aes(x=Group, y=Distance, fill=Group)) +
    geom_boxplot(alpha=0.7) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size=12),
      axis.title = element_text(size=14)
    ) +
    ggtitle("Distance to Group Center")
  
  # 合并两个图
  combined <- ggarrange(ellipse_plot, box_plot, ncol=2, widths=c(2,1))
  
  ggsave(filename = paste0(plot_dir, "/", microbe, "_beta_pcoa_permanova.png"), plot = combined, width=10, height=5)
  
  pcoa_list[[microbe]] <- pcoa_df
}

cat("所有分析完成，图已保存至：", plot_dir, "\n")
