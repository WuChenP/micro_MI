library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(vegan)
library(ggpubr)

# ================================
# 1. 文件路径（20%过滤后的新文件）
# ================================
abundance_files <- list(
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/古菌_filtered_20percent.csv",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/细菌_filtered_20percent.csv",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/真菌_filtered_20percent.csv",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/病毒_filtered_20percent.csv"
)

metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/sample_metadata.xlsx"
plot_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent_change/diversity_plots_OTU"

if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ================================
# 2. 读取样本元数据
# ================================
metadata <- readxl::read_excel(metadata_file)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$ID

# ================================
# 3. 读取 CSV 并转换为 样本 × 物种 数值矩阵
# ================================
process_abundance <- function(file){
  df <- readr::read_csv(file, show_col_types = FALSE)
  df <- as.data.frame(df)
  
  taxa <- df[[1]]           # 第一列为物种名
  df <- df[,-1, drop=FALSE] # 只保留样本列
  rownames(df) <- taxa
  
  mat <- t(as.matrix(df))
  mat <- apply(mat, 2, as.numeric)
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat)
  
  rownames(mat) <- colnames(df)
  colnames(mat) <- taxa
  
  return(mat)
}

# ================================
# 4. α/β 多样性循环
# ================================
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
  
  # Shannon 小提琴图 → PDF
  pdf(file = paste0(plot_dir, "/", microbe, "_alpha_shannon.pdf"), width = 6, height = 4)
  print(
    ggplot(alpha_df, aes(x=Group, y=Shannon, fill=Group)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1, fill="white", alpha=0.7) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14)
      ) +
      ggtitle(paste(microbe, "Shannon Diversity"))
  )
  dev.off()
  
  # Simpson 小提琴图 → PDF
  pdf(file = paste0(plot_dir, "/", microbe, "_alpha_simpson.pdf"), width = 6, height = 4)
  print(
    ggplot(alpha_df, aes(x=Group, y=Simpson, fill=Group)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1, fill="white", alpha=0.7) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14)
      ) +
      ggtitle(paste(microbe, "Simpson Diversity"))
  )
  dev.off()
  
  # β 多样性（Bray-Curtis + PCoA）
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
  
  # PCoA + Ellipse
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
  
  # 样本到组中心距离
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
  
  # 合并图 → PDF
  pdf(file = paste0(plot_dir, "/", microbe, "_beta_pcoa_permanova.pdf"), width = 10, height = 5)
  print(ggarrange(ellipse_plot, box_plot, ncol=2, widths=c(2,1)))
  dev.off()
  
  pcoa_list[[microbe]] <- pcoa_df
}

cat("所有分析完成，PDF 图已保存至：", plot_dir, "\n")
