# ============================================
# α & β 多样性分析（四类微生物，相对丰度表）
# α：Shannon、Simpson、InvSimpson + 组均值/中位数/组间p值（副标题自动换行）
# β：PCoA (Bray-Curtis) + 置信椭圆 + 样本数/ PERMANOVA p值
# 去掉背景网格线，图形对比明显
# ============================================

# 依赖包
if(!requireNamespace("phyloseq", quietly=TRUE)) install.packages("phyloseq")
if(!requireNamespace("vegan", quietly=TRUE)) install.packages("vegan")
if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if(!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
if(!requireNamespace("ape", quietly=TRUE)) install.packages("ape")
if(!requireNamespace("rlang", quietly=TRUE)) install.packages("rlang")
if(!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if(!requireNamespace("stringr", quietly=TRUE)) install.packages("stringr")

library(phyloseq); library(vegan); library(ggplot2)
library(readxl); library(ape); library(rlang); library(dplyr)
library(stringr)

# -------------------------
# 1. 读取样本元数据
# -------------------------
meta_file <- "E:\\Python\\MI_Analysis\\metagenome\\data_figures\\filtered_data_1percent\\sample_metadata.xlsx"
meta <- read_excel(meta_file)
meta <- as.data.frame(meta)
rownames(meta) <- meta$SampleID

# -------------------------
# 2. 丰度表路径列表
# -------------------------
microbe_files <- list(
  bacteria = "E:\\Python\\MI_Analysis\\metagenome\\data_figures\\filtered_data_1percent\\心梗组_细菌_filtered_1percent.csv",
  archaea  = "E:\\Python\\MI_Analysis\\metagenome\\data_figures\\filtered_data_1percent\\心梗组_古菌_filtered_1percent.csv",
  fungi    = "E:\\Python\\MI_Analysis\\metagenome\\data_figures\\filtered_data_1percent\\心梗组_真菌_filtered_1percent.csv",
  virus    = "E:\\Python\\MI_Analysis\\metagenome\\data_figures\\filtered_data_1percent\\心梗组_病毒(新)_filtered_1percent.csv"
)

# -------------------------
# 3. 输出目录
# -------------------------
out_dir <- "E:\\Python\\MI_Analysis\\metagenome\\data_figures\\filtered_data_1percent\\diversity_analysis_OTU"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# -------------------------
# 4. 批量分析
# -------------------------
for(name in names(microbe_files)){
  cat("\n===== Processing:", name, "=====\n")
  
  # 读取丰度表
  otu <- read.csv(microbe_files[[name]], header=TRUE, row.names=1, check.names = FALSE)
  
  # 样本交集及对齐
  shared_samples <- intersect(colnames(otu), meta$SampleID)
  if(length(shared_samples)==0) stop(paste("No matching samples for", name))
  otu <- otu[, shared_samples, drop=FALSE]
  meta_sub <- meta[shared_samples, , drop=FALSE]
  
  ps <- phyloseq(
    otu_table(as.matrix(otu), taxa_are_rows=TRUE),
    sample_data(meta_sub)
  )
  
  # -------------------------
  # α 多样性
  # -------------------------
  alpha_df <- estimate_richness(ps, measures=c("Shannon","Simpson","InvSimpson"))
  alpha_df$Group <- sample_data(ps)$Group
  
  write.csv(alpha_df, file.path(out_dir, paste0("alpha_diversity_", name, ".csv")), row.names=TRUE)
  
  metrics <- c("Shannon","Simpson","InvSimpson")
  for(metric in metrics){
    
    # 统计信息
    stat_df <- alpha_df %>%
      group_by(Group) %>%
      summarise(mean=mean(!!sym(metric)), median=median(!!sym(metric)))
    
    # 组间显著性
    if(length(unique(alpha_df$Group))==2){
      test <- wilcox.test(as.formula(paste(metric,"~ Group")), data=alpha_df)
      pval_text <- paste0("p = ", signif(test$p.value, 3))
    } else {
      test <- kruskal.test(as.formula(paste(metric,"~ Group")), data=alpha_df)
      pval_text <- paste0("p = ", signif(test$p.value, 3))
    }
    
    subtitle_text <- paste0(
      paste(stat_df$Group, ": mean=", round(stat_df$mean,2), ", median=", round(stat_df$median,2), collapse="; "),
      " | ", pval_text
    )
    # 自动换行副标题
    subtitle_text <- str_wrap(subtitle_text, width = 60)
    
    # 绘图，明显区分颜色，不透明，去网格线
    p <- ggplot(alpha_df, aes(x=Group, y=!!sym(metric), fill=Group)) +
      geom_violin(trim=FALSE, alpha=1) +
      geom_boxplot(width=0.12, outlier.size=0.5, fill="white") +
      scale_fill_manual(values=c("Control"="#1F78B4", "MI"="#E31A1C")) +
      theme_classic() +
      theme(legend.position="none") +
      ggtitle(paste("Alpha diversity -", metric, "(", name, ")", sep=""),
              subtitle=subtitle_text)
    
    ggsave(file.path(out_dir, paste0("alpha_", metric, "_", name, ".pdf")), p, width=6, height=5)
    cat("Alpha test (", metric, ") for", name, ":\n"); print(test)
  }
  
  # -------------------------
  # β PCoA (Bray-Curtis)
  # -------------------------
  ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
  otu_mat_rel <- as.matrix(otu_table(ps_rel))
  dist_bc <- vegdist(t(otu_mat_rel), method="bray")
  
  pcoa_res <- pcoa(dist_bc)
  pco_df <- data.frame(sample=rownames(pcoa_res$vectors),
                       PC1=pcoa_res$vectors[,1],
                       PC2=pcoa_res$vectors[,2])
  pco_df$Group <- sample_data(ps_rel)$Group
  
  # PERMANOVA
  if("Group" %in% colnames(meta_sub)){
    adonis_res <- adonis2(dist_bc ~ Group, data=data.frame(sample_data(ps_rel)))
    pval_text <- paste0("PERMANOVA p = ", signif(adonis_res$`Pr(>F)`[1], 3))
  } else {
    pval_text <- ""
  }
  
  # 样本数
  group_counts <- table(pco_df$Group)
  n_text <- paste(names(group_counts), "n=", group_counts, collapse="; ")
  
  # 绘图，置信椭圆，去网格线
  p_pcoa <- ggplot(pco_df, aes(x=PC1, y=PC2, color=Group)) +
    geom_point(size=3) +
    stat_ellipse(level=0.95, linetype=2) +
    scale_color_manual(values=c("Control"="#1F78B4", "MI"="#E31A1C")) +
    theme_classic() +
    labs(title=paste("PCoA (Bray-Curtis) -", name),
         subtitle=paste(n_text, "|", pval_text),
         x=paste0("PC1 (", round(pcoa_res$values$Relative_eig[1]*100,1), "%)"),
         y=paste0("PC2 (", round(pcoa_res$values$Relative_eig[2]*100,1), "%)"))
  
  ggsave(file.path(out_dir, paste0("beta_pcoa_bray_", name, ".pdf")), p_pcoa, width=6, height=5)
  cat("PERMANOVA for", name, ":\n"); print(adonis_res)
}
