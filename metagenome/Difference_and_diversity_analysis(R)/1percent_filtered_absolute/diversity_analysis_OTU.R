# ======================================================
# 四类微生物 α + β 多样性分析（含P值、R²、置信椭圆）
# 顶刊风格（P值放副标题）
# ======================================================

library(phyloseq)
library(vegan)
library(ggplot2)
library(openxlsx)
library(tidyverse)
library(reshape2)
library(ggpubr)

# ----------------------------
# 1. 输入路径
# ----------------------------
input_files <- c(
  "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_bacteria_species_level.csv",
  "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_archaea_species_level.csv",
  "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_fungi_species_level.csv",
  "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/species_level_data/心梗组_virus_species_level.csv"
)
names(input_files) <- c("Bacteria", "Archaea", "Fungi", "Virus")

meta_file <- "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/sample_metadata.xlsx"
output_dir <- "E:/Python/MI_Analysis/metagenome/Absolute_abundance_analysis/filtered_data_1percent/diversity_analysis_OTU"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------
# 2. 读取元数据
# ----------------------------
metadata <- read.xlsx(meta_file)
metadata$Group <- as.factor(metadata$Group)
rownames(metadata) <- metadata$SampleID

# ----------------------------
# 3. 创建保存统计结果的列表
# ----------------------------
alpha_results <- list()
beta_results <- list()

# ----------------------------
# 4. α 多样性分析函数（稳妥版）
# ----------------------------
calc_alpha_div <- function(file_path, species_name){
  cat("开始 α 多样性分析:", species_name, "\n")
  
  data <- read.csv(file_path, header=TRUE, check.names=FALSE)
  if("ID" %in% colnames(data)){
    rownames(data) <- data$ID
    data <- data[,-1, drop=FALSE]
  }
  
  common_samples <- intersect(colnames(data), metadata$SampleID)
  if(length(common_samples)==0) stop("⚠️ OTU 表与元数据没有匹配样本！")
  data <- data[, common_samples, drop=FALSE]
  meta <- metadata[common_samples, , drop=FALSE]
  
  ps <- phyloseq(otu_table(as.matrix(data), taxa_are_rows=TRUE),
                 sample_data(meta))
  
  alpha_metrics <- estimate_richness(ps,
                                     measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson"))
  alpha_metrics$SampleID <- rownames(alpha_metrics)
  alpha_metrics$Group <- meta[rownames(alpha_metrics),"Group"]
  alpha_metrics$Group <- droplevels(alpha_metrics$Group) # 去掉多余的 factor level
  
  # 打印分组情况
  cat("分组情况 (table):\n")
  print(table(alpha_metrics$Group))
  cat("Group levels:\n")
  print(levels(alpha_metrics$Group))
  
  # 计算每个指标的 Wilcoxon P值（稳妥处理 NA / 全零）
  metrics <- c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson")
  p_values <- sapply(metrics, function(m){
    x <- alpha_metrics[[m]][alpha_metrics$Group == levels(alpha_metrics$Group)[1]]
    y <- alpha_metrics[[m]][alpha_metrics$Group == levels(alpha_metrics$Group)[2]]
    # 只用非 NA 样本
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    if(length(x) < 1 | length(y) < 1){
      return(NA)
    } else {
      wilcox.test(x, y, exact=FALSE)$p.value
    }
  })
  
  # 保存α多样性结果
  alpha_results[[species_name]] <<- data.frame(
    Microorganism = species_name,
    Metric = metrics,
    P_value = p_values,
    stringsAsFactors = FALSE
  )
  
  subtitle_text <- paste0("Wilcoxon P-values: ", paste(metrics, ifelse(is.na(p_values), "NA", signif(p_values,3)), collapse=", "))
  
  # 输出 Excel（不含 Number）
  write.xlsx(alpha_metrics,
             file.path(output_dir,paste0(species_name,"_alpha_diversity.xlsx")),
             rowNames=TRUE)
  
  # 转长格式
  alpha_long <- melt(alpha_metrics,
                     id.vars=c("SampleID","Group"),
                     measure.vars=metrics)
  
  # 绘图
  p <- ggplot(alpha_long, aes(x=Group, y=value, fill=Group)) +
    geom_violin(trim=FALSE, alpha=0.6, color=NA) +
    geom_boxplot(width=0.12, color="black", alpha=0.9, outlier.shape=NA) +
    facet_wrap(~variable, scales="free_y", ncol=3) +
    scale_fill_manual(values=c("Control"="#4DAF4A","MI"="#E41A1C")) +
    theme_minimal(base_size=13) +
    theme(
      panel.grid=element_blank(),
      strip.background=element_rect(fill="#F0F0F0", color=NA),
      strip.text=element_text(size=11, face="bold"),
      legend.position="none",
      axis.title.x=element_blank(),
      axis.title.y=element_text(face="bold"),
      plot.title=element_text(hjust=0.5, face="bold", size=14),
      plot.subtitle=element_text(hjust=0.5, size=11)
    ) +
    labs(title=paste0("α Diversity - ", species_name),
         subtitle=subtitle_text,
         y="Index Value")
  
  ggsave(file.path(output_dir,paste0(species_name,"_alpha_diversity_violin_boxplot.pdf")),
         p, width=12, height=6)
  
  cat("α 多样性分析完成:", species_name, "\n\n")
}

# ----------------------------
# 5. β 多样性分析函数
# ----------------------------
calc_beta_div <- function(file_path, species_name){
  cat("开始 β 多样性分析:", species_name, "\n")
  
  data <- read.csv(file_path, header=TRUE, check.names=FALSE)
  if("ID" %in% colnames(data)){
    rownames(data) <- data$ID
    data <- data[,-1, drop=FALSE]
  }
  
  common_samples <- intersect(colnames(data), metadata$SampleID)
  if(length(common_samples)==0) stop("⚠️ OTU 表与元数据没有匹配样本！")
  data <- data[, common_samples, drop=FALSE]
  meta <- metadata[common_samples, , drop=FALSE]
  
  dist_bc <- vegdist(t(data), method="bray")
  pcoa_res <- cmdscale(dist_bc, eig=TRUE, k=2)
  pcoa_df <- data.frame(
    SampleID=rownames(pcoa_res$points),
    Axis1=pcoa_res$points[,1],
    Axis2=pcoa_res$points[,2],
    Group=meta$Group
  )
  
  eig1 <- round(100*pcoa_res$eig[1]/sum(pcoa_res$eig),2)
  eig2 <- round(100*pcoa_res$eig[2]/sum(pcoa_res$eig),2)
  
  adonis_res <- adonis2(dist_bc ~ Group, data=meta)
  R2 <- round(adonis_res$R2[1],3)
  pval <- adonis_res$`Pr(>F)`[1]
  
  # 保存β多样性结果
  beta_results[[species_name]] <<- data.frame(
    Microorganism = species_name,
    R_squared = R2,
    P_value = pval,
    stringsAsFactors = FALSE
  )
  
  subtitle_text <- paste0("PERMANOVA P = ", formatC(pval, format="e", digits=2), ", R² = ", R2)
  
  p <- ggplot(pcoa_df, aes(x=Axis1, y=Axis2, color=Group)) +
    geom_point(size=3.5, alpha=0.85) +
    stat_ellipse(level=0.95, linetype=1, linewidth=0.8) +
    scale_color_manual(values=c("Control"="#4DAF4A","MI"="#E41A1C")) +
    theme_minimal(base_size=13) +
    theme(
      panel.grid=element_blank(),
      axis.title=element_text(face="bold"),
      plot.title=element_text(hjust=0.5, face="bold", size=14),
      plot.subtitle=element_text(hjust=0.5, size=11)
    ) +
    labs(title=paste0("β Diversity (PCoA - Bray-Curtis) - ", species_name),
         subtitle=subtitle_text,
         x=paste0("PCoA1 (", eig1, "%)"),
         y=paste0("PCoA2 (", eig2, "%)"))
  
  ggsave(file.path(output_dir,paste0(species_name,"_beta_diversity_PCoA.pdf")),
         p, width=7, height=6)
  
  cat("β 多样性分析完成:", species_name, "\n\n")
}

# ----------------------------
# 6. 保存统计结果到txt文件的函数
# ----------------------------
save_statistical_results <- function() {
  # 合并α多样性结果
  alpha_df <- do.call(rbind, alpha_results)
  rownames(alpha_df) <- NULL
  
  # 合并β多样性结果
  beta_df <- do.call(rbind, beta_results)
  rownames(beta_df) <- NULL
  
  # 创建输出文件路径
  output_file <- file.path(output_dir, "diversity_statistical_results.txt")
  
  # 写入文件
  sink(output_file)
  
  cat("===============================================\n")
  cat("       微生物多样性分析统计结果\n")
  cat("===============================================\n\n")
  
  cat("α多样性分析 - Wilcoxon检验P值:\n")
  cat("-----------------------------------------------\n")
  print(alpha_df, row.names = FALSE)
  cat("\n\n")
  
  cat("β多样性分析 - PERMANOVA结果:\n")
  cat("-----------------------------------------------\n")
  print(beta_df, row.names = FALSE)
  cat("\n\n")
  
  cat("文件生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  sink()
  
  cat("✅ 统计结果已保存到:", output_file, "\n")
}

# ----------------------------
# 7. 主执行流程
# ----------------------------
for(sp in names(input_files)){
  calc_alpha_div(input_files[[sp]], sp)
  calc_beta_div(input_files[[sp]], sp)
}

# 保存统计结果
save_statistical_results()

cat("✅ 所有 α + β 多样性分析完成！结果输出到：", output_dir, "\n")