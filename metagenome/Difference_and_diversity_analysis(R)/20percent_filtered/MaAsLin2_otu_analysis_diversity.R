# ==========================================
# α & β 多样性分析 - 四类微生物（整合版，自适应PDF）
# ==========================================

library(tidyverse)
library(vegan)
library(readr)
library(readxl)
library(ggplot2)
library(stringr)

# -----------------------------
# 文件路径
# -----------------------------
files <- list(
  archaea = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/古菌_filtered_20percent.csv",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/细菌_filtered_20percent.csv",
  fungi   = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/真菌_filtered_20percent.csv",
  virus   = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/病毒_filtered_20percent.csv"
)

metadata_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/MaAsLin2/sample_metadata.xlsx"
output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data/MaAsLin2/diversity_analysis/"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# 读取样本元数据
# -----------------------------
metadata <- read_excel(metadata_file)
metadata$SampleID <- as.character(metadata$SampleID)

# -----------------------------
# 配色统一
# -----------------------------
group_colors <- c("Control" = "#1f77b4", "MI" = "#ff7f0e")

# -----------------------------
# α多样性计算函数
# -----------------------------
compute_alpha <- function(df, metadata){
  df_long <- df %>%
    column_to_rownames("Taxonomy") %>%
    t() %>%
    as.data.frame()
  
  alpha <- data.frame(
    SampleID = rownames(df_long),
    Shannon = diversity(df_long, index = "shannon"),
    Simpson = diversity(df_long, index = "simpson")
  )
  
  alpha <- alpha %>% left_join(metadata, by = "SampleID")
  return(alpha)
}

# -----------------------------
# β多样性计算函数
# -----------------------------
compute_beta <- function(df, metadata){
  df_long <- df %>%
    column_to_rownames("Taxonomy") %>%
    t() %>%
    as.data.frame()
  
  dist_mat <- vegdist(df_long, method = "bray")
  pcoa_res <- cmdscale(dist_mat, eig = TRUE, k = 2)
  
  beta <- data.frame(
    SampleID = rownames(pcoa_res$points),
    PCoA1 = pcoa_res$points[,1],
    PCoA2 = pcoa_res$points[,2]
  )
  
  beta <- beta %>% left_join(metadata, by = "SampleID")
  return(list(beta = beta, dist_mat = dist_mat))
}

# -----------------------------
# α多样性小提琴图函数（自适应PDF，副标题换行）
# -----------------------------
plot_alpha <- function(alpha_df, index_col, microbe, output_file){
  stats_df <- alpha_df %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      mean_val = mean(.data[[index_col]], na.rm = TRUE),
      sd_val = sd(.data[[index_col]], na.rm = TRUE)
    )
  
  subtitle_text <- paste0(
    paste0(stats_df$Group, ": n=", stats_df$n,
           ", mean=", round(stats_df$mean_val,3),
           ", sd=", round(stats_df$sd_val,3)),
    collapse = " | "
  )
  subtitle_text <- str_wrap(subtitle_text, width = 50)
  
  p <- ggplot(alpha_df, aes(x = Group, y = .data[[index_col]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 14) +
    ggtitle(paste0(microbe, " α diversity (", index_col, ")"),
            subtitle = subtitle_text) +
    ylab(index_col) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, face = "italic", lineheight = 1.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(c(2,2,2,2), "cm")
    )
  
  max_subtitle_len <- max(str_count(str_split(subtitle_text, "\n")[[1]], "\\S+"))
  n_groups <- length(unique(alpha_df$Group))
  plot_width <- max(6, n_groups*2, max_subtitle_len/5)
  
  y_range <- range(alpha_df[[index_col]], na.rm = TRUE)
  plot_height <- max(6, 4 + diff(y_range)*1.5)
  
  ggsave(filename = output_file, plot = p, width = plot_width, height = plot_height)
}

# -----------------------------
# β多样性 PCoA 图函数（简化副标题，只显示n和PERMANOVA p值）
# -----------------------------
plot_beta <- function(beta_list, microbe, output_file){
  beta_df <- beta_list$beta
  dist_mat <- beta_list$dist_mat
  
  # 每组样本数量
  stats_df <- beta_df %>%
    group_by(Group) %>%
    summarise(n = n())
  
  # PERMANOVA检验
  perm_res <- adonis2(dist_mat ~ Group, data = beta_df)
  perm_p <- perm_res$`Pr(>F)`[1]
  
  # 副标题
  subtitle_text <- paste0(
    paste0(stats_df$Group, ": n=", stats_df$n),
    collapse = " | "
  )
  subtitle_text <- paste0(subtitle_text, " | PERMANOVA p=", signif(perm_p,3))
  subtitle_text <- str_wrap(subtitle_text, width = 80)
  
  max_line_len <- max(str_count(str_split(subtitle_text, "\n")[[1]], "\\S+"))
  n_groups <- length(unique(beta_df$Group))
  plot_width <- max(6, n_groups*2, max_line_len/5)
  
  x_range <- range(beta_df$PCoA1, na.rm = TRUE)
  y_range <- range(beta_df$PCoA2, na.rm = TRUE)
  plot_height <- min(max(6, 6 + max(diff(x_range), diff(y_range))*1.5), 12)
  
  p <- ggplot(beta_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(aes(fill = Group), geom = "polygon", level = 0.95, alpha = 0.15, color = "black") +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 14) +
    ggtitle(paste0(microbe, " β diversity (Bray-Curtis PCoA)"),
            subtitle = subtitle_text) +
    theme(
      plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(size = 10, face = "italic", lineheight = 1.3),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.margin = unit(c(2,2,2,2), "cm")
    )
  
  ggsave(filename = output_file, plot = p, width = plot_width, height = plot_height)
}

# -----------------------------
# 循环处理四类微生物
# -----------------------------
for (microbe in names(files)){
  cat("Processing:", microbe, "\n")
  
  df <- read_csv(files[[microbe]], show_col_types = FALSE)
  
  # α多样性
  alpha <- compute_alpha(df, metadata)
  plot_alpha(alpha, "Shannon", microbe, paste0(output_dir, microbe, "_alpha_Shannon.pdf"))
  plot_alpha(alpha, "Simpson", microbe, paste0(output_dir, microbe, "_alpha_Simpson.pdf"))
  
  # β多样性
  beta <- compute_beta(df, metadata)
  plot_beta(beta, microbe, paste0(output_dir, microbe, "_beta_PCoA.pdf"))
}
