# ===============================
# Alpha 多样性分析
# ===============================

# 设置工作路径
setwd("D:/机器学习与生物信息学/北京朝阳医院项目/R")

# ---- 加载依赖包 ----
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

pkgs <- c("phyloseq", "ggplot2", "vegan", "microbiome")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ---- 读取数据 ----
df <- read.csv("D:/机器学习与生物信息学/北京朝阳医院项目/R/R_data.csv",
               header = TRUE, row.names = 2)

# 构建 OTU 矩阵（样本 × ASV）
otu_mat <- as.matrix(df[, grep("AMI|CON", colnames(df))])

# 分类学注释
tax_mat <- as.matrix(df[, c("Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus", "Species")])

# 分组信息
sample_ids <- colnames(otu_mat)
group <- ifelse(grepl("^AMI", sample_ids), "AMI", "Control")

sample_df <- data.frame(Group = group, row.names = sample_ids)

# 构建 phyloseq 对象
otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
tax <- tax_table(tax_mat)
samp <- sample_data(sample_df)
phy <- phyloseq(otu, tax, samp)

# ---- Alpha 多样性计算 ----
alpha_div <- estimate_richness(phy, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))

# 添加分组信息
alpha_div$Group <- sample_df$Group

# Good's coverage 计算
goods <- function(x) {
  1 - sum(x == 1) / sum(x)
}
alpha_div$Coverage <- apply(otu_mat, 2, goods)

# ---- 保存结果表 ----
write.csv(alpha_div, "Alpha_diversity_results.csv", row.names = TRUE)

# ---- 绘图函数 ----
plot_alpha <- function(metric, ylab_name) {
  p <- ggplot(alpha_div, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(title = paste0(metric, " Diversity"), x = "Group", y = ylab_name) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(metric, "_alpha_diversity.png"), p, width = 6, height = 5)
  return(p)
}

# ---- 绘制并保存图 ----
plot_alpha("Observed", "Observed ASVs")
plot_alpha("Chao1", "Chao1 Index")
plot_alpha("ACE", "ACE Index")
plot_alpha("Shannon", "Shannon Index")
plot_alpha("Simpson", "Simpson Index")
plot_alpha("Coverage", "Good's Coverage")

message("✅ Alpha 多样性结果表已保存到 Alpha_diversity_results.csv，图已保存到当前目录")
