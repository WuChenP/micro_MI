# ==================================
# Alpha diversity analysis (Richness, Shannon, Simpson, Chao1, ACE, Coverage)
# ==================================

# ---- 设置 CRAN 镜像 ----
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# ---- 检查并安装依赖包 ----
check_and_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

check_and_install("phyloseq")
check_and_install("vegan")
check_and_install("ggplot2")
check_and_install("ggpubr")
check_and_install("Rmisc")

# ---- 设置工作目录 ----
setwd("D:/机器学习与生物信息学/北京朝阳医院项目/R")   # 修改为你自己的保存路径

# ---- 读取数据 ----
df <- read.csv("R_data.csv", header = TRUE, row.names = 2)

# ---- 构建 OTU 矩阵 (ASV/OTU × 样本) ----
otu <- as.matrix(df[, grep("AMI|CON", colnames(df))])

# ---- 构建分组信息 ----
sample_ids <- colnames(otu)
group <- ifelse(grepl("^AMI", sample_ids), "AMI", "CON")
sample_df <- data.frame(
  Sample = sample_ids,
  Group = group
)

# ---- 计算 Alpha 多样性指数 ----
richness <- specnumber(t(otu))                  # Richness
shannon  <- diversity(t(otu), index = "shannon")# Shannon
simpson  <- diversity(t(otu), index = "simpson")# Simpson
chao1    <- estimateR(t(otu))["S.chao1", ]      # Chao1
ace      <- estimateR(t(otu))["S.ACE", ]        # ACE

# Coverage (Good's Coverage) 手动计算
coverage <- sapply(1:ncol(otu), function(i) {
  sample_counts <- otu[, i]
  N <- sum(sample_counts)
  F1 <- sum(sample_counts == 1)
  if (N == 0) return(NA)
  return(1 - F1 / N)
})

# ---- 合并结果 ----
alpha_div <- data.frame(
  Sample   = colnames(otu),
  Group    = group,
  Richness = richness,
  Shannon  = shannon,
  Simpson  = simpson,
  Chao1    = chao1,
  ACE      = ace,
  Coverage = coverage
)

# 保存结果表
write.csv(alpha_div, "alpha_diversity_results.csv", row.names = FALSE)

# ---- 绘图函数 ----
plot_alpha <- function(df, metric, ylab_text) {
  p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    labs(title = paste(ylab_text, "between groups"),
         x = "Group", y = ylab_text) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  return(p)
}

# ---- 生成并保存所有图 ----
metrics <- c("Richness", "Shannon", "Simpson", "Chao1", "ACE", "Coverage")

for (m in metrics) {
  p <- plot_alpha(alpha_div, m, m)
  ggsave(filename = paste0("alpha_", m, ".png"), plot = p,
         width = 6, height = 5, dpi = 300)
}

message("✅ Alpha diversity analysis done. Results saved in: ",
        getwd())
