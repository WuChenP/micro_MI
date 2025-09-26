# ================================
# Beta 多样性分析 (Bray-Curtis + PCoA + PERMANOVA)
# ================================

# ---- 加载包 ----
if (!require("vegan")) install.packages("vegan", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!require("ggplot2")) install.packages("ggplot2", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!require("ggpubr")) install.packages("ggpubr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

library(vegan)
library(ggplot2)
library(ggpubr)

# ---- 设置工作目录  ----
setwd("D:/机器学习与生物信息学/北京朝阳医院项目/R")

# ---- 读取数据 ----
df <- read.csv("R_data.csv", header = TRUE, row.names = 1, check.names = FALSE)

# ---- 提取数值型丰度矩阵 ----
otu <- df[, sapply(df, is.numeric)]
otu <- t(otu)  # 转置为样本 × 物种

# ---- 生成分组信息 (AMI vs Control) ----
sample_ids <- rownames(otu)
group <- ifelse(grepl("AMI", sample_ids, ignore.case = TRUE), "AMI", "Control")
metadata <- data.frame(SampleID = sample_ids, Group = factor(group, levels = c("Control", "AMI")))

# ---- 计算 Bray-Curtis 距离 ----
bray_dist <- vegdist(otu, method = "bray")

# ---- PCoA ----
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  SampleID = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2],
  Group = metadata$Group
)

# ---- 计算方差解释率 ----
eig <- pcoa$eig / sum(pcoa$eig)
xlab <- paste0("PC1 (", round(eig[1] * 100, 2), "%)")
ylab <- paste0("PC2 (", round(eig[2] * 100, 2), "%)")

# ---- 绘图 ----
p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(title = "PCoA (Bray-Curtis)", x = xlab, y = ylab) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# 保存图片
ggsave("Beta_diversity_PCoA.png", plot = p, width = 7, height = 6, dpi = 300)

# ---- PERMANOVA 检验 ----
permanova <- adonis2(bray_dist ~ Group, data = metadata)
print(permanova)

# 保存 PERMANOVA 结果
sink("Beta_diversity_PERMANOVA.txt")
print(permanova)
sink()
