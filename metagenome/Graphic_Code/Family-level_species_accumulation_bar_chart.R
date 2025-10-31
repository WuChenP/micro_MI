# ================================
# Family 水平物种组成堆叠柱状图（MI vs Control）
# Nature 顶刊风格最终版（渐进色 + 固定 Others 灰色 + 显示组样本数在同一行）
# ================================

library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(grid)  # unit() 需要

# -------------------------------
# 1. 设置路径
# -------------------------------
family_path <- "D:/PythonProject/micro_MI/metagenome/data_figures/filtered_data_1percent_change/family_data/virus_family_no_HF.xlsx"
group_path  <- "D:/PythonProject/micro_MI/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
output_path <- "D:/PythonProject/micro_MI/metagenome/Graphic/Family_Barplot/Family_Barplot.pdf"

# -------------------------------
# 2. 读取数据
# -------------------------------
family_df <- read_excel(family_path)
group_df <- read_excel(group_path)

# -------------------------------
# 3. 宽表转成长表
# -------------------------------
family_long <- family_df %>%
  pivot_longer(
    cols = -Family,
    names_to = "SampleID",
    values_to = "Abundance"
  )

# -------------------------------
# 4. 合并分组信息
# -------------------------------
data_merged <- family_long %>%
  left_join(group_df, by = "SampleID")

# -------------------------------
# 5. 计算相对丰度（百分比）
# -------------------------------
data_norm <- data_merged %>%
  group_by(SampleID) %>%
  mutate(RelAbund = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# -------------------------------
# 6. 选取平均丰度最高的前12个 Family
# -------------------------------
top12 <- data_norm %>%
  group_by(Family) %>%
  summarise(mean_abund = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 12) %>%
  pull(Family)

data_top <- data_norm %>%
  mutate(Family = ifelse(Family %in% top12, Family, "Others"))

# -------------------------------
# 7. 按组内丰度排序样本
# -------------------------------
sample_order <- data_top %>%
  group_by(Group, SampleID) %>%
  summarise(total = sum(RelAbund), .groups = "drop") %>%
  arrange(Group, desc(total)) %>%
  pull(SampleID)

data_top$SampleID <- factor(data_top$SampleID, levels = sample_order)

# -------------------------------
# 8. 设置渐进色 + 固定 Others 灰色
# -------------------------------
# 前12个Family名字
top12_names <- top12

# 创建渐进色
palette_top12 <- colorRampPalette(c("#4E79A7", "#76B7B2", "#59A14F", "#EDC948", "#F28E2B", "#E15759"))(length(top12_names))

# 命名向量，对应 Family 名称
names(palette_top12) <- top12_names

# 合并 Others 灰色
palette_family <- c(palette_top12, "Others" = "#B0B0B0")

# -------------------------------
# 9. 计算每组的样本数量
# -------------------------------
group_count <- data_top %>%
  group_by(Group) %>%
  summarise(group_size = n_distinct(SampleID))

# -------------------------------
# 10. 绘制堆叠柱状图（Nature 高级风格）
# -------------------------------
p <- ggplot(data_top, aes(x = SampleID, y = RelAbund, fill = Family)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.05) +  # 极细白边
  facet_wrap(~ Group, scales = "free_x", nrow = 1, labeller = labeller(Group = function(x) {
    # 查找组名对应的数量并显示在同一行
    group_label <- paste(x, "(n =", group_count$group_size[group_count$Group == x], ")")
    return(group_label)
  })) +
  scale_fill_manual(values = palette_family) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Family") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines"),                                # 最紧凑 facet 间距
    strip.background = element_rect(fill = "white", color = "gray30", linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold", family = "Helvetica", color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold", family = "Helvetica", size = 11),
    legend.text = element_text(family = "Helvetica", size = 10),
    legend.key.height = unit(0.35, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, family = "Helvetica")
  )

# -------------------------------
# 11. 保存图像
# -------------------------------
ggsave(output_path, plot = p, width = 12, height = 6, dpi = 300)

cat("✅ Nature 风格渐进色图像已成功保存到：", output_path, "\n")
