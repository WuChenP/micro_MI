# 加载必要的包
library(vegan)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

# 设置结果保存路径
output_dir <- "E:/Python/MI_Analysis/metagenome/Graphic/Virus_RDA_Biplot"

# 如果目录不存在则创建
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 读取病毒丰度数据
abundance_data <- read.csv("E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change_virus/心梗组_病毒_filtered_1percent.csv", 
                           header = TRUE, row.names = 1, check.names = FALSE)

# 转置数据，使样本为行，物种为列
abundance_t <- as.data.frame(t(abundance_data))

# 读取解释变量数据
env_data <- read_excel("E:/Python/MI_Analysis/origin_data/样本协变量数据.xlsx")

# 选择需要的列并重命名（去掉吸烟和饮酒）
env_selected <- env_data %>%
  select(分析名称, 分组, `1男2女`, 年龄, BMI) %>%
  rename(
    SampleID = 分析名称,
    Group = 分组,
    Gender = `1男2女`,
    Age = 年龄,
    BMI = BMI
  )

# 将分类变量转换为因子，修改参考水平
env_selected <- env_selected %>%
  mutate(
    # 修改参考水平：将CON设为参考水平，这样图中会显示MI的箭头
    Group = factor(Group, levels = c("CON", "MI")),
    # 将Female设为参考水平，这样图中会显示Male的箭头
    Gender = factor(Gender, levels = c(2, 1), labels = c("Female", "Male"))
  ) %>%
  as.data.frame()

# 确保样本ID匹配
rownames(env_selected) <- env_selected$SampleID

# 对齐样本
common_samples <- intersect(rownames(abundance_t), rownames(env_selected))
abundance_filtered <- abundance_t[common_samples, ]
env_filtered <- env_selected[common_samples, ]

# 检查数据
print(paste("Number of common samples:", length(common_samples)))
print("Environmental variables summary:")
print(summary(env_filtered))

# 对丰度数据进行Hellinger转换
abundance_hell <- decostand(abundance_filtered, "hellinger")

# 执行RDA分析（去掉Smoking和Drinking）
rda_result <- rda(abundance_hell ~ Group + Age + Gender + BMI, 
                  data = env_filtered)

# 显示RDA结果摘要
print(summary(rda_result))

# 计算RDA轴解释的方差比例
rda_summary <- summary(rda_result)
rda1_percent <- round(100 * rda_summary$cont$importance[2, "RDA1"], 2)
rda2_percent <- round(100 * rda_summary$cont$importance[2, "RDA2"], 2)

# 提取绘图数据
sites <- scores(rda_result, display = "sites")
species <- scores(rda_result, display = "species")
biplot <- scores(rda_result, display = "bp")

# 创建绘图数据框
plot_data <- data.frame(
  RDA1 = sites[, 1],
  RDA2 = sites[, 2],
  Group = env_filtered$Group
)

# 创建箭头数据框
biplot_data <- data.frame(
  Variable = rownames(biplot),
  RDA1 = biplot[, 1] * 3,
  RDA2 = biplot[, 2] * 3
)

# 简化变量名称显示（现在会显示MI而不是CON）
biplot_data$Variable <- recode(biplot_data$Variable,
                               "GroupMI" = "MI",  # 现在会显示MI
                               "Age" = "Age",
                               "GenderMale" = "Male",  # 现在会显示Male
                               "GenderFemale" = "Female",
                               "BMI" = "BMI")

# 创建带置信椭圆的双标图
rda_plot <- ggplot() +
  # 添加置信椭圆
  stat_ellipse(data = plot_data, 
               aes(x = RDA1, y = RDA2, color = Group, fill = Group),
               type = "norm", level = 0.95, alpha = 0.2, geom = "polygon") +
  
  # 添加样本点
  geom_point(data = plot_data, aes(x = RDA1, y = RDA2, color = Group, shape = Group), 
             size = 3, alpha = 0.7) +
  
  # 添加环境变量箭头
  geom_segment(data = biplot_data, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "red", linewidth = 1) +
  
  # 添加环境变量标签
  geom_text(data = biplot_data, 
            aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = Variable),
            color = "red", size = 4, fontface = "bold") +
  
  # 添加坐标轴标签
  labs(
    x = paste0("RDA1 (", rda1_percent, "%)"),
    y = paste0("RDA2 (", rda2_percent, "%)"),
    title = "RDA Biplot - Host Factors vs Viral Community",
    subtitle = paste("Arrow length = explanatory power | Total variance:", 
                     round(100 * rda_result$CCA$tot.chi / rda_result$tot.chi, 2), "%"),
    color = "Group", 
    shape = "Group",
    fill = "Group"  # 添加填充图例
  ) +
  
  # 主题设置
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  
  # 确保图形比例正确
  coord_fixed(ratio = 1) +
  
  # 设置颜色和形状
  scale_color_manual(values = c("MI" = "#E41A1C", "CON" = "#377EB8")) +
  scale_fill_manual(values = c("MI" = "#E41A1C", "CON" = "#377EB8")) +  # 椭圆填充色
  scale_shape_manual(values = c("MI" = 16, "CON" = 17))

# 显示图形
print(rda_plot)

# 保存图形到指定路径
ggsave(file.path(output_dir, "Virus_RDA_Biplot.png"), rda_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "Virus_RDA_Biplot.pdf"), rda_plot, width = 10, height = 8)

# 输出统计摘要到文件
sink(file.path(output_dir, "RDA_Statistical_Summary.txt"))
cat("=== RDA ANALYSIS STATISTICAL SUMMARY ===\n")
cat("Analysis date:", date(), "\n\n")
cat("Number of samples:", length(common_samples), "\n")
cat("Total variance:", round(rda_result$tot.chi, 4), "\n")
cat("Constrained variance:", round(rda_result$CCA$tot.chi, 4), "\n")
cat("Unconstrained variance:", round(rda_result$CA$tot.chi, 4), "\n")
cat("Proportion of variance explained:", round(rda_result$CCA$tot.chi / rda_result$tot.chi * 100, 2), "%\n\n")

cat("=== RDA AXIS VARIANCE EXPLANATION ===\n")
print(rda_summary$cont$importance)

cat("\n=== ANOVA RESULTS (Permutation test) ===\n")
anova_result <- anova(rda_result, by = "term", permutations = 999)
print(anova_result)

cat("\n=== ENVIRONMENTAL VARIABLE COEFFICIENTS ===\n")
print(coef(rda_result))

cat("\n=== VARIABLE LABELS EXPLANATION ===\n")
cat("MI: Myocardial Infarction group\n")
cat("CON: Control group\n")
cat("Male: Male (reference: Female)\n")
cat("Female: Female\n")
cat("Age: Age (continuous)\n")
cat("BMI: Body Mass Index (continuous)\n")
sink()

# 保存RDA结果对象
save(rda_result, file = file.path(output_dir, "RDA_result.RData"))

# 保存绘图数据
write.csv(plot_data, file.path(output_dir, "RDA_site_scores.csv"), row.names = TRUE)
write.csv(biplot_data, file.path(output_dir, "RDA_biplot_scores.csv"), row.names = FALSE)

# 输出完成信息
cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
cat("Results saved to:", output_dir, "\n")
cat("Generated files:\n")
cat("- Virus_RDA_Biplot.png/pdf: RDA biplot visualization\n")
cat("- RDA_Statistical_Summary.txt: Detailed statistical results\n")
cat("- RDA_result.RData: RDA analysis object\n")
cat("- RDA_site_scores.csv: Sample coordinates in RDA space\n")
cat("- RDA_biplot_scores.csv: Environmental variable arrow coordinates\n\n")

cat("RDA Summary:\n")
cat("- RDA1 explains", rda1_percent, "% of constrained variance\n")
cat("- RDA2 explains", rda2_percent, "% of constrained variance\n")
cat("- Total constrained variance:", round(100 * rda_result$CCA$tot.chi / rda_result$tot.chi, 2), "%\n")