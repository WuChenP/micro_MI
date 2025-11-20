# ================================
# R脚本：Wilcoxon 两组比较 + 顶刊风格绘图
# 适用于相对丰度乘以倍数值的数据
# ================================

library(tidyverse)
library(readxl)
library(ggpubr)
library(scales)

# ====================== 1. 读取数据和验证 ========================
file_path <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/virus_family_no_HF.xlsx"
raw <- read_xlsx(file_path)
df <- raw %>% rename(Family = 1)

# 验证数据格式
if(! "Family" %in% colnames(df)) {
  stop("Error: First column must be named 'Family'")
}

mi_cols  <- grep("^MI", colnames(df), value = TRUE)
con_cols <- grep("^CON", colnames(df), value = TRUE)

if(length(mi_cols) == 0 | length(con_cols) == 0) {
  stop("Error: No MI or CON samples found")
}

print(paste("MI samples:", length(mi_cols)))
print(paste("CON samples:", length(con_cols)))
print(paste("Total families:", nrow(df)))

# ====================== 2. 数据预处理 ========================
# 移除全为零的行
valid_rows <- rowSums(df[, c(mi_cols, con_cols)], na.rm = TRUE) > 0
df <- df[valid_rows, ]

if(nrow(df) == 0) {
  stop("Error: All family data are zero")
}

print(paste("Valid families:", nrow(df)))

# ====================== 3. Wilcoxon 检验 ========================
results <- map_df(1:nrow(df), function(i){
  fam <- df$Family[i]
  mi_vals  <- unlist(df[i, mi_cols])
  con_vals <- unlist(df[i, con_cols])
  
  # 数据质量检查
  if(any(is.na(mi_vals)) | any(is.na(con_vals))) {
    return(tibble(
      Family = fam,
      p_value = NA_real_,
      statistic = NA_real_,
      MI_mean = mean(mi_vals, na.rm = TRUE),
      CON_mean = mean(con_vals, na.rm = TRUE),
      Effect_size = mean(mi_vals, na.rm = TRUE) - mean(con_vals, na.rm = TRUE),
      note = "Contains NA values"
    ))
  }
  
  # 新增检查：如果CON组或MI组全为0，则跳过该Family
  if(all(con_vals == 0) | all(mi_vals == 0)) {
    return(tibble(
      Family = fam,
      p_value = NA_real_,
      statistic = NA_real_,
      MI_mean = mean(mi_vals),
      CON_mean = mean(con_vals),
      Effect_size = mean(mi_vals) - mean(con_vals),
      note = ifelse(all(con_vals == 0) & all(mi_vals == 0), "Both groups all zeros",
                    ifelse(all(con_vals == 0), "CON group all zeros", "MI group all zeros"))
    ))
  }
  
  # 检查方差
  if(var(mi_vals) == 0 & var(con_vals) == 0) {
    return(tibble(
      Family = fam,
      p_value = NA_real_,
      statistic = NA_real_,
      MI_mean = mean(mi_vals),
      CON_mean = mean(con_vals),
      Effect_size = 0,
      note = "Zero variance in both groups"
    ))
  }
  
  # 执行Wilcoxon检验
  test_result <- tryCatch({
    wt <- wilcox.test(mi_vals, con_vals, exact = FALSE)
    list(p_value = wt$p.value, statistic = wt$statistic)
  }, error = function(e) {
    list(p_value = NA_real_, statistic = NA_real_)
  })
  
  tibble(
    Family = fam,
    p_value = test_result$p_value,
    statistic = test_result$statistic,
    MI_mean = mean(mi_vals),
    CON_mean = mean(con_vals),
    Effect_size = mean(mi_vals) - mean(con_vals),
    note = NA_character_
  )
})

# ====================== 4. 结果整理 ========================
results <- results %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "Invalid",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**", 
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  )

# 统计结果
valid_tests <- sum(!is.na(results$p_value))
sig_results <- results %>% filter(p_value < 0.05 & !is.na(p_value))

print(paste("Successful tests:", valid_tests, "/", nrow(results)))
print(paste("Significant families (p < 0.05):", nrow(sig_results)))

# 输出被跳过的Family信息
skipped_families <- results %>% 
  filter(!is.na(note) & str_detect(note, "group all zeros"))
if(nrow(skipped_families) > 0) {
  cat("\n=== Skipped Families (one group all zeros) ===\n")
  print(skipped_families %>% select(Family, MI_mean, CON_mean, note))
}

# ====================== 5. 绘制顶刊风格差异Family柱状图 ========================
output_folder <- "E:/Python/MI_Analysis/metagenome/Graphic/Comparison_differences_two_groups_family_level"
if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# 保存完整结果
write.csv(results, file.path(output_folder, "virus_family_wilcoxon_results.csv"), 
          row.names = FALSE, na = "")

if(nrow(sig_results) > 0) {
  write.csv(sig_results, file.path(output_folder, "virus_family_significant_results.csv"), 
            row.names = FALSE, na = "")
  
  # 输出top结果
  cat("\n=== Top Significant Families (p < 0.05) ===\n")
  print(sig_results %>% 
          arrange(p_value) %>% 
          select(Family, Effect_size, p_value, significance) %>%
          head(10))
  
  # ====================== 顶刊风格绘图 ========================
  
  # 选择top显著差异的Family进行绘图（最多12个）
  top_sig_families <- sig_results %>% 
    arrange(p_value) %>% 
    head(12) %>% 
    pull(Family)
  
  # 准备绘图数据
  plot_data <- df %>% 
    filter(Family %in% top_sig_families) %>% 
    select(Family, all_of(c(mi_cols, con_cols))) %>% 
    pivot_longer(
      cols = -Family,
      names_to = "Sample",
      values_to = "Abundance"
    ) %>% 
    mutate(
      Group = ifelse(str_detect(Sample, "^MI"), "MI", "Control")
    )
  
  # 添加显著性标注
  sig_annotations <- sig_results %>% 
    filter(Family %in% top_sig_families) %>% 
    select(Family, p_value, significance)
  
  plot_data <- plot_data %>% 
    left_join(sig_annotations, by = "Family")
  
  # 计算每个Family的最大值用于标注位置
  annotation_data <- plot_data %>%
    group_by(Family) %>%
    summarise(
      max_abundance = max(Abundance),
      .groups = "drop"
    ) %>%
    left_join(sig_annotations, by = "Family")
  
  # 顶刊风格箱线图
  p1 <- ggplot(plot_data, aes(x = Family, y = Abundance, fill = Group)) +
    geom_boxplot(
      width = 0.6,
      alpha = 0.8,
      outlier.shape = 21,
      outlier.size = 1.5,
      outlier.fill = "white",
      outlier.color = "black"
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
      size = 1,
      alpha = 0.6,
      color = "black"
    ) +
    scale_fill_manual(
      values = c("MI" = "#D55E00", "Control" = "#0072B2"),  # Nature风格颜色
      name = "Group"
    ) +
    labs(
      title = "Differential Viral Family Abundance",
      subtitle = "Scaled relative abundance with Wilcoxon test",
      x = "Viral Family",
      y = "Scaled Relative Abundance"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10,
        face = "italic"
      ),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor.y = element_blank()
    ) +
    # 添加显著性标注
    geom_text(
      data = annotation_data,
      aes(x = Family, y = max_abundance * 1.15, label = significance),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    # 添加p值标注
    geom_text(
      data = annotation_data,
      aes(x = Family, y = max_abundance * 1.25, 
          label = paste0("p = ", format.pval(p_value, digits = 2))),
      inherit.aes = FALSE,
      size = 3,
      color = "grey40"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  # 顶刊风格均值柱状图
  mean_data <- plot_data %>% 
    group_by(Family, Group) %>% 
    summarise(
      Mean_Abundance = mean(Abundance),
      SE = sd(Abundance) / sqrt(n()),
      .groups = "drop"
    )
  
  p2 <- ggplot(mean_data, aes(x = Family, y = Mean_Abundance, fill = Group)) +
    geom_col(
      position = position_dodge(0.8),
      width = 0.7,
      color = "black",
      linewidth = 0.2
    ) +
    geom_errorbar(
      aes(ymin = Mean_Abundance - SE, ymax = Mean_Abundance + SE),
      position = position_dodge(0.8),
      width = 0.3,
      linewidth = 0.5
    ) +
    scale_fill_manual(
      values = c("MI" = "#D55E00", "Control" = "#0072B2"),
      name = "Group"
    ) +
    labs(
      title = "Mean Viral Family Abundance",
      subtitle = "Scaled relative abundance with standard error",
      x = "Viral Family",
      y = "Mean Scaled Abundance ± SE"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10,
        face = "italic"
      ),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # 保存图形
  ggsave(
    file.path(output_folder, "significant_families_nature_style.pdf"),
    p1,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(output_folder, "mean_abundance_nature_style.pdf"),
    p2,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  # 保存高分辨率PNG
  ggsave(
    file.path(output_folder, "significant_families_nature_style.png"),
    p1,
    width = 10,
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
  print(paste("Generated plots for", length(top_sig_families), "significant families"))
  
} else {
  cat("\n=== Top Near-Significant Families ===\n")
  top_results <- results %>% 
    filter(!is.na(p_value)) %>% 
    arrange(p_value) %>% 
    head(10)
  print(top_results %>% select(Family, Effect_size, p_value, significance))
}

cat("\n=== Analysis Complete ===\n")
cat("Output directory:", output_folder, "\n")