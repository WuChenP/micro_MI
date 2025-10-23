# =====================================================
# 四类微生物 Family 水平 α/β 多样性分析 (小提琴+箱线，β拼图)
# 包含所有多样性指标：Observed, Richness, Number, Shannon, Simpson, InvSimpson, Chao1, ACE
# 新增：保存所有样本的原始多样性指标数据
# =====================================================

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(patchwork) # 拼图

# ----------------------------
# 文件路径
# ----------------------------
microbe_files <- list(
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/fungi_family.xlsx",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/virus_family_no_HF.xlsx",
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/archaea_family.xlsx",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/family_data/bacteria_family.xlsx"
)

meta_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
output_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/diversity_plots_family/"

# 创建输出目录
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------
# 创建统计结果保存文件
# ----------------------------
stat_results_file <- file.path(output_dir, "family_diversity_statistical_results.txt")

# 初始化统计结果文件
cat("Family Level Diversity Analysis Statistical Results\n", file = stat_results_file)
cat("==================================================\n\n", file = stat_results_file, append = TRUE)
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = stat_results_file, append = TRUE)
cat("Diversity Indices: Observed, Richness, Number, Shannon, Simpson, InvSimpson, Chao1, ACE\n", file = stat_results_file, append = TRUE)
cat("Note: Chao1 and ACE based on converted data - interpret with caution\n\n", file = stat_results_file, append = TRUE)

# ----------------------------
# 创建原始数据保存文件
# ----------------------------
raw_data_file <- file.path(output_dir, "family_diversity_raw_data.xlsx")

# 初始化原始数据列表
raw_data_list <- list()

# ----------------------------
# 读取元数据
# ----------------------------
metadata <- read.xlsx(meta_file, rowNames = TRUE)
metadata$ID <- rownames(metadata)

# ----------------------------
# 函数：Excel 转 phyloseq
# ----------------------------
make_phyloseq_from_excel <- function(file, metadata){
  abundance <- read.xlsx(file, rowNames = TRUE)
  common_samples <- intersect(colnames(abundance), metadata$ID)
  abundance <- abundance[, common_samples, drop = FALSE]
  meta <- metadata[metadata$ID %in% common_samples, , drop = FALSE]
  rownames(meta) <- meta$ID
  otu <- otu_table(as.matrix(abundance), taxa_are_rows = TRUE)
  sam <- sample_data(meta)
  phyloseq(otu, sam)
}

# ----------------------------
# 改进的统计检验函数
# ----------------------------
perform_statistical_test <- function(data, metric, group_var = "Group") {
  groups <- unique(data[[group_var]])
  
  if(length(groups) == 2) {
    # 两组比较：Wilcoxon检验
    group1_data <- data[[metric]][data[[group_var]] == groups[1]]
    group2_data <- data[[metric]][data[[group_var]] == groups[2]]
    
    # 检查数据是否有效
    if(length(group1_data) < 3 || length(group2_data) < 3) {
      return(list(p.value = NA, method = "Insufficient data", note = "Need at least 3 samples per group"))
    }
    
    # 检查是否有方差
    if(all(group1_data == group1_data[1]) && all(group2_data == group2_data[1])) {
      return(list(p.value = NA, method = "No variance", note = "All values are identical within groups"))
    }
    
    # 执行Wilcoxon检验，抑制tied values警告
    test_result <- suppressWarnings({
      wilcox.test(as.formula(paste(metric, "~", group_var)), data = data)
    })
    
    note <- if(length(warnings()) > 0) "Tied values present" else ""
    
    return(list(
      p.value = test_result$p.value,
      method = "Wilcoxon rank sum test",
      note = note
    ))
    
  } else if(length(groups) > 2) {
    # 多组比较：Kruskal-Wallis检验
    if(nrow(data) >= length(groups) * 3) {
      test <- kruskal.test(as.formula(paste(metric, "~", group_var)), data = data)
      return(list(
        p.value = test$p.value,
        method = "Kruskal-Wallis test",
        note = ""
      ))
    } else {
      return(list(p.value = NA, method = "Insufficient data", note = "Need at least 3 samples per group"))
    }
  } else {
    return(list(p.value = NA, method = "Only one group", note = "Cannot perform group comparison with one group"))
  }
}

# ----------------------------
# 完整的α多样性计算函数（包含Number和ACE）
# ----------------------------
calculate_complete_alpha_diversity <- function(ps_obj) {
  otu_matrix <- as.matrix(otu_table(ps_obj))
  
  # 计算基础多样性指标
  alpha_df <- data.frame(
    Observed = apply(otu_matrix > 0, 2, sum),
    Richness = apply(otu_matrix > 0, 2, sum),
    Number = apply(otu_matrix > 0, 2, sum),  # 添加Number指标
    Shannon = diversity(t(otu_matrix), index = "shannon"),
    Simpson = diversity(t(otu_matrix), index = "simpson"),
    InvSimpson = 1 / diversity(t(otu_matrix), index = "simpson")
  )
  
  # 尝试计算Chao1和ACE（基于相对丰度数据，仅供参考）
  tryCatch({
    # 转换为模拟计数数据
    otu_int <- round(otu_matrix * 1000000)
    
    # 计算Chao1
    chao1_results <- estimateR(t(otu_int))[2, ]
    if(all(!is.na(chao1_results)) && all(chao1_results > 0)) {
      alpha_df$Chao1 <- chao1_results
    }
    
    # 计算ACE
    ace_results <- estimateR(t(otu_int))[4, ]
    if(all(!is.na(ace_results)) && all(ace_results > 0)) {
      alpha_df$ACE <- ace_results
    }
  }, error = function(e) {
    cat("Note: Chao1/ACE calculation may not be reliable with relative abundance data\n")
  })
  
  alpha_df$SampleID <- rownames(alpha_df)
  sample_data_df <- data.frame(sample_data(ps_obj))
  alpha_df <- merge(alpha_df, sample_data_df[, "Group", drop=FALSE],
                    by.x="SampleID", by.y="row.names")
  
  return(alpha_df)
}

# ----------------------------
# 函数：保存原始多样性数据
# ----------------------------
save_raw_diversity_data <- function(microbe, alpha_df) {
  # 重新排列列的顺序，使SampleID和Group在最后
  if("Chao1" %in% colnames(alpha_df) && "ACE" %in% colnames(alpha_df)) {
    raw_df <- alpha_df[, c("Observed", "Richness", "Number", "Shannon", "Simpson", 
                           "InvSimpson", "Chao1", "ACE", "Group", "SampleID")]
  } else if("Chao1" %in% colnames(alpha_df)) {
    raw_df <- alpha_df[, c("Observed", "Richness", "Number", "Shannon", "Simpson", 
                           "InvSimpson", "Chao1", "Group", "SampleID")]
  } else if("ACE" %in% colnames(alpha_df)) {
    raw_df <- alpha_df[, c("Observed", "Richness", "Number", "Shannon", "Simpson", 
                           "InvSimpson", "ACE", "Group", "SampleID")]
  } else {
    raw_df <- alpha_df[, c("Observed", "Richness", "Number", "Shannon", "Simpson", 
                           "InvSimpson", "Group", "SampleID")]
  }
  
  # 按SampleID排序
  raw_df <- raw_df[order(raw_df$SampleID), ]
  
  return(raw_df)
}

# ----------------------------
# 循环分析四类微生物
# ----------------------------
for(microbe in names(microbe_files)){
  message("开始分析: ", microbe)
  
  # 在统计结果文件中添加分类标题
  cat("\n", paste0("=== ", toupper(microbe), " (Family Level) ==="), "\n", file = stat_results_file, append = TRUE)
  cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = stat_results_file, append = TRUE)
  
  ps_obj <- make_phyloseq_from_excel(microbe_files[[microbe]], metadata)
  sample_data_df <- data.frame(sample_data(ps_obj))
  
  # ----------------------------
  # α 多样性 - 完整指标（包含Number和ACE）
  # ----------------------------
  alpha_df <- calculate_complete_alpha_diversity(ps_obj)
  
  # 动态确定可用的指标
  available_metrics <- c("Observed", "Richness", "Number", "Shannon", "Simpson", "InvSimpson")
  if("Chao1" %in% colnames(alpha_df) && !all(is.na(alpha_df$Chao1))) {
    available_metrics <- c(available_metrics, "Chao1")
    message("✓ Chao1 available for ", microbe)
  } else {
    message("✗ Chao1 not available for ", microbe)
  }
  
  if("ACE" %in% colnames(alpha_df) && !all(is.na(alpha_df$ACE))) {
    available_metrics <- c(available_metrics, "ACE")
    message("✓ ACE available for ", microbe)
  } else {
    message("✗ ACE not available for ", microbe)
  }
  
  message("Available alpha diversity metrics: ", paste(available_metrics, collapse = ", "))
  
  # 保存原始数据
  raw_data_list[[microbe]] <- save_raw_diversity_data(microbe, alpha_df)
  
  # 在统计结果文件中添加α多样性标题
  cat("ALPHA DIVERSITY RESULTS:\n", file = stat_results_file, append = TRUE)
  cat("------------------------\n", file = stat_results_file, append = TRUE)
  
  # 为每个α多样性指标绘图
  for(metric in available_metrics){
    
    # 检查该指标是否有有效数据
    if(!metric %in% colnames(alpha_df) || all(is.na(alpha_df[[metric]]))) {
      message("Warning: No valid data for ", metric, " in ", microbe)
      next
    }
    
    # 移除NA值
    alpha_clean <- alpha_df[!is.na(alpha_df[[metric]]), ]
    
    if(nrow(alpha_clean) == 0) {
      message("Warning: All values are NA for ", metric, " in ", microbe)
      next
    }
    
    # 统计信息
    stat_df <- alpha_clean %>%
      group_by(Group) %>%
      summarise(
        mean = mean(!!sym(metric), na.rm = TRUE),
        median = median(!!sym(metric), na.rm = TRUE),
        sd = sd(!!sym(metric), na.rm = TRUE)
      )
    
    # 组间显著性检验
    test_result <- perform_statistical_test(alpha_clean, metric)
    
    # 格式化p值文本
    if(!is.na(test_result$p.value)) {
      if(test_result$p.value < 0.001) {
        pval_text <- "p < 0.001"
        pval_numeric <- "< 0.001"
      } else if(test_result$p.value < 0.01) {
        pval_text <- "p < 0.01"
        pval_numeric <- sprintf("%.3f", test_result$p.value)
      } else if(test_result$p.value < 0.05) {
        pval_text <- paste0("p = ", sprintf("%.3f", test_result$p.value))
        pval_numeric <- sprintf("%.3f", test_result$p.value)
      } else {
        pval_text <- paste0("p = ", sprintf("%.3f", test_result$p.value))
        pval_numeric <- sprintf("%.3f", test_result$p.value)
      }
      pval_text <- paste0(pval_text, " (", test_result$method, ")")
      if(test_result$note != "") {
        pval_text <- paste0(pval_text, " [", test_result$note, "]")
      }
    } else {
      pval_text <- paste0("No test (", test_result$method, ")")
      pval_numeric <- "NA"
    }
    
    # 添加数据可靠性说明
    reliability_note <- ""
    if(metric == "Chao1") {
      reliability_note <- " [Based on converted data - interpret with caution]"
    } else if(metric == "ACE") {
      reliability_note <- " [Based on converted data - EXTREME caution advised]"
    } else {
      reliability_note <- " [Based on original data - high reliability]"
    }
    
    # 副标题文本
    subtitle_text <- paste0(
      paste(stat_df$Group, ": mean=", round(stat_df$mean, 2),
            "±", round(stat_df$sd, 2), collapse="; "),
      " | ", pval_text, reliability_note
    )
    
    # 将α多样性结果保存到文件
    cat(paste0(metric, ": ", pval_text, reliability_note, "\n"), file = stat_results_file, append = TRUE)
    if(!is.na(test_result$p.value)) {
      cat(paste0("  P-value: ", pval_numeric, "\n"), file = stat_results_file, append = TRUE)
    }
    cat(paste0("  Method: ", test_result$method, "\n"), file = stat_results_file, append = TRUE)
    if(test_result$note != "") {
      cat(paste0("  Note: ", test_result$note, "\n"), file = stat_results_file, append = TRUE)
    }
    cat(paste0("  Samples: Control (n=", sum(alpha_clean$Group == "Control"), "), MI (n=", sum(alpha_clean$Group == "MI"), ")\n"), file = stat_results_file, append = TRUE)
    cat("\n", file = stat_results_file, append = TRUE)
    
    # 绘图 - 根据可靠性使用不同颜色
    fill_colors <- if(metric == "ACE") {
      c("Control" = "#FF0000", "MI" = "#8B0000")  # 红色系 - 极度警告
    } else if(metric == "Chao1") {
      c("Control" = "#FF7F00", "MI" = "#B15928")  # 橙色系 - 警告色
    } else {
      c("Control" = "#1F78B4", "MI" = "#E31A1C")  # 蓝红色系 - 正常色
    }
    
    p_alpha <- ggplot(alpha_clean, aes(x=Group, y=!!sym(metric), fill=Group)) +
      geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
      geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
      geom_jitter(width = 0.1, size = 1, alpha = 0.6) +
      scale_fill_manual(values = fill_colors) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 8),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 10)
      ) +
      labs(
        title = paste(microbe, "Alpha Diversity -", metric),
        subtitle = subtitle_text,
        x = "Group",
        y = metric
      )
    
    # 保存 α-diversity 图
    ggsave(filename = file.path(output_dir, paste0(microbe, "_alpha_", metric, ".pdf")),
           p_alpha, width=8, height=6)
  }
  
  # ----------------------------
  # β 多样性 (PCoA Bray-Curtis)
  # ----------------------------
  dist <- phyloseq::distance(ps_obj, method = "bray")
  ord <- ordinate(ps_obj, method="PCoA", distance=dist)
  ord_df <- plot_ordination(ps_obj, ord, type="samples")$data
  if(!"Group" %in% colnames(ord_df)){
    ord_df$Group <- sample_data_df[rownames(ord_df), "Group"]
  }
  
  # PERMANOVA - 获取p值和R²值
  permanova_text <- "No PERMANOVA test"
  pval <- NA
  r_squared <- NA
  
  if("Group" %in% colnames(sample_data_df) && length(unique(ord_df$Group)) > 1){
    adonis_res <- adonis2(dist ~ Group, data = sample_data_df)
    pval <- adonis_res$`Pr(>F)`[1]
    r_squared <- adonis_res$R2[1]
    
    if(!is.na(pval) && !is.na(r_squared)) {
      # 格式化p值
      if(pval < 0.001) {
        pval_formatted <- "p < 0.001"
      } else if(pval < 0.01) {
        pval_formatted <- "p < 0.01"
      } else if(pval < 0.05) {
        pval_formatted <- paste0("p = ", sprintf("%.3f", pval))
      } else {
        pval_formatted <- paste0("p = ", sprintf("%.3f", pval))
      }
      
      # 格式化R²值
      r2_formatted <- sprintf("R² = %.3f", r_squared)
      
      # 组合文本
      permanova_text <- paste0("PERMANOVA: ", pval_formatted, ", ", r2_formatted)
    } else {
      permanova_text <- "PERMANOVA test failed"
    }
  }
  
  # 将β多样性结果保存到文件
  cat("BETA DIVERSITY RESULTS:\n", file = stat_results_file, append = TRUE)
  cat("-----------------------\n", file = stat_results_file, append = TRUE)
  cat(paste0("Method: PERMANOVA (Bray-Curtis distance)\n"), file = stat_results_file, append = TRUE)
  if(!is.na(pval) && !is.na(r_squared)) {
    cat(paste0("P-value: ", sprintf("%.6f", pval), "\n"), file = stat_results_file, append = TRUE)
    cat(paste0("R-squared: ", sprintf("%.6f", r_squared), "\n"), file = stat_results_file, append = TRUE)
    cat(paste0("Result: ", permanova_text, "\n"), file = stat_results_file, append = TRUE)
  } else {
    cat(paste0("Result: ", permanova_text, "\n"), file = stat_results_file, append = TRUE)
  }
  cat(paste0("Samples: Control (n=", sum(ord_df$Group == "Control"), "), MI (n=", sum(ord_df$Group == "MI"), ")\n"), file = stat_results_file, append = TRUE)
  cat("\n", file = stat_results_file, append = TRUE)
  
  # 样本数
  group_counts <- table(ord_df$Group)
  n_text <- paste(names(group_counts), "n =", group_counts, collapse = "; ")
  
  # PCoA图
  p_beta <- ggplot(ord_df, aes(x=Axis.1, y=Axis.2, color=Group, fill=Group)) +
    geom_point(size=3, alpha=0.8) +
    stat_ellipse(level=0.95, linetype=2, alpha=0.5) +
    scale_color_manual(values = c("Control" = "#1F78B4", "MI" = "#E31A1C")) +
    scale_fill_manual(values = c("Control" = "#1F78B4", "MI" = "#E31A1C")) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = paste(microbe, "Beta Diversity (PCoA-Bray)"),
      subtitle = paste(n_text, "|", permanova_text),
      x = paste0("PCoA1 (", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%)"),
      y = paste0("PCoA2 (", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%)")
    )
  
  # Bray-Curtis 组间距离箱线图
  dist_mat <- as.matrix(dist)
  n <- nrow(dist_mat)
  group_vec <- sample_data_df$Group
  beta_box <- data.frame()
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      beta_box <- rbind(beta_box, data.frame(
        Dist = dist_mat[i,j],
        Pair = paste(group_vec[i], group_vec[j], sep="-")
      ))
    }
  }
  
  p_beta_box <- ggplot(beta_box, aes(x=Pair, y=Dist, fill=Pair)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Control-Control" = "#1F78B4",
                                 "Control-MI" = "#B2DF8A",
                                 "MI-MI" = "#E31A1C")) +
    theme_classic() +
    labs(title=paste(microbe, "Bray-Curtis distances"), x="Group Pairs", y="Distance") +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.position = "none")
  
  # 拼接 β-diversity 图
  p_beta_combined <- p_beta / p_beta_box + plot_layout(heights = c(2,1))
  
  # 保存 β-diversity 拼图
  ggsave(filename = file.path(output_dir, paste0(microbe, "_beta_combined.pdf")),
         p_beta_combined, width=10, height=10)
  
  # 单独保存PCoA图
  ggsave(filename = file.path(output_dir, paste0(microbe, "_beta_pcoa.pdf")),
         p_beta, width=8, height=6)
  
  message("完成并保存: ", microbe)
}

# ----------------------------
# 保存所有原始数据到Excel文件
# ----------------------------
if(length(raw_data_list) > 0) {
  # 创建Excel工作簿
  wb <- createWorkbook()
  
  # 为每个微生物类型添加单独的工作表
  for(microbe in names(raw_data_list)) {
    addWorksheet(wb, microbe)
    writeData(wb, microbe, raw_data_list[[microbe]])
    
    # 设置列宽以便更好地查看数据
    setColWidths(wb, microbe, cols = 1:ncol(raw_data_list[[microbe]]), widths = "auto")
  }
  
  # 添加汇总表（所有数据合并）
  all_raw_data <- do.call(rbind, lapply(names(raw_data_list), function(microbe) {
    data <- raw_data_list[[microbe]]
    data$Microbe_Type <- microbe
    # 重新排列列，使Microbe_Type在第一列
    data <- data[, c("Microbe_Type", setdiff(colnames(data), "Microbe_Type"))]
    return(data)
  }))
  
  addWorksheet(wb, "All_Data_Combined")
  writeData(wb, "All_Data_Combined", all_raw_data)
  setColWidths(wb, "All_Data_Combined", cols = 1:ncol(all_raw_data), widths = "auto")
  
  # 保存Excel文件
  saveWorkbook(wb, raw_data_file, overwrite = TRUE)
  message("原始多样性数据已保存到: ", raw_data_file)
  
  # 显示数据预览
  message("\n数据预览:")
  for(microbe in names(raw_data_list)) {
    message("=== ", microbe, " ===")
    print(head(raw_data_list[[microbe]]))
    message("样本数量: ", nrow(raw_data_list[[microbe]]))
    message("Control组: ", sum(raw_data_list[[microbe]]$Group == "Control"))
    message("MI组: ", sum(raw_data_list[[microbe]]$Group == "MI"))
    message("")
  }
}

cat("\n===== Family Level Analysis Complete! =====\n")
cat("Statistical results saved to: 'family_diversity_statistical_results.txt'\n")
cat("Raw diversity data saved to: 'family_diversity_raw_data.xlsx'\n")
cat("Output directory: ", output_dir, "\n")