# ============================================
# α & β 多样性分析（四类微生物，相对丰度表）
# 专门修复ACE计算问题的版本 + 统计结果保存功能
# ============================================

# 设置编码和语言环境
Sys.setlocale(category = "LC_ALL", locale = "Chinese")
options(encoding = "UTF-8")

# 依赖包
suppressPackageStartupMessages({
  if(!requireNamespace("phyloseq", quietly=TRUE)) install.packages("phyloseq")
  if(!requireNamespace("vegan", quietly=TRUE)) install.packages("vegan")
  if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
  if(!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
  if(!requireNamespace("ape", quietly=TRUE)) install.packages("ape")
  if(!requireNamespace("rlang", quietly=TRUE)) install.packages("rlang")
  if(!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
  if(!requireNamespace("stringr", quietly=TRUE)) install.packages("stringr")
  
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(readxl)
  library(ape)
  library(rlang)
  library(dplyr)
  library(stringr)
})

# -------------------------
# 1. 读取样本元数据
# -------------------------
meta_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/sample_metadata.xlsx"
meta <- read_excel(meta_file)
meta <- as.data.frame(meta)
rownames(meta) <- meta$SampleID

# -------------------------
# 2. 丰度表路径列表
# -------------------------
microbe_files <- list(
  archaea  = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/古菌_filtered_20percent.csv",
  bacteria = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/细菌_filtered_20percent.csv",
  fungi    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/真菌_filtered_20percent.csv",
  virus    = "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/病毒_filtered_20percent.csv"
)

# -------------------------
# 3. 输出目录
# -------------------------
out_dir <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_20percent/diversity_plots_OTU"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# -------------------------
# 4. 创建统计结果保存文件
# -------------------------
stat_results_file <- file.path(out_dir, "diversity_statistical_results.txt")

# 初始化统计结果文件
cat("Diversity Analysis Statistical Results\n", file = stat_results_file)
cat("======================================\n\n", file = stat_results_file, append = TRUE)
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = stat_results_file, append = TRUE)

# -------------------------
# 5. 专门修复ACE计算的α多样性函数
# -------------------------
calculate_alpha_diversity_fixed <- function(otu_matrix) {
  
  cat("Note: Data conversion for diversity indices calculation\n")
  cat("Relative abundance × 1,000,000 → round → simulated count data\n")
  cat("Chao1 and ACE indices based on converted data, results for reference only\n\n")
  
  # 数据转换：相对丰度转换为模拟计数数据
  conversion_factor <- 1000000
  otu_int <- round(otu_matrix * conversion_factor)
  
  # 检查转换后的数据质量
  singletons <- sum(otu_int == 1)
  doubletons <- sum(otu_int == 2)
  cat(sprintf("Converted data stats: singletons = %d, doubletons = %d\n", singletons, doubletons))
  
  # 初始化结果数据框
  alpha_df <- data.frame(
    # 基于原始相对丰度的可靠指标
    Observed = apply(otu_matrix > 0, 2, sum),
    Richness = apply(otu_matrix > 0, 2, sum),
    Number = apply(otu_matrix > 0, 2, sum),
    Shannon = diversity(t(otu_matrix), index = "shannon"),
    Simpson = diversity(t(otu_matrix), index = "simpson"),
    InvSimpson = 1 / diversity(t(otu_matrix), index = "simpson")
  )
  
  # 计算Chao1
  tryCatch({
    chao1_results <- estimateR(t(otu_int))[2, ]
    if(all(!is.na(chao1_results)) && all(chao1_results > 0)) {
      alpha_df$Chao1 <- chao1_results
      cat("Chao1 calculation successful\n")
    } else {
      alpha_df$Chao1 <- NA
      cat("Warning: Chao1 calculation produced invalid results\n")
    }
  }, error = function(e) {
    cat("Error in Chao1 calculation:", e$message, "\n")
    alpha_df$Chao1 <<- NA
  })
  
  # 专门修复ACE计算
  cat("Attempting ACE calculation...\n")
  tryCatch({
    # 方法1: 直接使用estimateR
    ace_results <- estimateR(t(otu_int))[4, ]
    
    # 检查ACE结果是否有效
    if(all(!is.na(ace_results)) && all(ace_results > 0) && all(is.finite(ace_results))) {
      alpha_df$ACE <- ace_results
      cat("ACE calculation successful using estimateR\n")
    } else {
      # 方法2: 尝试手动计算ACE
      cat("Trying alternative ACE calculation method...\n")
      ace_manual <- apply(otu_int, 2, function(x) {
        # 手动实现ACE计算逻辑
        x <- x[x > 0]  # 移除零值
        if(length(x) == 0) return(NA)
        
        S_obs <- length(x)
        if(S_obs == 0) return(NA)
        
        # 计算稀有物种（abundance <= 10）
        rare_species <- x[x <= 10]
        S_rare <- length(rare_species)
        
        if(S_rare == 0) {
          # 如果没有稀有物种，返回观察到的物种数
          return(S_obs)
        }
        
        # 计算丰富物种（abundance > 10）
        abundant_species <- x[x > 10]
        S_abund <- length(abundant_species)
        
        # 计算ACE
        n_rare <- sum(rare_species)
        if(n_rare == 0) return(S_obs)
        
        # 计算变异系数
        f1 <- sum(rare_species == 1)  # singletons
        if(f1 == n_rare) return(S_obs)  # 避免除以零
        
        gamma_2 <- max(S_rare * sum(rare_species * (rare_species - 1)) / (n_rare * (n_rare - 1)) - 1, 0)
        ace <- S_abund + S_rare / (1 - f1 / n_rare) + (f1 * gamma_2) / (1 - f1 / n_rare)
        
        return(ace)
      })
      
      if(all(!is.na(ace_manual)) && all(ace_manual > 0)) {
        alpha_df$ACE <- ace_manual
        cat("ACE calculation successful using manual method\n")
      } else {
        alpha_df$ACE <- NA
        cat("Warning: ACE calculation failed with all methods\n")
      }
    }
  }, error = function(e) {
    cat("Error in ACE calculation:", e$message, "\n")
    alpha_df$ACE <<- NA
  })
  
  # 添加数据质量标记
  alpha_df$Data_Type <- "Relative_abundance"
  alpha_df$Conversion_Factor <- conversion_factor
  alpha_df$Singletons_Present <- singletons > 0
  
  # 根据ACE是否成功计算来标记可靠性
  if("ACE" %in% colnames(alpha_df) && !all(is.na(alpha_df$ACE))) {
    alpha_df$Reliability_Note <- "Chao1_ACE_based_on_converted_data_interpret_with_caution"
  } else {
    alpha_df$Reliability_Note <- "Chao1_only_ACE_calculation_failed"
  }
  
  return(alpha_df)
}

# -------------------------
# 6. 改进的统计检验函数
# -------------------------
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

# -------------------------
# 7. 批量分析
# -------------------------

# 创建分析说明文档
writeLines(
  c("Diversity Analysis Important Notes",
    "=================================",
    "Data Type: Relative abundance data",
    "Conversion Method: Relative abundance × 1,000,000 → round",
    "",
    "ACE Calculation Issue:",
    "- ACE index requires specific data conditions to calculate properly",
    "- May fail if there are no singletons or specific abundance patterns",
    "- Manual calculation method implemented as fallback",
    "",
    "Indicator Reliability Classification:",
    "[High Reliability]: Observed, Richness, Number, Shannon, Simpson, InvSimpson",
    "    - Based on original relative abundance data",
    "    - Results can be directly interpreted",
    "",
    "[Interpret with Caution]: Chao1",
    "    - Based on converted simulated count data",
    "    - Affected by conversion factor",
    "",
    "[If Available - Extra Caution]: ACE",
    "    - Most sensitive to data conversion artifacts",
    "    - Requires careful validation",
    "",
    "Statistical Notes:",
    "- Wilcoxon test used for two-group comparisons",
    "- Kruskal-Wallis test used for multiple groups",
    "- Tied values handled with normal approximation",
    "",
    "Recommendation: Focus on high reliability indicators for main results",
    "================================="),
  file.path(out_dir, "ANALYSIS_WARNING_README.txt")
)

for(name in names(microbe_files)){
  cat("\n===== Processing:", name, "=====\n")
  
  # 在统计结果文件中添加分类标题
  cat("\n", paste0("=== ", toupper(name), " ==="), "\n", file = stat_results_file, append = TRUE)
  cat("Microbial type:", name, "\n", file = stat_results_file, append = TRUE)
  cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = stat_results_file, append = TRUE)
  
  # 读取丰度表
  otu <- read.csv(microbe_files[[name]], header=TRUE, row.names=1, check.names = FALSE)
  
  # 样本交集及对齐
  shared_samples <- intersect(colnames(otu), meta$SampleID)
  if(length(shared_samples)==0) {
    cat("Warning: No matching samples for", name, ", skipping...\n")
    cat("No matching samples for", name, "\n\n", file = stat_results_file, append = TRUE)
    next
  }
  otu <- otu[, shared_samples, drop=FALSE]
  meta_sub <- meta[shared_samples, , drop=FALSE]
  
  # 检查数据是否为数值型
  if(!all(sapply(otu, is.numeric))) {
    cat("Converting data to numeric...\n")
    otu <- as.data.frame(apply(otu, 2, as.numeric))
    rownames(otu) <- rownames(otu)
  }
  
  # -------------------------
  # α 多样性 - 所有指数
  # -------------------------
  alpha_df <- calculate_alpha_diversity_fixed(as.matrix(otu))
  alpha_df$Group <- meta_sub$Group
  alpha_df$SampleID <- rownames(alpha_df)
  
  # 保存完整结果
  write.csv(alpha_df, file.path(out_dir, paste0("alpha_diversity_COMPLETE_", name, ".csv")), row.names=TRUE)
  
  # 动态确定可用的指标
  available_metrics <- c("Shannon", "Simpson", "InvSimpson", "Observed", "Richness", "Number")
  
  # 检查Chao1是否可用
  if("Chao1" %in% colnames(alpha_df) && !all(is.na(alpha_df$Chao1))) {
    available_metrics <- c(available_metrics, "Chao1")
    cat("✓ Chao1 available for", name, "\n")
  } else {
    cat("✗ Chao1 not available for", name, "\n")
  }
  
  # 检查ACE是否可用
  if("ACE" %in% colnames(alpha_df) && !all(is.na(alpha_df$ACE))) {
    available_metrics <- c(available_metrics, "ACE")
    cat("✓ ACE available for", name, "\n")
  } else {
    cat("✗ ACE not available for", name, "\n")
  }
  
  cat("Final available metrics for", name, ":", paste(available_metrics, collapse = ", "), "\n")
  
  # 在统计结果文件中添加α多样性标题
  cat("ALPHA DIVERSITY RESULTS:\n", file = stat_results_file, append = TRUE)
  cat("------------------------\n", file = stat_results_file, append = TRUE)
  
  for(metric in available_metrics){
    
    # 检查该指标是否有有效数据
    if(!metric %in% colnames(alpha_df) || all(is.na(alpha_df[[metric]]))) {
      cat("Warning: No valid data for", metric, "in", name, "\n")
      next
    }
    
    # 移除NA值
    alpha_clean <- alpha_df[!is.na(alpha_df[[metric]]), ]
    
    if(nrow(alpha_clean) == 0) {
      cat("Warning: All values are NA for", metric, "in", name, "\n")
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
    
    # 格式化p值文本 - 在α分析图中确认显示p值
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
    
    # 副标题文本 - 确保显示p值
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
    
    # 使用英文标题避免编码问题
    title_text <- paste("Alpha Diversity -", metric, "(", name, ")")
    
    p <- ggplot(alpha_clean, aes(x = Group, y = !!sym(metric), fill = Group)) +
      geom_violin(trim = FALSE, alpha = 1, width = 0.8) +
      geom_boxplot(width = 0.15, outlier.size = 0.8, fill = "white", alpha = 0.7) +
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
        title = title_text,
        subtitle = subtitle_text,
        x = "Group",
        y = metric
      )
    
    # 保存图片
    ggsave(
      file.path(out_dir, paste0("alpha_", metric, "_", name, ".pdf")),
      p, width = 8, height = 6
    )
    
    cat("Alpha test (", metric, ") for", name, ":", pval_text, reliability_note, "\n")
  }
  
  # -------------------------
  # β PCoA (Bray-Curtis) - 修改部分：显示p值和R²值
  # -------------------------
  # 使用相对丰度数据进行β多样性分析
  otu_rel <- apply(otu, 2, function(x) x/sum(x))
  
  # 检查是否有足够样本进行β多样性分析
  if(ncol(otu_rel) < 3) {
    cat("Warning: Not enough samples for beta diversity analysis in", name, "\n")
    cat("BETA DIVERSITY: Not enough samples for analysis\n\n", file = stat_results_file, append = TRUE)
    next
  }
  
  dist_bc <- vegdist(t(otu_rel), method = "bray")
  
  pcoa_res <- pcoa(dist_bc)
  pco_df <- data.frame(
    sample = rownames(pcoa_res$vectors),
    PC1 = pcoa_res$vectors[, 1],
    PC2 = pcoa_res$vectors[, 2]
  )
  pco_df$Group <- meta_sub$Group
  
  # PERMANOVA - 获取p值和R²值
  permanova_text <- "No PERMANOVA test"
  pval <- NA
  r_squared <- NA
  
  if("Group" %in% colnames(meta_sub) && length(unique(pco_df$Group)) > 1){
    adonis_res <- adonis2(dist_bc ~ Group, data = meta_sub)
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
      
      # 组合文本 - 在β分析图中同时显示p值和R²值
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
  cat(paste0("Samples: Control (n=", sum(pco_df$Group == "Control"), "), MI (n=", sum(pco_df$Group == "MI"), ")\n"), file = stat_results_file, append = TRUE)
  cat("\n", file = stat_results_file, append = TRUE)
  
  # 样本数
  group_counts <- table(pco_df$Group)
  n_text <- paste(names(group_counts), "n =", group_counts, collapse = "; ")
  
  # 绘图（使用英文标题）
  title_text <- paste("PCoA (Bray-Curtis) -", name)
  
  p_pcoa <- ggplot(pco_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 2, alpha = 0.5) +
    scale_color_manual(values = c("Control" = "#1F78B4", "MI" = "#E31A1C")) +
    scale_fill_manual(values = c("Control" = "#1F78B4", "MI" = "#E31A1C")) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = title_text,
      subtitle = paste(n_text, "|", permanova_text),  # 这里同时显示样本数和统计结果
      x = paste0("PC1 (", round(pcoa_res$values$Relative_eig[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(pcoa_res$values$Relative_eig[2] * 100, 1), "%)"),
      color = "Group", fill = "Group"
    )
  
  ggsave(
    file.path(out_dir, paste0("beta_pcoa_bray_", name, ".pdf")),
    p_pcoa, width = 7, height = 5
  )
  
  cat("PERMANOVA for", name, ":", permanova_text, "\n")
}

cat("\n===== Analysis Complete! =====\n")
cat("Important: Check 'ANALYSIS_WARNING_README.txt' in output directory\n")
cat("Statistical results saved to: 'diversity_statistical_results.txt'\n")
cat("Note: ACE calculation is challenging with converted relative abundance data\n")
cat("Recommendation: Focus on Shannon, Simpson, and Observed species for reliable results\n")