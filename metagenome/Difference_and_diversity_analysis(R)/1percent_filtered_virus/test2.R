# ======================================================
# 筛选 ANCOM-BC2 病毒结果
# - 根据富集方向和解释信息关键词分组
# - q值 < 0.05
# - |log2FC| > 1
# ======================================================

library(dplyr)
library(stringr)
library(openxlsx)

# ----------------------------
# 文件路径
# ----------------------------
result_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"

# ----------------------------
# 读取 virus_Result 工作表
# ----------------------------
df <- read.xlsx(result_file, sheet = "virus_Result")

# ----------------------------
# 关键词列表
# ----------------------------
keywords_control <- c(
  "O__Unclassifed","f__Uniclassified","g__Lambdavirus","g__Hubeivirus","o__Petitvirales",
  "d__Monodnaviria","p__Phixviicota","k__Sangervirae","c__Malgrandaviricetes",
  "f__Allomimiviridae","O__Algavirales","g__Bembunaguatrovius","g__Aguilavirus",
  "f__Vilmaviridae","f__Peduoviridae","g__Alegriavirus","S__Buchavirus_coli",
  "g__Lidleunavirus","g__Muvirus","g__Brunovirus","g__Delepguintavirus","g__Lederbergvirus",
  "g__Svunavirus","f__Filixviridae","g__Glaedevirus"
)

keywords_MI <- c(
  "g__Mushuvirus","g__Spinunavirus","S__Serangoonvirus_essarone",
  "g__Serangoonvirus","g__Punavirus","S__Salacisavirus_pssm2","g__Salacisavirus"
)

# ----------------------------
# 筛选条件
# ----------------------------
filtered_df <- df %>%
  filter(
    q_GroupMI < 0.05,
    abs(lfc_GroupMI) > 1,
  ) %>%
  filter(
    (lfc_GroupMI < 0 & str_detect(解释信息, str_c(keywords_control, collapse = "|"))) |
      (lfc_GroupMI > 0 & str_detect(解释信息, str_c(keywords_MI, collapse = "|")))
  )

cat("✅ 共筛选出", nrow(filtered_df), "行符合条件\n")

# ----------------------------
# 输出结果
# ----------------------------
out_file <- sub("\\.xlsx$", "_筛选结果_方向关键词.xlsx", result_file)
write.xlsx(filtered_df, out_file, overwrite = TRUE)

cat("🎉 筛选结果已保存到：", out_file, "\n")







library(dplyr)
library(openxlsx)

# ----------------------------
# 文件路径
# ----------------------------
result_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"

# 读取 virus_Result
df <- read.xlsx(result_file, sheet = "virus_Result")

# ----------------------------
# 只保留显著行
# ----------------------------
df_sig <- df %>%
  filter(q_GroupMI < 0.05)

# ----------------------------
# 按绝对 log2FC 排序，取前20个
# ----------------------------
top20 <- df_sig %>%
  arrange(desc(abs(lfc_GroupMI))) %>%
  slice_head(n = 20)

# ----------------------------
# 输出
# ----------------------------
out_file <- sub("\\.xlsx$", "_top20_diff.xlsx", result_file)
write.xlsx(top20, out_file, overwrite = TRUE)

cat("✅ 已筛选出差异最明显的20个微生物\n")
cat("输出文件：", out_file, "\n")






# ======================================================
# 筛选 ANCOM-BC2 病毒表显著差异行
# ======================================================

library(openxlsx)
library(dplyr)

# ----------------------------
# 文件路径
# ----------------------------
file_ancom <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/四类微生物_ANCOMBC2_results.xlsx"
sheet_name <- "virus_Result"  # 病毒对应 sheet

# ----------------------------
# 读取数据
# ----------------------------
df <- read.xlsx(file_ancom, sheet = sheet_name)

# ----------------------------
# 筛选显著行
# ----------------------------
# 这里假设 ANCOM-BC2 列名为 q_GroupMI 和 diff_GroupMI 或 lfc_GroupMI
df_sig <- df %>%
  filter(q_GroupMI < 0.05, abs(lfc_GroupMI) > 1)  # 可改为 lfc_GroupMI，如果你的 log2FC 在那列

cat("共筛选出", nrow(df_sig), "个显著病毒\n")

# ----------------------------
# 保存显著结果（可选）
# ----------------------------
out_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/virus_sig_filtered.xlsx"
write.xlsx(df_sig, out_file, overwrite = TRUE)
cat("✅ 已保存显著病毒表至：", out_file, "\n")


# ============================================
# α & β 多样性分析（病毒显著行）
# 仅使用显著行的丰度数据
# 输出图保存到显著病毒表所在目录
# ============================================

# ----------------------------
# 依赖包
# ----------------------------
if(!requireNamespace("phyloseq", quietly=TRUE)) install.packages("phyloseq")
if(!requireNamespace("vegan", quietly=TRUE)) install.packages("vegan")
if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if(!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
if(!requireNamespace("ape", quietly=TRUE)) install.packages("ape")
if(!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if(!requireNamespace("stringr", quietly=TRUE)) install.packages("stringr")

library(phyloseq); library(vegan); library(ggplot2)
library(readxl); library(ape); library(dplyr)
library(stringr)

# -------------------------
# 1. 读取样本元数据
# -------------------------
meta_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/sample_metadata.xlsx"
meta <- read_excel(meta_file)
meta <- as.data.frame(meta)
rownames(meta) <- meta$SampleID

# -------------------------
# 2. 筛选后的显著病毒表路径
# -------------------------
sig_virus_file <- "E:/Python/MI_Analysis/metagenome/data_figures/filtered_data_1percent_change/ancombc2_results_OTU/virus_sig_filtered.xlsx"
df_sig <- read.xlsx(sig_virus_file, sheet = 1)

# -------------------------
# 输出目录（和显著病毒表同目录）
# -------------------------
out_dir <- dirname(sig_virus_file)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------
# 3. 提取丰度矩阵
# -------------------------
abund_cols <- grep("^(CON|MI)", colnames(df_sig), value=TRUE)
otu <- df_sig[, abund_cols]
rownames(otu) <- df_sig$taxon

# -------------------------
# 4. 样本交集及对齐
# -------------------------
shared_samples <- intersect(colnames(otu), meta$SampleID)
if(length(shared_samples) == 0) stop("No matching samples found")
otu <- otu[, shared_samples, drop=FALSE]
meta_sub <- meta[shared_samples, , drop=FALSE]

# -------------------------
# 5. 构建 phyloseq 对象
# -------------------------
ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows=TRUE),
  sample_data(meta_sub)
)

# -------------------------
# 6. α 多样性分析（Shannon、Simpson、InvSimpson）
# -------------------------
alpha_df <- estimate_richness(ps, measures=c("Shannon","Simpson","InvSimpson"))
alpha_df$Group <- sample_data(ps)$Group

metrics <- c("Shannon","Simpson","InvSimpson")
for(metric in metrics){
  stat_df <- alpha_df %>% group_by(Group) %>%
    summarise(mean=mean(!!sym(metric)), median=median(!!sym(metric)))
  
  if(length(unique(alpha_df$Group))==2){
    test <- wilcox.test(as.formula(paste(metric,"~ Group")), data=alpha_df)
    pval_text <- paste0("p = ", signif(test$p.value,3))
  } else {
    test <- kruskal.test(as.formula(paste(metric,"~ Group")), data=alpha_df)
    pval_text <- paste0("p = ", signif(test$p.value,3))
  }
  
  subtitle_text <- str_wrap(
    paste(paste(stat_df$Group, ": mean=", round(stat_df$mean,2), ", median=", round(stat_df$median,2), collapse="; "),
          "|", pval_text),
    width=60
  )
  
  p <- ggplot(alpha_df, aes(x=Group, y=!!sym(metric), fill=Group)) +
    geom_violin(trim=FALSE, alpha=1) +
    geom_boxplot(width=0.12, outlier.size=0.5, fill="white") +
    scale_fill_manual(values=c("Control"="#1F78B4", "MI"="#E31A1C")) +
    theme_classic() +
    theme(legend.position="none") +
    ggtitle(paste("Alpha diversity -", metric, "(Virus Sig)"), subtitle=subtitle_text)
  
  ggsave(file.path(out_dir, paste0("alpha_", metric, "_virus_sig.pdf")), p, width=6, height=5)
}

# -------------------------
# 7. β 多样性分析 (Bray-Curtis PCoA)
# -------------------------
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
dist_bc <- vegdist(t(otu_table(ps_rel)), method="bray")
pcoa_res <- pcoa(dist_bc)

pco_df <- data.frame(
  sample=rownames(pcoa_res$vectors),
  PC1=pcoa_res$vectors[,1],
  PC2=pcoa_res$vectors[,2],
  Group=sample_data(ps_rel)$Group
)

adonis_res <- adonis2(dist_bc ~ Group, data=data.frame(sample_data(ps_rel)))
pval_text <- paste0("PERMANOVA p = ", signif(adonis_res$`Pr(>F)`[1],3))
group_counts <- table(pco_df$Group)
n_text <- paste(names(group_counts), "n=", group_counts, collapse="; ")

p_pcoa <- ggplot(pco_df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) +
  stat_ellipse(level=0.95, linetype=2) +
  scale_color_manual(values=c("Control"="#1F78B4", "MI"="#E31A1C")) +
  theme_classic() +
  labs(title="PCoA (Bray-Curtis) - Virus Sig",
       subtitle=paste(n_text, "|", pval_text),
       x=paste0("PC1 (", round(pcoa_res$values$Relative_eig[1]*100,1), "%)"),
       y=paste0("PC2 (", round(pcoa_res$values$Relative_eig[2]*100,1), "%)"))

ggsave(file.path(out_dir, "beta_pcoa_bray_virus_sig.pdf"), p_pcoa, width=6, height=5)

cat("✅ α & β 多样性分析完成，图已保存到：", out_dir, "\n")
