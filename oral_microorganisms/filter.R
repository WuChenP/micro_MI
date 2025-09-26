# ========== 设置工作路径 ==========
setwd("D:/机器学习与生物信息学/北京朝阳医院项目/R")  # 修改为你想保存的路径

# ========== 1. 读取数据 ==========
res <- read.csv("ANCOMBC2_results_sig.csv", header = TRUE, check.names = FALSE)

# ========== 2. 筛选 Species 不为空 ==========
res_species <- subset(res, Species != "" & !is.na(Species))

# ========== 3. 保存新表 ==========
write.csv(res_species, "ANCOMBC2_results_sig_with_species.csv", row.names = FALSE)

# ========== 4. 输出统计信息 ==========
cat("✅ 已生成 ANCOMBC2_results_sig_with_species.csv\n")
cat("原始数据行数:", nrow(res), "\n")
cat("筛选后行数:", nrow(res_species), "\n")
