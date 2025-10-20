#-----------------------------------------
# 🔹 WGCNA 完整分析脚本（增强版）
# 包含：网络构建 + 模块分析 + Hub代谢物识别 + 模块关系可视化 + Cytoscape网络
#-----------------------------------------

library(WGCNA)
library(tidyverse)
library(pheatmap)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

#-----------------------------------------
# 1️⃣ 数据读取
#-----------------------------------------
feature <- read.delim("abundance.tsv", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim("metadata.tsv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

if (!"Group" %in% colnames(metadata)) stop("❌ metadata.tsv 缺少 'Group' 列！")

datExpr <- t(feature)
datExpr <- datExpr[rownames(metadata), , drop = FALSE]
stopifnot(all(rownames(datExpr) == rownames(metadata)))

# 检查缺失比例
na_prop <- colMeans(is.na(datExpr))
cat("缺失比例超过20%的代谢物数量：", sum(na_prop > 0.2), "\n")

#-----------------------------------------
# 2️⃣ 去除低质量代谢物
#-----------------------------------------
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

#-----------------------------------------
# 3️⃣ 软阈值选择
#-----------------------------------------
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
cat("选择的软阈值 power 为：", softPower, "\n")

dir.create("WGCNA/结果2", recursive = TRUE, showWarnings = FALSE)

#-----------------------------------------
# 4️⃣ 网络构建 + 模块检测
#-----------------------------------------
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "metaboliteTOM",
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

#-----------------------------------------
# 5️⃣ 模块与性状相关性分析
#-----------------------------------------
metadata_num <- metadata
metadata_num$Group <- as.numeric(factor(metadata$Group, levels = c("CON", "AMI")))

moduleTraitCor <- cor(MEs, metadata_num, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

#-----------------------------------------
# 6️⃣ 计算 MM 与 GS
#-----------------------------------------
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
geneTraitSignificance <- as.data.frame(cor(datExpr, metadata_num$Group, use = "p"))
names(geneTraitSignificance) <- "Group"

#-----------------------------------------
# 7️⃣ Hub代谢物识别
#-----------------------------------------
hubThreshold_MM <- 0.8
hubThreshold_GS <- 0.2

dir.create("WGCNA/结果2/HubAnalysis", showWarnings = FALSE)

hubMetabolites <- c()
mod_colors <- unique(moduleColors)
for (i in seq_along(mod_colors)) {
  module <- mod_colors[i]
  modMetabolites <- names(datExpr)[moduleColors == module]
  ME_col <- paste0("ME", i)
  MM <- geneModuleMembership[modMetabolites, ME_col]
  GS <- geneTraitSignificance[modMetabolites, "Group"]
  
  hub <- modMetabolites[abs(MM) > hubThreshold_MM & abs(GS) > hubThreshold_GS]
  hubMetabolites <- c(hubMetabolites, hub)
  
  # 保存完整模块与 hub
  hubTable <- data.frame(Metabolite = modMetabolites, MM = MM, GS = GS)
  write.csv(hubTable, paste0("WGCNA/结果2/HubAnalysis/HubMetabolites_", module, ".csv"), row.names = FALSE)
  
  if (length(hub) > 0) {
    hubTableTop <- hubTable[hubTable$Metabolite %in% hub, ]
    write.csv(hubTableTop, paste0("WGCNA/结果2/HubAnalysis/HubMetabolites_Top_", module, ".csv"), row.names = FALSE)
  }
}

hubMetabolites <- unique(hubMetabolites)

#-----------------------------------------
# 8️⃣ Cytoscape 网络文件生成 (仅 hub 代谢物)
#-----------------------------------------
dir.create("WGCNA/结果2/Cytoscape", showWarnings = FALSE)

# 构建 TOM 网络矩阵
TOM <- TOMsimilarityFromExpr(datExpr[, hubMetabolites], power = softPower)

# 节点文件
nodeTable <- data.frame(
  nodeName = hubMetabolites,
  altName = hubMetabolites,
  module = moduleColors[match(hubMetabolites, names(datExpr))],
  MM = geneModuleMembership[hubMetabolites, match(paste0("ME", 1:length(MEs)), names(MEs))],
  GS = geneTraitSignificance[hubMetabolites, "Group"]
)
write.table(nodeTable, "WGCNA/结果2/Cytoscape/CytoscapeNodes.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# 边文件，只保留 TOM > 0.1
edges <- which(TOM > 0.1, arr.ind = TRUE)
edges <- edges[edges[,1] < edges[,2], ]  # 去掉重复
edgeTable <- data.frame(
  fromNode = hubMetabolites[edges[,1]],
  toNode = hubMetabolites[edges[,2]],
  weight = TOM[edges],
  direction = "undirected",
  fromAltName = hubMetabolites[edges[,1]],
  toAltName = hubMetabolites[edges[,2]]
)
write.table(edgeTable, "WGCNA/结果2/Cytoscape/CytoscapeEdges.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ Cytoscape 网络文件已生成 (仅 hub 代谢物)。\n")
cat("节点文件：WGCNA/结果2/Cytoscape/CytoscapeNodes.txt\n")
cat("边文件：WGCNA/结果2/Cytoscape/CytoscapeEdges.txt\n")
