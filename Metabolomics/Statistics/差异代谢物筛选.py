import pandas as pd

# 读取数据
df = pd.read_excel("E:/代谢组学/AMI_vs_CON_lev1_Statistics.xlsx")

# 设置阈值
fdr_threshold = 0.05
fc_up = 2
fc_down = 0.5

# 新增 type 列
df["type"] = "insignificant"  # 先默认所有为 insignificant

# 仅对 FDR 满足条件的进行分类
df.loc[(df["FDR"] < fdr_threshold) & (df["FoldChange"] > fc_up), "type"] = "up"
df.loc[(df["FDR"] < fdr_threshold) & (df["FoldChange"] < fc_down), "type"] = "down"

# 打印结果统计
print(df["type"].value_counts())

# 预览前几行
print(df[["Compounds", "FDR", "FoldChange", "type"]].head())

# 保存结果
df.to_excel("E:/代谢组学/传统统计学/AMI_vs_CON_lev1_Statistics_filter.xlsx", index=False)
print("结果已保存到E:/代谢组学/传统统计学/AMI_vs_CON_lev1_Statistics_filter.xlsx")
