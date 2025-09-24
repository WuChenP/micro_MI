import pandas as pd
import numpy as np

# 读取数据
df = pd.read_excel("../AMI_vs_CON_base(里面有数据对照).xlsx")

# 提取实验组和对照组列
ami_cols = [col for col in df.columns if col.startswith("AMI")]
con_cols = [col for col in df.columns if col.startswith("CON")]
# 计算均值
df["AMI_mean"] = df[ami_cols].mean(axis=1)
df["CON_mean"] = df[con_cols].mean(axis=1)

# Fold Change
df["FoldChange"] = df["AMI_mean"] / df["CON_mean"]
df["Log2FC"] = np.log2(df["FoldChange"].replace(0, np.nan))  # 避免除零
df.to_excel("../AMI_vs_CON_base(里面有数据对照).xlsx", index=False)
