import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
# 读取数据
df = pd.read_excel("../AMI_vs_CON_pre.xlsx")

# 提取实验组和对照组列
ami_cols = [col for col in df.columns if col.startswith("AMI")]
con_cols = [col for col in df.columns if col.startswith("CON")]

# 计算 P-value
#独立样本 t 检验（两组均值比较）
p_values = []
for i in range(len(df)):
    ami_vals = df.loc[i, ami_cols].astype(float).values
    con_vals = df.loc[i, con_cols].astype(float).values
    stat, p = ttest_ind(ami_vals, con_vals, equal_var=False)  # Welch's t-test
    p_values.append(p)

df["P_value"] = p_values

# FDR 校正
_, fdrs, _, _ = multipletests(df["P_value"], method='fdr_bh')
df["FDR"] = fdrs


df.to_excel("../AMI_vs_CON_pre.xlsx",index=False)



