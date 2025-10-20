
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.cross_decomposition import PLSRegression
# 读取数据
df = pd.read_excel("E:/代谢组学/AMI_vs_CON_lev1.xlsx")

# 提取实验组和对照组列
ami_cols = [col for col in df.columns if col.startswith("AMI")]
con_cols = [col for col in df.columns if col.startswith("CON")]

#########Step1:  计算FC和log2FC
# 计算均值
df["AMI_mean"] = df[ami_cols].mean(axis=1)
df["CON_mean"] = df[con_cols].mean(axis=1)

# Fold Change
df["FoldChange"] = df["AMI_mean"] / df["CON_mean"]
df["Log2FC"] = np.log2(df["FoldChange"].replace(0, np.nan))  # 避免除零


#########Step2:  计算p值和FDR(q值)
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


#########Step3:  计算VIP值
# 构建 X 和 y
X = df[ami_cols + con_cols].values
y = np.array([1]*len(ami_cols) + [0]*len(con_cols))  # 1=AMI, 0=CON
##########################
#y = np.tile(y, (len(df), 1)).T  # 每个化合物一行，需要转置

# PLS 模型
pls = PLSRegression(n_components=2)
pls.fit(X.T, y)

# VIP 计算函数
def calculate_vip(pls, X, y):
    t = pls.x_scores_
    w = pls.x_weights_
    q = pls.y_loadings_
    p, h = w.shape
    vip = np.zeros((p,))
    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    total_s = np.sum(s)
    for i in range(p):
        weight = np.array([(w[i,j]**2) * s[j,0] for j in range(h)])
        vip[i] = np.sqrt(p * (weight.sum())/total_s)
    return vip
# 计算 VIP
vip_scores = calculate_vip(pls, X.T, y)
df["VIP"] = vip_scores


df.to_excel("E:/代谢组学/AMI_vs_CON_lev1_Statistics.xlsx", index=False)
