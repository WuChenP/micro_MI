import pandas as pd
import numpy as np

from sklearn.cross_decomposition import PLSRegression
# 读取数据
df = pd.read_excel("../AMI_vs_CON_base_data.xlsx")


# 提取实验组和对照组列
ami_cols = [col for col in df.columns if col.startswith("AMI")]
con_cols = [col for col in df.columns if col.startswith("CON")]


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
df["VIP2"] = vip_scores

df.to_excel("../AMI_vs_CON_base_data.xlsx",index=False)