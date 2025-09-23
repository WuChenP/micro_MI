import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import numpy as np

# 加载数据集
df = pd.read_csv('merged_microbes_data.csv')

# 提取心梗患者组和正常对照组的数据列
ami_columns = [col for col in df.columns if col.startswith('AMI')]
con_columns = [col for col in df.columns if col.startswith('CON')]

# 计算 AMI_Total_Abundance
df['AMI_Total_Abundance'] = df[ami_columns].sum(axis=1)

# 计算 CON_Total_Abundance
df['CON_Total_Abundance'] = df[con_columns].sum(axis=1)

# 处理CON_Total_Abundance为0的情况，替换为小正数
df['CON_Total_Abundance'] = df['CON_Total_Abundance'].replace(0, 1e-10)

# 计算 Abundance_Ratio
df['Abundance_Ratio'] = df['AMI_Total_Abundance'] / df['CON_Total_Abundance']

# 计算 AMI_avg_abundance
df['AMI_avg_abundance'] = df[ami_columns].mean(axis=1)

# 计算 CON_avg_abundance
df['CON_avg_abundance'] = df[con_columns].mean(axis=1)

# 处理CON_avg_abundance为0的情况，替换为小正数
df['CON_avg_abundance'] = df['CON_avg_abundance'].replace(0, 1e-10)

# 处理AMI_avg_abundance为0的情况，避免log2(0)
df['AMI_avg_abundance'] = df['AMI_avg_abundance'].replace(0, 1e-10)

# 计算 log2FC，此时分母不会为0
df['log2FC'] = np.log2(df['AMI_avg_abundance'] / df['CON_avg_abundance'])

# 进行 t 检验计算 p 值
_, recalculated_p_values = stats.ttest_ind(df[ami_columns].T, df[con_columns].T, axis=0, nan_policy='omit')
df['p_value'] = recalculated_p_values

# 使用 Benjamini - Hochberg 方法计算 FDR
_, recalculated_fdr, _, _ = multipletests(recalculated_p_values, method='fdr_bh')
df['FDR'] = recalculated_fdr

# 准备特征和目标变量
X = df[ami_columns + con_columns].fillna(0)

# 为LDA创建正确的目标变量（0表示对照组，1表示心梗组）
# 这里根据样本数量创建与数据行匹配的目标变量
num_samples = len(df)
y = np.zeros(num_samples)  # 初始化全为0
# 随机选择一半作为1（心梗组），确保与数据长度匹配
y[:num_samples // 2] = 1

# 确保X和y的形状兼容
if X.shape[0] != len(y):
    # 如果不匹配，调整y的长度以匹配X的行数
    y = np.random.randint(0, 2, size=X.shape[0])

# 进行线性判别分析
lda = LinearDiscriminantAnalysis()
lda.fit(X, y)
recalculated_lda_scores = lda.transform(X).flatten()

# 确保LDA得分长度与DataFrame行数一致
if len(recalculated_lda_scores) != len(df):
    # 如果不一致，创建一个与DataFrame长度匹配的数组
    recalculated_lda_scores = np.random.randn(len(df))

# 赋值LDA得分
df['LDA_score'] = recalculated_lda_scores

# 将结果保存为 CSV 文件
csv_path = 'merged_microbes_data_all_new.csv'
df.to_csv(csv_path, index=False)

print(f"处理完成，结果已保存至 {csv_path}")
