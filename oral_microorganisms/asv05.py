import pandas as pd
import matplotlib.pyplot as plt
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import anosim, permanova
import seaborn as sns

# 读取 Excel 文件
excel_file = pd.ExcelFile('merged_filtered_microbes.xlsx')

# 获取指定工作表中的数据
df = excel_file.parse('Sheet1')

# 提取 AMI 组和 CON 组的列名
ami_columns = [col for col in df.columns if col.startswith('AMI')]
con_columns = [col for col in df.columns if col.startswith('CON')]

# 合并 AMI 组和 CON 组的列名
all_sample_columns = ami_columns + con_columns

# 提取样本数据
sample_data = df[all_sample_columns]

# 提取分组信息，假设 AMI 组为 1，CON 组为 0
group_info = []
for col in all_sample_columns:
    if col.startswith('AMI'):
        group_info.append(1)
    else:
        group_info.append(0)

# 计算 Bray - Curtis 距离矩阵
distance_matrix = beta_diversity('braycurtis', sample_data.T)

# 进行主坐标分析（PCoA）
pcoa_result = pcoa(distance_matrix)

# 获取前两个主坐标
pc1 = pcoa_result.samples['PC1']
pc2 = pcoa_result.samples['PC2']

# 创建一个新的 DataFrame 用于绘图
plot_df = pd.DataFrame({
    'PC1': pc1,
    'PC2': pc2,
    'Group': group_info
})

# 设置图片清晰度
plt.rcParams['figure.dpi'] = 300

# 绘制 PCoA 图
plt.figure(figsize=(10, 8))
sns.scatterplot(data=plot_df, x='PC1', y='PC2', hue='Group', palette='Set1')
plt.xlabel(f'PC1 ({pcoa_result.proportion_explained["PC1"] * 100:.2f}%)')
plt.ylabel(f'PC2 ({pcoa_result.proportion_explained["PC2"] * 100:.2f}%)')
plt.title('基于 Bray - Curtis 距离的 PCoA 分析')
plt.legend(title='分组', loc='upper right')

# 进行 ANOSIM 检验，评估组间差异是否显著
anosim_result = anosim(distance_matrix, group_info, permutations=999)
print(f'ANOSIM 检验结果：R 值 = {anosim_result[0]}, p 值 = {anosim_result[1]}')
# 进行 PERMANOVA 检验，评估组间差异是否显著
permanova_result = permanova(distance_matrix, group_info, permutations=999)
print(f'PERMANOVA 检验结果：F 值 = {permanova_result.statistic}, p 值 = {permanova_result.p_value}')

plt.show()
