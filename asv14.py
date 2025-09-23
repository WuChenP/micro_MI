import pandas as pd

# 加载数据
from statsmodels.stats.multitest import multipletests

df = pd.read_csv('merged_microbes_with_metrics_fixed.csv')

print('数据基本信息：')
df.info()

# 查看数据集行数和列数
rows, columns = df.shape

if rows < 100 and columns < 20:
    # 短表数据（行数少于100且列数少于20）查看全量数据信息
    print('数据全部内容信息：')
    print(df.to_csv(sep='\t', na_rep='nan'))
else:
    # 长表数据查看数据前几行信息
    print('数据前几行内容信息：')
    print(df.head().to_csv(sep='\t', na_rep='nan'))

# 检查 p_value 列的数据类型
print('p_value 列数据类型：', df['p_value'].dtype)

# 检查 p_value 列是否存在缺失值
print('p_value 列缺失值数量：', df['p_value'].isnull().sum())

# 检查 p_value 列是否存在非数值数据
non_numeric = pd.to_numeric(df['p_value'], errors='coerce').isnull()
print('p_value 列非数值数据数量：', non_numeric.sum())

# 若存在缺失值，用 1 填充
if df['p_value'].isnull().sum() > 0:
    df['p_value'] = df['p_value'].fillna(1)

# 重新应用 Benjamini - Hochberg 方法计算 FDR
df['FDR'] = multipletests(df['p_value'], method='fdr_bh')[1]

# 将结果保存为 CSV 文件
csv_path = 'merged_microbes_with_metrics_fixed_FDR.csv'
df.to_csv(csv_path)