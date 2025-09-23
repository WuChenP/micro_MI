import pandas as pd
from statsmodels.stats.multitest import multipletests

# 读取数据
df = pd.read_csv('merged_microbes_with_metrics_fixed.csv')

# 检查p值列是否存在
if 'p_value' not in df.columns:
    raise ValueError("数据中未找到'p_value'列，请检查列名是否正确")

# 提取p值
p_values = df['p_value'].values

# 使用Benjamini-Hochberg方法计算FDR（最常用的FDR校正方法）
# 返回值包括：是否拒绝原假设、校正后的p值（即FDR）、校正前的临界值、校正后的临界值
rejected, fdr_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

# 将FDR值添加到数据框
df['fdr'] = fdr_values

# 保存结果
output_file = 'merged_microbes_with_metrics_fixed_with_fdr.csv'
df.to_csv(output_file, index=False)

print(f"FDR计算完成，结果已保存至 {output_file}")
print(f"前5行FDR结果预览：")
print(df[['p_value', 'fdr']].head())
