import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

# 读取数据
df = pd.read_csv('merged_microbes_with_metrics_fixed_FDR.csv')

# 提取AMI和CON组的样本列
ami_cols = [col for col in df.columns if col.startswith('AMI') and not col.startswith('AMI_')]
con_cols = [col for col in df.columns if col.startswith('CON') and not col.startswith('CON_')]

# 重新计算p_value
new_p_values = []
for idx, row in df.iterrows():
    ami_data = row[ami_cols].values
    con_data = row[con_cols].values

    # 确保数据有效
    if len(ami_data) > 0 and len(con_data) > 0:
        try:
            # 使用Mann-Whitney U检验
            stat, p_value = mannwhitneyu(ami_data, con_data, alternative='two-sided')
            new_p_values.append(p_value)
        except:
            # 如果检验失败，使用原值
            new_p_values.append(row['p_value'])
    else:
        new_p_values.append(row['p_value'])

# 更新p_value列
df['p_value'] = new_p_values

# 由于p值改变，可能需要重新计算FDR（错误发现率）
from statsmodels.stats.multitest import multipletests

rejected, p_adjusted, _, _ = multipletests(df['p_value'], alpha=0.05, method='fdr_bh')
df['FDR'] = p_adjusted

# 保存更新后的数据
df.to_csv('merged_microbes_with_metrics_recalculated.csv', index=False)