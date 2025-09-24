import pandas as pd
import numpy as np
# 读取数据
df = pd.read_excel("AMI_vs_CON_base_data.xlsx")

'''
# 基于统计学指标筛选 &(df['FDR'] < 0.05)

df_filtered= df[(df['P_value'] < 0.05) &(df['VIP'] >= 1) &
              ((df['FoldChange'] >= 1.2) | (df['FoldChange'] <= 0.83))]'''
# 首先创建筛选条件
significant_condition = (df['VIP'] >= 1) & (df['P_value'] < 0.05)
up_condition = significant_condition & (df['FoldChange'] >2)
down_condition = significant_condition & (df['FoldChange'] <0.5)

# 创建type列
conditions = [
    up_condition,
    down_condition,
    ~significant_condition  # 不显著的都是insig
]

choices = ['up', 'down', 'insig']

df['type'] = np.select(conditions, choices, default='insig')
print((df['type'] == 'up').sum(),(df['type'] == 'down').sum(),(df['type'] == 'insig').sum())



#保存文档
df.to_excel("AMI_vs_CON_base_data.xlsx", index=False)