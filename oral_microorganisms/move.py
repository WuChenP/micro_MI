import pandas as pd

# 加载数据集
df = pd.read_csv('merged_microbes_data_all_new.csv')

# 要移动的列
columns_to_move = ['AMI_Total_Abundance', 'CON_Total_Abundance', 'Abundance_Ratio']

# 目标列
target_column = 'AMI_avg_abundance'

# 获取目标列的索引
target_index = df.columns.get_loc(target_column)

# 新的列顺序：在目标列索引前插入要移动的列，然后拼接其他列
new_columns = df.columns.tolist()
for col in reversed(columns_to_move):
    new_columns.insert(target_index, new_columns.pop(new_columns.index(col)))

# 重新排列 DataFrame 的列
df = df[new_columns]

# 将结果保存为 CSV 文件
csv_path = 'merged_microbes_data_all_new1.csv'
df.to_csv(csv_path, index=False)