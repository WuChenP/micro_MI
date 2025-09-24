import pandas as pd

# 加载数据
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

# 获取要移动的列的数据
columns_to_move = ['AMI_Total_Abundance', 'CON_Total_Abundance', 'Abundance_Ratio']
moved_data = df[columns_to_move]

# 删除原位置的列
df = df.drop(columns=columns_to_move)

# 获取 con9 列的索引
con9_index = df.columns.get_loc('CON9')

# 在 con9 列后插入移动的列
for col in reversed(columns_to_move):
    df.insert(con9_index + 1, col, moved_data[col])

# 将结果保存为csv文件
csv_path = 'merged_microbes_with_metrics_fixed.csv'
df.to_csv(csv_path)