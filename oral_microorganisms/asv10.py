import pandas as pd

# 读取文件
excel_file = pd.ExcelFile('ASV_table_filter.xlsx')

# 获取所有表名
sheet_names = excel_file.sheet_names
sheet_names

# 获取指定工作表中的数据
df1 = excel_file.parse('Sheet1')

# 查看数据的基本信息
print('数据基本信息：')
df1.info()

# 查看数据集行数和列数
rows, columns = df1.shape

if rows < 100 and columns < 20:
    # 短表数据（行数少于100且列数少于20）查看全量数据信息
    print('数据全部内容信息：')
    print(df1.to_csv(sep='\t', na_rep='nan'))
else:
    # 长表数据查看数据前几行信息
    print('数据前几行内容信息：')
    print(df1.head().to_csv(sep='\t', na_rep='nan'))

# 读取数据
df2 = pd.read_csv('all_microbes.csv')

print('数据基本信息：')
df2.info()

# 查看数据集行数和列数
rows, columns = df2.shape

if rows < 100 and columns < 20:
    # 短表数据（行数少于100且列数少于20）查看全量数据信息
    print('数据全部内容信息：')
    print(df2.to_csv(sep='\t', na_rep='nan'))
else:
    # 长表数据查看数据前几行信息
    print('数据前几行内容信息：')
    print(df2.head().to_csv(sep='\t', na_rep='nan'))

# 将 df1 的 #OTU ID 列设置为索引
df1.set_index('#OTU ID', inplace=True)

# 将 df2 的 Microbe_ID 列设置为索引
df2.set_index('Microbe_ID', inplace=True)

# 根据索引合并数据
merged_df = pd.concat([df2, df1], axis=1, join='inner')

# # 保存合并后的数据
# csv_path = 'merged_microbes_data.csv'
# merged_df.to_csv(csv_path)

import pandas as pd

# 读取文件
df = pd.read_csv('merged_microbes_data.csv')

# 删除最后一个字段
df = df.iloc[:, :-1]

# 保存处理后的数据
csv_path = 'merged_microbes_data.csv'
df.to_csv(csv_path, index=False)