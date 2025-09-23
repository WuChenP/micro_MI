import pandas as pd

# 1. 读取数据文件
# 读取表一（包含需要筛选的微生物ID及分类信息）
df1 = pd.read_csv('significant_microbes_filtered.csv')

# 读取表二（包含所有微生物数据）
excel_file = pd.ExcelFile('ASV_table_filter.xlsx')
df2 = excel_file.parse('Sheet1')  # 假设只有一个工作表'Sheet1'

# 2. 从表二中筛选出表一中存在的微生物
filtered_df = df2[df2['#OTU ID'].isin(df1['Microbe_ID'])]

# 3. 按照表一的样式整合数据
# 合并两个数据框，基于ID列进行匹配
merged_df = pd.merge(df1, filtered_df, left_on='Microbe_ID', right_on='#OTU ID', how='outer')

# 从taxonomy列中提取分类信息，补充到对应的分类列中
taxonomy_levels = {
    'Kingdom': 'k__',
    'Phylum': 'p__',
    'Class': 'c__',
    'Order': 'o__',
    'Family': 'f__',
    'Genus': 'g__',
    'Species': 's__'
}

for col, prefix in taxonomy_levels.items():
    merged_df[col] = merged_df.apply(
        lambda row: row[col] if pd.notnull(row[col])
        else (row['taxonomy'].split('; ')[[p for p in taxonomy_levels.values()].index(prefix)].replace(prefix, '')
              if pd.notnull(row['taxonomy']) and prefix in row['taxonomy']
              else None),
        axis=1
    )

# 4. 清理不需要的列
merged_df = merged_df.drop(columns=['#OTU ID', 'taxonomy'])

# 5. 保存结果
merged_df.to_excel('merged_filtered_microbes.xlsx', index=False)

print(f"整合完成，共处理 {len(merged_df)} 条微生物数据")
