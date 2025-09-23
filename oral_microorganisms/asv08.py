# 心梗相关口腔微生物筛选整合代码
# 功能：从口腔微生物数据中分析可能与心梗相关的微生物，输出所有微生物的丰度比值结果

# 1. 导入所需库
import pandas as pd

# 2. 数据导入与初步探查
# 2.1 读取数据文件
excel_file = pd.ExcelFile('ASV_table_filter.xlsx')  # 数据文件路径，需根据实际路径调整
sheet_names = excel_file.sheet_names
print(f"文件中工作表名称：{sheet_names}")  # 查看工作表名称，确认数据存储位置

# 提取Sheet1数据（根据实际工作表名称调整，此处默认Sheet1）
df = excel_file.parse('Sheet1')

# 2.2 查看数据基本信息
print("\n=== 数据基本信息 ===")
df.info()  # 查看数据类型、缺失值等
rows, columns = df.shape
print(f"数据规模：{rows}行（微生物样本）× {columns}列（字段）")

# 查看前5行数据，了解数据结构
print("\n=== 数据前5行内容 ===")
print(df.head().to_csv(sep='\t', na_rep='nan'))

# 3. 数据预处理
# 3.1 列名重命名：将"#OTU ID"改为"Microbe_ID"，增强可读性
df = df.rename(columns={'#OTU ID': 'Microbe_ID'})

# 3.2 提取微生物分类信息（从taxonomy列拆分界、门、纲、目、科、属、种）
# 按"; "分隔taxonomy列，拆分后生成7个分类层级列
taxonomy_df = df['taxonomy'].str.split('; ', expand=True)
taxonomy_df.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# 合并分类信息到原数据，删除原始taxonomy列
df = pd.concat([df, taxonomy_df], axis=1)
df = df.drop('taxonomy', axis=1)

print("\n=== 预处理后数据基本信息 ===")
print(df.info())  # 确认预处理后的数据结构
print(f"预处理后数据规模：{df.shape[0]}行 × {df.shape[1]}列")

# 4. 微生物丰度计算
# 4.1 筛选分组列：区分AMI组（心梗患者）和CON组（正常人）
ami_columns = [col for col in df.columns if col.startswith('AMI')]  # AMI组样本列
con_columns = [col for col in df.columns if col.startswith('CON')]  # CON组样本列
print(f"\n=== 分组样本数量 ===")
print(f"AMI组（心梗患者）样本数：{len(ami_columns)}")
print(f"CON组（正常人）样本数：{len(con_columns)}")

# 4.2 计算每组总丰度：按行求和，得到每个微生物在两组中的总丰度
df['AMI_Total_Abundance'] = df[ami_columns].sum(axis=1)  # AMI组总丰度
df['CON_Total_Abundance'] = df[con_columns].sum(axis=1)  # CON组总丰度

print("\n=== 丰度计算结果示例（前5行） ===")
print(df[['Microbe_ID', 'AMI_Total_Abundance', 'CON_Total_Abundance']].head())

# 5. 微生物丰度比值计算（不进行阈值筛选）
# 5.1 计算丰度比值（AMI组总丰度 / CON组总丰度），避免分母为0的情况
def calculate_abundance_ratio(ami_abd, con_abd):
    if con_abd == 0:
        return 0  # 若CON组丰度为0，比值设为0（表示仅在AMI组存在）
    else:
        return round(ami_abd / con_abd, 4)  # 保留4位小数，便于查看

df['Abundance_Ratio'] = df.apply(
    lambda row: calculate_abundance_ratio(row['AMI_Total_Abundance'], row['CON_Total_Abundance']),
    axis=1
)

# 5.2 不进行阈值筛选，保留所有微生物数据
all_microbes = df.copy()  # 复制所有结果

# 5.3 整理结果：保留关键字段（微生物标识、分类信息、丰度比值）
result_columns = ['Microbe_ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
                  'AMI_Total_Abundance', 'CON_Total_Abundance', 'Abundance_Ratio']
all_microbes = all_microbes[result_columns]

print(f"\n=== 结果统计 ===")
print(f"分析的微生物总种类数：{len(all_microbes)}")
print("\n=== 结果示例（前10行） ===")
print(all_microbes.head(10).to_csv(sep='\t', na_rep='nan'))

# 6. 结果输出：保存所有微生物的分析结果到CSV文件
output_path = 'all_microbes.csv'  # 输出文件路径，可根据需求调整
all_microbes.to_csv(output_path, index=False, encoding='utf-8-sig')
print(f"\n=== 结果输出完成 ===")
print(f"分析结果已保存至：{output_path}")
print(f"文件包含字段：{', '.join(result_columns)}")
