# 心梗相关口腔微生物筛选整合代码
# 功能：从口腔微生物数据中筛选可能与心梗相关的微生物，输出筛选结果文件（包含log2FC、p-value和LDA score）

# 1. 导入所需库
import pandas as pd
import numpy as np
from scipy import stats
from scipy.sparse import csr_matrix
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from statsmodels.stats.multitest import multipletests
import warnings

warnings.filterwarnings('ignore')  # 忽略警告信息，使输出更清晰

# 2. 数据导入与初步探查
# 2.1 读取数据文件
excel_file = pd.ExcelFile('/mnt/ASV_table_filter.xlsx')  # 数据文件路径，需根据实际路径调整
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
# 4.1 筛选分组列：区分AMI组（心梗患者）和CON组（正常人