import pandas as pd

# 读入你的原始表格
df = pd.read_excel("E:/代谢组学/AMI_vs_CON_lev1.xlsx")  # 修改为你的文件名

# 1. 提取样本列
ami_cols = [c for c in df.columns if c.startswith('AMI')]
con_cols = [c for c in df.columns if c.startswith('CON')]
sample_cols = ami_cols + con_cols

# 2. 构建 abundance 表（代谢物×样本）
abundance = df[sample_cols].copy()
abundance.index = df['Compounds']  # 把代谢物名作为行索引

# 🚩在这里加过滤步骤
zero_cut = 0.9   # ≥90%零值就删
lib_cut = 1000   # 总丰度<1000的样本就删

# 先按零值比例过滤代谢物
keep_rows = (abundance == 0).sum(axis=1) / abundance.shape[1] < zero_cut
abundance = abundance.loc[keep_rows]

# 再按总丰度过滤样本
keep_samples = abundance.sum(axis=0) > lib_cut
abundance = abundance.loc[:, keep_samples]

# 3. 构建 metadata 表（样本×分组）
# 只保留还存在的样本
all_samples = abundance.columns
groups = ['AMI' if s.startswith('AMI') else 'CON' for s in all_samples]

metadata = pd.DataFrame({
    'SampleID': all_samples,
    'Group': groups
})

# 4. 保存文件
abundance.to_csv('abundance.tsv', sep='\t')
metadata.to_csv('metadata.tsv', sep='\t', index=False)

print('丰度表和分组信息已保存：abundance.tsv 和 metadata.tsv')
