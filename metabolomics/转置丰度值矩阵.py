import pandas as pd

# 1. 读取 tsv 文件
df = pd.read_csv("abundance.tsv", sep="\t")

# 2. 查看前几行（可选）
print("原始数据：")
print(df.head())

# 3. 检查行列格式（WGCNA要求行为代谢物，列为样本）
# 如果你当前是“样本为行，代谢物为列”，就要转置
df_t = df.set_index(df.columns[0]).T

# 4. 输出转置后的结果
print("\n转置后数据：")
print(df_t.head())

# 5. 保存成新的文件
df_t.to_csv("abundance_WGCNA.tsv", sep="\t", index=True)
print("\n✅ abundance_WGCNA.tsv 文件")
