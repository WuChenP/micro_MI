import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# 读取 Excel 文件
excel_file = pd.ExcelFile('merged_filtered_microbes.xlsx')
# 获取指定工作表中的数据
df = excel_file.parse('Sheet1')

# 提取 AMI 组和 CON 组的列名
ami_columns = [col for col in df.columns if col.startswith('AMI')]
con_columns = [col for col in df.columns if col.startswith('CON')]

# 计算每个样本的观测物种数（Observed OTUs）
df['AMI_Observed_OTUs'] = df[ami_columns].astype(bool).sum(axis=1)
df['CON_Observed_OTUs'] = df[con_columns].astype(bool).sum(axis=1)

# 定义计算 Chao1 指数的函数
def chao1(s_obs, f1, f2):
    if f2 == 0:
        return s_obs
    return s_obs + (f1 ** 2) / (2 * f2)

# 计算每个样本的 Chao1 指数
df['AMI_Chao1'] = df[ami_columns].apply(lambda x: chao1(x.astype(bool).sum(), (x == 1).sum(), (x == 2).sum()), axis=1)
df['CON_Chao1'] = df[con_columns].apply(lambda x: chao1(x.astype(bool).sum(), (x == 1).sum(), (x == 2).sum()), axis=1)

# 定义计算 Shannon 指数的函数
def shannon(abundances):
    total = np.sum(abundances)
    if total == 0:
        return 0
    p_i = abundances / total
    p_i = p_i[p_i > 0]  # 避免 log(0) 的情况
    return -np.sum(p_i * np.log(p_i))

# 计算每个样本的 Shannon 指数
df['AMI_Shannon'] = df[ami_columns].apply(shannon, axis=1)
df['CON_Shannon'] = df[con_columns].apply(shannon, axis=1)

# 定义计算 Simpson 指数的函数
def simpson(abundances):
    total = np.sum(abundances)
    if total == 0:
        return 0
    p_i = abundances / total
    return 1 - np.sum(p_i ** 2)

# 计算每个样本的 Simpson 指数
df['AMI_Simpson'] = df[ami_columns].apply(simpson, axis=1)
df['CON_Simpson'] = df[con_columns].apply(simpson, axis=1)

# 整理数据以便进行统计检验和绘图
alpha_indices = ['Observed_OTUs', 'Chao1', 'Shannon', 'Simpson']
results = []
# 用于存储统计检验结果
stats_results = {}

for index in alpha_indices:
    ami_values = df[f'AMI_{index}']
    con_values = df[f'CON_{index}']

    # 正态性检验
    _, p_ami = stats.shapiro(ami_values)
    _, p_con = stats.shapiro(con_values)

    # 根据正态性选择检验方法
    if p_ami > 0.05 and p_con > 0.05:
        # 数据符合正态分布，使用独立样本 t 检验
        stat, p_val = stats.ttest_ind(ami_values, con_values)
        test_type = 't-test'
    else:
        # 数据不符合正态分布，使用 Wilcoxon 秩和检验
        stat, p_val = stats.mannwhitneyu(ami_values, con_values)
        test_type = 'Wilcoxon test'

    stats_results[index] = {
        'test_type': test_type,
        'statistic': stat,
        'p_value': p_val
    }

    # 整理数据用于绘图
    for value, group in zip(ami_values, ['AMI'] * len(ami_values)):
        results.append({
            'Alpha_Index': index,
            'Group': group,
            'Value': value
        })
    for value, group in zip(con_values, ['CON'] * len(con_values)):
        results.append({
            'Alpha_Index': index,
            'Group': group,
            'Value': value
        })

# 将结果转换为 DataFrame
results_df = pd.DataFrame(results)

# 绘制箱线图
plt.figure(figsize=(12, 8))
sns.boxplot(x='Alpha_Index', y='Value', hue='Group', data=results_df)
sns.stripplot(x='Alpha_Index', y='Value', hue='Group', data=results_df, dodge=True, alpha=0.5, jitter=True, color='black')
plt.title('α多样性指数在 AMI 组和 CON 组的差异')
plt.xlabel('α多样性指数')
plt.ylabel('指数值')
plt.legend(title='分组')

# 在图上标注统计检验结果
for i, index in enumerate(alpha_indices):
    ami_mean = df[f'AMI_{index}'].mean()
    con_mean = df[f'CON_{index}'].mean()
    test_res = stats_results[index]
    y_pos = max(df[f'AMI_{index}'].max(), df[f'CON_{index}'].max()) * 1.1
    plt.text(i, y_pos, f'{test_res["test_type"]}\np={test_res["p_value"]:.3f}',
             ha='center', va='center', fontsize=9)

plt.tight_layout()
plt.savefig('alpha_diversity_plot.png', dpi=300)
plt.show()

# 打印统计检验结果
for index, res in stats_results.items():
    print(f'{index} 指数统计检验结果：')
    print(f'检验方法：{res["test_type"]}')
    print(f'统计量：{res["statistic"]:.3f}')
    print(f'p 值：{res["p_value"]:.3f}')
    print('---')