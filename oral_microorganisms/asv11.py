import pandas as pd
from scipy.stats import mannwhitneyu  # 非参数检验，适配稀疏数据
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import numpy as np
import warnings

# 抑制无关警告
warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in log2')
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in log10')


# 1. 读取数据并强制数据类型校验
def load_and_validate_data(file_path):
    """读取数据并确保样本列为数值类型，非数值转为0（微生物数据中代表未检出）"""
    df = pd.read_csv(file_path)

    # 识别样本列（AMI/CON开头）
    ami_columns = [col for col in df.columns if col.startswith('AMI')]
    con_columns = [col for col in df.columns if col.startswith('CON')]

    print(f"识别到AMI组样本列: {len(ami_columns)}个")
    print(f"识别到CON组样本列: {len(con_columns)}个")

    # 强制样本列为数值类型，无法转换的设为0（而非NaN，避免后续删除有效数据）
    for col in ami_columns + con_columns:
        # 先尝试转换为数值，错误值设为0
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        # 确保是float类型，避免后续计算报错
        df[col] = df[col].astype(float)

    # 检查是否有样本列（防止识别失败）
    if len(ami_columns) == 0 or len(con_columns) == 0:
        raise ValueError("未识别到AMI或CON样本列，请检查列名格式")

    return df, ami_columns, con_columns


# 2. 计算p_value（核心修复：数据类型校验+宽松样本量+异常捕获）
def calculate_p_value(row, ami_cols, con_cols):
    """
    计算p值：
    1. 强制数据类型为float
    2. 样本量门槛：每组至少1个检出样本（非零），总样本量≥3
    3. 对数转换前加极小值，避免log10(0)错误
    """
    # 提取样本值并强制为float（解决类型异常）
    ami_values = row[ami_cols].astype(float).values
    con_values = row[con_cols].astype(float).values

    # 统计"检出样本数"（非零值，代表该微生物在样本中存在）
    ami_detected = np.sum(ami_values > 0)  # AMI组检出样本数
    con_detected = np.sum(con_values > 0)  # CON组检出样本数
    total_detected = ami_detected + con_detected  # 总检出样本数

    # 样本量门槛：适配微生物稀疏数据（每组至少1个检出，总检出≥3）
    if ami_detected < 1 or con_detected < 1 or total_detected < 3:
        return np.nan

    try:
        # 对数转换（加极小值避免log10(0)，同时减少极端值影响）
        min_val = 1e-10  # 不影响原始数据量级的极小值
        ami_log = np.log10(ami_values + min_val)
        con_log = np.log10(con_values + min_val)

        # Mann-Whitney U检验（非参数，无需正态分布，适配微生物数据）
        stat, p_value = mannwhitneyu(ami_log, con_log, alternative='two-sided')
        return p_value

    except Exception as e:
        # 详细异常信息，便于排查个别异常行
        asv_id = row.get('Unnamed: 0', '未知ASV')
        print(f"ASV {asv_id} 计算p值失败: {str(e)[:50]}")  # 只打印前50字符，避免日志过长
        return np.nan


# 3. 计算LDA_score（修复特征维度问题）
def calculate_lda_score(df, ami_cols, con_cols):
    """计算LDA分数，过滤近零方差特征，避免计算错误"""
    try:
        # 准备LDA输入：样本为行，ASV为列（需转置）
        X = df[ami_cols + con_cols].values.T  # 形状：(样本数, ASV数)
        y = np.array([1] * len(ami_cols) + [0] * len(con_cols))  # 标签：1=AMI，0=CON

        # 过滤近零方差特征（避免LDA无法收敛）
        var_threshold = 1e-8  # 极低方差阈值，保留多数有效特征
        variances = np.var(X, axis=0)  # 计算每个ASV的方差
        valid_features = variances > var_threshold  # 筛选有效ASV

        if np.sum(valid_features) == 0:
            print("警告：所有ASV方差接近零，LDA无法计算，返回NaN")
            return np.full(len(df), np.nan)

        # 仅用有效特征训练LDA
        X_filtered = X[:, valid_features]
        lda = LinearDiscriminantAnalysis()
        lda.fit(X_filtered, y)

        # 为所有ASV分配LDA分数（无效特征分数设为0）
        lda_scores = np.zeros(len(df))
        lda_scores[valid_features] = np.abs(lda.coef_[0])  # 取系数绝对值作为LDA分数
        return lda_scores

    except Exception as e:
        print(f"LDA计算警告: {str(e)}")
        return np.full(len(df), np.nan)


# 主流程
if __name__ == "__main__":
    # 读取并校验数据（替换为你的数据路径）
    df, ami_columns, con_columns = load_and_validate_data('merged_microbes_data.csv')

    # 初始化新列存储结果
    new_columns = pd.DataFrame(index=df.index)

    # 计算平均丰度
    new_columns['AMI_avg_abundance'] = df[ami_columns].mean(axis=1)
    new_columns['CON_avg_abundance'] = df[con_columns].mean(axis=1)

    # 处理平均丰度的零值（避免log2错误）
    min_value = 1e-10
    new_columns['AMI_avg_abundance'] = new_columns['AMI_avg_abundance'].clip(lower=min_value)
    new_columns['CON_avg_abundance'] = new_columns['CON_avg_abundance'].clip(lower=min_value)

    # 计算log2FC
    new_columns['log2FC'] = np.log2(new_columns['AMI_avg_abundance'] / new_columns['CON_avg_abundance'])
    new_columns['log2FC'] = new_columns['log2FC'].replace([np.inf, -np.inf], np.nan)

    # 计算p_value（关键修复）
    print("\n开始计算p值...")
    new_columns['p_value'] = df.apply(
        calculate_p_value,
        args=(ami_columns, con_columns),  # 传递样本列名
        axis=1
    )

    # 计算LDA_score
    print("\n开始计算LDA分数...")
    new_columns['LDA_score'] = calculate_lda_score(df, ami_columns, con_columns)

    # 合并结果并保存
    df_result = pd.concat([df, new_columns], axis=1)
    output_path = 'merged_microbes_with_metrics_fixed.csv'
    df_result.to_csv(output_path, index=False)

    # 输出结果统计（验证修复效果）
    valid_p = new_columns['p_value'].notna().sum()
    total_p = len(new_columns)
    print(f"\n处理完成！结果保存至: {output_path}")
    print(f"p值统计：有效p值 {valid_p}/{total_p} ({valid_p / total_p * 100:.2f}%)")
    if valid_p > 0:
        print(f"p值范围：{new_columns['p_value'].min():.6f} ~ {new_columns['p_value'].max():.6f}")
    else:
        print("警告：仍无有效p值，建议检查原始数据的样本量和数值分布")