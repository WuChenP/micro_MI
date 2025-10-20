import pandas as pd
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.stats import kruskal, mannwhitneyu

def preprocess_abundance(abundance_df):
    """log10+1 转换 + 样本内归一化"""
    df_t = abundance_df.T.copy()
    # log 转换避免极差
    df_t = np.log10(df_t + 1)
    # 样本内归一化
    df_t = df_t.div(df_t.sum(axis=1), axis=0)
    return df_t.T

def run_simple_lefse_with_results(abundance_df, metadata_df, class_col,
                                  p_threshold=0.05, lda_threshold=2.0, lda_scale=1.0):
    """
    在 abundance_df 上增加 KW_p、MW_p、LDA_score 列，并返回筛选结果
    """
    abundance_df_proc = preprocess_abundance(abundance_df)

    # 合并分组信息
    df = abundance_df_proc.T.copy()
    df = df.reset_index().rename(columns={'index': 'SampleID'})
    merged = pd.merge(metadata_df, df, on='SampleID', how='inner')

    features = abundance_df_proc.index.tolist()
    groups = merged[class_col].unique().tolist()

    results_df = pd.DataFrame(index=features, columns=['KW_p','MW_p','LDA_score'], dtype=float)

    for feat in features:
        # 取各组的值
        group_vals = [merged.loc[merged[class_col] == g, feat].dropna().values for g in groups]
        # KW检验
        try:
            kw_stat, kw_p = kruskal(*group_vals)
        except Exception:
            kw_p = np.nan
        # MW检验（两组时）
        if len(groups) == 2:
            try:
                mw_stat, mw_p = mannwhitneyu(group_vals[0], group_vals[1], alternative='two-sided')
            except Exception:
                mw_p = np.nan
        else:
            mw_p = np.nan

        # LDA 每个特征单独做判别
        X = merged[[feat]].values
        y = merged[class_col].astype('category').cat.codes.values
        try:
            lda = LinearDiscriminantAnalysis(n_components=1)
            lda.fit(X, y)
            coef = abs(lda.coef_[0][0])
            # 更接近官方LEfSe：取log10再乘放大倍数
            lda_score = np.log10(coef + 1e-6) * lda_scale
        except Exception:
            lda_score = np.nan

        results_df.loc[feat] = [kw_p, mw_p, lda_score]

    # 合并原始丰度与统计指标
    full_df = pd.concat([abundance_df, results_df], axis=1)

    # 筛选显著差异
    filtered = full_df[(full_df['KW_p'] < p_threshold) &
                       (full_df['LDA_score'] > lda_threshold)]

    return full_df, filtered

def prepare_abundance_metadata(excel_file):
    """读取原始 Excel，返回 abundance_df 和 metadata_df"""
    df = pd.read_excel(excel_file)
    ami_cols = [c for c in df.columns if c.startswith('AMI')]
    con_cols = [c for c in df.columns if c.startswith('CON')]
    sample_cols = ami_cols + con_cols

    abundance_df = df[sample_cols].copy()
    abundance_df.index = df['Compounds']

    groups = ['AMI' if s.startswith('AMI') else 'CON' for s in sample_cols]
    metadata_df = pd.DataFrame({
        'SampleID': sample_cols,
        'Group': groups
    })
    return abundance_df, metadata_df

if __name__ == '__main__':
    excel_file = 'E:/代谢组学/AMI_vs_CON_lev1.xlsx'  # 改成你的文件名
    abundance_df, metadata_df = prepare_abundance_metadata(excel_file)

    full_df, lefse_res = run_simple_lefse_with_results(
        abundance_df, metadata_df,
        class_col='Group',
        p_threshold=0.01,
        lda_threshold=3.0,
        lda_scale=1.0  # 放大10倍可调整
    )

    print("完整结果：")
    print(full_df.head())
    print("\n显著差异代谢物：")
    print(lefse_res)

    full_df.to_excel('abundance_with_results_lefse.xlsx')
    lefse_res.to_excel('E:/代谢组学/AMI_vs_CON_LEfSe.xlsx')
    print('结果已保存到 abundance_with_results.xlsx 和 lefse_result.xlsx')
