import pandas as pd
import os
from statsmodels.stats.multitest import multipletests


def process_data(input_file_path, output_file_path):
    try:
        # 读取 CSV 文件
        df = pd.read_csv(input_file_path)

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

        # 计算 fdr
        rejected, fdr_values, _, _ = multipletests(df['p_value'], method='fdr_bh')

        # 将 fdr 值添加到数据框中
        df['fdr'] = fdr_values

        # 将结果保存为 csv 文件
        df.to_csv(output_file_path)
        print(f'数据已成功处理并保存到 {output_file_path}')

        # 筛选 fdr 小于 0.05 的数据
        filtered_df = df[df['fdr'] < 0.05]
        filtered_output_file = os.path.join(os.path.dirname(output_file_path),
                                            'significant_microbes_filtered_fdr_less_005.csv')
        filtered_df.to_csv(filtered_output_file)
        print(f'fdr 小于 0.05 的数据已成功保存到 {filtered_output_file}')

    except FileNotFoundError:
        print(f'文件 {input_file_path} 未找到，请检查文件路径是否正确。')
    except KeyError as e:
        print(f'数据框中不存在列 {e}，请检查数据结构。')
    except Exception as e:
        print(f'发生未知错误：{e}')


if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir,'significant_microbes_filtered.csv')
    output_file = os.path.join(current_dir,'significant_microbes_filtered_fdr.csv')

    process_data(input_file, output_file)