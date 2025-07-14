#!/usr/bin/env python3
"""
行列の条件数を計算するプログラム
data/Aフォルダにある行列ファイルを読み込み、各行列の条件数を計算する
"""

import numpy as np
import pandas as pd
import os
import glob
from pathlib import Path

def calculate_condition_number(matrix):
    """
    行列の条件数を計算する

    Args:
        matrix: numpy配列の行列

    Returns:
        condition_number: 条件数（最大特異値 / 最小特異値）
    """
    try:
        # 特異値分解を使用して条件数を計算
        singular_values = np.linalg.svd(matrix, compute_uv=False)
        condition_number = singular_values[0] / singular_values[-1]
        return condition_number
    except np.linalg.LinAlgError:
        # 行列が特異な場合
        return np.inf

def load_matrix_from_csv(file_path):
    """
    CSVファイルから行列を読み込む

    Args:
        file_path: CSVファイルのパス

    Returns:
        matrix: numpy配列の行列
    """
    try:
        # CSVファイルを読み込み
        data = pd.read_csv(file_path, header=None)
        matrix = data.values
        return matrix
    except Exception as e:
        print(f"エラー: {file_path} の読み込みに失敗しました: {e}")
        return None

def main():
    """
    メイン関数
    """
    # data/Aフォルダのパス
    data_dir = Path("data/A")

    # 結果を格納するリスト
    results = []

    # CSVファイルを取得（数字順にソート）
    csv_files = sorted(data_dir.glob("*.csv"), key=lambda x: int(x.stem))

    print("行列の条件数を計算中...")
    print("-" * 50)

    for csv_file in csv_files:
        matrix_size = int(csv_file.stem)

        # 行列を読み込み
        matrix = load_matrix_from_csv(csv_file)

        if matrix is not None:
            # 条件数を計算
            condition_number = calculate_condition_number(matrix)

            # 結果を保存
            results.append({
                'size': matrix_size,
                'condition_number': condition_number,
                'log10_condition': np.log10(condition_number) if condition_number != np.inf else np.inf
            })

            # 進捗を表示
            print(f"サイズ {matrix_size:3d}: 条件数 = {condition_number:.2e}")

    # 結果をDataFrameに変換
    df = pd.DataFrame(results)

    # 結果をCSVファイルに保存
    output_file = "condition_numbers.csv"
    df.to_csv(output_file, index=False)
    print(f"\n結果を {output_file} に保存しました")

    # 統計情報を表示
    print("\n統計情報:")
    print(f"計算した行列数: {len(df)}")
    print(f"条件数の最小値: {df['condition_number'].min():.2e}")
    print(f"条件数の最大値: {df['condition_number'].max():.2e}")
    print(f"条件数の平均値: {df['condition_number'].mean():.2e}")
    print(f"条件数の中央値: {df['condition_number'].median():.2e}")

    # 条件数が非常に大きい行列を特定
    large_condition = df[df['condition_number'] > 1e10]
    if not large_condition.empty:
        print(f"\n条件数が非常に大きい行列 (>1e10):")
        for _, row in large_condition.iterrows():
            print(f"  サイズ {row['size']:3d}: 条件数 = {row['condition_number']:.2e}")

    # 特異な行列を特定
    singular_matrices = df[df['condition_number'] == np.inf]
    if not singular_matrices.empty:
        print(f"\n特異な行列:")
        for _, row in singular_matrices.iterrows():
            print(f"  サイズ {row['size']:3d}")

    return df

if __name__ == "__main__":
    df = main()
