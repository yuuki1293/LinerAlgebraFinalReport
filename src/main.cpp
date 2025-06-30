#include "../include/linear_algebra.hpp"
#include <print>
#include <vector>

int main() {
    std::print("=== 応用線形代数 最終レポート ===\n");
    std::print("C++26 std::printを使用した実装\n\n");

    // テスト行列の定義
    std::vector<std::vector<double>> matrix1 = {
        {2.0, 1.0, 0.0},
        {1.0, 3.0, 1.0},
        {0.0, 1.0, 2.0}
    };

    std::vector<std::vector<double>> matrix2 = {
        {4.0, 1.0, 0.0},
        {1.0, 4.0, 1.0},
        {0.0, 1.0, 4.0}
    };

    // 対角化可能な行列（異なる固有値を持つ）
    std::vector<std::vector<double>> diagonalizable_matrix = {
        {3.0, 1.0, 0.0},
        {0.0, 2.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    // 対角化できない行列（ジョルダン標準形）
    std::vector<std::vector<double>> non_diagonalizable_matrix = {
        {2.0, 1.0, 0.0},
        {0.0, 2.0, 1.0},
        {0.0, 0.0, 2.0}
    };

    // 1. 基本行列操作のテスト
    std::print("1. 基本行列操作のテスト\n");
    std::print("------------------------\n");

    LinearAlgebra::MatrixOperations::printMatrix(matrix1, "行列A");

    // 行列式の計算
    double det1 = LinearAlgebra::MatrixOperations::determinant(matrix1);
    std::print("行列式: {}\n", det1);

    // 条件数
    double cond = LinearAlgebra::MatrixOperations::conditionNumber(matrix1);
    std::print("条件数: {}\n", cond);

    // ランク
    int rank = LinearAlgebra::MatrixOperations::rank(matrix1);
    std::print("ランク: {}\n", rank);

    // 転置
    auto transposed = LinearAlgebra::MatrixOperations::transpose(matrix1);
    LinearAlgebra::MatrixOperations::printMatrix(transposed, "転置行列");

    // 2. 線形方程式の解法テスト
    std::print("\n2. 線形方程式の解法テスト\n");
    std::print("------------------------\n");

    std::vector<double> b = {1.0, 2.0, 3.0};
    std::print("右辺ベクトル b:\n");
    for (size_t i = 0; i < b.size(); i++) {
        std::print("b[{}] = {}\n", i, b[i]);
    }
    std::print("\n");

    // LU分解による解法
    auto x_lu = LinearAlgebra::LinearSolver::solveLU(matrix1, b);
    std::print("LU分解による解:\n");
    for (size_t i = 0; i < x_lu.size(); i++) {
        std::print("x[{}] = {}\n", i, x_lu[i]);
    }
    std::print("\n");

    // 3. 固有値・固有ベクトル計算テスト
    std::print("\n3. 固有値・固有ベクトル計算テスト\n");
    std::print("--------------------------------\n");

    // 対角化可能な行列のテスト
    std::print("対角化可能な行列のテスト:\n");
    LinearAlgebra::MatrixOperations::printMatrix(diagonalizable_matrix, "対角化可能な行列");
    LinearAlgebra::EigenvalueAnalysis::eigenvalueDecomposition(diagonalizable_matrix);

    // 対角化できない行列のテスト
    std::print("対角化できない行列のテスト:\n");
    LinearAlgebra::MatrixOperations::printMatrix(non_diagonalizable_matrix, "対角化できない行列");
    LinearAlgebra::EigenvalueAnalysis::eigenvalueDecomposition(non_diagonalizable_matrix);

    // 4. 誤差解析テスト
    std::print("\n4. 誤差解析テスト\n");
    std::print("----------------\n");

    // 真の解（簡単な例）
    std::vector<double> x_exact = {0.5, 0.5, 1.0};
    std::print("真の解:\n");
    for (size_t i = 0; i < x_exact.size(); i++) {
        std::print("x_exact[{}] = {}\n", i, x_exact[i]);
    }
    std::print("\n");

    LinearAlgebra::NumericalAnalysis::errorAnalysis(matrix1, b, x_exact, x_lu);

    // 5. 行列の保存・読み込みテスト
    std::print("\n5. 行列の保存・読み込みテスト\n");
    std::print("----------------------------\n");

    LinearAlgebra::MatrixOperations::saveMatrix(matrix1, "data/test_matrix.txt");
    std::print("行列を data/test_matrix.txt に保存しました。\n");

    auto loaded_matrix = LinearAlgebra::MatrixOperations::loadMatrix("data/test_matrix.txt");
    LinearAlgebra::MatrixOperations::printMatrix(loaded_matrix, "読み込んだ行列");

    // 等価性チェック
    bool isEqual = LinearAlgebra::MatrixOperations::isEqual(matrix1, loaded_matrix);
    std::print("保存・読み込みの等価性: {}\n", isEqual ? "成功" : "失敗");

    // 6. ランダム行列テスト
    std::print("\n6. ランダム行列テスト\n");
    std::print("--------------------\n");

    // 小規模テスト（n=1~10, 各サイズ1回）
    std::print("小規模ランダム行列テスト (n=1~10) を実行します...\n");
    LinearAlgebra::RandomMatrixAnalysis::runRandomMatrixTest(10, 1);

    // 中規模テスト（n=1~50, 各サイズ1回）
    std::print("\n中規模ランダム行列テスト (n=1~50) を実行します...\n");
    LinearAlgebra::RandomMatrixAnalysis::runRandomMatrixTest(50, 1);

    // 大規模テスト（n=1~100, 各サイズ1回）- 時間がかかる場合があるため注意
    std::print("\n大規模ランダム行列テスト (n=1~100) を実行します...\n");
    std::print("注意: このテストは時間がかかる場合があります。\n");
    LinearAlgebra::RandomMatrixAnalysis::runRandomMatrixTest(100, 1);

    std::print("\n=== テスト完了 ===\n");

    return 0;
}
