#include "../include/linear_algebra.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cout << "=== 応用線形代数 最終レポート ===" << std::endl;
    std::cout << "C++20 std::coutを使用した実装" << std::endl << std::endl;

    using namespace LinearAlgebra;

    // テスト用の行列
    std::vector<std::vector<double>> A = {
        {2, 1, 0},
        {1, 3, 1},
        {0, 1, 2}
    };

    // 1. 基本行列操作のテスト
    std::cout << "1. 基本行列操作のテスト" << std::endl;
    std::cout << "------------------------" << std::endl;
    MatrixOperations::printMatrix(A, "行列A");

    double det1 = MatrixOperations::determinant(A);
    std::cout << "行列式: " << det1 << std::endl;

    double cond = MatrixOperations::conditionNumber(A);
    std::cout << "条件数: " << cond << std::endl;

    int rank = MatrixOperations::rank(A);
    std::cout << "ランク: " << rank << std::endl;

    auto A_transpose = MatrixOperations::transpose(A);
    MatrixOperations::printMatrix(A_transpose, "転置行列");

    // 2. 線形方程式の解法テスト
    std::cout << std::endl << "2. 線形方程式の解法テスト" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::vector<double> b = {1, 2, 3};
    std::cout << "右辺ベクトル b:" << std::endl;
    for (size_t i = 0; i < b.size(); i++) {
        std::cout << "b[" << i << "] = " << b[i] << std::endl;
    }
    std::cout << std::endl;

    auto x_lu = LinearSolver::solveLU(A, b);
    std::cout << "LU分解による解:" << std::endl;
    for (size_t i = 0; i < x_lu.size(); i++) {
        std::cout << "x[" << i << "] = " << x_lu[i] << std::endl;
    }
    std::cout << std::endl;

    // 3. 固有値・固有ベクトル計算テスト
    std::cout << std::endl << "3. 固有値・固有ベクトル計算テスト" << std::endl;
    std::cout << "--------------------------------" << std::endl;

    std::cout << "対角化可能な行列のテスト:" << std::endl;
    std::vector<std::vector<double>> diagonalizable = {
        {3, 1, 0},
        {0, 2, 0},
        {0, 0, 1}
    };
    MatrixOperations::printMatrix(diagonalizable, "対角化可能な行列");
    EigenvalueAnalysis::eigenvalueDecomposition(diagonalizable);

    std::cout << "対角化できない行列のテスト:" << std::endl;
    std::vector<std::vector<double>> nonDiagonalizable = {
        {2, 1, 0},
        {0, 2, 1},
        {0, 0, 2}
    };
    MatrixOperations::printMatrix(nonDiagonalizable, "対角化できない行列");
    EigenvalueAnalysis::eigenvalueDecomposition(nonDiagonalizable);

    // 4. 誤差解析テスト
    std::cout << std::endl << "4. 誤差解析テスト" << std::endl;
    std::cout << "----------------" << std::endl;
    std::vector<double> x_exact = {0.5, 0.5, 1.0};
    std::vector<double> x_computed = {1.0, 0.0, 0.5};

    std::cout << "真の解:" << std::endl;
    for (size_t i = 0; i < x_exact.size(); i++) {
        std::cout << "x_exact[" << i << "] = " << x_exact[i] << std::endl;
    }
    std::cout << std::endl;

    NumericalAnalysis::errorAnalysis(A, b, x_exact, x_computed);

    // 5. 行列の保存・読み込みテスト
    std::cout << std::endl << "5. 行列の保存・読み込みテスト" << std::endl;
    std::cout << "----------------------------" << std::endl;
    MatrixOperations::saveMatrix(A, "data/test_matrix.txt");
    std::cout << "行列を data/test_matrix.txt に保存しました。" << std::endl;

    auto loaded_matrix = MatrixOperations::loadMatrix("data/test_matrix.txt");
    MatrixOperations::printMatrix(loaded_matrix, "読み込んだ行列");

    bool isEqual = MatrixOperations::isEqual(A, loaded_matrix);
    std::cout << "保存・読み込みの等価性: " << (isEqual ? "成功" : "失敗") << std::endl;

    // 6. ランダム行列テスト
    std::cout << std::endl << "6. ランダム行列テスト" << std::endl;
    std::cout << "--------------------" << std::endl;

#ifdef DEBUG_MODE
    int maxSize = 10;
    std::cout << "デバッグモード: ランダム行列テスト (n=1~" << maxSize << ") を実行します..." << std::endl;
#else
    int maxSize = 100;
    std::cout << "ランダム行列テスト (n=1~" << maxSize << ") を実行します..." << std::endl;
#endif
    std::cout << "注意: このテストは時間がかかる場合があります。" << std::endl;
    RandomMatrixAnalysis::runRandomMatrixTest(maxSize, 1);

    std::cout << std::endl << "=== テスト完了 ===" << std::endl;
    return 0;
}
