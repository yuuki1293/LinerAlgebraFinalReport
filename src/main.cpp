// QR法の反復回数設定（必要に応じて変更可能）
#define QR_MAX_ITERATIONS 50000
#define QR_TOLERANCE 1e-10

#include "../include/linear_algebra.hpp"
#include <iostream>
#include <vector>
#include <iomanip>

using namespace LinearAlgebra;

// 行列式計算テスト関数
void runDeterminantTest() {
    std::cout << "1. 行列式計算テスト" << std::endl;
    std::cout << "------------------" << std::endl;

    std::vector<std::vector<double>> test1 = {
        {1, 0, 2},
        {3, 4, 5},
        {5, 6, 7}
    };
    MatrixOperations::printMatrix(test1, "テスト行列1");
    double det_test1 = MatrixOperations::determinant(test1);
    std::cout << "行列式: " << det_test1 << std::endl;
    MatrixOperations::saveDeterminantToCSV(test1, 3, det_test1);

    std::vector<std::vector<double>> test2 = {
        {1, 0, 0},
        {2, 3, 5},
        {4, 1, 3}
    };
    MatrixOperations::printMatrix(test2, "テスト行列2");
    double det_test2 = MatrixOperations::determinant(test2);
    std::cout << "行列式: " << det_test2 << std::endl;
    MatrixOperations::saveDeterminantToCSV(test2, 3, det_test2);
}

// 連立方程式解法テスト関数
void runLinearSolverTest() {
    std::cout << std::endl << "2. 連立方程式解法テスト" << std::endl;
    std::cout << "------------------------" << std::endl;

    std::vector<std::vector<double>> A = {
        {2, 3, 1},
        {1, 2, 3},
        {3, 1, 2}
    };
    MatrixOperations::printMatrix(A, "係数行列A");

    std::vector<double> b = {9, 6, 8};
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
}

// 固有値・固有ベクトル計算テスト関数
void runEigenvalueTest() {
    std::cout << std::endl << "3. 固有値・固有ベクトル計算テスト" << std::endl;
    std::cout << "--------------------------------" << std::endl;

    std::vector<std::vector<double>> matrix = {
        {1, -2, 0},
        {-2, 2, -2},
        {0, -2, 3}
    };
    MatrixOperations::printMatrix(matrix, "テスト行列");
    auto [eigenvalues, eigenvectors] = EigenvalueAnalysis::qrEigenDecomposition(matrix);
    EigenvalueAnalysis::printEigenvalues(eigenvalues, "QR法による固有値");

    // 固有ベクトルの表示
    std::cout << "固有ベクトル:" << std::endl;
    for (size_t i = 0; i < eigenvectors.size(); i++) {
        std::cout << "v[" << i << "] = [";
        for (size_t j = 0; j < eigenvectors[i].size(); j++) {
            std::cout << std::fixed << std::setprecision(4) << eigenvectors[i][j];
            if (j < eigenvectors[i].size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
}

// ランダム行列テスト関数
void runRandomMatrixTest() {
    std::cout << std::endl << "4. ランダム行列テスト" << std::endl;
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
}

int main() {
    std::cout << "=== 応用線形代数 最終レポート ===" << std::endl;
    std::cout << "C++20 std::coutを使用した実装" << std::endl << std::endl;

    // 各テスト関数を実行
    runDeterminantTest();
    runLinearSolverTest();
    runEigenvalueTest();
    runRandomMatrixTest();

    std::cout << std::endl << "=== テスト完了 ===" << std::endl;
    return 0;
}
