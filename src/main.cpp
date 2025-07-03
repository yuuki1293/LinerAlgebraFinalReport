// QR法の反復回数設定（必要に応じて変更可能）
#define QR_MAX_ITERATIONS 1000
#define QR_TOLERANCE 1e-10

#include "../include/linear_algebra.hpp"
#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>

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
    MatrixOperations::saveDeterminantToCSV(3, det_test1);

    std::vector<std::vector<double>> test2 = {
        {1, 0, 0},
        {2, 3, 5},
        {4, 1, 3}
    };
    MatrixOperations::printMatrix(test2, "テスト行列2");
    double det_test2 = MatrixOperations::determinant(test2);
    std::cout << "行列式: " << det_test2 << std::endl;
    MatrixOperations::saveDeterminantToCSV(3, det_test2);
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
    MatrixOperations::printMatrix(matrix, "テスト行列（対称行列）");

    // 通常のQR法による計算
    auto [eigenvalues, eigenvectors] = EigenvalueAnalysis::qrEigenDecomposition(matrix);
    EigenvalueAnalysis::printEigenvalues(eigenvalues, "QR法による固有値");

    // 固有ベクトルの表示
    std::cout << "QR法による固有ベクトル:" << std::endl;
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

// パフォーマンス比較テスト関数
void runPerformanceComparisonTest() {
    std::cout << std::endl << "5. パフォーマンス比較テスト" << std::endl;
    std::cout << "------------------------" << std::endl;

    // テスト用の対称行列を生成
    std::vector<std::vector<double>> matrix = {
        {4, 1, 0, 0},
        {1, 3, 1, 0},
        {0, 1, 2, 1},
        {0, 0, 1, 1}
    };

    std::cout << "4x4対称行列での比較:" << std::endl;
    MatrixOperations::printMatrix(matrix, "テスト行列");

    // 理論値の概算（対称行列のため実数固有値）
    std::cout << "理論固有値の概算: λ ≈ 4.5, 3.0, 1.5, 0.0" << std::endl;
    std::cout << std::endl;

    // 通常QR法の時間測定
    auto start = std::chrono::high_resolution_clock::now();
    auto [qr_eigenvalues, qr_eigenvectors] = EigenvalueAnalysis::qrEigenDecomposition(matrix);
    auto end = std::chrono::high_resolution_clock::now();
    auto qr_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // シフト付きQR法の時間測定
    start = std::chrono::high_resolution_clock::now();
    auto [shifted_eigenvalues, shifted_eigenvectors] = EigenvalueAnalysis::shiftedQREigenDecomposition(matrix);
    end = std::chrono::high_resolution_clock::now();
    auto shifted_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "通常QR法: " << qr_time.count() << " μs" << std::endl;
    std::cout << "シフト付きQR法: " << shifted_time.count() << " μs" << std::endl;
    std::cout << "高速化率: " << std::fixed << std::setprecision(1)
              << (double)qr_time.count() / shifted_time.count() << "倍" << std::endl;

    std::cout << std::endl << "結果比較:" << std::endl;
    std::cout << "通常QR法 - 固有値: ";
    for (const auto& val : qr_eigenvalues) {
        std::cout << std::fixed << std::setprecision(4) << val.real() << " ";
    }
    std::cout << std::endl;
    std::cout << "シフト付きQR法 - 固有値: ";
    for (const auto& val : shifted_eigenvalues) {
        std::cout << std::fixed << std::setprecision(4) << val.real() << " ";
    }
    std::cout << std::endl;
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
