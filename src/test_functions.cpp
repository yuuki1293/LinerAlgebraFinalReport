#include "../include/linear_algebra.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

using namespace LinearAlgebra;

// 行列式計算のテスト
void testDeterminant() {
    std::cout << "=== 行列式計算テスト ===" << std::endl;

    // テストケース1: 2x2行列
    std::vector<std::vector<double>> A2 = {
        {2, 1},
        {1, 3}
    };
    double det2 = MatrixOperations::determinant(A2);
    double expected2 = 2*3 - 1*1; // 5
    std::cout << "2x2行列の行列式: " << det2 << " (期待値: " << expected2 << ")" << std::endl;
    assert(std::abs(det2 - expected2) < 1e-10);

    // テストケース2: 3x3行列
    std::vector<std::vector<double>> A3 = {
        {1, 0, 2},
        {3, 4, 5},
        {5, 6, 7}
    };
    double det3 = MatrixOperations::determinant(A3);
    double expected3 = -6; // 手計算で確認済み
    std::cout << "3x3行列の行列式: " << det3 << " (期待値: " << expected3 << ")" << std::endl;
    assert(std::abs(det3 - expected3) < 1e-10);

    // テストケース3: 単位行列
    auto I = MatrixOperations::identity(3);
    double detI = MatrixOperations::determinant(I);
    std::cout << "単位行列の行列式: " << detI << " (期待値: 1)" << std::endl;
    assert(std::abs(detI - 1.0) < 1e-10);

    std::cout << "行列式計算テスト: 成功" << std::endl << std::endl;
}

// 連立方程式解法のテスト
void testLinearSolver() {
    std::cout << "=== 連立方程式解法テスト ===" << std::endl;

    // テストケース1: 簡単な2x2システム
    std::vector<std::vector<double>> A = {
        {2, 1},
        {1, 3}
    };
    std::vector<double> b = {5, 6};
    auto x = LinearSolver::solveLU(A, b);

    // 解の検証: Ax = b
    std::vector<std::vector<double>> x_col(x.size(), std::vector<double>(1));
    for (size_t i = 0; i < x.size(); ++i) x_col[i][0] = x[i];
    auto Ax = MatrixOperations::matrixMultiply(A, x_col);
    std::cout << "2x2システム:" << std::endl;
    std::cout << "解: [" << x[0] << ", " << x[1] << "]" << std::endl;
    std::cout << "Ax: [" << Ax[0][0] << ", " << Ax[1][0] << "]" << std::endl;
    std::cout << "b: [" << b[0] << ", " << b[1] << "]" << std::endl;

    assert(std::abs(Ax[0][0] - b[0]) < 1e-10);
    assert(std::abs(Ax[1][0] - b[1]) < 1e-10);

    // テストケース2: 3x3システム
    std::vector<std::vector<double>> A3 = {
        {2, 3, 1},
        {1, 2, 3},
        {3, 1, 2}
    };
    std::vector<double> b3 = {9, 6, 8};
        auto x3 = LinearSolver::solveLU(A3, b3);

    std::vector<std::vector<double>> x3_col(x3.size(), std::vector<double>(1));
    for (size_t i = 0; i < x3.size(); ++i) x3_col[i][0] = x3[i];
    auto Ax3 = MatrixOperations::matrixMultiply(A3, x3_col);
    std::cout << "3x3システム:" << std::endl;
    std::cout << "解: [" << x3[0] << ", " << x3[1] << ", " << x3[2] << "]" << std::endl;

    for (size_t i = 0; i < b3.size(); i++) {
        assert(std::abs(Ax3[i][0] - b3[i]) < 1e-10);
    }

    std::cout << "連立方程式解法テスト: 成功" << std::endl << std::endl;
}

// QR分解のテスト
void testQRDecomposition() {
    std::cout << "=== QR分解テスト ===" << std::endl;

    // テストケース: 3x3行列
    std::vector<std::vector<double>> A = {
        {1, 2, 3},
        {0, 1, -3},
        {0, -3, 1}
    };

    auto [Q, R] = EigenvalueAnalysis::qrDecomposition(A);

    std::cout << "元の行列A:" << std::endl;
    MatrixOperations::printMatrix(A, "A");

    std::cout << "Q行列:" << std::endl;
    MatrixOperations::printMatrix(Q, "Q");

    std::cout << "R行列:" << std::endl;
    MatrixOperations::printMatrix(R, "R");

    // 検証1: A = QR
    auto QR = MatrixOperations::matrixMultiply(Q, R);
    std::cout << "QR積:" << std::endl;
    MatrixOperations::printMatrix(QR, "QR");

    // 検証2: Q^T Q = I (Qの直交性)
    auto QT = MatrixOperations::transpose(Q);
    auto QTQ = MatrixOperations::matrixMultiply(QT, Q);
    std::cout << "Q^T Q (単位行列になるはず):" << std::endl;
    MatrixOperations::printMatrix(QTQ, "Q^T Q");

    // 検証3: Rが上三角行列
    bool isUpperTriangular = true;
    for (size_t i = 0; i < R.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            if (std::abs(R[i][j]) > 1e-10) {
                isUpperTriangular = false;
                break;
            }
        }
    }

    std::cout << "Rが上三角行列: " << (isUpperTriangular ? "Yes" : "No") << std::endl;
    assert(isUpperTriangular);

    // A = QRの検証
    bool AequalsQR = true;
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++) {
            if (std::abs(A[i][j] - QR[i][j]) > 1e-10) {
                AequalsQR = false;
                break;
            }
        }
    }

    std::cout << "A = QR: " << (AequalsQR ? "Yes" : "No") << std::endl;
    assert(AequalsQR);

    std::cout << "QR分解テスト: 成功" << std::endl << std::endl;
}

// 固有値・固有ベクトル計算のテスト
void testEigenvalueDecomposition() {
    std::cout << "=== 固有値・固有ベクトル計算テスト ===" << std::endl;

    // テストケース1: 対角行列（固有値が明確）
    std::vector<std::vector<double>> diagonal = {
        {3, 0, 0},
        {0, 2, 0},
        {0, 0, 1}
    };

    auto [eigenvalues1, eigenvectors1] = EigenvalueAnalysis::qrEigenDecomposition(diagonal, 100, 1e-8);

    std::cout << "対角行列の固有値:" << std::endl;
    for (size_t i = 0; i < eigenvalues1.size(); i++) {
        std::cout << "λ[" << i << "] = " << eigenvalues1[i].real() << std::endl;
    }

    // 期待値: 3, 2, 1
    std::vector<double> expected1 = {3.0, 2.0, 1.0};
    for (size_t i = 0; i < eigenvalues1.size(); i++) {
        assert(std::abs(eigenvalues1[i].real() - expected1[i]) < 1e-6);
    }

    // テストケース2: 元のテスト行列
    std::vector<std::vector<double>> matrix = {
        {1, 2, 3},
        {0, 1, -3},
        {0, -3, 1}
    };

    auto [eigenvalues2, eigenvectors2] = EigenvalueAnalysis::qrEigenDecomposition(matrix, 100, 1e-8);

    std::cout << "テスト行列の固有値:" << std::endl;
    for (size_t i = 0; i < eigenvalues2.size(); i++) {
        std::cout << "λ[" << i << "] = " << eigenvalues2[i].real();
        if (std::abs(eigenvalues2[i].imag()) > 1e-10) {
            std::cout << " + " << eigenvalues2[i].imag() << "i";
        }
        std::cout << std::endl;
    }

    // 固有値がNaNでないことを確認
    for (size_t i = 0; i < eigenvalues2.size(); i++) {
        assert(!std::isnan(eigenvalues2[i].real()));
        assert(!std::isnan(eigenvalues2[i].imag()));
    }

        // 固有ベクトルの検証: Av = λv
    std::cout << "固有ベクトルの検証:" << std::endl;
    for (size_t i = 0; i < eigenvectors2.size(); i++) {
        std::vector<double> v = eigenvectors2[i];
        std::vector<std::vector<double>> v_col(v.size(), std::vector<double>(1));
        for (size_t j = 0; j < v.size(); ++j) v_col[j][0] = v[j];
        auto Av = MatrixOperations::matrixMultiply(matrix, v_col);
        auto lambda_v = std::vector<double>{v[0] * eigenvalues2[i].real(),
                                          v[1] * eigenvalues2[i].real(),
                                          v[2] * eigenvalues2[i].real()};

        std::cout << "固有値 " << eigenvalues2[i].real() << " の固有ベクトル:" << std::endl;
        std::cout << "v = [" << v[0] << ", " << v[1] << ", " << v[2] << "]" << std::endl;
        std::cout << "Av = [" << Av[0][0] << ", " << Av[1][0] << ", " << Av[2][0] << "]" << std::endl;
        std::cout << "λv = [" << lambda_v[0] << ", " << lambda_v[1] << ", " << lambda_v[2] << "]" << std::endl;

        // Av ≈ λv の検証
        for (size_t j = 0; j < v.size(); j++) {
            assert(std::abs(Av[j][0] - lambda_v[j]) < 1e-6);
        }
        std::cout << std::endl;
    }

    std::cout << "固有値・固有ベクトル計算テスト: 成功" << std::endl << std::endl;
}

// 行列演算のテスト
void testMatrixOperations() {
    std::cout << "=== 行列演算テスト ===" << std::endl;

    // 行列の積
    std::vector<std::vector<double>> A = {
        {1, 2},
        {3, 4}
    };
    std::vector<std::vector<double>> B = {
        {2, 0},
        {1, 2}
    };

    auto AB = MatrixOperations::matrixMultiply(A, B);
    std::cout << "行列積 A*B:" << std::endl;
    MatrixOperations::printMatrix(AB, "AB");

    // 期待値: [[4, 4], [10, 8]]
    assert(std::abs(AB[0][0] - 4) < 1e-10);
    assert(std::abs(AB[0][1] - 4) < 1e-10);
    assert(std::abs(AB[1][0] - 10) < 1e-10);
    assert(std::abs(AB[1][1] - 8) < 1e-10);

    // 転置
    auto AT = MatrixOperations::transpose(A);
    std::cout << "Aの転置:" << std::endl;
    MatrixOperations::printMatrix(AT, "A^T");

    assert(AT[0][0] == A[0][0]);
    assert(AT[0][1] == A[1][0]);
    assert(AT[1][0] == A[0][1]);
    assert(AT[1][1] == A[1][1]);

    // 単位行列
    auto I = MatrixOperations::identity(3);
    std::cout << "3x3単位行列:" << std::endl;
    MatrixOperations::printMatrix(I, "I");

    for (size_t i = 0; i < I.size(); i++) {
        for (size_t j = 0; j < I[i].size(); j++) {
            if (i == j) {
                assert(std::abs(I[i][j] - 1.0) < 1e-10);
            } else {
                assert(std::abs(I[i][j]) < 1e-10);
            }
        }
    }

    std::cout << "行列演算テスト: 成功" << std::endl << std::endl;
}

int main() {
    std::cout << "=== 単体テスト開始 ===" << std::endl << std::endl;

    try {
        testMatrixOperations();
        testDeterminant();
        testLinearSolver();
        testQRDecomposition();
        testEigenvalueDecomposition();

        std::cout << "=== 全てのテストが成功しました ===" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "テスト失敗: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
