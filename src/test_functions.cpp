#include "../include/linear_algebra.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace LinearAlgebra;

// 行列式計算のテスト
void testDeterminant() {
    // テストケース2: 3x3行列
    std::vector<std::vector<double>> A = {
        {1, 0, 2},
        {3, 4, 5},
        {5, 6, 7}
    };
    double det = MatrixOperations::determinant(A);
    double expected = -6; // 手計算で確認済み
    std::cout << "行列式: " << det << " (真値: " << expected << ")" << std::endl;
    assert(std::abs(det - expected) < 1e-10);

    std::cout << "行列式計算テスト: 成功" << std::endl << std::endl;
}

// 連立方程式解法のテスト
void testLinearSolver() {
    std::cout << "=== 連立方程式解法テスト ===" << std::endl;
    // テストケース2: 3x3システム
    std::vector<std::vector<double>> A3 = {
        {2, -3, 5},
        {1, 1, -1},
        {-3, -6, 2}
    };
    std::vector<double> b3 = {-3, 0, -7};
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

// 固有値・固有ベクトル計算のテスト
void testEigenvalueDecomposition() {
    std::cout << "=== 固有値・固有ベクトル計算テスト ===" << std::endl;

    // テストケース1: 対角行列（固有値が明確）
    std::vector<std::vector<double>> diagonal = {
        {2, -2},
        {-2, -1}
    };

    auto [eigenvalues1, eigenvectors1] = EigenvalueAnalysis::qrEigenDecomposition(diagonal);

    std::cout << "対角行列の固有値:" << std::endl;
    for (size_t i = 0; i < eigenvalues1.size(); i++) {
        std::cout << "λ[" << i << "] = " << eigenvalues1[i].real() << std::endl;
    }
    std::cout << "対角行列の固有ベクトル:" << std::endl;
    for (size_t i = 0; i < eigenvectors1.size(); i++) {
        std::cout << "v[" << i << "] = [" << eigenvectors1[i][0] << ", " << eigenvectors1[i][1] << "]" << std::endl;
    }

    // 期待値: (3 ± √17)/2 ≈ 3.56155, -0.561553
    std::vector<double> expected1 = {3, -2};
    for (size_t i = 0; i < eigenvalues1.size(); i++) {
        assert(std::abs(eigenvalues1[i].real() - expected1[i]) < 1e-5);
    }

    std::cout << "固有値・固有ベクトル計算テスト: 成功" << std::endl << std::endl;
}

int main() {
    std::cout << "=== テスト開始 ===" << std::endl << std::endl;

    try {
        testDeterminant();
        testLinearSolver();
        testEigenvalueDecomposition();

        std::cout << "=== 全てのテストが成功しました ===" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "テスト失敗: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
