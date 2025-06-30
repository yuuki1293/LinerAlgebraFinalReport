#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <random>
#include <chrono>

// Eigenライブラリのインクルード（インストールが必要）
#ifdef EIGEN_AVAILABLE
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif

namespace LinearAlgebra {

// 行列の基本操作
class MatrixOperations {
public:
    // 行列の表示
    static void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name = "Matrix");

    // 行列の保存
    static void saveMatrix(const std::vector<std::vector<double>>& matrix, const std::string& filename);

    // 行列の読み込み
    static std::vector<std::vector<double>> loadMatrix(const std::string& filename);

    // 行列の条件数計算
    static double conditionNumber(const std::vector<std::vector<double>>& matrix);

    // 行列のランク計算
    static int rank(const std::vector<std::vector<double>>& matrix);

    // 行列式の計算（LU分解によるアプローチ）
    static double determinant(const std::vector<std::vector<double>>& matrix);

    // 小行列の取得（行列式計算用）
    static std::vector<std::vector<double>> getMinor(const std::vector<std::vector<double>>& matrix, int row, int col);

    // 行列の転置
    static std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix);

    // 行列の積
    static std::vector<std::vector<double>> matrixMultiply(const std::vector<std::vector<double>>& A,
                                                           const std::vector<std::vector<double>>& B);

    // 単位行列の作成
    static std::vector<std::vector<double>> identity(int n);

    // 行列の等価性チェック
    static bool isEqual(const std::vector<std::vector<double>>& A,
                       const std::vector<std::vector<double>>& B, double tolerance = 1e-10);

    // 行列のコピー
    static std::vector<std::vector<double>> copyMatrix(const std::vector<std::vector<double>>& matrix);
};

// 固有値・固有ベクトル計算
class EigenvalueAnalysis {
public:
    // 固有値分解（対角化による）
    static void eigenvalueDecomposition(const std::vector<std::vector<double>>& matrix);

    // 対角化による固有値計算
    static std::vector<std::complex<double>> diagonalizationEigenvalues(const std::vector<std::vector<double>>& matrix);

    // 対角化可能性の判定
    static bool isDiagonalizable(const std::vector<std::vector<double>>& matrix);

    // 固有値の重複度をチェック
    static std::vector<int> checkEigenvalueMultiplicity(const std::vector<std::complex<double>>& eigenvalues,
                                                       double tolerance = 1e-10);

    // べき乗法による最大固有値と固有ベクトルの計算
    static std::pair<double, std::vector<double>> powerMethod(const std::vector<std::vector<double>>& matrix,
                                                             int maxIterations = 1000, double tolerance = 1e-6);

    // QR法による固有値・固有ベクトル計算
    // 戻り値: pair(固有値リスト, 固有ベクトル行列)
    static std::pair<std::vector<std::complex<double>>, std::vector<std::vector<double>>>
    qrEigenDecomposition(const std::vector<std::vector<double>>& matrix, int maxIterations = 1000, double tolerance = 1e-8);

    // QR分解
    static std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> qrDecomposition(
        const std::vector<std::vector<double>>& matrix);

    // グラム・シュミット直交化
    static std::vector<std::vector<double>> gramSchmidt(const std::vector<std::vector<double>>& vectors);

    // 固有値の表示
    static void printEigenvalues(const std::vector<std::complex<double>>& eigenvalues, const std::string& name = "固有値");

    // 固有ベクトルの表示
    static void printEigenvector(const std::vector<double>& eigenvector, double eigenvalue, const std::string& name = "固有ベクトル");

    // NaNチェック関数
    static bool isNaN(const std::complex<double>& value);
    static bool isNaN(const std::vector<double>& vector);
    static bool isNaN(const std::vector<std::complex<double>>& eigenvalues);

    // 特異値分解
    static void singularValueDecomposition(const std::vector<std::vector<double>>& matrix);

    // 主成分分析
    static void principalComponentAnalysis(const std::vector<std::vector<double>>& data);
};

// 線形方程式の解法
class LinearSolver {
public:
    // LU分解による解法
    static std::vector<double> solveLU(const std::vector<std::vector<double>>& A, const std::vector<double>& b);

    // QR分解による解法
    static std::vector<double> solveQR(const std::vector<std::vector<double>>& A, const std::vector<double>& b);

    // 最小二乗法
    static std::vector<double> leastSquares(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
};

// 数値解析
class NumericalAnalysis {
public:
    // 行列のノルム計算
    static double matrixNorm(const std::vector<std::vector<double>>& matrix, const std::string& normType = "frobenius");

    // ベクトルのノルム計算
    static double vectorNorm(const std::vector<double>& vector, const std::string& normType = "euclidean");

    // 誤差解析
    static void errorAnalysis(const std::vector<std::vector<double>>& A, const std::vector<double>& b,
                             const std::vector<double>& x_exact, const std::vector<double>& x_computed);

    // 安定性解析
    static void stabilityAnalysis(const std::vector<std::vector<double>>& matrix);
};

// ランダム行列生成とCSV保存
class RandomMatrixAnalysis {
public:
    // ランダム行列の生成
    static std::vector<std::vector<double>> generateRandomMatrix(int n, double minVal = -10.0, double maxVal = 10.0);

    // ランダムベクトルの生成
    static std::vector<double> generateRandomVector(int n, double minVal = -10.0, double maxVal = 10.0);

    // 行列のCSV保存
    static void saveMatrixToCSV(const std::vector<std::vector<double>>& matrix, const std::string& filename);

    // ベクトルのCSV保存
    static void saveVectorToCSV(const std::vector<double>& vector, const std::string& filename);

    // 結果のCSV保存
    static void saveResultsToCSV(const std::string& filename,
                                const std::vector<int>& sizes,
                                const std::vector<double>& determinants,
                                const std::vector<double>& conditionNumbers,
                                const std::vector<int>& ranks,
                                const std::vector<std::vector<std::complex<double>>>& allEigenvalues,
                                const std::vector<double>& computationTimes);

    // 計算精度・時間のCSV保存
    static void savePerformanceToCSV(const std::string& filename,
                                    const std::vector<int>& sizes,
                                    const std::vector<double>& computationTimes,
                                    const std::vector<double>& conditionNumbers);

    // 行列特性のCSV保存
    static void saveMatrixPropertiesToCSV(const std::string& filename,
                                         const std::vector<int>& sizes,
                                         const std::vector<double>& determinants,
                                         const std::vector<int>& ranks,
                                         const std::vector<std::vector<std::complex<double>>>& allEigenvalues);

    // n=1~100のランダム行列テスト実行
    static void runRandomMatrixTest(int maxSize = 100, int numTests = 1);

    // 単一サイズのテスト実行
    static void runSingleSizeTest(int n, int testIndex = 0);
};

} // namespace LinearAlgebra

#endif // LINEAR_ALGEBRA_HPP
