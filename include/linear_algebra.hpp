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

// QR法の反復回数設定
#ifndef QR_MAX_ITERATIONS
#define QR_MAX_ITERATIONS 1000
#endif

#ifndef QR_TOLERANCE
#define QR_TOLERANCE 1e-8
#endif

// 計算時間を各項目ごとに記録する構造体
struct ComputationTimes {
    double determinantTime = 0.0;      // 行列式計算時間
    double eigenvalueTime = 0.0;       // 固有値・固有ベクトル計算時間
    double linearSolverTime = 0.0;     // 線形方程式解法時間
    double totalTime = 0.0;            // 総計算時間

    // 合計時間を計算
    void calculateTotal() {
        totalTime = determinantTime + eigenvalueTime + linearSolverTime;
    }
};

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

    // 行列式の保存（data/determinants.csv）
    static void saveDeterminantToCSV(int n, double determinant);

    // 行列式の保存（data/det/<N>）
    static void saveDeterminantToFile(int n, double determinant);
};

// 固有値・固有ベクトル計算
class EigenvalueAnalysis {
public:
    // QR法による固有値・固有ベクトル計算
    // 戻り値: pair(固有値リスト, 固有ベクトル行列)
    static std::pair<std::vector<std::complex<double>>, std::vector<std::vector<double>>>
    qrEigenDecomposition(const std::vector<std::vector<double>>& matrix, int maxIterations = QR_MAX_ITERATIONS, double tolerance = QR_TOLERANCE);

    // QR分解
    static std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> qrDecomposition(
        const std::vector<std::vector<double>>& matrix);

    // 固有値の表示
    static void printEigenvalues(const std::vector<std::complex<double>>& eigenvalues, const std::string& name = "固有値");

    // シフト付きQR法による固有値・固有ベクトル計算（高速版）
    static std::pair<std::vector<std::complex<double>>, std::vector<std::vector<double>>>
    shiftedQREigenDecomposition(
        const std::vector<std::vector<double>>& matrix,
        int maxIterations = 200,
        double tolerance = 1e-8
    );
};

// 線形方程式の解法
class LinearSolver {
public:
    // LU分解による解法
    static std::vector<double> solveLU(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
};

// 数値解析
class NumericalAnalysis {
public:
    // ベクトルのノルム計算
    static double vectorNorm(const std::vector<double>& vector, const std::string& normType = "euclidean");
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



    // 詳細計算時間のCSV保存
    static void saveDetailedTimesToCSV(const std::string& filename,
                                       const std::vector<int>& sizes,
                                       const std::vector<ComputationTimes>& times);

    // 各計算時間を別々のファイルに保存
    static void saveDeterminantTimesToCSV(const std::string& filename,
                                         const std::vector<int>& sizes,
                                         const std::vector<ComputationTimes>& times);

    static void saveEigenvalueTimesToCSV(const std::string& filename,
                                        const std::vector<int>& sizes,
                                        const std::vector<ComputationTimes>& times);

    static void saveLinearSolverTimesToCSV(const std::string& filename,
                                          const std::vector<int>& sizes,
                                          const std::vector<ComputationTimes>& times);

    // 行列特性のCSV保存
    static void saveMatrixPropertiesToCSV(const std::string& filename,
                                         const std::vector<int>& sizes,
                                         const std::vector<double>& determinants,
                                         const std::vector<int>& ranks,
                                         const std::vector<double>& conditionNumbers,
                                         const std::vector<std::vector<std::complex<double>>>& allEigenvalues);

    // n=1~100のランダム行列テスト実行
    static void runRandomMatrixTest(int maxSize = 100, int numTests = 1);

    // 単一サイズのテスト実行
    static void runSingleSizeTest(int n, int testIndex = 0);
};

} // namespace LinearAlgebra

#endif // LINEAR_ALGEBRA_HPP
