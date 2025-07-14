#include "../include/linear_algebra.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <limits>
#include <random>
#include <chrono>
#include <fstream>
#include <complex>
#include <cmath>
#include <set>

namespace LinearAlgebra {

// MatrixOperations クラスの実装
void MatrixOperations::printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << ":" << std::endl;
    for (const auto& row : matrix) {
        for (double val : row) {
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << val;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void MatrixOperations::saveMatrix(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& row : matrix) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i < row.size() - 1) file << " ";
            }
            file << std::endl;
        }
        file.close();
    }
}

std::vector<std::vector<double>> MatrixOperations::loadMatrix(const std::string& filename) {
    std::vector<std::vector<double>> matrix;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream iss(line);
        double val;
        while (iss >> val) {
            row.push_back(val);
        }
        if (!row.empty()) {
            matrix.push_back(row);
        }
    }
    file.close();
    return matrix;
}

double MatrixOperations::conditionNumber(const std::vector<std::vector<double>>& matrix) {
    // 簡易な条件数の計算（最大値 / 最小値）
    // より正確な計算は固有値から行うため、ここでは簡易計算のみ
    double maxVal = 0.0, minVal = std::numeric_limits<double>::max();

    for (const auto& row : matrix) {
        for (double val : row) {
            maxVal = std::max(maxVal, std::abs(val));
            minVal = std::min(minVal, std::abs(val));
        }
    }

    return (minVal > 1e-10) ? maxVal / minVal : std::numeric_limits<double>::infinity();
}

int MatrixOperations::rank(const std::vector<std::vector<double>>& matrix) {
    // 簡単なランク計算（非ゼロ行の数）
    int rank = 0;
    for (const auto& row : matrix) {
        bool hasNonZero = false;
        for (double val : row) {
            if (std::abs(val) > 1e-10) {
                hasNonZero = true;
                break;
            }
        }
        if (hasNonZero) rank++;
    }
    return rank;
}

// 小行列の取得（行列式計算用）
std::vector<std::vector<double>> MatrixOperations::getMinor(const std::vector<std::vector<double>>& matrix, int row, int col) {
    int n = matrix.size();
    std::vector<std::vector<double>> minor(n-1, std::vector<double>(n-1));

    int minorRow = 0;
    for (int i = 0; i < n; i++) {
        if (i == row) continue;
        int minorCol = 0;
        for (int j = 0; j < n; j++) {
            if (j == col) continue;
            minor[minorRow][minorCol] = matrix[i][j];
            minorCol++;
        }
        minorRow++;
    }

    return minor;
}

// 行列式の計算（LU分解によるアプローチ）
double MatrixOperations::determinant(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));

    // LU分解
    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
        for (int j = i; j < n; j++) {
            U[i][j] = matrix[i][j];
            for (int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
        for (int j = i + 1; j < n; j++) {
            L[j][i] = matrix[j][i];
            for (int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }

    // 行列式 = det(L) * det(U) = 1 * (対角要素の積)
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= U[i][i];
    }

    return det;
}

// 行列の転置
std::vector<std::vector<double>> MatrixOperations::transpose(const std::vector<std::vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

// 行列の積
std::vector<std::vector<double>> MatrixOperations::matrixMultiply(const std::vector<std::vector<double>>& A,
                                                                  const std::vector<std::vector<double>>& B) {
    int rowsA = A.size();
    int colsA = A[0].size();
    int colsB = B[0].size();

    std::vector<std::vector<double>> result(rowsA, std::vector<double>(colsB, 0.0));

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
            for (int k = 0; k < colsA; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

// 単位行列の作成
std::vector<std::vector<double>> MatrixOperations::identity(int n) {
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        I[i][i] = 1.0;
    }
    return I;
}

// 行列の等価性チェック
bool MatrixOperations::isEqual(const std::vector<std::vector<double>>& A,
                              const std::vector<std::vector<double>>& B, double tolerance) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        return false;
    }

    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++) {
            if (std::abs(A[i][j] - B[i][j]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

// 行列のコピー
std::vector<std::vector<double>> MatrixOperations::copyMatrix(const std::vector<std::vector<double>>& matrix) {
    return matrix;
}

// 行列式の保存（data/det/<N>）
void MatrixOperations::saveDeterminantToCSV(int n, double determinant) {
    std::string filename = "data/det/" + std::to_string(n);
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    if (file.is_open()) {
        // 行列式の値を科学記数法で保存
        file << std::scientific << std::setprecision(10) << determinant;
        file.close();
    }
}

void MatrixOperations::saveDeterminantToFile(int n, double determinant) {
    // data/detディレクトリの作成
    std::string mkdir_cmd = "mkdir -p data/det";
    system(mkdir_cmd.c_str());

    std::string filename = "data/det/" + std::to_string(n);
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    if (file.is_open()) {
        // 行列式の値を科学記数法で上書き保存（最新のみ）
        file << std::scientific << std::setprecision(10) << determinant;
        file.close();
    }
}

// EigenvalueAnalysis クラスの実装

// ベクトルのノルム計算
double NumericalAnalysis::vectorNorm(const std::vector<double>& vector, const std::string& normType) {
    if (normType == "euclidean") {
        double sum = 0.0;
        for (double val : vector) {
            sum += val * val;
        }
        return std::sqrt(sum);
    }
    return 0.0;
}

// LinearSolver クラスの実装
std::vector<double> LinearSolver::solveLU(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    int n = A.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    std::vector<double> x(n, 0.0);
    std::vector<double> y(n, 0.0);

    // LU分解
    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
        for (int j = i; j < n; j++) {
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
        for (int j = i + 1; j < n; j++) {
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }

    // 前進代入 (Ly = b)
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }

    // 後退代入 (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

// NumericalAnalysis クラスの実装

// EigenvalueAnalysis クラスの実装

// 改良されたQR分解による固有値・固有ベクトル計算
std::pair<std::vector<std::complex<double>>, std::vector<std::vector<double>>>
EigenvalueAnalysis::qrEigenDecomposition(const std::vector<std::vector<double>>& matrix, int maxIterations, double tolerance) {
    int n = matrix.size();
    std::vector<std::vector<double>> A = MatrixOperations::copyMatrix(matrix);
    std::vector<std::vector<double>> Q_total = MatrixOperations::identity(n);
    std::vector<double> prev_diagonal(n);

    for (int iter = 0; iter < maxIterations; iter++) {
        // 前回の対角成分を保存
        for (int i = 0; i < n; i++) {
            prev_diagonal[i] = A[i][i];
        }

        // 改良されたHouseholder変換によるQR分解
        auto [Q, R] = EigenvalueAnalysis::qrDecomposition(A);

        // R*Qの計算（行列乗算を最適化）
        std::vector<std::vector<double>> RQ(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    RQ[i][j] += R[i][k] * Q[k][j];
                }
            }
        }
        A = RQ;

        // Q_totalの更新（行列乗算を最適化）
        std::vector<std::vector<double>> Q_new(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    Q_new[i][j] += Q_total[i][k] * Q[k][j];
                }
            }
        }
        Q_total = Q_new;

        // 収束判定（対角成分の変化をチェック）
        if (iter > 0) {
            bool converged = true;
            double max_change = 0.0;
            for (int i = 0; i < n; i++) {
                double change = std::abs(A[i][i] - prev_diagonal[i]);
                max_change = std::max(max_change, change);
                if (change > tolerance) {
                    converged = false;
                }
            }
            // 収束した場合は反復を終了
            if (converged) {
                break;
            }
        }
    }

    // 固有値の抽出（対角成分）
    std::vector<std::complex<double>> eigenvalues(n);
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = std::complex<double>(A[i][i], 0.0);
    }

    // 固有ベクトルは累積Q行列の列ベクトル
    std::vector<std::vector<double>> eigenvectors = Q_total;

    return {eigenvalues, eigenvectors};
}

// 固有値計算のみ（QR法）
std::vector<std::complex<double>>
EigenvalueAnalysis::qrEigenvalueComputation(const std::vector<std::vector<double>>& matrix, int maxIterations, double tolerance) {
    int n = matrix.size();
    std::vector<std::vector<double>> A = MatrixOperations::copyMatrix(matrix);
    std::vector<double> prev_diagonal(n);

    for (int iter = 0; iter < maxIterations; iter++) {
        // 前回の対角成分を保存
        for (int i = 0; i < n; i++) {
            prev_diagonal[i] = A[i][i];
        }

        // QR分解
        auto [Q, R] = EigenvalueAnalysis::qrDecomposition(A);

        // R*Qの計算（行列乗算を最適化）
        std::vector<std::vector<double>> RQ(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    RQ[i][j] += R[i][k] * Q[k][j];
                }
            }
        }
        A = RQ;

        // 収束判定（対角成分の変化をチェック）
        if (iter > 0) {
            bool converged = true;
            double max_change = 0.0;
            for (int i = 0; i < n; i++) {
                double change = std::abs(A[i][i] - prev_diagonal[i]);
                max_change = std::max(max_change, change);
                if (change > tolerance) {
                    converged = false;
                }
            }
            // 収束した場合は反復を終了
            if (converged) {
                break;
            }
        }
    }

    // 固有値の抽出（対角成分）
    std::vector<std::complex<double>> eigenvalues(n);
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = std::complex<double>(A[i][i], 0.0);
    }

    return eigenvalues;
}

// 固有ベクトル計算のみ（累積Q行列から抽出）
std::vector<std::vector<double>>
EigenvalueAnalysis::qrEigenvectorComputation(const std::vector<std::vector<double>>& matrix, int maxIterations, double tolerance) {
    int n = matrix.size();
    std::vector<std::vector<double>> A = MatrixOperations::copyMatrix(matrix);
    std::vector<std::vector<double>> Q_total = MatrixOperations::identity(n);
    std::vector<double> prev_diagonal(n);

    for (int iter = 0; iter < maxIterations; iter++) {
        // 前回の対角成分を保存
        for (int i = 0; i < n; i++) {
            prev_diagonal[i] = A[i][i];
        }

        // QR分解
        auto [Q, R] = EigenvalueAnalysis::qrDecomposition(A);

        // R*Qの計算（行列乗算を最適化）
        std::vector<std::vector<double>> RQ(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    RQ[i][j] += R[i][k] * Q[k][j];
                }
            }
        }
        A = RQ;

        // Q_totalの更新（行列乗算を最適化）
        std::vector<std::vector<double>> Q_new(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    Q_new[i][j] += Q_total[i][k] * Q[k][j];
                }
            }
        }
        Q_total = Q_new;

        // 収束判定（対角成分の変化をチェック）
        if (iter > 0) {
            bool converged = true;
            double max_change = 0.0;
            for (int i = 0; i < n; i++) {
                double change = std::abs(A[i][i] - prev_diagonal[i]);
                max_change = std::max(max_change, change);
                if (change > tolerance) {
                    converged = false;
                }
            }
            // 収束した場合は反復を終了
            if (converged) {
                break;
            }
        }
    }

    // 固有ベクトルは累積Q行列の列ベクトル
    std::vector<std::vector<double>> eigenvectors = Q_total;

    return eigenvectors;
}

// 固有値の表示
void EigenvalueAnalysis::printEigenvalues(const std::vector<std::complex<double>>& eigenvalues, const std::string& name) {
    std::cout << name << ":" << std::endl;
    for (size_t i = 0; i < eigenvalues.size(); i++) {
        std::cout << "λ[" << i << "] = " << eigenvalues[i].real();
        if (std::abs(eigenvalues[i].imag()) > 1e-10) {
            std::cout << " + " << eigenvalues[i].imag() << "i";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Householder変換によるQR分解
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
EigenvalueAnalysis::qrDecomposition(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> R = matrix;
    std::vector<std::vector<double>> Q = MatrixOperations::identity(n);

    for (int k = 0; k < n - 1; ++k) {
        // xベクトルの抽出
        std::vector<double> x(n - k);
        for (int i = k; i < n; ++i) x[i - k] = R[i][k];
        // ノルム計算
        double norm_x = 0.0;
        for (double xi : x) norm_x += xi * xi;
        norm_x = std::sqrt(norm_x);
        if (norm_x < 1e-12) continue;
        // vベクトルの構築
        std::vector<double> v = x;
        v[0] += (x[0] >= 0 ? norm_x : -norm_x);
        // vの正規化
        double norm_v = 0.0;
        for (double vi : v) norm_v += vi * vi;
        norm_v = std::sqrt(norm_v);
        if (norm_v < 1e-12) continue;
        for (double& vi : v) vi /= norm_v;
        // R = H_k * R
        for (int j = k; j < n; ++j) {
            double dot = 0.0;
            for (int i = k; i < n; ++i) dot += v[i - k] * R[i][j];
            for (int i = k; i < n; ++i) R[i][j] -= 2.0 * v[i - k] * dot;
        }
        // Q = Q * H_k
        for (int i = 0; i < n; ++i) {
            double dot = 0.0;
            for (int j = k; j < n; ++j) dot += Q[i][j] * v[j - k];
            for (int j = k; j < n; ++j) Q[i][j] -= 2.0 * v[j - k] * dot;
        }
    }
    return {Q, R};
}

// RandomMatrixAnalysis クラスの実装

// ランダム行列の生成
std::vector<std::vector<double>> RandomMatrixAnalysis::generateRandomMatrix(int n, double minVal, double maxVal) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(minVal, maxVal);

    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = dis(gen);
        }
    }
    return matrix;
}

// 対称ランダム行列の生成（固有値計算用）
std::vector<std::vector<double>> RandomMatrixAnalysis::generateSymmetricMatrix(int n, double minVal, double maxVal) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(minVal, maxVal);

    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));

    // 上三角部分を生成
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            matrix[i][j] = dis(gen);
        }
    }

    // 対称性を保つために下三角部分をコピー
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            matrix[i][j] = matrix[j][i];
        }
    }

    return matrix;
}

// ランダムベクトルの生成
std::vector<double> RandomMatrixAnalysis::generateRandomVector(int n, double minVal, double maxVal) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(minVal, maxVal);

    std::vector<double> vector(n);
    for (int i = 0; i < n; i++) {
        vector[i] = dis(gen);
    }

    return vector;
}

// 行列のCSV保存
void RandomMatrixAnalysis::saveMatrixToCSV(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[i].size(); j++) {
                file << matrix[i][j];
                if (j < matrix[i].size() - 1) file << ",";
            }
            file << std::endl;
        }
        file.close();
    }
}

// ベクトルのCSV保存
void RandomMatrixAnalysis::saveVectorToCSV(const std::vector<double>& vector, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (size_t i = 0; i < vector.size(); i++) {
            file << vector[i];
            if (i < vector.size() - 1) file << ",";
        }
        file << std::endl;
        file.close();
    }
}

// 詳細計算時間のCSV保存
void RandomMatrixAnalysis::saveDetailedTimesToCSV(const std::string& filename,
                                                  const std::vector<int>& sizes,
                                                  const std::vector<ComputationTimes>& times) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // ヘッダー行
    file << "Size,DeterminantTime,EigenTime,LinearSolverTime\n";

    // データ行
    for (size_t i = 0; i < sizes.size(); i++) {
        file << sizes[i] << ","
             << times[i].determinantTime << ","
             << times[i].eigenvalueComputationTime << ","
             << times[i].linearSolverTime << "\n";
    }

    file.close();
}

// 行列式計算時間のCSV保存
void RandomMatrixAnalysis::saveDeterminantTimesToCSV(const std::string& filename,
                                                    const std::vector<int>& sizes,
                                                    const std::vector<ComputationTimes>& times) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // ヘッダー行
        file << "Size,DeterminantTime" << std::endl;

        // データ行
        for (size_t i = 0; i < sizes.size(); i++) {
            file << sizes[i] << ","
                 << std::fixed << std::setprecision(6) << times[i].determinantTime << std::endl;
        }
        file.close();
    }
}

// 固有値・固有ベクトル計算時間のCSV保存（統合）
void RandomMatrixAnalysis::saveEigenTimesToCSV(const std::string& filename,
                                              const std::vector<int>& sizes,
                                              const std::vector<ComputationTimes>& times) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // ヘッダー行
    file << "Size,EigenTime\n";

    // データ行
    for (size_t i = 0; i < sizes.size(); i++) {
        file << sizes[i] << ","
             << times[i].eigenvalueComputationTime << "\n";
    }

    file.close();
}

// 線形方程式解法時間のCSV保存
void RandomMatrixAnalysis::saveLinearSolverTimesToCSV(const std::string& filename,
                                                     const std::vector<int>& sizes,
                                                     const std::vector<ComputationTimes>& times) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // ヘッダー行
        file << "Size,LinearSolverTime" << std::endl;

        // データ行
        for (size_t i = 0; i < sizes.size(); i++) {
            file << sizes[i] << ","
                 << std::fixed << std::setprecision(6) << times[i].linearSolverTime << std::endl;
        }
        file.close();
    }
}

// 行列特性のCSV保存
void RandomMatrixAnalysis::saveMatrixPropertiesToCSV(const std::string& filename,
                                                    const std::vector<int>& sizes,
                                                    const std::vector<double>& determinants,
                                                    const std::vector<int>& ranks,
                                                    const std::vector<double>& conditionNumbers,
                                                    const std::vector<std::vector<std::complex<double>>>& allEigenvalues) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // ヘッダー行
        file << "Size,Determinant,Rank,ConditionNumber,Eigenvalues" << std::endl;

        // データ行
        for (size_t i = 0; i < sizes.size(); i++) {
            file << sizes[i] << ","
                 << determinants[i] << ","
                 << ranks[i] << ","
                 << std::fixed << std::setprecision(6) << conditionNumbers[i] << ",";

            // 固有値を文字列として保存
            std::string eigenStr = "";
            for (size_t j = 0; j < allEigenvalues[i].size(); j++) {
                if (j > 0) eigenStr += ";";
                eigenStr += std::to_string(allEigenvalues[i][j].real());
            }
            file << eigenStr << std::endl;
        }
        file.close();
    }
}

// 単一サイズのテスト実行
void RandomMatrixAnalysis::runSingleSizeTest(int n, int testIndex) {
    std::cout << "サイズ " << n << " のランダム行列テスト " << testIndex << " を実行中..." << std::endl;

    // runSingleSizeTestではディレクトリ作成は行わない

    // 対称行列の生成（固有値計算用）
    auto symmetricMatrix = generateSymmetricMatrix(n);

    // 右辺ベクトルの生成
    auto b = generateRandomVector(n);

    // 連立1次方程式の解を計算（時間測定）
    auto linearSolverStart = std::chrono::high_resolution_clock::now();
    auto x = LinearSolver::solveLU(symmetricMatrix, b);
    auto linearSolverEnd = std::chrono::high_resolution_clock::now();
    auto linearSolverDuration = std::chrono::duration_cast<std::chrono::microseconds>(linearSolverEnd - linearSolverStart);
    double linearSolverTime = linearSolverDuration.count() / 1000.0;

    // runSingleSizeTestではファイル保存は行わない

    // 行列式計算（時間測定）
    auto detStart = std::chrono::high_resolution_clock::now();
    double determinant = MatrixOperations::determinant(symmetricMatrix);
    auto detEnd = std::chrono::high_resolution_clock::now();
    auto detDuration = std::chrono::duration_cast<std::chrono::microseconds>(detEnd - detStart);
    double determinantTime = detDuration.count() / 1000.0;

    // runSingleSizeTestではファイル保存は行わない

    // 条件数とランク計算
    int rank = MatrixOperations::rank(symmetricMatrix);

    // 固有値・固有ベクトル計算（統合して時間測定）
    auto eigenStart = std::chrono::high_resolution_clock::now();
    auto [eigenvalues, eigenvectors] = EigenvalueAnalysis::qrEigenDecomposition(symmetricMatrix);
    auto eigenEnd = std::chrono::high_resolution_clock::now();
    auto eigenDuration = std::chrono::duration_cast<std::chrono::microseconds>(eigenEnd - eigenStart);
    double eigenTime = eigenDuration.count() / 1000.0;

    // 条件数計算（固有値から計算）
    double conditionNumber = 0.0;
    if (!eigenvalues.empty()) {
        double maxEigenvalue = 0.0;
        double minEigenvalue = std::numeric_limits<double>::max();

        for (const auto& eigenval : eigenvalues) {
            double absVal = std::abs(eigenval.real());
            if (absVal > 1e-10) {  // ゼロに近い値は除外
                maxEigenvalue = std::max(maxEigenvalue, absVal);
                minEigenvalue = std::min(minEigenvalue, absVal);
            }
        }

        if (minEigenvalue > 1e-10) {
            conditionNumber = maxEigenvalue / minEigenvalue;
        } else {
            conditionNumber = std::numeric_limits<double>::infinity();
        }
    }

    // 詳細時間を記録
    ComputationTimes times;
    times.determinantTime = determinantTime;
    times.eigenvalueComputationTime = eigenTime; // 固有値・固有ベクトル計算時間
    times.eigenvectorComputationTime = 0.0; // 統合されたため0
    times.linearSolverTime = linearSolverTime;

    // 結果の表示
    std::cout << "サイズ " << n << ": 行列式=" << std::fixed << std::setprecision(6) << determinant << ", 条件数=" << std::fixed << std::setprecision(6) << conditionNumber << ", ランク=" << rank << ", 計算時間=" << (determinantTime + eigenTime + linearSolverTime) << "ms" << std::endl;

    // 固有値の表示（最初の5個まで）
    std::cout << "固有値（最初の5個）:" << std::endl;
    for (size_t i = 0; i < std::min(eigenvalues.size(), size_t(5)); i++) {
        std::cout << "  λ[" << i << "] = " << std::fixed << std::setprecision(6) << eigenvalues[i].real() << std::endl;
    }
    if (eigenvalues.size() > 5) {
        std::cout << "  ... (他 " << eigenvalues.size() - 5 << " 個の固有値)" << std::endl;
    }
    std::cout << std::endl;

    // runSingleSizeTestではファイル保存は行わない
}

// n=1~100のランダム行列テスト実行
void RandomMatrixAnalysis::runRandomMatrixTest(int maxSize, int numTests) {
    std::cout << "=== ランダム行列テスト (n=1~" << maxSize << " 各サイズ" << numTests << "回) ===" << std::endl;

    // データディレクトリの作成
    system("mkdir -p data");
    system("mkdir -p data/eigen");
    system("mkdir -p data/A");
    system("mkdir -p data/B");
    system("mkdir -p data/x");

    std::vector<int> sizes;
    std::vector<double> determinants;
    std::vector<double> conditionNumbers;
    std::vector<int> ranks;
    std::vector<std::vector<std::complex<double>>> allEigenvalues;
    std::vector<double> computationTimes;
    std::vector<ComputationTimes> detailedTimes;

    // 各サイズでテスト実行
    for (int n = 1; n <= maxSize; n++) {
        for (int test = 0; test < numTests; test++) {
            // 対称行列の生成（固有値計算用）
            auto symmetricMatrix = generateSymmetricMatrix(n);

            // 右辺ベクトルの生成
            auto b = generateRandomVector(n);

            // 連立1次方程式の解を計算（時間測定）
            auto linearSolverStart = std::chrono::high_resolution_clock::now();
            auto x = LinearSolver::solveLU(symmetricMatrix, b);
            auto linearSolverEnd = std::chrono::high_resolution_clock::now();
            auto linearSolverDuration = std::chrono::duration_cast<std::chrono::microseconds>(linearSolverEnd - linearSolverStart);
            double linearSolverTime = linearSolverDuration.count() / 1000.0;

            // 対称行列と右辺ベクトルを保存
            std::string matrixAFilename = "data/A/" + std::to_string(n) + ".csv";
            std::string vectorBFilename = "data/B/" + std::to_string(n) + ".csv";
            std::string solutionXFilename = "data/x/" + std::to_string(n) + ".csv";
            saveMatrixToCSV(symmetricMatrix, matrixAFilename);
            saveVectorToCSV(b, vectorBFilename);
            saveVectorToCSV(x, solutionXFilename);

            // 行列式計算（時間測定）
            auto detStart = std::chrono::high_resolution_clock::now();
            double determinant = MatrixOperations::determinant(symmetricMatrix);
            auto detEnd = std::chrono::high_resolution_clock::now();
            auto detDuration = std::chrono::duration_cast<std::chrono::microseconds>(detEnd - detStart);
            double determinantTime = detDuration.count() / 1000.0;

            // 行列式をdata/det/<N>ファイルに保存
            MatrixOperations::saveDeterminantToFile(n, determinant);

            // 条件数とランク計算
            int rank = MatrixOperations::rank(symmetricMatrix);

            // 固有値・固有ベクトル計算（統合して時間測定）
            auto eigenStart = std::chrono::high_resolution_clock::now();
            auto [eigenvalues, eigenvectors] = EigenvalueAnalysis::qrEigenDecomposition(symmetricMatrix);
            auto eigenEnd = std::chrono::high_resolution_clock::now();
            auto eigenDuration = std::chrono::duration_cast<std::chrono::microseconds>(eigenEnd - eigenStart);
            double eigenTime = eigenDuration.count() / 1000.0;

            // 条件数計算（固有値から計算）
            double conditionNumber = 0.0;
            if (!eigenvalues.empty()) {
                double maxEigenvalue = 0.0;
                double minEigenvalue = std::numeric_limits<double>::max();

                for (const auto& eigenval : eigenvalues) {
                    double absVal = std::abs(eigenval.real());
                    if (absVal > 1e-10) {  // ゼロに近い値は除外
                        maxEigenvalue = std::max(maxEigenvalue, absVal);
                        minEigenvalue = std::min(minEigenvalue, absVal);
                    }
                }

                if (minEigenvalue > 1e-10) {
                    conditionNumber = maxEigenvalue / minEigenvalue;
                } else {
                    conditionNumber = std::numeric_limits<double>::infinity();
                }
            }

            // 詳細時間を記録
            ComputationTimes times;
            times.determinantTime = determinantTime;
            times.eigenvalueComputationTime = eigenTime; // 固有値・固有ベクトル計算時間
            times.eigenvectorComputationTime = 0.0; // 統合されたため0
            times.linearSolverTime = linearSolverTime;

            // 結果を保存
            sizes.push_back(n);
            determinants.push_back(determinant);
            conditionNumbers.push_back(conditionNumber);
            ranks.push_back(rank);
            allEigenvalues.push_back(eigenvalues);
            computationTimes.push_back(eigenTime);
            detailedTimes.push_back(times);

            // 固有値と固有ベクトルをCSV保存（横に並べて）
            std::string eigenFilename = "data/eigen/" + std::to_string(n) + ".csv";
            std::ofstream eigenFile(eigenFilename);
            if (eigenFile.is_open()) {
                // ヘッダー行
                eigenFile << "Index,Eigenvalue";
                for (int i = 0; i < n; i++) {
                    eigenFile << ",Eigenvector_" << i;
                }
                eigenFile << std::endl;

                // データ行（各固有値と対応する固有ベクトル）
                for (int i = 0; i < n; i++) {
                    eigenFile << i << "," << eigenvalues[i].real();
                    for (int j = 0; j < n; j++) {
                        eigenFile << "," << eigenvectors[j][i]; // 列ベクトルとして保存
                    }
                    eigenFile << std::endl;
                }
                eigenFile.close();
            }

            // 進捗表示（10サイズごと）
            if (n % 10 == 0) {
                std::cout << "サイズ " << n << " 完了" << std::endl;
            }
        }
    }

    // 結果をCSVファイルに保存
    std::string propertiesFilename = "data/matrix_properties.csv";
    std::string detailedTimesFilename = "data/detailed_computation_times.csv";
    std::string determinantTimesFilename = "data/determinant_times.csv";
    std::string eigenTimesFilename = "data/eigen_times.csv";
    std::string linearSolverTimesFilename = "data/linear_solver_times.csv";

    saveMatrixPropertiesToCSV(propertiesFilename, sizes, determinants, ranks, conditionNumbers, allEigenvalues);
    saveDetailedTimesToCSV(detailedTimesFilename, sizes, detailedTimes);
    saveDeterminantTimesToCSV(determinantTimesFilename, sizes, detailedTimes);
    saveEigenTimesToCSV(eigenTimesFilename, sizes, detailedTimes);
    saveLinearSolverTimesToCSV(linearSolverTimesFilename, sizes, detailedTimes);

    std::cout << "\n=== テスト完了 ===" << std::endl;
    std::cout << "結果を以下のファイルに保存しました:" << std::endl;
    std::cout << "  行列特性: " << propertiesFilename << std::endl;
    std::cout << "  詳細計算時間: " << detailedTimesFilename << std::endl;
    std::cout << "  行列式計算時間: " << determinantTimesFilename << std::endl;
    std::cout << "  固有値・固有ベクトル計算時間: " << eigenTimesFilename << std::endl;
    std::cout << "  線形方程式解法時間: " << linearSolverTimesFilename << std::endl;

    // 統計情報の表示
    if (!computationTimes.empty()) {
        double avgTime = std::accumulate(computationTimes.begin(), computationTimes.end(), 0.0) / computationTimes.size();
        double maxTime = *std::max_element(computationTimes.begin(), computationTimes.end());
        double minTime = *std::min_element(computationTimes.begin(), computationTimes.end());

        std::cout << "計算時間統計:" << std::endl;
        std::cout << "  平均: " << avgTime << "ms" << std::endl;
        std::cout << "  最大: " << maxTime << "ms" << std::endl;
        std::cout << "  最小: " << minTime << "ms" << std::endl;
    }
}

} // namespace LinearAlgebra

