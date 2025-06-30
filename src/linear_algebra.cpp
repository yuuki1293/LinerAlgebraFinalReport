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
    // 簡単な条件数の計算（最大固有値 / 最小固有値）
    // 実際の実装ではより複雑な計算が必要
    double maxVal = 0.0, minVal = std::numeric_limits<double>::max();

    for (const auto& row : matrix) {
        for (double val : row) {
            maxVal = std::max(maxVal, std::abs(val));
            minVal = std::min(minVal, std::abs(val));
        }
    }

    return (minVal > 0) ? maxVal / minVal : std::numeric_limits<double>::infinity();
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

// 固有値の重複度をチェック
std::vector<int> EigenvalueAnalysis::checkEigenvalueMultiplicity(const std::vector<std::complex<double>>& eigenvalues,
                                                                double tolerance) {
    std::vector<int> multiplicities;
    std::vector<bool> counted(eigenvalues.size(), false);

    for (size_t i = 0; i < eigenvalues.size(); i++) {
        if (counted[i]) continue;

        int multiplicity = 1;
        for (size_t j = i + 1; j < eigenvalues.size(); j++) {
            if (!counted[j] && std::abs(eigenvalues[i] - eigenvalues[j]) < tolerance) {
                multiplicity++;
                counted[j] = true;
            }
        }
        multiplicities.push_back(multiplicity);
        counted[i] = true;
    }

    return multiplicities;
}

// 対角化可能性の判定
bool EigenvalueAnalysis::isDiagonalizable(const std::vector<std::vector<double>>& matrix) {
    // まず固有値を計算
    auto [eigenvalues, eigenvectors] = qrEigenDecomposition(matrix);

    // 固有値の重複度をチェック
    std::vector<int> multiplicities = checkEigenvalueMultiplicity(eigenvalues);

    // 重複度が1の固有値の数をカウント
    int distinctEigenvalues = 0;
    for (int multiplicity : multiplicities) {
        if (multiplicity == 1) {
            distinctEigenvalues++;
        }
    }

    // 行列のサイズと異なる固有値の数が一致すれば対角化可能
    bool hasDistinctEigenvalues = (distinctEigenvalues == static_cast<int>(matrix.size()));

    // 追加のチェック: ジョルダン標準形の判定
    // 上三角行列で対角要素が同じ場合、対角化不可能
    bool isUpperTriangular = true;
    bool hasRepeatedDiagonal = false;

    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            if (std::abs(matrix[i][j]) > 1e-10) {
                isUpperTriangular = false;
                break;
            }
        }
        if (!isUpperTriangular) break;
    }

    if (isUpperTriangular) {
        for (size_t i = 0; i < matrix.size() - 1; i++) {
            if (std::abs(matrix[i][i] - matrix[i+1][i+1]) < 1e-10) {
                hasRepeatedDiagonal = true;
                break;
            }
        }
    }

    // ジョルダン標準形の場合は対角化不可能
    if (isUpperTriangular && hasRepeatedDiagonal) {
        return false;
    }

    return hasDistinctEigenvalues;
}

// 対角化による固有値計算
std::vector<std::complex<double>> EigenvalueAnalysis::diagonalizationEigenvalues(const std::vector<std::vector<double>>& matrix) {
    std::cout << "対角化可能性をチェック中..." << std::endl;

    if (!isDiagonalizable(matrix)) {
        std::cout << "エラー: この行列は対角化できません。" << std::endl;
        std::cout << "理由: 固有値の重複度が高すぎるか、固有ベクトルが線形独立でない可能性があります。" << std::endl;

        // NaNの固有値を返す
        int n = matrix.size();
        return std::vector<std::complex<double>>(n, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), 0.0));
    }

    std::cout << "行列は対角化可能です。対角化による固有値計算を実行します。" << std::endl;

    // QR法で固有値を計算（対角化可能な場合）
    auto [eigenvalues, eigenvectors] = qrEigenDecomposition(matrix);
    return eigenvalues;
}

// べき乗法による最大固有値と固有ベクトルの計算
std::pair<double, std::vector<double>> EigenvalueAnalysis::powerMethod(const std::vector<std::vector<double>>& matrix,
                                                                       int maxIterations, double tolerance) {
    int n = matrix.size();

    // 初期ベクトル（ランダム）
    std::vector<double> x(n, 1.0);
    double norm = NumericalAnalysis::vectorNorm(x, "euclidean");
    for (int i = 0; i < n; i++) {
        x[i] /= norm;
    }

    double eigenvalue = 0.0;
    std::vector<double> eigenvector = x;

    for (int iter = 0; iter < maxIterations; iter++) {
        // Ax の計算
        std::vector<double> Ax(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Ax[i] += matrix[i][j] * x[j];
            }
        }

        // 新しい固有値の推定
        double newEigenvalue = 0.0;
        for (int i = 0; i < n; i++) {
            newEigenvalue += x[i] * Ax[i];
        }

        // 固有ベクトルの正規化
        double norm = NumericalAnalysis::vectorNorm(Ax, "euclidean");
        for (int i = 0; i < n; i++) {
            x[i] = Ax[i] / norm;
        }

        // 収束判定
        if (std::abs(newEigenvalue - eigenvalue) < tolerance) {
            eigenvalue = newEigenvalue;
            eigenvector = x;
            std::cout << "べき乗法が収束しました。反復回数: " << iter + 1 << std::endl;
            break;
        }

        eigenvalue = newEigenvalue;
        eigenvector = x;
    }

    return std::make_pair(eigenvalue, eigenvector);
}

// グラム・シュミット直交化
std::vector<std::vector<double>> EigenvalueAnalysis::gramSchmidt(const std::vector<std::vector<double>>& vectors) {
    int n = vectors.size();
    std::vector<std::vector<double>> orthonormal = vectors;

    for (int i = 0; i < n; i++) {
        // 正規化
        double norm = NumericalAnalysis::vectorNorm(orthonormal[i], "euclidean");
        for (size_t j = 0; j < orthonormal[i].size(); j++) {
            orthonormal[i][j] /= norm;
        }

        // 後続のベクトルから直交成分を除去
        for (int j = i + 1; j < n; j++) {
            double dotProduct = 0.0;
            for (size_t k = 0; k < orthonormal[i].size(); k++) {
                dotProduct += orthonormal[i][k] * orthonormal[j][k];
            }

            for (size_t k = 0; k < orthonormal[j].size(); k++) {
                orthonormal[j][k] -= dotProduct * orthonormal[i][k];
            }
        }
    }

    return orthonormal;
}

// QR分解
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> EigenvalueAnalysis::qrDecomposition(
    const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();

    // Q行列の初期化（Aの列ベクトル）
    std::vector<std::vector<double>> Q(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q[j][i] = matrix[j][i];
        }
    }

    // グラム・シュミット直交化
    Q = gramSchmidt(Q);

    // R行列の計算: R = Q^T * A
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                R[i][j] += Q[k][i] * matrix[k][j];  // Q^T[i][k] * A[k][j]
            }
        }
    }

    return std::make_pair(Q, R);
}

// QR法による固有値・固有ベクトル計算
std::pair<std::vector<std::complex<double>>, std::vector<std::vector<double>>>
EigenvalueAnalysis::qrEigenDecomposition(const std::vector<std::vector<double>>& matrix, int maxIterations, double tolerance) {
    int n = matrix.size();
    std::vector<std::vector<double>> A = MatrixOperations::copyMatrix(matrix);
    std::vector<std::vector<double>> V = MatrixOperations::identity(n); // 固有ベクトル格納用

    for (int iter = 0; iter < maxIterations; iter++) {
        // QR分解
        auto [Q, R] = qrDecomposition(A);
        A = MatrixOperations::matrixMultiply(R, Q);
        V = MatrixOperations::matrixMultiply(V, Q); // 固有ベクトルの更新

        // 収束判定
        bool converged = true;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && std::abs(A[i][j]) > tolerance) {
                    converged = false;
                    break;
                }
            }
            if (!converged) break;
        }
        if (converged) break;
    }
    std::vector<std::complex<double>> eigenvalues;
    for (int i = 0; i < n; i++) {
        eigenvalues.push_back(std::complex<double>(A[i][i], 0.0));
    }
    // Vの各列が固有ベクトル
    return {eigenvalues, V};
}

// 固有値分解（対角化による）
void EigenvalueAnalysis::eigenvalueDecomposition(const std::vector<std::vector<double>>& matrix) {
    std::cout << "=== 固有値分解（対角化による）===" << std::endl;

    // 対角化による固有値計算
    std::vector<std::complex<double>> eigenvalues = diagonalizationEigenvalues(matrix);

    if (eigenvalues.empty() || isNaN(eigenvalues)) {
        std::cout << "対角化できないため、固有値と固有ベクトルをNaNに設定します。" << std::endl;
        std::cout << "理由: 固有値の重複度が高すぎるか、固有ベクトルが線形独立でない可能性があります。" << std::endl;

        // NaNの固有値と固有ベクトルを作成
        int n = matrix.size();
        std::vector<std::complex<double>> nanEigenvalues(n, std::complex<double>(std::numeric_limits<double>::quiet_NaN(), 0.0));
        std::vector<double> nanEigenvector(n, std::numeric_limits<double>::quiet_NaN());

        printEigenvalues(nanEigenvalues, "対角化による固有値（NaN）");
        printEigenvector(nanEigenvector, std::numeric_limits<double>::quiet_NaN(), "対応する固有ベクトル（NaN）");
        return;
    }

    printEigenvalues(eigenvalues, "対角化による固有値");

    // べき乗法による最大固有値の検証
    std::cout << "\n=== べき乗法による検証 ===" << std::endl;
    auto [maxEigenvalue, eigenvector] = powerMethod(matrix);
    std::cout << "べき乗法による最大固有値: " << maxEigenvalue << std::endl;
    printEigenvector(eigenvector, maxEigenvalue, "対応する固有ベクトル");
}

// 固有値の表示
void EigenvalueAnalysis::printEigenvalues(const std::vector<std::complex<double>>& eigenvalues, const std::string& name) {
    std::cout << name << ":" << std::endl;
    for (size_t i = 0; i < eigenvalues.size(); i++) {
        if (isNaN(eigenvalues[i])) {
            std::cout << "λ[" << i << "] = NaN" << std::endl;
        } else if (std::abs(eigenvalues[i].imag()) < 1e-10) {
            std::cout << "λ[" << i << "] = " << std::fixed << std::setprecision(6) << eigenvalues[i].real() << std::endl;
        } else {
            std::cout << "λ[" << i << "] = " << std::fixed << std::setprecision(6) << eigenvalues[i].real() << " + " << std::fixed << std::setprecision(6) << eigenvalues[i].imag() << "i" << std::endl;
        }
    }
    std::cout << std::endl;
}

// 固有ベクトルの表示
void EigenvalueAnalysis::printEigenvector(const std::vector<double>& eigenvector, double eigenvalue, const std::string& name) {
    std::cout << name << " (λ = " << eigenvalue << "):" << std::endl;
    for (size_t i = 0; i < eigenvector.size(); i++) {
        if (std::isnan(eigenvector[i])) {
            std::cout << "v[" << i << "] = NaN" << std::endl;
        } else {
            std::cout << "v[" << i << "] = " << std::fixed << std::setprecision(6) << eigenvector[i] << std::endl;
        }
    }
    std::cout << std::endl;
}

// NaNチェック関数
bool EigenvalueAnalysis::isNaN(const std::complex<double>& value) {
    return std::isnan(value.real()) || std::isnan(value.imag());
}

bool EigenvalueAnalysis::isNaN(const std::vector<double>& vector) {
    for (double val : vector) {
        if (std::isnan(val)) return true;
    }
    return false;
}

bool EigenvalueAnalysis::isNaN(const std::vector<std::complex<double>>& eigenvalues) {
    for (const auto& eigenval : eigenvalues) {
        if (isNaN(eigenval)) return true;
    }
    return false;
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
double NumericalAnalysis::matrixNorm(const std::vector<std::vector<double>>& matrix, const std::string& normType) {
    if (normType == "frobenius") {
        double sum = 0.0;
        for (const auto& row : matrix) {
            for (double val : row) {
                sum += val * val;
            }
        }
        return std::sqrt(sum);
    }
    return 0.0;
}

void NumericalAnalysis::errorAnalysis(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& b,
                                     const std::vector<double>& x_exact,
                                     const std::vector<double>& x_computed) {
    int n = x_exact.size();

    // 誤差の計算
    double maxError = 0.0;
    double relativeError = 0.0;
    double exactNorm = 0.0;

    for (int i = 0; i < n; i++) {
        double error = std::abs(x_computed[i] - x_exact[i]);
        maxError = std::max(maxError, error);
        exactNorm += x_exact[i] * x_exact[i];
    }

    exactNorm = std::sqrt(exactNorm);
    relativeError = maxError / exactNorm;

    std::cout << "誤差解析結果:" << std::endl;
    std::cout << "最大絶対誤差: " << maxError << std::endl;
    std::cout << "相対誤差: " << relativeError << std::endl;
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
    if (file.is_open()) {
        // ヘッダー行
        file << "Size,DeterminantTime,EigenvalueTime,LinearSolverTime,TotalTime" << std::endl;

        // データ行
        for (size_t i = 0; i < sizes.size(); i++) {
            file << sizes[i] << ","
                 << std::fixed << std::setprecision(6) << times[i].determinantTime << ","
                 << std::fixed << std::setprecision(6) << times[i].eigenvalueTime << ","
                 << std::fixed << std::setprecision(6) << times[i].linearSolverTime << ","
                 << std::fixed << std::setprecision(6) << times[i].totalTime << std::endl;
        }
        file.close();
    }
}

// 行列特性のCSV保存
void RandomMatrixAnalysis::saveMatrixPropertiesToCSV(const std::string& filename,
                                                    const std::vector<int>& sizes,
                                                    const std::vector<double>& determinants,
                                                    const std::vector<int>& ranks,
                                                    const std::vector<std::vector<std::complex<double>>>& allEigenvalues) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // ヘッダー行
        file << "Size,Determinant,Rank,Eigenvalues" << std::endl;

        // データ行
        for (size_t i = 0; i < sizes.size(); i++) {
            file << sizes[i] << ","
                 << determinants[i] << ","
                 << ranks[i] << ",";

            // 固有値を文字列として保存
            std::string eigenStr = "";
            for (size_t j = 0; j < allEigenvalues[i].size(); j++) {
                if (j > 0) eigenStr += ";";
                if (std::abs(allEigenvalues[i][j].imag()) < 1e-10) {
                    eigenStr += std::to_string(allEigenvalues[i][j].real());
                } else {
                    eigenStr += std::to_string(allEigenvalues[i][j].real()) + "+" +
                               std::to_string(allEigenvalues[i][j].imag()) + "i";
                }
            }
            file << eigenStr << std::endl;
        }
        file.close();
    }
}

// 単一サイズのテスト実行
void RandomMatrixAnalysis::runSingleSizeTest(int n, int testIndex) {
    std::cout << "サイズ " << n << " のランダム行列テスト " << testIndex << " を実行中..." << std::endl;

    // データディレクトリの作成
    system("mkdir -p data");
    system("mkdir -p data/eigen");
    system("mkdir -p data/A");
    system("mkdir -p data/B");
    system("mkdir -p data/x");

    // ランダム行列の生成
    auto matrix = generateRandomMatrix(n);

    // 右辺ベクトルの生成
    auto b = generateRandomVector(n);

    // 連立1次方程式の解を計算（時間測定）
    auto linearSolverStart = std::chrono::high_resolution_clock::now();
    auto x = LinearSolver::solveLU(matrix, b);
    auto linearSolverEnd = std::chrono::high_resolution_clock::now();
    auto linearSolverDuration = std::chrono::duration_cast<std::chrono::microseconds>(linearSolverEnd - linearSolverStart);
    double linearSolverTime = linearSolverDuration.count() / 1000.0;

    // 元の行列と右辺ベクトルを保存
    std::string matrixAFilename = "data/A/" + std::to_string(n) + ".csv";
    std::string vectorBFilename = "data/B/" + std::to_string(n) + ".csv";
    std::string solutionXFilename = "data/x/" + std::to_string(n) + ".csv";
    saveMatrixToCSV(matrix, matrixAFilename);
    saveVectorToCSV(b, vectorBFilename);
    saveVectorToCSV(x, solutionXFilename);

    // 行列式計算（時間測定）
    auto detStart = std::chrono::high_resolution_clock::now();
    double determinant = MatrixOperations::determinant(matrix);
    auto detEnd = std::chrono::high_resolution_clock::now();
    auto detDuration = std::chrono::duration_cast<std::chrono::microseconds>(detEnd - detStart);
    double determinantTime = detDuration.count() / 1000.0;

    // 条件数とランク計算
    double conditionNumber = MatrixOperations::conditionNumber(matrix);
    int rank = MatrixOperations::rank(matrix);

    // 固有値・固有ベクトル計算（時間測定）
    auto eigenStart = std::chrono::high_resolution_clock::now();
    auto [eigenvalues, eigenvectors] = EigenvalueAnalysis::qrEigenDecomposition(matrix);
    auto eigenEnd = std::chrono::high_resolution_clock::now();
    auto eigenDuration = std::chrono::duration_cast<std::chrono::microseconds>(eigenEnd - eigenStart);
    double eigenvalueTime = eigenDuration.count() / 1000.0;

    // 総計算時間
    double totalTime = determinantTime + eigenvalueTime + linearSolverTime;

    // 結果の表示
    std::cout << "サイズ " << n << ": 行列式=" << std::fixed << std::setprecision(6) << determinant << ", 条件数=" << std::fixed << std::setprecision(6) << conditionNumber << ", ランク=" << rank << ", 計算時間=" << totalTime << "ms" << std::endl;

    // 固有値の表示（最初の5個まで）
    std::cout << "固有値（最初の5個）:" << std::endl;
    for (size_t i = 0; i < std::min(eigenvalues.size(), size_t(5)); i++) {
        if (std::abs(eigenvalues[i].imag()) < 1e-10) {
            std::cout << "  λ[" << i << "] = " << std::fixed << std::setprecision(6) << eigenvalues[i].real() << std::endl;
        } else {
            std::cout << "  λ[" << i << "] = " << std::fixed << std::setprecision(6) << eigenvalues[i].real() << " + " << std::fixed << std::setprecision(6) << eigenvalues[i].imag() << "i" << std::endl;
        }
    }
    if (eigenvalues.size() > 5) {
        std::cout << "  ... (他 " << eigenvalues.size() - 5 << " 個の固有値)" << std::endl;
    }
    std::cout << std::endl;

    // 行列と固有値・固有ベクトルの保存
    std::string matrixFilename = "data/random_matrix_" + std::to_string(n) + "_" + std::to_string(testIndex) + ".csv";
    std::string eigenFilename = "data/eigen/eigen_" + std::to_string(n) + "_" + std::to_string(testIndex) + ".csv";

    saveMatrixToCSV(matrix, matrixFilename);

    // 固有値と固有ベクトルをCSV保存（横に並べて）
    std::ofstream eigenFile(eigenFilename);
    if (eigenFile.is_open()) {
        // ヘッダー行
        eigenFile << "Index,Eigenvalue_Real,Eigenvalue_Imaginary";
        for (int i = 0; i < n; i++) {
            eigenFile << ",Eigenvector_" << i;
        }
        eigenFile << std::endl;

        // データ行（各固有値と対応する固有ベクトル）
        for (int i = 0; i < n; i++) {
            eigenFile << i << "," << eigenvalues[i].real() << "," << eigenvalues[i].imag();
            for (int j = 0; j < n; j++) {
                eigenFile << "," << eigenvectors[j][i]; // 列ベクトルとして保存
            }
            eigenFile << std::endl;
        }
        eigenFile.close();
    }
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
            // ランダム行列の生成
            auto matrix = generateRandomMatrix(n);

            // 右辺ベクトルの生成
            auto b = generateRandomVector(n);

            // 連立1次方程式の解を計算（時間測定）
            auto linearSolverStart = std::chrono::high_resolution_clock::now();
            auto x = LinearSolver::solveLU(matrix, b);
            auto linearSolverEnd = std::chrono::high_resolution_clock::now();
            auto linearSolverDuration = std::chrono::duration_cast<std::chrono::microseconds>(linearSolverEnd - linearSolverStart);
            double linearSolverTime = linearSolverDuration.count() / 1000.0;

            // 元の行列と右辺ベクトルを保存
            std::string matrixAFilename = "data/A/" + std::to_string(n) + ".csv";
            std::string vectorBFilename = "data/B/" + std::to_string(n) + ".csv";
            std::string solutionXFilename = "data/x/" + std::to_string(n) + ".csv";
            saveMatrixToCSV(matrix, matrixAFilename);
            saveVectorToCSV(b, vectorBFilename);
            saveVectorToCSV(x, solutionXFilename);

            // 行列式計算（時間測定）
            auto detStart = std::chrono::high_resolution_clock::now();
            double determinant = MatrixOperations::determinant(matrix);
            auto detEnd = std::chrono::high_resolution_clock::now();
            auto detDuration = std::chrono::duration_cast<std::chrono::microseconds>(detEnd - detStart);
            double determinantTime = detDuration.count() / 1000.0;

            // 条件数とランク計算
            double conditionNumber = MatrixOperations::conditionNumber(matrix);
            int rank = MatrixOperations::rank(matrix);

            // 固有値・固有ベクトル計算（時間測定）
            auto eigenStart = std::chrono::high_resolution_clock::now();
            auto [eigenvalues, eigenvectors] = EigenvalueAnalysis::qrEigenDecomposition(matrix);
            auto eigenEnd = std::chrono::high_resolution_clock::now();
            auto eigenDuration = std::chrono::duration_cast<std::chrono::microseconds>(eigenEnd - eigenStart);
            double eigenvalueTime = eigenDuration.count() / 1000.0;

            // 総計算時間
            double totalTime = determinantTime + eigenvalueTime + linearSolverTime;

            // 詳細時間を記録
            ComputationTimes times;
            times.determinantTime = determinantTime;
            times.eigenvalueTime = eigenvalueTime;
            times.linearSolverTime = linearSolverTime;
            times.totalTime = totalTime;

            // 結果を保存
            sizes.push_back(n);
            determinants.push_back(determinant);
            conditionNumbers.push_back(conditionNumber);
            ranks.push_back(rank);
            allEigenvalues.push_back(eigenvalues);
            computationTimes.push_back(totalTime);
            detailedTimes.push_back(times);

            // 固有値と固有ベクトルをCSV保存（横に並べて）
            std::string eigenFilename = "data/eigen/eigen_" + std::to_string(n) + "_" + std::to_string(test) + ".csv";
            std::ofstream eigenFile(eigenFilename);
            if (eigenFile.is_open()) {
                // ヘッダー行
                eigenFile << "Index,Eigenvalue_Real,Eigenvalue_Imaginary";
                for (int i = 0; i < n; i++) {
                    eigenFile << ",Eigenvector_" << i;
                }
                eigenFile << std::endl;

                // データ行（各固有値と対応する固有ベクトル）
                for (int i = 0; i < n; i++) {
                    eigenFile << i << "," << eigenvalues[i].real() << "," << eigenvalues[i].imag();
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

    saveMatrixPropertiesToCSV(propertiesFilename, sizes, determinants, ranks, allEigenvalues);
    saveDetailedTimesToCSV(detailedTimesFilename, sizes, detailedTimes);

    std::cout << "\n=== テスト完了 ===" << std::endl;
    std::cout << "結果を以下のファイルに保存しました:" << std::endl;
    std::cout << "  行列特性: " << propertiesFilename << std::endl;
    std::cout << "  詳細計算時間: " << detailedTimesFilename << std::endl;

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

