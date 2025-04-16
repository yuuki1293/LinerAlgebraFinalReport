#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_lib.h"

void print_matrix(int n, const char *name, double mat[n][n])
{
    printf("%s =\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%8.4f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(int n, const char *name, double vec[n])
{
    printf("%s = [", name);
    for (int i = 0; i < n; i++)
    {
        printf("%8.4f", vec[i]);
        if (i < n - 1)
            printf(", ");
    }
    printf("]\n\n");
}

void swap_rows(int n, double mat[n][n], int row1, int row2)
{
    for (int j = 0; j < n; j++)
    {
        double tmp = mat[row1][j];
        mat[row1][j] = mat[row2][j];
        mat[row2][j] = tmp;
    }
}

void swap_int(int *a, int *b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

void PLU(int n, double A[n][n], double L[n][n], double U[n][n], int P[n])
{
    // 初期化
    for (int i = 0; i < n; i++)
    {
        P[i] = i;
        for (int j = 0; j < n; j++)
        {
            L[i][j] = (i == j) ? 1.0 : 0.0;
            U[i][j] = A[i][j];
        }
    }

    for (int k = 0; k < n; k++)
    {
        // ピボット選択（最大絶対値の行を探す）
        int pivot = k;
        double max = fabs(U[k][k]);
        for (int i = k + 1; i < n; i++)
        {
            if (fabs(U[i][k]) > max)
            {
                max = fabs(U[i][k]);
                pivot = i;
            }
        }

        // 行を入れ替える（U, P, Lの一部）
        if (pivot != k)
        {
            swap_rows(n, U, k, pivot);
            swap_int(&P[k], &P[pivot]);

            // Lの左側（すでに埋まった部分）も入れ替える
            for (int j = 0; j < k; j++)
            {
                double tmp = L[k][j];
                L[k][j] = L[pivot][j];
                L[pivot][j] = tmp;
            }
        }

        // Uの下を消去して、Lに係数を入れる
        for (int i = k + 1; i < n; i++)
        {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; j++)
            {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }
}

void print_permutation(int n, const int P[n])
{
    printf("P =\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d ", (P[i] == j) ? 1 : 0);
        }
        printf("\n");
    }
    printf("\n");
}

double det(int n, double U[n][n], int P[n])
{
    double detU = 1.0;
    for (int i = 0; i < n; i++)
        detU *= U[i][i];

    // 置換のパリティ（偶奇）を数える
    int sign = 1;
    int visited[n];
    for (int i = 0; i < n; i++)
    {
        visited[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        if (visited[i])
            continue;
        int cycle_size = 0;
        int j = i;
        while (!visited[j])
        {
            visited[j] = 1;
            j = P[j];
            cycle_size++;
        }
        if (cycle_size % 2 == 0)
            sign *= -1;
    }

    return sign * detU;
}

void forward_assign(int n, double L[n][n], double b[n], double y[n])
{
    for (int i = 0; i < 3; i++)
    {
        y[i] = b[i];
        for (int j = 0; j < i; j++)
        {
            y[i] -= L[i][j] * y[j];
        }
    }
}

void backward_assign(int n, double U[n][n], double y[n], double x[n])
{
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++)
        {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

void apply_permutation(int n, int P[n], double b[n], double pb[n])
{
    for (int i = 0; i < n; i++)
    {
        pb[i] = b[P[i]];
    }
}

void power_iteration(int n, double A[n][n], double eigenvector[n], double *eigenvalue, int max_iter, double tol)
{
    double b_k[n];
    for (int i = 0; i < n; i++)
        b_k[i] = 1.0; // 初期ベクトルを1で初期化

    for (int iter = 0; iter < max_iter; iter++)
    {
        // A * b_k を計算
        double b_k1[n];
        for (int i = 0; i < n; i++)
        {
            b_k1[i] = 0.0;
            for (int j = 0; j < n; j++)
            {
                b_k1[i] += A[i][j] * b_k[j];
            }
        }

        // 正規化
        double norm = 0.0;
        for (int i = 0; i < n; i++)
            norm += b_k1[i] * b_k1[i];
        norm = sqrt(norm);

        for (int i = 0; i < n; i++)
            b_k1[i] /= norm;

        // 収束判定
        double diff = 0.0;
        for (int i = 0; i < n; i++)
            diff += fabs(b_k1[i] - b_k[i]);

        if (diff < tol)
        {
            for (int i = 0; i < n; i++)
                eigenvector[i] = b_k1[i];

            *eigenvalue = 0.0;
            for (int i = 0; i < n; i++)
            {
                double temp = 0.0;
                for (int j = 0; j < n; j++)
                    temp += A[i][j] * b_k1[j];
                *eigenvalue += temp * b_k1[i];
            }
            return;
        }

        for (int i = 0; i < n; i++)
            b_k[i] = b_k1[i];
    }

    printf("Power iteration did not converge within the maximum number of iterations.\n");
}