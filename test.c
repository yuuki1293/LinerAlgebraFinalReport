#include <stdio.h>
#include "matrix_lib.h"

#define N 3

int main(){
    double A[N][N] = {
        {2, 1, 1},
        {4, -6, 4},
        {-2, 7, 2}};
    double b[N] = {5, -2, 9};

    double L[N][N], U[N][N];
    int P[N];

    PLU(N, A, L, U, P);
    double d = det(N, U, P);

    print_matrix(N, "L", L);
    print_matrix(N, "U", U);
    print_permutation(N, P);
    printf("det = %f\n", d);

    double pb[N], y[N], x[N];
    apply_permutation(N, P, b, pb);
    forward_assign(N, L, pb, y);
    backward_assign(N, U, y, x);

    printf("解 x =\n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.4f\n", i, x[i]);
    }

    double eigenvector[N], eigenvalue;
    power_iteration(N, A, eigenvector, &eigenvalue, 10000, 1);

    printf("固有値 = %f\n", eigenvalue);
    print_vector(N, "固有ベクトル", eigenvector);
    return 0;
}
