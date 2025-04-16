#ifndef MATRIX_LIB_H
#define MATRIX_LIB_H

void print_matrix(int n, const char *name, double mat[n][n]);

void print_vector(int n, const char *name, double vec[n]);

void swap_rows(int n, double mat[n][n], int row1, int row2);

void swap_int(int *a, int *b);

void PLU(int n, double A[n][n], double L[n][n], double U[n][n], int P[n]);

void print_permutation(int n, const int P[n]);

double det(int n, double U[n][n], int P[n]);

void forward_assign(int n, double L[n][n], double b[n], double y[n]);

void backward_assign(int n, double U[n][n], double y[n], double x[n]);

void apply_permutation(int n, int P[n], double b[n], double pb[n]);

void power_iteration(int n, double A[n][n], double eigenvector[n], double *eigenvalue, int max_iter, double tol);

#endif // MATRIX_LIB_H