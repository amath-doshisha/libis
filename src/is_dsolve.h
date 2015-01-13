#ifndef IS_DSOLVE_H
#define IS_DSOLVE_H

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int dsolve(int n, int NRHS, double *B, int LDB, double *A, int LDA);

// R=B-A*X
void dsolve_residual(int n, int NRHS, double *R, int LDR, const double *A, int LDA, const double *X, int LDX, const double *B, int LDB);

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int dsolve_gauss_LU(int n, int NRHS, double *B, int LDB, double *A, int LDA);

// Sove A*X=B by Gauss sweeper
// X: input=B, output=X
// A: input=A, output=destroyed
int dsolve_gauss_sweeper(int n, int NRHS, double *B, int LDB, double *A, int LDA);

#endif
