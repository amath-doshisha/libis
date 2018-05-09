#ifndef IS_ZSOLVE_H
#define IS_ZSOLVE_H

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int zsolve(int n, int NRHS, dcomplex *B, int LDB, dcomplex *A, int LDA);

// R=B-A*X
void zsolve_residual(int n, int NRHS, dcomplex *R, int LDR, dcomplex *A, int LDA, dcomplex *X, int LDX, dcomplex *B, int LDB);

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int zsolve_gauss_LU(int n, int NRHS, dcomplex *B, int LDB, dcomplex *A, int LDA);

// Sove A*X=B by Gauss sweeper
// X: input=B, output=X
// A: input=A, output=destroyed
int zsolve_gauss_sweeper(int n, int NRHS, dcomplex *B, int LDB, dcomplex *A, int LDA);

#endif
