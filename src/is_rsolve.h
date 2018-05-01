#ifndef IS_RSOLVE_H
#define IS_RSOLVE_H

#include<is_rmulti.h>

/**
 @file  is_rsolve.h
 @brief 多倍長精度実数型rmultiによる線形方程式A*X=Bの解法に関する関数の宣言.
 */


// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
void rsolve(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *info);

// R=B-A*X
void rsolve_residual(int n, int NRHS, rmulti **R, int LDR, rmulti **A, int LDA, rmulti **X, int LDX, rmulti **B, int LDB);

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
void rsolve_lu(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *info);
void rsolve_lu_decomp(int n, rmulti **A, int LDA, int *p, int *info);
void rsolve_lu_backsubs(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *p);

// Sove A*X=B by Gauss sweeper
// X: input=B, output=X
// A: input=A, output=destroyed
void rsolve_gauss_sweeper(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *info);

#endif
