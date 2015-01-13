#ifndef IS_CSOLVE_H
#define IS_CSOLVE_H

#include<is_cmulti.h>

/**
 @file  is_csolve.h
 @brief 多倍長精度実数型cmultiによる線形方程式A*X=Bの解法に関する関数の宣言.
 */

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int csolve(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *info);

// R=B-A*X
int csolve_residual(int n, int NRHS, cmulti **R, int LDR, cmulti **A, int LDA, cmulti **X, int LDX, cmulti **B, int LDB);

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int csolve_lu(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *info);
int csolve_lu_decomp(int n, cmulti **A, int LDA, int *p, int *info);
int csolve_lu_backsubs(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *p);

// Sove A*X=B by Gauss sweeper
// X: input=B, output=X
// A: input=A, output=destroyed
int csolve_gauss_sweeper(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *info);

#endif
