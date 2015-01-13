#ifndef IS_REIG_H
#define IS_REIG_H

#include<is_rmulti.h>

/**
 @file  is_reig.h
 @brief 多倍長精度実数型rmultiによる固有値計算に関する関数の宣言.
 */


// F=A*x-lambda*x
void reig_residual(int n, rmulti **F, rmulti **A, int LDA, rmulti **x, rmulti *lambda);
// E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))), k=1,2,..,n
void reig_max_abs_residuals(int n, rmulti **E, rmulti **A, int LDA, rmulti **X, int LDX, rmulti **lambda);
// sort by eigenvalues
void reig_sort(int n, int k, rmulti **lambda, rmulti **X, int LDX);
// sort by indecies
void reig_sort_index(int n, int k, rmulti **lambda, rmulti **X, int LDX, const int *index);
// sort by another vectors
void reig_sort_vector_guide(int n, int k, rmulti **lambda, rmulti **X, int LDX, rmulti **X0, int LDX0);
// sort by another values
void reig_sort_value_guide(int n, int k, rmulti **lambda, rmulti **X, int LDX, rmulti **lambda0);

#endif
