#ifndef IS_CEIG_H
#define IS_CEIG_H

#include<is_cmulti.h>

/**
 @file  is_ceig.h
 @brief 多倍長精度実数型cmultiによる固有値計算に関する関数の宣言.
 */

// F=A*X-lambda*X
void ceig_residual(int n, cmulti **F, cmulti **A, int LDA, cmulti **x, cmulti *lambda);
// E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))), k=1,2,..,n
void ceig_residual_norm_max(int n, rmulti **E, cmulti **A, int LDA, cmulti **X, int LDX, cmulti **lambda);
// sort by eigenvalues
void ceig_sort(int n, int k, cmulti **lambda, cmulti **X, int LDX);
// sort by indecies
void ceig_sort_index(int n, int k, cmulti **lambda, cmulti **X, int LDX, const int *index);
// sort by another vectors
void ceig_sort_vector_guide(int n, int k, cmulti **lambda, cmulti **X, int LDX, cmulti **X0, int LDX0);
// sort by another values
void ceig_sort_value_guide(int n, int k, cmulti **lambda, cmulti **X, int LDX, cmulti **lambda0);

#endif
