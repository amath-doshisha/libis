#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_rvec.h"
#include"is_cvec.h"
#include"is_cmat.h"
#include"is_ceig.h"

/**
 @file  ceig.c
 @brief 多倍長精度実数型cmultiによる固有値計算に関する関数の定義.
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
 */

/**
 @brief 固有値問題の残差 F=A*x-lambda*x.
 */
void ceig_residual(int n, cmulti **F, cmulti **A, int LDA, cmulti **x, cmulti *lambda)
{
  cvec_mul_cscalar(n,F,x,lambda);      // F=x*lambda
  cvec_neg(n,F,F);               // F=-x*lambda
  cvec_add_lintr(n,n,F,A,LDA,x); // F=A*x-lambda*x
}

/**
 @brief 固有値問題の全ての残差の最大値ノルム E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))).
 */
void ceig_residual_norm_max(int n, rmulti **E, cmulti **A, int LDA, cmulti **X, int LDX, cmulti **lambda)
{
  int j;
  cmulti **F=NULL;
  F=cvec_allocate_prec(n,rvec_get_prec_max(n,E));
  for(j=0; j<n; j++){
    ceig_residual(n,F,A,LDA,&COL(X,j,LDX),lambda[j]);
    cvec_max_abs(E[j],n,F);
  }
  F=cvec_free(n,F); // free  
}

/**
 @brief 固有対を固有値の大きい順に並び替える.
 */
void ceig_sort(int n, int k, cmulti **lambda, cmulti **X, int LDX)
{
  int *index=NULL;
  // allocate
  index=ivec_allocate(k);
  // sort
  cvec_sort(k,lambda,index);
  ivec_reverse(k,index);
  cvec_reverse(k,lambda);
  cmat_cols_swap_index(n,k,X,LDX,index);
  // free
  index=ivec_free(index);
}

/**
 @brief 固有対を配列Iに従って並び替える.
 */
void ceig_sort_index(int n, int k, cmulti **lambda, cmulti **X, int LDX, const int *index)
{
  cvec_swap_index(k,lambda,index);
  cmat_cols_swap_index(n,k,X,LDX,index);
}

/**
 @brief 固有対をベクトルの組X0に従って並び替える.
 */
void ceig_sort_vector_guide(int n, int k, cmulti **lambda, cmulti **X, int LDX, cmulti **X0, int LDX0)
{
  int i,j,prec;
  rmulti *value,**a;
  prec=cmat_get_prec_max(n,k,X0,LDX0);
  value=rallocate_prec(prec);
  a=rvec_allocate_prec(k,prec);
  for(j=0; j<k; j++){
    cmat_cols_max_abs_sub_cvec(a,n,k,X,LDX,&COL(X0,j,LDX0));
    rvec_min_abs_index(value,k,a,&i);
    if(i!=j){
      cvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX));
      if(lambda!=NULL){ cswap(lambda[i],lambda[j]); }
    }
  }
  value=rfree(value);
  a=rvec_free(n,a);
}

/**
 @brief 固有対を配列lambda0に従って並び替える.
 */
void ceig_sort_value_guide(int n, int k, cmulti **lambda, cmulti **X, int LDX, cmulti **lambda0)
{
  int i,j,prec;
  rmulti *value,**a;
  prec=cvec_get_prec_max(k,lambda0);
  value=rallocate_prec(prec);
  a=rvec_allocate_prec(k,prec);
  for(j=0; j<k; j++){
    cvec_abs_sub_c(k,a,lambda,lambda0[j]);
    rvec_min_abs_index(value,k,a,&i);
    if(i!=j){
      if(X!=NULL){ cvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX)); }
      cswap(lambda[i],lambda[j]);
    }
  }
  value=rfree(value);
  a=rvec_free(n,a);
}

//EOF
