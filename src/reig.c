#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_ivec.h"
#include"is_rvec.h"
#include"is_rmat.h"
#include"is_reig.h"


/**
 @file  reig.c
 @brief 多倍長精度実数型rmultiによる固有値計算に関する関数の定義.
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */

/**
 @brief 固有値問題の残差 F=A*x-lambda*x.
 */
void reig_residual(int n, rmulti **F, rmulti **A, int LDA, rmulti **x, rmulti *lambda)
{
  rvec_mul_r(n,F,x,lambda);      // F=x*lambda
  rvec_neg(n,F,F);               // F=-x*lambda
  rvec_add_lintr(n,n,F,A,LDA,x); // F=A*x-lambda*x
}

/**
 @brief 固有値問題の全ての残差の最大値ノルム E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))).
 */
void reig_max_abs_residuals(int n, rmulti **E, rmulti **A, int LDA, rmulti **X, int LDX, rmulti **lambda)
{
  int j;
  rmulti **F=NULL;
  F=rvec_allocate_prec(n,rvec_get_prec_max(n,E)); // allocate
  for(j=0; j<n; j++){
    reig_residual(n,F,A,LDA,&COL(X,j,LDX),lambda[j]);
    rvec_max_abs(E[j],n,F);
  }
  F=rvec_free(n,F); // free 
}

/**
 @brief 固有対を固有値の大きい順に並び替える.
 */
void reig_sort(int n, int k, rmulti **lambda, rmulti **X, int LDX)
{
  int *index=NULL;
  // allocate
  index=ivec_allocate(k);
  // sort
  rvec_sort(k,lambda,index);
  ivec_reverse(k,index);
  rvec_reverse(k,lambda);
  rmat_cols_swap_index(n,k,X,LDX,index);
  // free
  index=ivec_free(index);
}

/**
 @brief 固有対を配列Iに従って並び替える.
 */
void reig_sort_index(int n, int k, rmulti **lambda, rmulti **X, int LDX, const int *index)
{
  rvec_swap_index(k,lambda,index);
  rmat_cols_swap_index(n,k,X,LDX,index);
}

/**
 @brief 固有対をベクトルの組X0に従って並び替える.
 */
void reig_sort_vector_guide(int n, int k, rmulti **lambda, rmulti **X, int LDX, rmulti **X0, int LDX0)
{
  int i,j,prec;
  rmulti *value,**a;
  prec=rmat_get_prec_max(n,k,X0,LDX0);
  value=rallocate_prec(prec);
  a=rvec_allocate_prec(k,prec);
  for(j=0; j<k; j++){
    rmat_cols_max_abs_sub_rvec(a,n,k,X,LDX,&COL(X0,j,LDX0));
    rvec_min_abs_index(value,k,a,&i);
    if(i!=j){
      rvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX));
      if(lambda!=NULL){ rswap(lambda[i],lambda[j]); }
    }
  }
  value=rfree(value);
  a=rvec_free(n,a);
}

/**
 @brief 固有対を配列lambda0に従って並び替える.
 */
void reig_sort_value_guide(int n, int k, rmulti **lambda, rmulti **X, int LDX, rmulti **lambda0)
{
  int i,j,prec;
  rmulti *value,**a;
  prec=rvec_get_prec_max(k,lambda0);
  value=rallocate_prec(prec);
  a=rvec_allocate_prec(k,prec);
  for(j=0; j<k; j++){
    rvec_abs_sub_r(k,a,lambda,lambda0[j]);
    rvec_min_abs_index(value,k,a,&i);
    if(i!=j){
      if(X!=NULL){ rvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX)); }
      rswap(lambda[i],lambda[j]);
    }
  }
  value=rfree(value);
  a=rvec_free(n,a);
}

//EOF
