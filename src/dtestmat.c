#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_dmat.h"
#include"is_dtestmat.h"

/**
 @file  dmat.c
 @brief 倍精度の行列のテスト行列に関する関数の定義
 */

/** @name 倍精度の行列のテスト行列に関する関数 */
/** @{ */

/**
 @brief  倍精度のToeplitz型行列の新規生成.
 @param[in]  m      行列Aの行の個数.
 @param[in]  n      行列Aの列の個数.
 @param[in]  A      初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA    Aの第1次元.
 @param[in]  k      配列aのサイズ.
 @param[in]  offset 配置場所.
 @param[out] A      Toeplitz型行列.
*/
void dmat_toeplitz(int m, int n, double *A, int LDA, int k, double *a, int offset)
{
  int i;
  dmat_set_zeros(m,n,A,LDA);
  for(i=0; i<k; i++){
    dmat_diag_set_scalar(m,n,A,LDA,a[i],i-offset);
  }
}


/**
 @brief  倍精度のFrank行列の新規生成.
 @param[in]  m      行列Aの行の個数.
 @param[in]  n      行列Aの列の個数.
 @param[in]  A      初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA    Aの第1次元.
 @param[in]  param  パラメータ.
 @param[out] A      Frank型行列.
*/
void dmat_frank(int m, int n, double *A, int LDA, int param)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(param){
	if(i==j+1)    { MAT(A,i,j,LDA)=i; }
	else if(i<j+1){ MAT(A,i,j,LDA)=i+1; }
	else          { MAT(A,i,j,LDA)=0; }
      }else{
	if(i==j+1)    { MAT(A,i,j,LDA)=m-j-1; }
	else if(i<j+1){ MAT(A,i,j,LDA)=m-j; }
	else          { MAT(A,i,j,LDA)=0; }
      }
    }
  }
}

double dfact(double x)
{
  if(x<=1){ return 1.0; }
  else { return x*dfact(x-1); }
}

/**
 @brief  倍精度のipjfact行列の新規生成.
 @param[in]  m      行列Aの行の個数.
 @param[in]  n      行列Aの列の個数.
 @param[in]  A      初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA    Aの第1次元.
 @param[in]  param  パラメータ.
 @param[out] A      ipjfact型行列.
*/
void dmat_ipjfact(int m, int n, double *A, int LDA, int param)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=dfact((i+1)+(j+1));
      if(param){ MAT(A,i,j,LDA)=1.0/MAT(A,i,j,LDA); }
    }
  }
}


/** @} */

//EOF
