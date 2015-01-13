#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<isys.h>

/**
 @file  chsnbrg.c
 @brief 多倍長精度実数型cmultiによるヘッセンバーグ型行列に関する関数の定義.
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
 */

/**
 @brief    密行列AをHessenberg型行列Bに相似変換する.
 @param[in]  n    正方行列A,Bのサイズ
 @param[in]  A    行列
 @param[in]  LDA  Aの第1次元
 @param[out] B    相似変換により得られたHessenberg型行列
 @param[in]  LDB  Bの第1次元
 */
void chsnbrg_simtr(int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,prec=53;
  rmulti *alpha=NULL;
  cmulti **h=NULL;
  // allocate
  prec=cmat_get_prec_max(n,n,B,LDB);
  alpha=rallocate_prec(prec);
  h=cvec_allocate_prec(n,prec);
  // copy
  cmat_copy(n,n,B,LDB,A,LDA);
  // compute
  for(i=0; i<n-2; i++){
    // Householder vector
    chouseholder_vec(n,i+1,h,alpha,&COL(B,i,LDB));
    // similarity transformation
    chouseholder_right(n,n,B,LDB,B,LDB,i+1,h,alpha);
    chouseholder_left(n,n,B,LDB,B,LDB,i+1,h,alpha);
    // set lower part as zeros
    cvec_set_zeros(n-2-i,&MAT(B,i+2,i,LDB));
  }
  // done
  alpha=rfree(alpha);
  h=cvec_free(n,h);
}

/**
 @brief    多倍長Hessenbergに相似変換する(mステップで終了)
 @details  Aは書き換え,cmat型
 @param[in]  n    行列サイズ
 @param[in]  A    行列
 @param[in]  LDA 
 @param[in]  step 相似変換ステップ数(m-2が上限)
 @param[out] A    hessenberg行列
 */
void chsnbrg_simtr_step(int n, cmulti **A, int LDA, int step)
{
  int prec,i;
  rmulti *alpha=NULL;
  cmulti **h=NULL;
  if(step<=0 || step>(n-2)){ printf("ERROR! Necessary step =< n-2 \n"); exit(0); }
  // allocate
  prec=cmat_get_prec_max(n,n,A,LDA);
  alpha=rallocate_prec(prec);
  h=cvec_allocate_prec(n,prec);
  // compute
  for(i=0; i<step; i++){
    chouseholder_vec(n,i+1,h,alpha,&COL(A,i,LDA));
    chouseholder_right(n,n,A,LDA,A,LDA,i+1,h,alpha);
    chouseholder_left(n,n,A,LDA,A,LDA,i+1,h,alpha);
  }
  // done
  alpha=rfree(alpha);
  h=cvec_free(n,h);
}

// EOF
