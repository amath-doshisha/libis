#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<isys.h>

/**
 @file  rhsnbrg.c
 @brief 多倍長精度実数型rmultiによるヘッセンバーグ型行列に関する関数の定義.
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */


/**
 @brief    密行列AをHessenberg型行列Bに相似変換する.
 @param[in]  n    正方行列A,Bのサイズ
 @param[in]  A    行列
 @param[in]  LDA  Aの第1次元
 @param[out] B    相似変換により得られたHessenberg型行列
 @param[in]  LDB  Bの第1次元
 */
void rhsnbrg_simtr(int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,prec=53;
  rmulti *alpha=NULL,**h=NULL;
  // allocate
  prec=rmat_get_prec_max(n,n,B,LDB);
  alpha=rallocate_prec(prec);
  h=rvec_allocate_prec(n,prec);
  // copy
  rmat_copy(n,n,B,LDB,A,LDA);
  // compute
  for(i=0; i<n-2; i++){
    // Householder vector
    rhouseholder_vec(n,i+1,h,alpha,&COL(B,i,LDB));
    // similarity transformation
    rhouseholder_right(n,n,B,LDB,B,LDB,i+1,h,alpha);
    rhouseholder_left(n,n,B,LDB,B,LDB,i+1,h,alpha);
    // set lower part as zeros
    rvec_set_zeros(n-2-i,&MAT(B,i+2,i,LDB));
  }
  // done
  alpha=rfree(alpha);
  h=rvec_free(n,h);
}



//EOF
