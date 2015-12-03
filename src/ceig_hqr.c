#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<isys.h>

/**
 @file  ceig_hqr.c
 @brief 多倍長精度実数型cmultiによるHessenberg型行列に前処理をした後にQR法で固有値を計算.
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
          cmulti型の固有値問題に関する定義は@link ceig.c@endlinkを参照のこと.
 */


/** @name 固有値計算に関する関数 */
/** @{ */

/**
 @brief   Hessenberg型行列に前処理をした後にQR法で固有値を計算.
 @param[in]  A      固有値計算をする(n,n)型行列.
 @param[in]  n      行列サイズ.
 @param[in]  LDA    Aの第1次元.
 @param[in]  debug  デバッグ.
 @param[out] lambda 計算された固有値を返す.
 */
void ceig_hqr(int n, cmulti **lambda, cmulti **A, int LDA, int debug)
{
  rmulti **B=NULL;
  B=rmat_allocate_prec(n,n,cmat_get_prec_max(n,n,A,LDA));
  cmat_real_clone(n,n,B,n,A,LDA);
  reig_hqr(n,lambda,B,n,debug);
  B=rmat_free(n,n,B);
}

/** @} */

/** @name 固有値計算に関する関数 */
/** @{ */

/**
 @brief   Hessenberg型行列に前処理をした後にQR法で固有値を計算.
 @param[in]  A      固有値計算をする(n,n)型行列.
 @param[in]  n      行列サイズ.
 @param[in]  LDA    Aの第1次元.
 @param[in]  debug  デバッグ.
 @param[out] B      QR法の相似変換後の行列.
 @param[out] lambda 計算された固有値を返す.
 */
void ceig_hqr_mt(int m, int n, cmulti **B, int LDB, cmulti **lambda, cmulti **A, int LDA, int debug)
{
  rmulti **C=NULL;
  C=rmat_allocate_prec(m,n,cmat_get_prec_max(m,n,B,LDB));
  cmat_real_clone(m,n,C,LDB,A,LDA);
  reig_hqr_mt(m,n,C,LDB,lambda,C,LDB,debug);
  cmat_copy_rmat(m,n,B,LDB,C,LDB);
  C=rmat_free(m,n,C);
}

/** @} */


//EOF
