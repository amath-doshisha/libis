#include<stdio.h>
#include<stdlib.h>
#include<isys.h>

/**
 @file  reig_ii.c
 @brief 多倍長精度実数型rmultiにおける逆反復法を用いて固有値ベクトルを計算.
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
          rmulti型の固有値問題に関する定義は@link reig.c@endlinkを参照のこと.
 */

/** @name 逆反復法に関する関数 */
/** @{ */

#define NAME_REIGII       "reig_ii"
#define NAME_REIGII_1PAIR "reig_ii_1pair"

/**
 @brief      逆反復法で行列Aの固有ベクトルXをすべて計算.
 @param[in]      n      行列サイズ.
 @param[in]      A      固有値計算をする(n,n)型行列.
 @param[in]      LDA    Aの第1次元.
 @param[in]      lambda 計算済みの固有値.
 @param[in,out]  X      [in]初期化済みの(n,n)型行列.[out]計算されたすべての固有ベクトル.
 @param[in]      LDX    Xの第1次元.
 @param[in]      debug  デバッグ.
 */
void reig_ii(int n, rmulti **X, int LDX, rmulti **A, int LDA, rmulti **lambda, int debug)
{
  int i;
  for(i=0; i<n; i++){
    if(debug>=1){ printf("[%s] %04d/%04d\n",NAME_REIGII,i,n); }
    reig_ii_1pair(n,&COL(X,i,LDX),A,LDA,lambda[i],debug-1);
  }
}


/**
 @brief      逆反復法で行列Aの固有ベクトルXを1個のみ計算.
 @param[in]      n      行列サイズ.
 @param[in]      A      固有値計算をする(n,n)型行列.
 @param[in]      LDA    Aの第1次元.
 @param[in]      lambda 計算済みの固有値.
 @param[in,out]  x      [in]初期化済みのn次ベクトル.[out]計算された固有ベクトル.
 @param[in]      debug  デバッグ.
 */
void reig_ii_1pair(int n, rmulti **x, rmulti **A, int LDA, rmulti *lambda, int debug)
{
  int prec=53,*p=NULL,info,i;
  rmulti **B=NULL,*eps=NULL;
  // allocate
  p=ivec_allocate(n);
  prec=rvec_get_prec_max(n,x);
  B=rmat_allocate_prec(n,n,prec);
  // LU decomposition
  info=1; i=0;
  while(info){
    // B=A-diag(lambda)
    rmat_diag_sub_r(n,n,B,n,A,LDA,lambda);
    // B=B+diag(i*eps)
    eps=rmepsilon(prec);
    rmul_rd(eps,eps,i);
    rmat_diag_add_r(n,n,B,n,B,n,eps);
    eps=rfree(eps);
    // B=L*U
    rsolve_lu_decomp(n,B,n,p,&info);
    // next
    if(debug>=1){ printf("[%s] LU: step=%02d\n",NAME_REIGII_1PAIR,i); }
    i++;
  }
  // compute
  rvec_set_rand(n,x,2,-1);
  rvec_normalize_sgn_rvec(n,x,x);
  for(i=0; i<2; i++){
    rsolve_lu_backsubs(n,1,x,n,B,n,p);
    rvec_normalize_sgn_rvec(n,x,x);
    if(debug>=1){ printf("[%s] II: step=%02d\n",NAME_REIGII_1PAIR,i); }
  }
  // done
  p=ivec_free(p);
  B=rmat_free(n,n,B);
}

/** @} */


//EOF
