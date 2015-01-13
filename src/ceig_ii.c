#include<stdio.h>
#include<stdlib.h>
#include<isys.h>

/**
 @file  ceig_ii.c
 @brief 多倍長精度実数型cmultiにおける逆反復法を用いて固有値ベクトルを計算.
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
          cmulti型の固有値問題に関する定義は@link ceig.c@endlinkを参照のこと.
 */

/** @name 逆反復法に関する関数 */
/** @{ */

#define NAME_CEIGII       "ceig_ii"
#define NAME_CEIGII_1PAIR "ceig_ii_1pair"

/**
 @brief      逆反復法で行列Aの固有ベクトルXをすべて計算.
 @param[in]  n      行列サイズ.
 @param[in]  A      固有値計算をする(n,n)型行列.
 @param[in]  LDA    Aの第1次元.
 @param[in]  lambda 計算済みの固有値.
 @param[in]  X      初期化済みの(n,n)型行列.
 @param[in]  LDX    Xの第1次元.
 @param[out] X      計算されたすべての固有ベクトル.
 @param[in]  debug  デバッグ.
 */
void ceig_ii(int n, cmulti **X, int LDX, cmulti **A, int LDA, cmulti **lambda, int debug)
{
  int i;
  for(i=0; i<n; i++){
    if(debug>=1){ printf("[%s] %04d/%04d\n",NAME_CEIGII,i,n); }
    ceig_ii_1pair(n,&COL(X,i,LDX),A,LDA,lambda[i],debug-1);
  }
}


/**
 @brief      逆反復法で行列Aの固有ベクトルXを1個のみ計算.
 @param[in]  n      行列サイズ.
 @param[in]  A      固有値計算をする(n,n)型行列.
 @param[in]  LDA    Aの第1次元.
 @param[in]  lambda 計算済みの固有値.
 @param[in]  x      初期化済みのn次ベクトル.
 @param[out] x      計算された固有ベクトル.
 @param[in]  debug  デバッグ.
 */
void ceig_ii_1pair(int n, cmulti **x, cmulti **A, int LDA, cmulti *lambda, int debug)
{
  int prec=53,*p=NULL,info,i;
  cmulti **B=NULL;
  rmulti *eps=NULL;
  // allocate
  p=ivec_allocate(n);
  prec=cvec_get_prec_max(n,x);
  B=cmat_allocate_prec(n,n,prec);
  // LU decomposition
  info=1; i=0;
  while(info){
    // B=A-diag(lambda)
    cmat_diag_sub_c(n,n,B,n,A,LDA,lambda);
    // B=B+diag(i*eps)
    eps=rmepsilon(prec);
    rmul_si(eps,eps,i);
    cmat_diag_add_r(n,n,B,n,B,n,eps);
    eps=rfree(eps);
    // B=L*U
    csolve_lu_decomp(n,B,n,p,&info);
    // next
    if(debug>=1){ printf("[%s] LU: step=%02d\n",NAME_CEIGII_1PAIR,i); }
    i++;
  };
  // compute
  cvec_set_rand(n,x,2,-1);
  cvec_normalize_sgn(n,x,x);
  for(i=0; i<2; i++){
    csolve_lu_backsubs(n,1,x,n,B,n,p);
    cvec_normalize_sgn(n,x,x);
    if(debug>=1){ printf("[%s] II: step=%02d\n",NAME_CEIGII_1PAIR,i); }
  }
  // done
  p=ivec_free(p);
  B=cmat_free(n,n,B);
}

/** @} */


//EOF
