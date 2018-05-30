#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_ivec.h"
#include"is_cvec.h"
#include"is_cmat.h"
#include"is_chshldr.h"

/**
 @file  chshldr.c
 @brief 多倍長精度実数型cmultiによるハウスホルダー変換に関する関数の定義.
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
 */

/**
 @brief ハウスホルダー・ベクトルへの変換.
 @details h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
 @param[in,out]  h [in] 初期化済みのベクトル.サイズはn. [out] ハウスホルダー・ベクトル.
 @param[in]  x 初期化済みのベクトル.サイズはn.
 @param[in]  n ベクトルのサイズ.
 @param[in]  k 第k要素が基準.
 */
void chouseholder_vec(int n, int k, cmulti **h, rmulti *alpha, cmulti **x)
{
  int p0,p1,prec;
  rmulti *eta=NULL,*zeta=NULL,*xi=NULL,*axk=NULL;
  // allocate
  p0=rget_prec(alpha);
  p1=cvec_get_prec_max(n,h);
  prec=MAX2(p0,p1);
  eta=rallocate_prec(prec);
  zeta=rallocate_prec(prec);
  xi=rallocate_prec(prec);
  axk=rallocate_prec(prec);
  //----------- norm
  rsum_pow2_abs_cvec(xi,n-k-1,&x[k+1]);     // xi=sum(abs(x((k+1):end)).^2);
  rabs2_c(axk,x[k]);                     // axk=|x[k]|^2
  radd_rr(eta,axk,xi);                    // eta=|x[k]|^2+...
  rsqrt_r(axk,axk);                      // axk=|x[k]|
  rsqrt_r(eta,eta);                      // eta=sqrt(|x[k]|^2+...)
  if(eq_rd(eta,0)){rsub_rr(xi,eta,axk);}  // xi=eta-|x(k)|
  else{                                // xi=xi/(|x(k)|+eta)
    radd_rr(zeta,axk,eta);
    rdiv_rr(xi,xi,zeta);
  }
  //----------- h
  cvec_set_zeros(k,h);
  cvec_copy_cvec(n-k-1,&h[k+1],&x[k+1]);      // h((k+1):end)=x((k+1):end);
  if(cis_zero(x[k])){
    cset_r(h[k],xi); cneg_c(h[k],h[k]);        //h[k]=-xi
  }else{
    rdiv_rr(zeta,xi,axk); rneg_r(zeta,zeta);    // zeta=-xi/axk;
    cmul_cr(h[k],x[k],zeta);                // h[k]=zeta*x[k];
  }
  //----------- alpha
  if(eq_rd(xi,0) || eq_rd(eta,0)){
    rset_d(alpha,0);
  }else{
    rmul_rr(alpha,xi,eta);                  // alpha=1/(xi*eta)
    rinv_r(alpha,alpha);
  }
  // free
  eta=rmfree(eta);
  zeta=rmfree(zeta);
  xi=rmfree(xi);
  axk=rmfree(axk);
}

/**
 @brief ハウスホルダー変換を右から作用 B=A*H, H=I-alpha*h*h'
 @details h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
 @param[in]  m     行列Aの行の個数.
 @param[in]  n     行列Aの列の個数.
 @param[in]  k     第k要素が基準.
 @param[in]  h     ハウスホルダー・ベクトル.
 @param[in]  alpha ハウスホルダー・ベクトルの規格化定数.
 @param[out] B     変換された行列B.
*/
void chouseholder_right(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int k, cmulti **h, rmulti *alpha)
{
  int prec;
  cmulti *a=NULL,**p=NULL;
  // allocate
  prec=cmat_get_prec_max(m,n,B,LDB);
  a=callocate_prec(prec);
  p=cvec_allocate_prec(m,prec);
  // p=A*h
  cvec_mul_cmat_cvec(m,n-k,p,&COL(A,k,LDA),LDA,&h[k]);
  // B=A*H=A-alpha*(A*h)*h'=A-alpha*p*h'
  cset_r(a,alpha); cneg_c(a,a);
  cmat_rank1op(m,n-k,&COL(B,k,LDB),LDB,&COL(A,k,LDA),LDA,a,p,&h[k]);
  // done
  a=cmfree(a);
  p=cvec_free(m,p);
}

/**
 @brief ハウスホルダー変換を左から作用 B=H*A, H=I-alpha*h*h'
 @details h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
 @param[in]  m     行列Aの行の個数.
 @param[in]  n     行列Aの列の個数.
 @param[in]  k     第k要素が基準.
 @param[in]  h     ハウスホルダー・ベクトル.
 @param[in]  alpha ハウスホルダー・ベクトルの規格化定数.
 @param[out] B     変換された行列B.
*/
void chouseholder_left(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int k, cmulti **h, rmulti *alpha)
{
  int prec;
  cmulti *a=NULL,**p=NULL;
  // allocate
  prec=cmat_get_prec_max(m,n,B,LDB);
  a=callocate_prec(prec);
  p=cvec_allocate_prec(m,prec);
  // p=A'*h
  cvec_lintr_ct(m-k,n,p,&MAT(A,k,0,LDA),LDA,&h[k]);
  // B=H*A=A-alpha*h*(A'*h)'=A-alpha*h*p'
  cset_r(a,alpha); cneg_c(a,a);
  cmat_rank1op(m-k,n,&MAT(B,k,0,LDB),LDB,&MAT(A,k,0,LDA),LDA,a,&h[k],p);
  // done
  a=cmfree(a);
  p=cvec_free(m,p);
}

/**
 @brief ハウスホルダー変換によるQR分解の変則版.
 @param[in]  n     ベクトルxのサイズ.ベクトルhのサイズ.
 @param[in]  k0    Hの最初のハウスホルダー・ベクトルの基準は第k0要素.
 @param[in]  nH    配列Alphaのサイズ.行列Hの列の個数.
 @param[in]  k     第k要素が基準.
 @param[in]  H     ハウスホルダー・ベクトルの格納用の行列.サイズは(n,nH).
 @param[in]  Alpha ハウスホルダー・ベクトルの規格化定数の格納用の配列.サイズはnH.
 @param[in]  x     変換されるx.
 @param[out] h     生成されたハウスホルダー・ベクトル.
 @param[out] alpha 生成されたハウスホルダー・ベクトルの規格化定数.
*/
void chouseholder(int n, int k0, int nH, int k, cmulti **h, rmulti *alpha, cmulti **H, int LDH, rmulti **Alpha, cmulti **x)
{
  int j,l,p0,p1,prec;
  cmulti *value=NULL,**R=NULL; 
  // allocate
  p0=rget_prec(alpha);
  p1=cvec_get_prec_max(n,h);
  prec=MAX2(p0,p1);
  value=callocate_prec(prec);
  R=cvec_allocate_prec(n,prec);
  cvec_copy_cvec(n,R,x);
  // R=(I-alpha*h*h')*R=R-alpha*h*(h'*R)    
  for(j=0; j<nH; j++){
    // R=R-alpha[j]*(H(:,j)'*R)*H(:,j)
    l=k0+j;
    cinnprod_cvec_cvec(value,n-l,&MAT(H,l,j,LDH),&R[l]);
    cmul_cr(value,value,Alpha[j]);
    cvec_sub_mul_cvec_cscalar(n-l,&R[l],&MAT(H,l,j,LDH),value);
  }
  chouseholder_vec(n,k,h,alpha,R);
  // free
  value=cmfree(value);
  R=cvec_free(n,R);
}


//////////////////////////////////////////////////

/** @name QR分解に関する関数 */
/** @{ */


/**
 @brief   A=Q*R分解でRを生成.
 @param[in]  n      Aの行列サイズ
 @param[in]  A      行列
 @param[in]  LDA    行列の第１次元．
 @param[in,out]  R      [in]初期化済み行列、サイズ(n,n)．[out]QR分解後の行列R．
 @param[in,out]  H      [out]初期化済み行列、サイズ(n,n-1)．[out]householderベクトル(行列の形)
 @param[in,out]  alpha  [in]初期化済みベクトル、サイズ(m)．[out]householderベクトルの規格化定数．
 */
void chouseholder_qr(int n, cmulti **R, int LDR, cmulti **H, rmulti **alpha, cmulti **A, int LDA)
{
  int i;
  cmat_copy(n,n,R,LDR,A,LDA);
  for(i=0; i<n-1; i++){
    chouseholder_vec(n,i,&COL(H,i,n),alpha[i],&COL(R,i,LDR));
    chouseholder_left(n,n,R,LDR,R,LDR,i,&COL(H,i,n),alpha[i]);
  }
}

/**
 @brief   QR分解で相似変換 (A2=Q'A1 Q)
 @details Aは書き換え、nステップごとにストップ可能、householder変換後行列を確認できる
 @param[in]     n      Aの行列サイズ
 @param[in,out] A      [in]行列．[out]相似変換後の行列．
 @param[in,out] LDA    行列Aの第1次元.
 @param[in]     H      サイズ(n,step)の行列. Householderベクトル(行列の形)
 @param[in]     alpha  サイズがstepのベクトル. Householderベクトルの規格化定数
 */
void chouseholder_qr_nstep(int n, cmulti **A, int LDA, cmulti **H, rmulti **alpha, int step)
{
  int i;
  for(i=0; i<step; i++){
    chouseholder_right(n,n,A,LDA,A,LDA,i,&COL(H,i,n),alpha[i]);
    chouseholder_left(n,n,A,LDA,A,LDA,i,&COL(H,i,n),alpha[i]);
  }
}

/** @} */


//EOF
