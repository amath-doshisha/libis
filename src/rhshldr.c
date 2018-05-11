#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_ivec.h"
#include"is_rvec.h"
#include"is_rmat.h"
#include"is_rhshldr.h"

/**
 @file  rhshldr.c
 @brief 多倍長精度実数型rmultiによるハウスホルダー変換に関する関数の定義.
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */

/**
 @brief ハウスホルダー・ベクトルへの変換.
 @details h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
 @param[in,out]  h   [in]初期化済みのベクトル.サイズはn.[out]ハウスホルダー・ベクトル.
 @param[in]      x   初期化済みのベクトル.サイズはn.
 @param[in]      n   ベクトルのサイズ.
 @param[in]      k   第k要素が基準.
 */
void rhouseholder_vec(int n, int k, rmulti **h, rmulti *alpha, rmulti **x)
{
  int p0,p1,prec;
  rmulti *eta=NULL,*zeta=NULL,*xi=NULL,*axk=NULL;
  // allocate
  p0=rget_prec(alpha);
  p1=rvec_get_prec_max(n,h);
  prec=MAX2(p0,p1);
  eta=rallocate_prec(prec);
  zeta=rallocate_prec(prec);
  xi=rallocate_prec(prec);
  axk=rallocate_prec(prec);
  //----------- norm
  rvec_sum_pow2(xi,n-k-1,&x[k+1]);    // xi=sum(abs(x((k+1):end)).^2);
  rmul_r(axk,x[k],x[k]);                // axk=|x[k]|^2
  radd_r(eta,axk,xi);                   // eta=|x[k]|^2+...
  rsqrt(axk,axk);                     // axk=|x[k]|
  rsqrt(eta,eta);                     // eta=sqrt(|x[k]|^2+...)
  if(req_d(eta,0)){rsub_r(xi,eta,axk);} // xi=eta-|x(k)|
  else{                               // xi=xi/(|x(k)|+eta)
    radd_r(zeta,axk,eta);
    rdiv_r(xi,xi,zeta);
  }
  //----------- h
  rvec_set_zeros(k,h);
  rvec_copy(n-k-1,&h[k+1],&x[k+1]);     // h((k+1):end)=x((k+1):end);
  if(ris_zero(x[k])){
    rcopy(h[k],xi); rneg(h[k],h[k]);    // h[k]=-xi
  }else{
    rdiv_r(zeta,xi,axk); rneg(zeta,zeta); // zeta=-xi/axk;
    rmul_r(h[k],x[k],zeta);               // h[k]=zeta*x[k];
  }
  //----------- alpha
  if(req_d(xi,0) || req_d(eta,0)){
    rset_d(alpha,0);
  }else{
    rmul_r(alpha,xi,eta);                 // alpha=1/(xi*eta)
    rinv(alpha,alpha);
  }
  // free
  eta=rfree(eta);
  zeta=rfree(zeta);
  xi=rfree(xi);
  axk=rfree(axk);
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
void rhouseholder_right(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int k, rmulti **h, rmulti *alpha)
{
  int prec;
  rmulti *a=NULL,**p=NULL;
  // allocate
  prec=rmat_get_prec_max(m,n,B,LDB);
  a=rallocate_prec(prec);
  p=rvec_allocate_prec(m,prec);
  // p=A*h
  rvec_lintr(m,n-k,p,&COL(A,k,LDA),LDA,&h[k]);
  // B=A*H=A-alpha*(A*h)*h'=A-alpha*p*h'
  rneg(a,alpha);
  rmat_rank1op(m,n-k,&COL(B,k,LDB),LDB,&COL(A,k,LDA),LDA,a,p,&h[k]);
  // done
  a=rfree(a);
  p=rvec_free(m,p);
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
void rhouseholder_left(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int k, rmulti **h, rmulti *alpha)
{
  int prec;
  rmulti *a=NULL,**p=NULL;
  // allocate
  prec=rmat_get_prec_max(m,n,B,LDB);
  a=rallocate_prec(prec);
  p=rvec_allocate_prec(m,prec);
  // p=A'*h
  rvec_lintr_t(m-k,n,p,&MAT(A,k,0,LDA),LDA,&h[k]);
  // B=H*A=A-alpha*h*(A'*h)'=A-alpha*h*p'
  rneg(a,alpha);
  rmat_rank1op(m-k,n,&MAT(B,k,0,LDB),LDB,&MAT(A,k,0,LDA),LDA,a,&h[k],p);
  // done
  a=rfree(a);
  p=rvec_free(m,p);
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
void rhouseholder(int n, int k0, int nH, int k, rmulti **h, rmulti *alpha, rmulti **H, int LDH, rmulti **Alpha, rmulti **x)
{
  int j,l,p0,p1,prec;
  rmulti *value=NULL,**R=NULL; 
  // allocate
  p0=rget_prec(alpha);
  p1=rvec_get_prec_max(n,h);
  prec=MAX2(p0,p1);
  value=rallocate_prec(prec);
  R=rvec_allocate_prec(n,prec);
  rvec_copy(n,R,x);
  // R=(I-alpha*h*h')*R=R-alpha*h*(h'*R)    
  for(j=0; j<nH; j++){
    // R=R-alpha[j]*(H(:,j)'*R)*H(:,j)
    l=k0+j;
    rvec_sum_mul(value,n-l,&MAT(H,l,j,LDH),&R[l]);
    rmul_r(value,value,Alpha[j]);
    rvec_sub_mul_r(n-l,&R[l],&MAT(H,l,j,LDH),value);
  }
  rhouseholder_vec(n,k,h,alpha,R);
  // free
  value=rfree(value);
  R=rvec_free(n,R);
}


//EOF
