#include"is_rmulti.h"
#include"is_rmat.h"
#include"is_rsolve.h"
#include"is_reig.h"
#include"is_irmulti.h"
#include"is_irvec.h"
#include"is_irmat.h"
#include"is_ireig.h"
#include"is_print.h"

/**
 @file  ireig.c
 @brief 多倍長精度実数型rmultiの機械区間演算の固有値計算に関する関数の定義.
 @details スカラーに関しては@link irmulti.c@endlinkを参照のこと.
          ベクトルに関しては@link irvec.c@endlinkを参照のこと.
          行列に関しては@link irmat.c@endlinkを参照のこと.
 */

#define RMA(A,M,N,P) { A=rmat_allocate_prec(M,N,P); }
#define RMF(A,M,N)   { A=rmat_free(M,N,A); }
#define RVA(X,N,P)   { X=rvec_allocate_prec(N,P); }
#define RVF(X,N)     { X=rvec_free(N,X); }
#define RA(X,P)      { X=rallocate_prec(P); }
#define RF(X)        { X=rfree(X); }

///////////////////////////////////////

/** @name 固有値問題 */
/** @{ */

/**
 @brief 固有値問題の残差 F=A*X-lambda*X
*/
void ireig_residual(int n, rmulti **F0, rmulti **F1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **X0, rmulti **X1, rmulti *lambda0, rmulti *lambda1)
{
  irvec_mul_r(n,F0,F1,X0,X1,lambda0,lambda1);       // F=x*lambda
  irvec_neg(n,F0,F1,F0,F1);                         // F=-x*lambda
  irvec_add_lintr(n,n,F0,F1,A0,LDA0,A1,LDA1,X0,X1); // F=A*X-lambda*X
}

/** @} */

///////////////////////////////////////

/** @name クラフチック法 */
/** @{ */

#define NAME_KRAWCZYK "ireig_krawczyk"

/**
 @brief クラフチック法による固有値分解の誤差評価
 */
int ireig_krawczyk(int n, int k, rmulti **E_vec, int LDE, rmulti **E_val, rmulti **A, int LDA, rmulti **X, int LDX, rmulti **lambda, int debug)
{
  int p0,p1,prec,i,ret=0,r;
  rmulti **E=NULL;
  p0=rmat_get_prec_max(n,n,E_vec,LDE);
  p1=rvec_get_prec_max(n,E_val);
  prec=MAX2(p0,p1);
  E=rvec_allocate_prec(n+1,prec);
  for(i=0; i<k; i++){
    r=ireig_1pair_krawczyk(n,E,A,LDA,&COL(X,i,LDX),lambda[i],debug-1);
    rvec_copy(n,&COL(E_vec,i,LDE),E);
    rcopy(E_val[i],E[n]);
    if(r){ ret=1; }
    if(debug>=1){
      if(r){ print_red();   printf("[%s] %d/%d failed",NAME_KRAWCZYK,i,n);  print_reset(); printf("\n"); }
      else { print_green(); printf("[%s] %d/%d success",NAME_KRAWCZYK,i,n); print_reset(); printf("\n"); }
    }
  }
  if(debug>=1){ printf("[%s] result=%s\n",NAME_KRAWCZYK,(ret?"failed":"success"));  }
  E=rvec_free(n+1,E);
  return ret;
}


/**
 @brief クラフチック法による固有値対の誤差評価
 @param[in]  x      固有ベクトル
 @param[in]  lambda 固有値
 @param[in]  A      行列
 @param[in]  n      サイズ
 @param[out] e      誤差
 @return     偽値は成功，真値は失敗.
 */
int ireig_1pair_krawczyk(int n, rmulti **e, rmulti **A, int LDA, rmulti **x, rmulti *lambda, int debug)
{
  int prec,info,ret=1,m;
  double alpha=4.0; // e=alpha*|R*F|
  rmulti **J=NULL,**J0=NULL,**J1=NULL;
  rmulti **R=NULL;
  rmulti **F=NULL,**F0=NULL,**F1=NULL;
  rmulti **X0=NULL,**X1=NULL;
  rmulti **M0=NULL,**M1=NULL,**H0=NULL,**H1=NULL;
  // begin
  debug=0;
  m=n+1;
  prec=rvec_get_prec_max(m,e);
  RMA(J,m,m,prec); RMA(J0,m,m,prec); RMA(J1,m,m,prec); RMA(R,m,m,prec); RMA(M0,m,m,prec); RMA(M1,m,m,prec);
  RVA(F,m,prec); RVA(F0,m,prec); RVA(F1,m,prec); RVA(X0,m,prec); RVA(X1,m,prec); RVA(H0,m,prec); RVA(H1,m,prec);
  // J=[A -x; x' 0]
  rmat_set_zeros(m,m,J,m);              // J=zeors(m,m)
  rmat_copy(n,n,J,m,A,LDA);             // J(1:n,1:m)=A
  rvec_neg(n,&COL(J,n,m),x);            // J(:,m)=-x
  rmat_copy_t(n,1,&MAT(J,n,0,m),m,x,n); // J(n,:)=x'
  // [J0,J1]=[A -x; x' 0]
  irmat_set_zeros(m,m,J0,m,J1,m);                              // J=zeors(m,m)
  irmat_copy(n,n,J0,m,J1,m,A,LDA,A,LDA);                       // J(1:n,1:m)=A
  irvec_neg(n,&COL(J0,n,m),&COL(J1,n,m),x,x);                  // J(:,m)=-x
  irmat_copy_t(n,1,&MAT(J0,n,0,m),m,&MAT(J1,n,0,m),m,x,n,x,n); // J(n,:)=x'
  // R=inv(J)
  rmat_set_eye(m,m,R,m);
  rsolve(m,m,R,m,J,m,&info);
  if(info){ printf("Error info=%d, rsolve() in ireig_1pair_krawczyk()\n",info); exit(0); }
  // F=[A*x-lambda*x; (sum(x.^2)-1)/2]
  reig_residual(n,F,A,LDA,x,lambda); // F(1:n)=A*x-lambda*x
  rvec_sum_pow2(F[n],n,x);           // F(m)=(sum(x.^2)-1)/2
  rsub_d2(F[n],F[n],1);
  rmul_d(F[n],F[n],0.5);
  // [F0,F1]=[A*x-lambda*x; (sum(x.^2)-1)/2]
  ireig_residual(n,F0,F1,A,LDA,A,LDA,x,x,lambda,lambda); // F(1:n)=A*x-lambda*x
  irvec_sum_pow2(F0[n],F1[n],n,x,x);                     // F(m)=(sum(x.^2)-1)/2
  irsub_d2(F0[n],F1[n],F0[n],F1[n],1);
  irmul_d(F0[n],F1[n],F0[n],F1[n],0.5);
  // e=abs(R*F)
  rvec_lintr(m,m,e,R,m,F);
  rvec_abs(m,e,e);
  rvec_mul_d(m,e,e,alpha);
  // [X0,X1]=[-e,e]
  irvec_pm(m,X0,X1,e);
  // [M0,M1]=I-R*[J0,J1]
  rmat_set_eye(m,m,M0,m); rmat_set_eye(m,m,M1,m);
  irmat_sub_prod(m,m,m,M0,m,M1,m,R,m,R,m,J0,m,J1,m);
  // [H0,H1]=-R*[F0,F1]+[M0,M1]*[X0,X1]
  rvec_set_zeros(m,H0); rvec_set_zeros(m,H1);
  irvec_sub_lintr(m,m,H0,H1,R,m,R,m,F0,F1);
  irvec_add_lintr(m,m,H0,H1,M0,m,M1,m,X0,X1);
  // check X0<=H0 and H1<=X1
  if(!rvec_has_nan(m,X0) && !rvec_has_nan(m,X1) &&
     !rvec_has_nan(m,H0) && !rvec_has_nan(m,H1) &&
     rvec_le(m,X0,H0) && rvec_le(m,H1,X1)){ ret=0; }
  else                                    { ret=1; }
  // done
  RMF(J,m,m); RMF(J0,m,m); RMF(J1,m,m); RMF(R,m,m); RMF(M0,m,m); RMF(M1,m,m);
  RVF(F,m); RVF(F0,m); RVF(F1,m); RVF(X0,m); RVF(X1,m); RVF(H0,m); RVF(H1,m);
  return ret;
}

/** @} */

//EOF
