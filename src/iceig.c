#include"is_cmulti.h"
#include"is_cmat.h"
#include"is_cvec.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_icvec.h"
#include"is_icmat.h"
#include"is_iceig.h"
#include"is_csolve.h"
#include"is_ceig.h"
#include"is_print.h"

/**
 @file  iceig.c
 @brief 多倍長精度実数型cmultiの機械区間演算の固有値計算に関する関数の定義.
 @details スカラーに関しては@link icmulti.c@endlinkを参照のこと.
          ベクトルに関しては@link icvec.c@endlinkを参照のこと.
          行列に関しては@link icmat.c@endlinkを参照のこと.
 */


#define CMA(A,M,N,P)     { A=cmat_allocate_prec(M,N,P); }
#define CMF(A,M,N)       { A=cmat_free(M,N,A); }
#define CMR(M,N,A,LDA,B) { cmat_round((M),(N),(A),(LDA),(B)); }
#define CVA(X,N,P)       { X=cvec_allocate_prec(N,P); }
#define CVF(X,N)         { X=cvec_free(N,X); }
#define CVR(N,X,B)       { cvec_round((N),(X),(B)); }


///////////////////////////////////////

/** @name 固有値問題 */
/** @{ */

/**
 @brief 固有値問題の残差 F=A*X-lambda*X
*/
void iceig_residual(int n, cmulti **F0, cmulti **F1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **X0, cmulti **X1, cmulti *lambda0, cmulti *lambda1)
{
  icvec_mul_c(n,F0,F1,X0,X1,lambda0,lambda1);       // F=lambda*X
  icvec_neg(n,F0,F1,F0,F1);                         // F=-lambda*X
  icvec_add_lintr(n,n,F0,F1,A0,LDA0,A1,LDA1,X0,X1); // F=A*X-lambda*X
}

/** @} */

///////////////////////////////////////

/** @name クラフチック法 */
/** @{ */

#define NAME_KRAWCZYK "iceig_krawczyk"

/**
 @brief クラフチック法による固有値分解の誤差評価
 */
int iceig_krawczyk(int n, int k, cmulti **E_vec, int LDE, cmulti **E_val, cmulti **A, int LDA, cmulti **X, int LDX, cmulti **lambda, int debug)
{
  int p0,p1,prec,i,ret=0,r;
  cmulti **E=NULL;
  p0=cmat_get_prec_max(n,n,E_vec,LDE);
  p1=cvec_get_prec_max(n,E_val);
  prec=MAX2(p0,p1);
  E=cvec_allocate_prec(n+1,prec);
  for(i=0; i<k; i++){
    if(debug>=1){  }
    r=iceig_1pair_krawczyk(n,E,A,LDA,&COL(X,i,LDX),lambda[i],debug-1);
    cvec_copy(n,&COL(E_vec,i,LDE),E);
    ccopy(E_val[i],E[n]);
    if(r){ ret=1; }
    if(debug>=1){
      if(r){ print_red();   printf("[%s] %d/%d failed ",NAME_KRAWCZYK,i,n);  print_reset(); printf("\n"); }
      else { print_green(); printf("[%s] %d/%d success",NAME_KRAWCZYK,i,n); print_reset(); printf("\n"); }
    }
  }
  if(debug>=1){ printf("[%s] result=%s\n",NAME_KRAWCZYK,(ret?"failed":"success"));  }
  E=cvec_free(n+1,E);
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
int iceig_1pair_krawczyk(int n, cmulti **e, cmulti **A, int LDA, cmulti **x, cmulti *lambda, int debug)
{
  int prec,info,ret=1,m;
  double alpha=2; // e=alpha*|R*F|
  cmulti **J=NULL,**J0=NULL,**J1=NULL;
  cmulti **R=NULL;
  cmulti **F=NULL,**F0=NULL,**F1=NULL;
  cmulti **X0=NULL,**X1=NULL;
  cmulti **M0=NULL,**M1=NULL,**H0=NULL,**H1=NULL;
  // begin
  m=n+1;
  prec=cvec_get_prec_max(m,e);
  CMA(J,m,m,prec); CMA(J0,m,m,prec); CMA(J1,m,m,prec); CMA(R,m,m,prec); CMA(M0,m,m,prec); CMA(M1,m,m,prec);
  CVA(F,m,prec); CVA(F0,m,prec); CVA(F1,m,prec); CVA(X0,m,prec); CVA(X1,m,prec); CVA(H0,m,prec); CVA(H1,m,prec);
  debug=0;
  // compute
  // J=[A -x; x' 0]
  cmat_set_zeros(m,m,J,m);              // J=zeors(m,m)
  cmat_copy(n,n,J,m,A,LDA);             // J(1:n,1:m)=A
  cvec_neg(n,&COL(J,n,m),x);            // J(:,m)=-x
  cmat_copy_ct(n,1,&MAT(J,n,0,m),m,x,n); // J(n,:)=conj(x)'
  // [J0,J1]=[A -x; x' 0]
  icmat_set_zeros(m,m,J0,m,J1,m);                              // J=zeors(m,m)
  icmat_copy(n,n,J0,m,J1,m,A,LDA,A,LDA);                       // J(1:n,1:m)=A
  icvec_neg(n,&COL(J0,n,m),&COL(J1,n,m),x,x);                  // J(:,m)=-x
  icmat_copy_ct(n,1,&MAT(J0,n,0,m),m,&MAT(J1,n,0,m),m,x,n,x,n); // J(n,:)=conj(x)'
  // R=inv(J)
  cmat_set_eye(m,m,R,m); csolve(m,m,R,m,J,m,&info);
  // F=[A*x-lambda*x; (x'.*x-1)]
  ceig_residual(n,F,A,LDA,x,lambda); // F(1:n)=A*x-lambda*x
  cvec_sum_abs2(C_R(F[n]),n,x);      // F(m)=(sum(abs(x).^2)-1)/2
  rset_d(C_I(F[n]),0);
  csub_d(F[n],F[n],1);
  // [F0,F1]=[A*x-lambda*x; (sum(x.^2)-1)/2]
  iceig_residual(n,F0,F1,A,LDA,A,LDA,x,x,lambda,lambda); // F(1:n)=A*x-lambda*x
  icset_d(F0[n],F1[n],0,0);                              // F(m)=0
  icvec_sum_abs2(C_R(F0[n]),C_R(F1[n]),n,x,x);           // F(m).r=sum(x.^2)
  icsub_d(F0[n],F1[n],F0[n],F1[n],1,1);                  // F(m)=sum(x.^2)-1
  // e=abs(R*F)
  cvec_lintr(m,m,e,R,m,F);
  cvec_absc(m,e,e);
  cvec_mul_d(m,e,e,alpha);
  // [X0,X1]=[-e,e]
  icvec_pm(m,X0,X1,e);
  // [M0,M1]=I-R*[J0,J1]
  cmat_set_eye(m,m,M0,m); cmat_set_eye(m,m,M1,m);
  icmat_sub_prod(m,m,m,M0,m,M1,m,R,m,R,m,J0,m,J1,m);
  // H=-R*[F0,F1]+[M0,M1]*X
  cvec_set_zeros(m,H0); cvec_set_zeros(m,H1);
  icvec_sub_lintr(m,m,H0,H1,R,m,R,m,F0,F1);
  icvec_add_lintr(m,m,H0,H1,M0,m,M1,m,X0,X1);  
  // check X0<=H0 and H1<=X1
  if(!cvec_has_nan(m,X0) && !cvec_has_nan(m,X1) &&
     !cvec_has_nan(m,H0) && !cvec_has_nan(m,H1) &&
     cvec_lec(m,X0,H0) && cvec_lec(m,H1,X1)){ ret=0; }
  else                                      { ret=1; }
  // done
  CMF(J,m,m); CMF(J0,m,m); CMF(J1,m,m); CMF(R,m,m); CMF(M0,m,m); CMF(M1,m,m);
  CVF(F,m); CVF(F0,m); CVF(F1,m); CVF(X0,m); CVF(X1,m);CVF(H0,m); CVF(H1,m);
  return ret;
}

/** @} */

//EOF
