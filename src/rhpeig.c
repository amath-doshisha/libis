#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_print.h"
#include"is_rmulti.h"
#include"is_ivec.h"
#include"is_rvec.h"
#include"is_rmat.h"
#include"is_rsolve.h"
#include"is_reig.h"
#include"is_rhshldr.h"
#include"is_rhpeig.h"
#include"mt19937ar.h"

// enum { RHPEIG_NONE=0, RHPEIG_CONVERGENT, RHPEIG_DIVERGENT, RHPEIG_SINGULAR, RHPEIG_NUM };
const char *rhpeig_str[]={ "    ", "conv", "divr", "sing", NULL };

/**
 @brief 固有値問題 F=A*x-lambda*x, lambda=w'*x/Cのヤコビ行列M=A-lambda*I-x*w'/Cの計算.
 @param[in]  n      サイズ.
 @param[in]  A      サイズ(n,n)の行列.
 @param[in]  LDA    行列Aの第1次元.
 @param[in]  x      サイズnのベクトル.
 @param[in]  w      サイズnのベクトル.
 @param[in]  lambda 固有値.
 @param[in]  C      定数.
 @param[in]  LDJM   行列JMの第1次元.
 @param[in]  JM     サイズ(n,n)の初期化済み行列.
 @param[out] JM     ヤコビ行列.
*/
void rhpeig_jacobi_mat(int n, rmulti **JM, int LDJM, rmulti **A, int LDA, rmulti **x, rmulti **w, rmulti *lambda, rmulti *C)
{
  rmulti *a=NULL;
  a=rallocate_prec(rmat_get_prec_max(n,n,JM,LDJM));
  rmat_diag_sub_r(n,n,JM,LDJM,A,LDA,lambda); // JM=A-lambda*I
  rinv(a,C); rneg(a,a);                      // a=-1/C
  rmat_rank1op(n,n,JM,LDJM,JM,LDJM,a,x,w);   // JM=A-lambda*I-(1/C)*x*w'
  a=rfree(a);
}

enum { RHPEIG_STATUS_PRE=0, RHPEIG_STATUS_CONV, RHPEIG_STATUS_END };
const char *rhpeig_status[]={" ", "o", "." };

// hyperplane constrained method for a pair
int rhpeig_1pair(int n, rmulti **x, rmulti *lambda, rmulti *E, int *Step, rmulti **A, int LDA, rmulti **z, int debug)
{
#define PUT_1PAIR { mpfr_printf("[%s] status=%s done=%s step=%03d lambda=%+.16Re  E=%.0Re eta=%.0Re mu=%.0Re\n",fname,rhpeig_status[status],rhpeig_str[done],step,lambda,E,eta,mu); }
  char fname[]="rhpeig_1pair";
  int p0,p1,prec,status=RHPEIG_STATUS_PRE,step=0,step_max=0,info=-1,done=0,LDJM;
  rmulti **w=NULL,**F=NULL,**JM=NULL,*C=NULL;
  rmulti *eps=NULL,*eps_half=NULL,*eta=NULL,*eta0=NULL,*mu=NULL;
  // precision
  p0=rvec_get_prec_max(n,x);
  p1=rget_prec(lambda);
  prec=MAX2(p0,p1);
  // allocation
  w=rvec_allocate_prec(n,prec);
  F=rvec_allocate_prec(n,prec);
  LDJM=n; JM=rmat_allocate_prec(LDJM,n,prec);
  C=rallocate_prec(prec);
  eta=rallocate_prec(prec);
  eta0=rallocate_prec(prec);
  mu=rallocate_prec(prec);
  eps=rmepsilon(prec);
  eps_half=rmepsilon(prec/2);
  // constants
  rset_nan(eta); rset_nan(eta0); rset_nan(mu);
  step_max=RHPEIG_STEP_MAX;  
  rvec_lintr_t(n,n,w,A,LDA,z); // w=A'*z
  if(debug>=1){
    mpfr_printf("[%s] prec=%d\n",fname,prec);
    mpfr_printf("[%s] eps=%.1Re=2^%d\n",fname,eps,(int)(log(fabs(rget_d(eps)))/log(2)));
    mpfr_printf("[%s] eps_half=%.1Re=2^%d\n",fname,eps_half,(int)(log(fabs(rget_d(eps_half)))/log(2)));
    mpfr_printf("[%s] step_max=%d\n",fname,step_max);
  }
  // initial vector
  step=0;
  rvec_normalize_sgn(n,x,x);          // x=x/sqrt(x'*x)
  rvec_sum_mul(C,n,z,x);              // C=z'*x
  rvec_sum_mul(lambda,n,w,x);         // lambda=w'*x
  rdiv(lambda,lambda,C);              // lambda=(w'*x)/C
  reig_residual(n,F,A,LDA,x,lambda);  // F=A*x-lambda*x
  rvec_max_abs(E,n,F);                // E=max(abs(F))
  if(ris_zero(E)) { done=RHPEIG_CONVERGENT; } else { done=RHPEIG_NONE; } // convergent or not?
  if(debug>0) { PUT_1PAIR; }
  // compute
  step=1;
  while(done==RHPEIG_NONE){
    rhpeig_jacobi_mat(n,JM,LDJM,A,LDA,x,w,lambda,C); // JM=A-lambda*I-x*w'/C
    rsolve_lu(n,1,F,n,JM,LDJM,&info);                // F=JM\F;
    if(info<0)      { ERROR_EXIT("Error in rhpeig_1pair(), rsolve, info=%d.\n",info); }
    else if(info>0) { done=RHPEIG_SINGULAR; }
    else{
      rvec_max_abs(eta,n,F);              // eta=max(abs(F))
      rvec_sub(n,x,x,F);                  // x=x-F
      rvec_normalize_sgn(n,x,x);          // x=x/sqrt(x'*x)
      rvec_sum_mul(C,n,z,x);              // C=z'*x
      rvec_sum_mul(lambda,n,w,x);         // lambda=w'*x
      rdiv(lambda,lambda,C);              // lambda=(w'*x)/C
      reig_residual(n,F,A,LDA,x,lambda);  // F=A*x-lambda*x
      rvec_max_abs(E,n,F);                // E=max(abs(F))
    }
    rdiv(mu,eta,eta0); rmul_d(mu,mu,2);   // mu=2*eta/eta0
    if     (status==RHPEIG_STATUS_PRE && rlt(eta,eps_half)){ status=RHPEIG_STATUS_CONV; } // start to convergent
    else if(status==RHPEIG_STATUS_CONV && rgt_d2(mu,1.0))  { status=RHPEIG_STATUS_END; }  // end to convergent
    else if(status==RHPEIG_STATUS_CONV && rget_exp(E)-rget_exp(lambda)<prec){ status=RHPEIG_STATUS_END; }  // end to convergent
    if     (status==RHPEIG_STATUS_END)                     { done=RHPEIG_CONVERGENT; }    // convergent or not?
    else if(status==RHPEIG_STATUS_PRE && step>=step_max)   { done=RHPEIG_DIVERGENT; }     // divergent or not?
    if(debug>0) { PUT_1PAIR; }
    if(done==RHPEIG_NONE) { step++; rcopy(eta0,eta); }  // next step
  }
  // free
  w=rvec_free(n,w);
  F=rvec_free(n,F);
  JM=rmat_free(LDJM,n,JM);
  C=rfree(C);
  eta=rfree(eta);
  eta0=rfree(eta0);
  mu=rfree(mu);
  eps=rfree(eps);
  eps_half=rfree(eps_half);
  // done
  (*Step)=step;
  return done;
}

/**
 @brief 超平面制約法による固有値分解.
 @param[in] debug   デバグ用.
 @param[in] A       (n,n)型行列.
 @param[in] LDA     Aの第1次元.
 @param[in] X       (n,n)型行列.
 @param[in] LDX     Xの第1次元.
 @param[out] X      計算された固有ベクトル.
 @param[out] Lambda 計算された固有値.
 @return            計算された固有対の個数.
 */
int rhpeig(int n, rmulti **X, int LDX, rmulti **Lambda, rmulti **A, int LDA, int debug)
{
#define Qk COL(Q,k,LDQ)
#define Xk COL(X,k,LDX)
#define Hk COL(H,k,LDH)
  char msg,fname[]="rhpeig";
  int p0=0,p1=0,prec=53,ret=0,conv,t,k,step,LDQ,LDH,fail,info=0,done,fail_max,itr=0;
  rmulti **z=NULL,**Q=NULL,**H=NULL;
  rmulti **alpha=NULL,*E=NULL;
  // precision
  p0=rvec_get_prec_max(n,Lambda);
  p1=rmat_get_prec_max(n,n,X,LDX);
  prec=MAX2(p0,p1);
  // allocation
  LDQ=n; Q=rmat_allocate_prec(LDQ,n,prec);
  LDH=n; H=rmat_allocate_prec(LDH,n,prec);
  alpha=rvec_allocate_prec(n,prec);
  z=rvec_allocate_prec(n,prec);
  E=rallocate_prec(prec);
  // constants
  fail_max=n;
  if(debug>1){
    printf("[%s] matrix size is %d.\n",fname,n);
    printf("[%s] max of fails is %d.\n",fname,fail_max);
  }
  // init 
  if(rmat_has_nan(n,n,X,LDX)){ rmat_set_nan(n,n,X,LDX); }
  rmat_set_zeros(n,n,H,LDH);
  rvec_set_zeros(n,Lambda);
  rvec_set_zeros(n,alpha);
  rmat_set_eye(n,n,Q,LDQ);
  // initialize MT
  init_genrand(0);
  // loop
  t=0; done=0;
  for(k=0; !done && k<n; k++){
    rvec_copy(n,z,&Qk);                                                    // set normal vector z as Q(:,k).
    if(rvec_has_nan(n,&Xk)){
      if(k==0){ rvec_set_rand(n,&Xk,2,-1); rvec_normalize_sgn(n,&Xk,&Xk); } // set initial vector X(:,k) as unit random vecot if k=0.
      else    { rvec_copy(n,&Xk,z); }                                       // set initial vector X(:,k) as z if k!=0.
    }
    fail=0; conv=0; rset_one(E);
    for(fail=0; !conv && fail<fail_max; fail++){
      info=rhpeig_1pair(n,&Xk,Lambda[k],E,&step,A,LDA,z,debug-1); // compute one of eigenpairs
      t++;
      itr+=step;
      if(info==RHPEIG_CONVERGENT) { conv=1; msg='o'; } // convergent
      else                        { conv=0; msg='x'; } // divergent
      if(debug>0){
	if(msg=='x' && fail==fail_max-1) print_light_red(); else if(msg=='x') print_red(); else print_green();
	mpfr_printf("[%s] trial=%03d %03d/%03d [%c] ",fname,t,k+1,n,msg);
	if(msg=='x') mpfr_printf("fail=%03d ",fail+1); else printf("         ");
	mpfr_printf("step=%03d log2(E)=%+6.1f lambda=%+.16Re\n",step,log(fabs(rget_d(E)))/log(2),Lambda[k]);
	print_reset();
      }
      if(!conv) { rvec_set_rand(n,&Xk,2,-1); rvec_normalize_sgn(n,&Xk,&Xk); } // reset initial vector X(:,k) as unit random vecot
    }
    if(!conv){ done=1; }
    if(!done){
      rhouseholder(n,0,k,k,&Hk,alpha[k],H,LDH,alpha,&Xk); // create R and create Householder vector and normalize constant
      rhouseholder_right(n,n,Q,LDH,Q,LDH,k,&Hk,alpha[k]); // Q=Q*H; H=I-alpha*h*h'
    }
  }
  if(done){
    k--;
    rmat_set_zeros(n,n-k,&Xk,LDX);
    rvec_set_zeros(n-k,&Lambda[k]);
  }
  ret=k;
  if(debug>1){
    if(done) printf("[%s] failed.\n",fname); else printf("[%s] succeeded.\n",fname);
    printf("[%s] computed vetros is %d.\n",fname,k);
    printf("[%s] total trials is %d.\n",fname,t);
    printf("[%s] sum of iterations is %d.\n",fname,itr);
    printf("[%s] average of iterations is %.2f.\n",fname,(double)itr/n);
  }
  // free
  Q=rmat_free(LDQ,n,Q);
  H=rmat_free(LDH,n,H);
  z=rvec_free(n,z);   
  alpha=rvec_free(n,alpha);  
  E=rfree(E);
  // done 
  return ret;
}


//EOF
