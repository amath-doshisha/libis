#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_print.h"
#include"is_rmulti.h"
#include"is_ivec.h"
#include"is_rvec.h"
#include"is_cvec.h"
#include"is_cmat.h"
#include"is_csolve.h"
#include"is_ceig.h"
#include"is_chshldr.h"
#include"is_chpeig.h"
#include"is_icmulti.h"
#include"is_iceig.h"
#include"mt19937ar.h"


/////////////////////////////////////////////////////////

#define RA(X,P) ((X)=rallocate_prec(prec))
#define RF(X)   ((X)=rfree(X))
#define CA(X,P) ((X)=callocate_prec(prec))
#define CF(X)   ((X)=cfree(X))
#define RVA(X,N,P){ X=rvec_allocate_prec(N,P); }
#define RVF(X,N)  { X=rvec_free(N,X); }
#define CVA(X,N,P){ X=cvec_allocate_prec(N,P); }
#define CVF(X,N)  { X=cvec_free(N,X); }

/////////////////////////////////////////////////////////

// enum { CHPEIG_NONE=0, CHPEIG_CONVERGENT, CHPEIG_DIVERGENT, CHPEIG_SINGULAR, CHPEIG_RECOMPUTE, CHPEIG_NUM };
const char *chpeig_str[]={ "    ", "conv", "divr", "sing", NULL };

/**
 @brief 固有値問題 F=A*x-lambda*x, lambda=w'*x/Cのヤコビ行列M=A-lambda*I-x*w'/Cの計算.
 @param[in]      n      サイズ.
 @param[in]      A      サイズ(n,n)の行列.
 @param[in]      LDA    行列Aの第1次元.
 @param[in]      x      サイズnのベクトル.
 @param[in]      w      サイズnのベクトル.
 @param[in]      lambda 固有値.
 @param[in]      C      定数.
 @param[in]      LDJM   行列JMの第1次元.
 @param[in,out]  JM     [in]サイズ(n,n)の初期化済み行列.[out]ヤコビ行列.
*/
void chpeig_jacobi_mat(int n, cmulti **JM, int LDJM, cmulti **A, int LDA, cmulti **x, cmulti **w, cmulti *lambda, cmulti *C)
{
  cmulti *a=NULL;
  a=callocate_prec(cmat_get_prec_max(n,n,JM,LDJM));
  cmat_diag_sub_c(n,n,JM,LDJM,A,LDA,lambda); // JM=A-lambda*I
  cinv(a,C); cneg(a,a);                      // a=-1/C
  cmat_rank1op(n,n,JM,LDJM,JM,LDJM,a,x,w);   // JM=A-lambda*I-(1/C)*x*w'
  a=cfree(a);
}

enum { CHPEIG_STATUS_PRE=0, CHPEIG_STATUS_CONV, CHPEIG_STATUS_END };
const char *chpeig_status[]={" ", "o", "." };

// hyperplane constrained method for a pair
int chpeig_1pair(int n, cmulti **x, cmulti *lambda, rmulti *E, int *Step, cmulti **A, int LDA, cmulti **z, int debug)
{
#define PUT_1PAIR { mpfr_printf("[%s] status=%s done=%s step=%03d lambda=%+.16Re %+.16Re  E=%.0Re eta=%.0Re mu=%.0Re\n",fname,chpeig_status[status],chpeig_str[done],step,C_R(lambda),C_I(lambda),E,eta,mu); }
  char fname[]="chpeig_1pair";
  int p0,p1,prec,status=CHPEIG_STATUS_PRE,step=0,step_max=0,info=-1,done=0,LDJM;
  cmulti **w=NULL,**F=NULL,**JM=NULL,*C=NULL;
  rmulti *eps=NULL,*eps_half=NULL,*eta=NULL,*eta0=NULL,*mu=NULL;
  // precision
  p0=cvec_get_prec_max(n,x);
  p1=cget_prec(lambda);
  prec=MAX2(p0,p1);
  // allocation
  LDJM=n; JM=cmat_allocate_prec(LDJM,n,prec);
  CVA(w,n,prec); CVA(F,n,prec); CA(C,prec); RA(eta,prec); RA(eta0,prec); RA(mu,prec);
  eps=rmepsilon(prec); eps_half=rmepsilon(prec*0.4);
  // constants
  rset_nan(eta); rset_nan(eta0); rset_nan(mu);
  step_max=CHPEIG_STEP_MAX;
  cvec_lintr_ct(n,n,w,A,LDA,z); // w=A'*z
  if(debug>=1){
    mpfr_printf("[%s] prec=%d\n",fname,prec);
    mpfr_printf("[%s] eps=%.1Re=2^%d\n",fname,eps,(int)(log(fabs(rget_d(eps)))/log(2)));
    mpfr_printf("[%s] eps_half=%.1Re=2^%d\n",fname,eps_half,(int)(log(fabs(rget_d(eps_half)))/log(2)));
    mpfr_printf("[%s] step_max=%d\n",fname,step_max);
  }
  // initial vector
  step=0;
  cvec_normalize_sgn(n,x,x);           // x=x/sqrt(x'*x)
  cvec_sum_dot(C,n,z,x);               // C=z'*x
  cvec_sum_dot(lambda,n,w,x);          // lambda=w'*x
  cdiv(lambda,lambda,C);               // lambda=(w'*x)/C
  ceig_residual(n,F,A,LDA,x,lambda);   // F=A*x-lambda*x
  cvec_max_abs(E,n,F);                 // E=max(abs(F))
  if(ris_zero(E)){ done=CHPEIG_CONVERGENT; } // convergent?
  else { done=CHPEIG_NONE; }                 // not convergent?
  if(debug>=1) { PUT_1PAIR; }
  // compute
  step=1;
  while(done==CHPEIG_NONE){
    chpeig_jacobi_mat(n,JM,LDJM,A,LDA,x,w,lambda,C); // JM=A-lambda*I-x*w'/C
    csolve_lu(n,1,F,n,JM,LDJM,&info);                // F=JM\F;
    if(info<0)      { ERROR_EXIT("Error in chpeig_1pair(), csolve, info=%d.\n",info); }
    else if(info>0) { done=CHPEIG_SINGULAR; }
    else{
      cvec_max_abs(eta,n,F);             // eta=max(abs(F))
      cvec_sub_cvec(n,x,x,F);                 // x=x-F
      cvec_normalize_sgn(n,x,x);         // x=x/sqrt(x'*x)
      cvec_sum_dot(C,n,z,x);             // C=z'*x
      cvec_sum_dot(lambda,n,w,x);        // lambda=w'*x
      cdiv(lambda,lambda,C);             // lambda=(w'*x)/C
      ceig_residual(n,F,A,LDA,x,lambda); // F=A*x-lambda*x
      cvec_max_abs(E,n,F);               // E=max(abs(F))
    }
    rdiv(mu,eta,eta0); rmul_d(mu,mu,2);  // mu=2*eta/eta0
    if     (status==CHPEIG_STATUS_PRE && rlt(eta,eps_half)){ status=CHPEIG_STATUS_CONV; } // start to convergent
    else if(status==CHPEIG_STATUS_CONV && rgt_d2(mu,1.0))  { status=CHPEIG_STATUS_END; }  // end to convergent
    else if(status==CHPEIG_STATUS_CONV && rget_exp(E)-cget_exp(lambda)<prec){ status=CHPEIG_STATUS_END; }  // end to convergent
    if     (status==CHPEIG_STATUS_END){ done=CHPEIG_CONVERGENT; }    // convergent or not?
    else if(status==CHPEIG_STATUS_PRE && step>=step_max)   { done=CHPEIG_DIVERGENT; }     // divergent or not?
    if(debug>0) { PUT_1PAIR; }
    if(done==CHPEIG_NONE) { step++; rcopy(eta0,eta); } // next step
  }
  // free
  w=cvec_free(n,w);  
  F=cvec_free(n,F); 
  JM=cmat_free(LDJM,n,JM);
  CF(C); RF(eta); RF(eta0); RF(mu); RF(eps); RF(eps_half);
  // done
  (*Step)=step;
  return done;
}

/**
 @brief 超平面制約法による固有値分解.
 @param[in]     debug   デバグ用.
 @param[in]     A       (n,n)型行列.
 @param[in]     LDA     Aの第1次元.
 @param[in,out] X       [in](n,n)型行列.[out]計算された固有ベクトル.
 @param[in]     LDX     Xの第1次元.
 @param[out]    Lambda 計算された固有値.
 @return               計算された固有対の個数.
 */
int chpeig(int n, cmulti **X, int LDX, cmulti **Lambda, cmulti **A, int LDA, int debug)
{
#define Qk COL(Q,k,LDQ)
#define Xk COL(X,k,LDX)
#define Hk COL(H,k,LDH)
  char msg,fname[]="chpeig";
  int p0,p1,prec,ret=0,conv,t,k,step,LDQ,LDH,fail,info=0,done,fail_max,itr=0;
  cmulti **z=NULL,**Q=NULL,**H=NULL;
  rmulti **alpha=NULL,*E=NULL;
  // precision
  p0=cvec_get_prec_max(n,Lambda);
  p1=cmat_get_prec_max(n,n,X,LDX);
  prec=MAX2(p0,p1);
  // allocation
  LDQ=n; Q=cmat_allocate_prec(LDQ,n,prec);
  LDH=n; H=cmat_allocate_prec(LDH,n,prec);
  alpha=rvec_allocate_prec(n,prec);
  z=cvec_allocate_prec(n,prec);
  E=rallocate_prec(prec);
  // constants
  fail_max=n;
  if(debug>1){
    printf("[%s] matrix size is %d.\n",fname,n);
    printf("[%s] max of fails is %d.\n",fname,fail_max);
  }
  // init 
  if(cmat_has_nan(n,n,X,LDX)){ cmat_set_nan(n,n,X,LDX); }
  cmat_set_zeros(n,n,H,LDH);
  cvec_set_zeros(n,Lambda);
  rvec_set_zeros(n,alpha);
  cmat_set_eye(n,n,Q,LDQ);
  // initialize MT
  init_genrand(0);
  // loop
  t=0; done=0;
  for(k=0; !done && k<n; k++){
    cvec_copy(n,z,&Qk);                                                      // set normal vector z as Q(:,k).
    if(cvec_has_nan(n,&Xk)){
      if(k==0) { cvec_set_rand(n,&Xk,2,-1); cvec_normalize_sgn(n,&Xk,&Xk); } // set initial vector X(:,k) as unit random vecot if k=0.
      else     { cvec_copy(n,&Xk,z); }                                       // set initial vector X(:,k) as z if k!=0.
    }
    fail=0; conv=0; rset_d(E,1);
    for(fail=0; !conv && fail<fail_max; fail++){
      info=chpeig_1pair(n,&Xk,Lambda[k],E,&step,A,LDA,z,debug-1); // compute one of eigenpairs
      t++;
      itr+=step;
      if(info==CHPEIG_CONVERGENT) { conv=1; msg='o'; } // convergent
      else                        { conv=0; msg='x'; } // divergent
      if(debug>0){
	if(msg=='x' && fail==fail_max-1) print_light_red(); else if(msg=='x') print_red(); else print_green();
	mpfr_printf("[%s] trial=%03d %03d/%03d [%c] ",fname,t,k+1,n,msg);
	if(msg=='x') mpfr_printf("fail=%03d ",fail+1); else printf("         ");
	mpfr_printf("step=%03d log2(E)=%+6.1f lambda=%+.16Re %+.16Re\n",step,log(fabs(rget_d(E)))/log(2),C_R(Lambda[k]),C_I(Lambda[k]));
	print_reset();
      }
      if(!conv) { cvec_set_rand(n,&Xk,2,-1); cvec_normalize_sgn(n,&Xk,&Xk); } // reset initial vector X(:,k) as unit random vecot
    }
    if(!conv){ done=1; }
    if(!done){
      chouseholder(n,0,k,k,&Hk,alpha[k],H,LDH,alpha,&Xk); // create R and create Householder vector and normalize constant
      chouseholder_right(n,n,Q,LDH,Q,LDH,k,&Hk,alpha[k]); // Q=Q*H; H=I-alpha*h*h'
    }
  }
  if(done){
    k--;
    cmat_set_zeros(n,n-k,&Xk,LDX);
    cvec_set_zeros(n-k,&Lambda[k]);
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
  Q=cmat_free(LDQ,n,Q);
  H=cmat_free(LDH,n,H);
  z=cvec_free(n,z); 
  alpha=rvec_free(n,alpha);
  E=rfree(E);
  // done 
  return ret;
}

int get_next_prec(int prec)
{
  if     (prec<=24){ return 53; }
  else if(prec<=53){ return 128; }
  else             { return prec*2; }
}

/**
 @brief 超平面制約法による固有値分解.
 @param[in]     debug   デバグ用.
 @param[in]     A       (n,n)型行列.
 @param[in]     LDA     Aの第1次元.
 @param[in,out] X       [in](n,n)型行列.[out]計算された固有ベクトル.
 @param[in]     LDX     Xの第1次元.
 @param[out]    Lambda  計算された固有値.
 @return                計算された固有対の個数.
 */
int chpeig_verify(int n, cmulti **X, int LDX, cmulti **Lambda, cmulti **XE, int LDXE, cmulti **LE, cmulti **A, int LDA, int prec_verify, int *prec, int *kprec, int debug)
{
#define Qk COL(Q,k,LDQ)
#define Xk COL(X,k,LDX)
#define Hk COL(H,k,LDH)
  char msg,fname[]="chpeig-verify";
  int p0,p1,ret=0,conv,t,k,step,LDQ,LDH,fail,info=0,done,fail_max,itr=0,kfalse=0;
  cmulti **z=NULL,**Q=NULL,**H=NULL,**E_K=NULL;
  rmulti **alpha=NULL,*E=NULL,*XEmax=NULL,*LEmax=NULL,*LEmaxRel=NULL,*L=NULL;
  // precision
  p0=cvec_get_prec_max(n,Lambda);
  p1=cmat_get_prec_max(n,n,X,LDX);
  (*prec)=MAX2(p0,p1);
  (*kprec)=get_next_prec(*prec);
  // allocation
  LDQ=n; Q=cmat_allocate_prec(LDQ,n,(*prec));
  LDH=n; H=cmat_allocate_prec(LDH,n,(*prec));
  alpha=rvec_allocate_prec(n,(*prec));
  z=cvec_allocate_prec(n,(*prec));
  E=rallocate_prec(*prec);
  L=rallocate_prec((*prec));
  E_K=cvec_allocate_prec(n+1,(*kprec));
  XEmax=rallocate_prec((*kprec));
  LEmax=rallocate_prec((*kprec));
  LEmaxRel=rallocate_prec((*kprec));
  // constants
  fail_max=n;
  if(debug>1){
    printf("[%s] matrix size is %d.\n",fname,n);
    printf("[%s] max of fails is %d.\n",fname,fail_max);
  }
  // init 
  if(cmat_has_nan(n,n,X,LDX)){ cmat_set_nan(n,n,X,LDX); }
  cmat_set_zeros(n,n,H,LDH);
  cvec_set_zeros(n,Lambda);
  rvec_set_zeros(n,alpha);
  cmat_set_eye(n,n,Q,LDQ);
  // initialize MT
  init_genrand(0);
  // loop
  t=0; done=0;
  for(k=0; !done && k<n; k++){
    cvec_copy(n,z,&Qk);                                                      // set normal vector z as Q(:,k).
    if(cvec_has_nan(n,&Xk)){
      if(k==0) { cvec_set_rand(n,&Xk,2,-1); cvec_normalize_sgn(n,&Xk,&Xk); } // set initial vector X(:,k) as unit random vecot if k=0.
      else     { cvec_copy(n,&Xk,z); }                                       // set initial vector X(:,k) as z if k!=0.
    }
    fail=0; conv=CHPEIG_DIVERGENT; rset_d(E,1);
    for(fail=0; conv!=CHPEIG_CONVERGENT && fail<fail_max; fail++){
      info=chpeig_1pair(n,&Xk,Lambda[k],E,&step,A,LDA,z,debug-1); // compute one of eigenpairs
      t++;
      itr+=step;
      if(info==CHPEIG_CONVERGENT){
	// convergent
	conv=CHPEIG_CONVERGENT;
	msg='o';
      }else{
	// divergent
	conv=CHPEIG_DIVERGENT;
	msg='x';
      }
      if(debug>0){
	if(msg=='x' && fail==fail_max-1) print_light_red(); else if(msg=='x') print_red(); else print_green();
	mpfr_printf("[%s] trial=%03d %03d/%03d [%c] prec=%4d ",fname,t,k+1,n,msg,(*prec));
	if(msg=='x') mpfr_printf("fail=%03d ",fail+1); else printf("         ");
	mpfr_printf("step=%03d F=%8.1Re L=%+18.10Re %+18.10Re",step,E,C_R(Lambda[k]),C_I(Lambda[k]));
	print_reset(); if(msg=='x'){ printf("\n"); }
      }
      if(conv==CHPEIG_DIVERGENT){
	// gain prec
	(*prec)=get_next_prec(*prec);
	(*kprec)=get_next_prec(*prec);
	cmat_round(n,LDX,X,LDX,(*prec));
	cmat_round(n,LDXE,XE,LDXE,(*kprec));
	cvec_round(n,Lambda,(*prec));
	cvec_round(n,LE,(*kprec));
	cmat_round(LDQ,n,Q,LDQ,(*prec));
	cmat_round(LDH,n,H,LDH,(*prec));
	rvec_round(n,alpha,(*prec));
	cvec_round(n,z,(*prec));
	rround(E,*prec);
	rround(L,(*prec));
	cvec_round(n+1,E_K,(*kprec));
	rround(XEmax,(*kprec));
	rround(LEmax,(*kprec)); 
	rround(LEmaxRel,(*kprec)); 
	// reset initial vector X(:,k) as unit random vecot
	cvec_set_rand(n,&Xk,2,-1); cvec_normalize_sgn(n,&Xk,&Xk); 
      }
    }
    if(conv==CHPEIG_DIVERGENT){
      // terminate
      done=1;
    }else if(conv==CHPEIG_CONVERGENT){
      // convergent
      kfalse=iceig_1pair_krawczyk(n,E_K,A,LDA,&Xk,Lambda[k],debug-1);
      cvec_clone(n,&COL(XE,k,LDXE),E_K);
      cclone(LE[k],E_K[n]);
      cvec_max_absc(XEmax,n,&COL(XE,k,LDXE)); // absolute error
      cmax_absc(LEmax,LE[k]); // absolute error
      cmax_absc(L,Lambda[k]); rdiv_2exp(LEmaxRel,LEmax,rget_exp(L)); // relative error
      if(!kfalse && rget_exp(XEmax)<-prec_verify && rget_exp(LEmaxRel)<-prec_verify){
	conv=CHPEIG_CONVERGENT; msg='o';
      }else if(!kfalse && rget_exp(XEmax)<-prec_verify && rget_exp(LEmax)<-prec_verify && icin_pm(Lambda[k],Lambda[k],LE[k])){
	conv=CHPEIG_CONVERGENT; msg='0';
      }else{
	conv=CHPEIG_RECOMPUTE; msg='x';
      }
      if(conv==CHPEIG_RECOMPUTE){
	(*prec)=get_next_prec(*prec);
	(*kprec)=get_next_prec(*prec);
	cmat_round(n,LDX,X,LDX,(*prec));
	cmat_round(n,LDXE,XE,LDXE,(*kprec));
	cvec_round(n,Lambda,(*prec));
	cvec_round(n,LE,(*kprec));
	cmat_round(LDQ,n,Q,LDQ,(*prec));
	cmat_round(LDH,n,H,LDH,(*prec));
	rvec_round(n,alpha,(*prec));
	cvec_round(n,z,(*prec));
	rround(E,*prec);
	rround(L,(*prec));
	cvec_round(n+1,E_K,(*kprec));
	rround(XEmax,(*kprec));
	rround(LEmax,(*kprec)); 
	rround(LEmaxRel,(*kprec)); 
      }
      if(debug>=1){
	if(msg=='x'){ print_red();   mpfr_printf(" failed  LE=%8.1Re LEr=%8.1Re XE=%8.1Re [%c]",LEmax,LEmaxRel,XEmax,msg); print_reset(); printf("\n"); }
	else        { print_green(); mpfr_printf(" success LE=%8.1Re LEr=%8.1Re XE=%8.1Re [%c]",LEmax,LEmaxRel,XEmax,msg); print_reset(); printf("\n"); }
      }
    }else{ ERROR_AT; exit(0); }
    // continue if succeeded and not done
    if(conv==CHPEIG_CONVERGENT){
      chouseholder(n,0,k,k,&Hk,alpha[k],H,LDH,alpha,&Xk); // create R and create Householder vector and normalize constant
      chouseholder_right(n,n,Q,LDH,Q,LDH,k,&Hk,alpha[k]); // Q=Q*H; H=I-alpha*h*h'
    }else if(conv==CHPEIG_RECOMPUTE){
      // restart
      k=-1;
      cmat_set_eye(n,n,Q,LDQ); // Q=eye(n)
    }else{ ERROR_AT; exit(0); }
  }
  // failed?
  if(done){
    k--;
    cmat_set_zeros(n,n-k,&Xk,LDX);
    cvec_set_zeros(n-k,&Lambda[k]);
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
  Q=cmat_free(LDQ,n,Q);
  H=cmat_free(LDH,n,H);
  E_K=cvec_free(n+1,E_K);
  z=cvec_free(n,z); 
  alpha=rvec_free(n,alpha);
  E=rfree(E);
  L=rfree(L);
  XEmax=rfree(XEmax);
  LEmax=rfree(LEmax);
  LEmaxRel=rfree(LEmaxRel);
  // done 
  return ret;
}

//EOF
