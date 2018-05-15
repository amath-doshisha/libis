#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"mt19937ar.h"
#include"is_macros.h"
#include"is_print.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_zvec.h"
#include"is_zmat.h"
#include"is_zsolve.h"
#include"is_zeig.h"
#include"is_zhshldr.h"
#include"is_zhpeig.h"

// enum { ZHPEIG_NONE=0, ZHPEIG_CONVERGENT, ZHPEIG_DIVERGENT, ZHPEIG_SINGULAR, ZHPEIG_NUM };
char *zhpeig_str[]={ "    ", "conv", "divr", "sing", NULL };

// JM=A-lambda*I-x*w'/C
void zhpeig_jacobi_mat(int n, dcomplex *JM, int LDJM, dcomplex *A, int LDA, dcomplex *X, dcomplex *W, dcomplex lambda, dcomplex C)
{
  dcomplex a;
  zmat_copy(n,n,JM,LDJM,A,LDA);            // JM=A
  zmat_diag_sub_scalar(n,JM,LDJM,lambda);  // JM=A-lambda*I
  Z_SET_INV(a,C); Z_NEG(a);                // a=-1/C
  zmat_rank1op(n,n,JM,LDJM,a,X,W);         // JM=A-lambda*I-(1/C)*x*w'
}

enum { ZHPEIG_STATUS_PRE=0, ZHPEIG_STATUS_CONV, ZHPEIG_STATUS_END };
char *zhpeig_status[]={" ", "o", "." };

// hyperplane constrained method for a pair
int zhpeig_1pair(int n, dcomplex *A, int LDA, dcomplex *z, dcomplex *x, dcomplex *lambda, double *E, int *Step, int debug)
{
#define PUT_1PAIR { printf("[%s] status=%s done=%s step=%03d lambda=%+.16e %+.16e  E=2^%+04d eta=2^%+04d=%.0e mu=%.0e\n",fname,zhpeig_status[status],zhpeig_str[done],step,Z_R(*lambda),Z_I(*lambda),(int)(log(fabs(*E))/log(2)),(int)(log(fabs(eta))/log(2)),eta,mu); }
  char fname[]="zhpeig_1pair";
  int status=ZHPEIG_STATUS_PRE,step=0,step_max=0,info=-1,done=0,LDJM;
  dcomplex *w=NULL,*F=NULL,*JM=NULL,C,dot;
  double eps,eps_half,eta,eta0,mu;
  // allocation
  w=zvec_allocate(n);
  F=zvec_allocate(n);
  LDJM=n; JM=zmat_allocate(LDJM,n);
  // constants
  eps=pow(2,-53);
  eps_half=pow(2,-26);
  eta=1.0/eps; eta0=0; mu=eta/eta0;
  step_max=ZHPEIG_STEP_MAX;  
  zvec_lintr_ct(n,n,w,A,LDA,z); // w=A'*z
  if(debug>1){
    printf("[%s] eps=%.1e=2^%d\n",fname,eps,(int)(log(eps)/log(2)));
    printf("[%s] eps_half=%.1e=2^%d\n",fname,eps_half,(int)(log(eps_half)/log(2)));
    printf("[%s] max of itr steps is %d.\n",fname,step_max);
  }
  // initial vector
  step=0;
  zvec_normalize_sgn_zvec(n,x,x);                       // x=x/sqrt(x'*x)
  C=zvec_dot(n,z,x);                               // C=z'*x
  dot=zvec_dot(n,w,x);                             // lambda=w'*x
  Z_SET_DIV((*lambda),dot,C);                      // lambda=(w'*x)/C
  zeig_residual(n,F,A,LDA,x,(*lambda));            // F=A*x-lambda*x;
  (*E)=dnorm_max_zvec(n,F);                         // E=max(abs(F));
  if((*E)==0) { done=ZHPEIG_CONVERGENT; } else { done=ZHPEIG_NONE; } // convergent or not?
  if(debug>0) { PUT_1PAIR; }
  // compute
  step=1;
  while(done==ZHPEIG_NONE){
    zhpeig_jacobi_mat(n,JM,LDJM,A,LDA,x,w,(*lambda),C); // JM=A-lambda*I-x*w'/C
    info=zsolve(n,1,F,n,JM,LDJM);                       // F=JM\F;
    if(info<0)      { ERROR_EXIT("Error in zhpeig_1pair(), zsolve, info=%d.\n",info); }
    else if(info>0) { done=ZHPEIG_SINGULAR; }
    else{
      eta=dnorm_max_zvec(n,F);                        // eta=norm_max(F)
      zvec_sub_zvec(n,x,x,F);                        // x=x-F
      zvec_normalize_sgn_zvec(n,x,x);                     // x=x/sqrt(x'*x)
      C=zvec_dot(n,z,x);                             // C=z'*x
      dot=zvec_dot(n,w,x);                           // lambda=w'*x
      Z_SET_DIV((*lambda),dot,C);                    // lambda=(w'*x)/C
      zeig_residual(n,F,A,LDA,x,(*lambda));          // F=A*x-lambda*x;
      (*E)=dnorm_max_zvec(n,F);                       // E=max(abs(F));
    }
    mu=2*eta/eta0; // mu=2*eta/eta0
    if(status==ZHPEIG_STATUS_PRE && eta<eps_half)        { status=ZHPEIG_STATUS_CONV; } // start to convergent
    else if(status==ZHPEIG_STATUS_CONV && mu>1.0)        { status=ZHPEIG_STATUS_END; }  // end to convergent
    if(status==ZHPEIG_STATUS_END)                        { done=ZHPEIG_CONVERGENT; }    // convergent or not?
    else if(status==ZHPEIG_STATUS_PRE && step>=step_max) { done=ZHPEIG_DIVERGENT; }     // divergent or not?
    if(debug>0) { PUT_1PAIR; }
    if(done==ZHPEIG_NONE) { step++; eta0=eta; } // next step
  }
  // free
  w=zvec_free(w);
  F=zvec_free(F);
  JM=zmat_free(JM);
  // done
  (*Step)=step;
  return done;
}

int zhpeig(int n, dcomplex *A, int LDA, dcomplex *X, int LDX, dcomplex *Lambda, int debug)
{
#define Qk COL(Q,k,LDQ)
#define Xk COL(X,k,LDX)
#define Hk COL(H,k,LDH)
  char msg,fname[]="zhpeig";
  int ret=0,conv,t,k,step,LDQ,LDH,fail,info=0,done,fail_max,itr=0;
  dcomplex *z=NULL,*Q=NULL,*r=NULL,*H=NULL;
  double E,*alpha=NULL;
  // allocation
  LDQ=n; Q=zmat_allocate(LDQ,n);
  LDH=n; H=zmat_allocate(LDH,n);
  alpha=dvec_allocate(n);
  z=zvec_allocate(n);
  r=zvec_allocate(n);
  // constants
  fail_max=n;
  if(debug>1){
    printf("[%s] matrix size is %d.\n",fname,n);
    printf("[%s] max of fails is %d.\n",fname,fail_max);
  }
  // init 
  zmat_set_zeros(n,n,X,LDX);
  zmat_set_zeros(n,n,H,LDH);
  zvec_set_zeros(n,Lambda);
  dvec_set_zeros(n,alpha); 
  zmat_set_eye(n,n,Q,LDQ);
  // initialize MT
  init_genrand(0);
  // loop
  t=0; done=0;
  for(k=0; !done && k<n; k++){
    zvec_copy_zvec(n,z,&Qk);                                                  // set normal vector z as Q(:,k).
    if(k==0) { zvec_set_rand(n,&Xk,2,-1); zvec_normalize_sgn_zvec(n,&Xk,&Xk); } // set initial vector X(:,k) as unit random vecot if k=0.
    else     { zvec_copy_zvec(n,&Xk,z); }                                     // set initial vector X(:,k) as z if k!=0.
    fail=0; conv=0; E=1;
    for(fail=0; !conv && fail<fail_max; fail++){
      info=zhpeig_1pair(n,A,LDA,z,&Xk,&Lambda[k],&E,&step,debug-1);
      t++;
      itr+=step;
      if(info==ZHPEIG_CONVERGENT) { conv=1; msg='o'; } // convergent
      else                        { conv=0; msg='x'; } // divergent
      if(debug>0){
	if(msg=='x' && fail==fail_max-1) print_light_red(); else if(msg=='x') print_red(); else print_green();
	printf("[%s] trial=%03d %03d/%03d [%c] ",fname,t,k+1,n,msg);
	if(msg=='x') printf("fail=%03d ",fail+1); else printf("         ");
	printf("step=%03d log2(E)=%+6.1f lambda=%+.16e %+.16e\n",step,log(fabs(E))/log(2),Z_R(Lambda[k]),Z_I(Lambda[k]));
	print_reset();
      }
      if(!conv) { zvec_set_rand(n,&Xk,2,-1); zvec_normalize_sgn_zvec(n,&Xk,&Xk); } // reset initial vector X(:,k) as unit random vecot
    }
    if(!conv){ done=1; }
    if(!done){
      zhouseholder(n,0,k,k,&Hk,&alpha[k],H,LDH,alpha,&Xk); // create R and create Householder vector and normalize constant
      zhouseholder_right(n,n,Q,LDH,k,&Hk,alpha[k]);        // Q=Q*H; H=I-alpha*h*h'
    }
  }
  if(done){
    k--;
    zmat_set_zeros(n,n-k,&Xk,LDX);
    zvec_set_zeros(n-k,&Lambda[k]);
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
  Q=zmat_free(Q);
  H=zmat_free(H);
  z=zvec_free(z);
  alpha=dvec_free(alpha);
  r=zvec_free(r);
  // done 
  return ret;
}

