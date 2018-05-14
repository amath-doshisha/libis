#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"mt19937ar.h"
#include"is_macros.h"
#include"is_print.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_dmat.h"
#include"is_dsolve.h"
#include"is_deig.h"
#include"is_dhshldr.h"
#include"is_dhpeig.h"

// enum { DHPEIG_NONE=0, DHPEIG_CONVERGENT, DHPEIG_DIVERGENT, DHPEIG_SINGULAR, DHPEIG_NUM };
char *dhpeig_str[]={ "    ", "conv", "divr", "sing", NULL };

// JM=A-lambda*I-x*w'/C
void dhpeig_jacobi_mat(int n, double *JM, int LDJM, double *A, int LDA, double *X, double *W, double lambda, double C)
{
  dmat_copy(n,n,JM,LDJM,A,LDA);
  dmat_diag_sub_scalar(n,JM,LDJM,lambda);
  dmat_rank1op(n,n,JM,LDJM,(-1.0/C),X,W);
}

enum { DHPEIG_STATUS_PRE=0, DHPEIG_STATUS_CONV, DHPEIG_STATUS_END };
char *dhpeig_status[]={" ", "o", "." };

// hyperplane constrained method for a pair
int dhpeig_1pair(int n, double *A, int LDA, double *z, double *x, double *lambda, double *E, int *Step, int debug)
{
#define PUT_1PAIR { printf("[%s] status=%s done=%s step=%03d lambda=%+.16e  E=2^%+04d eta=2^%+04d=%.0e mu=%.0e\n",fname,dhpeig_status[status],dhpeig_str[done],step,(*lambda),(int)(log(fabs(*E))/log(2)),(int)(log(fabs(eta))/log(2)),eta,mu); }
  char fname[]="zhpeig_1pair";
  int status=DHPEIG_STATUS_PRE,step=0,step_max=0,info=-1,done=0,LDJM;
  double *w=NULL,*F=NULL,*JM=NULL,C;
  double eps,eps_half,eta=0,eta0=0,mu=0;
  // allocation
  w=dvec_allocate(n);
  F=dvec_allocate(n);
  LDJM=n; JM=dmat_allocate(LDJM,n);
  // constants
  eps=pow(2,-53);
  eps_half=pow(2,-26);
  step_max=DHPEIG_STEP_MAX;  
  dvec_lintr_t(n,n,w,A,LDA,z); // w=A'*z
  if(debug>1){
    printf("[%s] eps=%.1e=2^%d\n",fname,eps,(int)(log(eps)/log(2)));
    printf("[%s] eps_half=%.1e=2^%d\n",fname,eps_half,(int)(log(eps_half)/log(2)));
    printf("[%s] max of itr steps is %d.\n",fname,step_max);
  }
  // initial vector
  step=0;
  dvec_normalize_sgn(n,x,x);            // x=x/sqrt(x'*x)
  C=dvec_dot(n,z,x);                    // C=z'*x
  (*lambda)=dvec_dot(n,w,x)/C;          // lambda=(w'*x)/C
  deig_residual(n,F,A,LDA,x,(*lambda)); // F=A*x-lambda*x;
  (*E)=dvec_norm_max(n,F);              // E=max(abs(F));
  if((*E)==0) { done=DHPEIG_CONVERGENT; } else { done=DHPEIG_NONE; } // convergent or not?
  if(debug>0) { PUT_1PAIR; }
  // compute
  step=1;
  while(done==DHPEIG_NONE){
    dhpeig_jacobi_mat(n,JM,LDJM,A,LDA,x,w,(*lambda),C); // JM=A-lambda*I-x*w'/C
    info=dsolve(n,1,F,n,JM,LDJM);                       // F=JM\F;
    if(info<0)      { ERROR_EXIT("Error in dhpeig_1pair(), dsolve, info=%d.\n",info); }
    else if(info>0) { done=DHPEIG_SINGULAR; }
    else{
      eta=dvec_norm_max(n,F);               // h=norm_max(F);
      dvec_sub_dvec(n,x,x,F);               // x=x-F
      dvec_normalize_sgn(n,x,x);            // x=x/sqrt(x'*x)
      C=dvec_dot(n,z,x);                    // C=z'*x
      (*lambda)=dvec_dot(n,w,x)/C;          // lambda=(w'*x)/C
      deig_residual(n,F,A,LDA,x,(*lambda)); // F=A*x-lambda*x;
      (*E)=dvec_norm_max(n,F);              // E=max(abs(F));
    }
    mu=2*eta/eta0; // mu=2*eta/eta0
    if(status==DHPEIG_STATUS_PRE && eta<eps_half)        { status=DHPEIG_STATUS_CONV; } // start to convergent
    else if(status==DHPEIG_STATUS_CONV && mu>1.0)        { status=DHPEIG_STATUS_END; }  // end to convergent
    if(status==DHPEIG_STATUS_END)                        { done=DHPEIG_CONVERGENT; }    // convergent or not?
    else if(status==DHPEIG_STATUS_PRE && step>=step_max) { done=DHPEIG_DIVERGENT; }     // divergent or not?
    if(debug>0) { PUT_1PAIR; }
    if(done==DHPEIG_NONE) { step++; eta0=eta; } // next step
  }
  // free
  dvec_free(w);
  dvec_free(F);
  dmat_free(JM);
  // done
  (*Step)=step;
  return done;
}

int dhpeig(int n, double *A, int LDA, double *X, int LDX, double *Lambda, int debug)
{
#define Qk COL(Q,k,LDQ)
#define Xk COL(X,k,LDX)
#define Hk COL(H,k,LDH)
  char msg,fname[]="zhpeig";
  int ret=0,conv,t,k,step,LDQ,LDH,fail,info=0,done,fail_max,itr=0;
  double *z=NULL,*Q=NULL,*r=NULL,*H=NULL,E,*alpha=NULL;
  // allocation
  LDQ=n; Q=dmat_allocate(LDQ,n);
  LDH=n; H=dmat_allocate(LDH,n);
  alpha=dvec_allocate(n);
  z=dvec_allocate(n);
  r=dvec_allocate(n);
  // constants
  fail_max=n;
  if(debug>1){
    printf("[%s] matrix size is %d.\n",fname,n);
    printf("[%s] max of fails is %d.\n",fname,fail_max);
  }
  // init 
  dmat_set_zeros(n,n,X,LDX);
  dmat_set_zeros(n,n,H,LDH);
  dvec_set_zeros(n,Lambda);
  dvec_set_zeros(n,alpha);
  dmat_set_eye(n,n,Q,LDQ);
  // initialize MT
  init_genrand(0);
  // loop
  t=0; done=0;
  for(k=0; !done && k<n; k++){
    dvec_copy(n,z,&Qk);                                              // set normal vector z as Q(:,k).
    if(k==0) { dvec_set_rand(n,&Xk,2,-1); dvec_normalize(n,&Xk,&Xk); } // set initial vector X(:,k) as unit random vecot if k=0.
    else     { dvec_copy(n,&Xk,z); }                                 // set initial vector X(:,k) as z if k!=0.
    fail=0; conv=0; E=1;
    for(fail=0; !conv && fail<fail_max; fail++){
      info=dhpeig_1pair(n,A,LDA,z,&Xk,&Lambda[k],&E,&step,debug-1);
      t++;
      itr+=step;
      if(info==DHPEIG_CONVERGENT) { conv=1; msg='o'; } // convergent
      else                        { conv=0; msg='x'; } // divergent
      if(debug>0){
	if(msg=='x' && fail==fail_max-1) print_light_red(); else if(msg=='x') print_red(); else print_green();
	printf("[%s] trial=%03d %03d/%03d [%c] ",fname,t,k+1,n,msg);
	if(msg=='x') printf("fail=%03d ",fail+1); else printf("         ");
	printf("step=%03d log2(E)=%+6.1f lambda=%+.16e\n",step,log(fabs(E))/log(2),Lambda[k]);
	print_reset();
      }
      if(!conv) { dvec_set_rand(n,&Xk,2,-1); dvec_normalize_sgn(n,&Xk,&Xk); } // reset initial vector X(:,k) as unit random vecot
    }
    if(!conv){ done=1; }
    if(!done){
      dhouseholder(n,0,k,k,&Hk,&alpha[k],H,LDH,alpha,&Xk); // create R and create Householder vector and normalize constant
      dhouseholder_right(n,n,Q,LDH,k,&Hk,alpha[k]);        // Q=Q*H; H=I-alpha*h*h'
    }
  }
  if(done){
    k--;
    dmat_set_zeros(n,n-k,&Xk,LDX);
    dvec_set_zeros(n-k,&Lambda[k]);
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
  dmat_free(Q);
  dmat_free(H);
  dvec_free(z);
  dvec_free(alpha);
  dvec_free(r);
  // done
  return ret;
}

//EOF
