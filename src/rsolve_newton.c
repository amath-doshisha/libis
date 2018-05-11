#include<stdio.h>
#include<stdlib.h>
#include<isys.h>

#define RMA(A,M,N,P) { A=rmat_allocate_prec(M,N,P); }
#define RMF(A,M,N)   { A=rmat_free(M,N,A); }
#define RVA(X,N,P)   { X=rvec_allocate_prec(N,P); }
#define RVF(X,N)     { X=rvec_free(N,X); }
#define RA(X,P)      { X=rallocate_prec(P); }
#define RF(X)        { X=rfree(X); }

// Input: x,dF
// Output: e
// Return: 0 if success, 1 if fail
int rsolve_krawczyk(int m, rmulti **e, rmulti **x, func_t *fF, int debug)
{
  int info,ret=1,prec;
  func_t *fJ=NULL;
  rmulti **L=NULL,**L0=NULL,**L1=NULL;
  rmulti **R=NULL;
  rmulti **F=NULL,**F0=NULL,**F1=NULL;
  rmulti **X0=NULL,**X1=NULL,**T0=NULL,**T1=NULL;
  rmulti **M0=NULL,**M1=NULL,**r0=NULL,**r1=NULL,**H0=NULL,**H1=NULL;
  // begin
  prec=rvec_get_prec_max(m,e);
  fJ=func_grad(func_retain(fF),func_var1_list(m));
  if(debug>1){ printf("fJ="); func_print(fJ); printf("\n"); }
  RMA(L,m,m,prec); RMA(L0,m,m,prec); RMA(L1,m,m,prec); RMA(R,m,m,prec); RMA(M0,m,m,prec); RMA(M1,m,m,prec);
  RVA(F,m,prec); RVA(F0,m,prec); RVA(F1,m,prec); RVA(X0,m,prec); RVA(X1,m,prec);
  RVA(T0,m,prec); RVA(T1,m,prec); RVA(r0,m,prec); RVA(r1,m,prec); RVA(H0,m,prec); RVA(H1,m,prec);
  // compute
  // L=F'(x)
  rmat_func_list2(m,m,L,m,fJ,m,x);
  // R=inv(L)
  rmat_set_eye(m,m,R,m); rsolve(m,m,R,m,L,m,&info);
  // F=F(x)
  rvec_func_list(m,F,fF,m,x); if(debug>0){ rvec_print(m,F,"F=",'e',1); }
  // [F0,F1]=F(x)
  irvec_func_list(m,F0,F1,fF,m,x,x); if(debug>0){ irvec_print(m,F0,F1,"[F]=",'e',1); }
  // e=abs(R*F)
  rvec_lintr(m,m,e,R,m,F); rvec_abs(m,e,e); rvec_mul_dscalar(m,e,e,2);
  // [X0,X1]=[-e,e]
  irvec_pm(m,X0,X1,e); if(debug>0){ irvec_print(m,X0,X1,"[X]=",'e',1); }
  // [T0,T1]=[x-e,x+e]
  irvec_add_rvec(m,T0,T1,X0,X1,x,x); if(debug>1){ irvec_print(m,T0,T1,"[T]=",'f',20); }
  // [L0,L1]=F'([T0,T1])
  irmat_func_list2(m,m,L0,L1,m,fJ,m,T0,T1);
  // [M0,M1]=I-R*[L0,L1]
  rmat_set_eye(m,m,M0,m); rmat_set_eye(m,m,M1,m);
  irmat_sub_prod(m,m,m,M0,m,M1,m,R,m,R,m,L0,m,L1,m);
  if(debug>1){ irmat_print(m,m,M0,m,M1,m,"[M]=",'e',1); }
  // H=-R*[F0,F1]+[M0,M1]*X
  rvec_set_zeros(m,H0); rvec_set_zeros(m,H1);
  irvec_sub_lintr(m,m,H0,H1,R,m,R,m,F0,F1);
  irvec_add_lintr(m,m,H0,H1,M0,m,M1,m,X0,X1);  
  if(debug>0){ irvec_print(m,H0,H1,"[H]=",'e',1); }
  // check X0<=H0 and H1<=X1
  if(rvec_le(m,X0,H0) && rvec_le(m,H1,X1)){ ret=0; }
  else                                    { ret=1; }
  // done
  RMF(L,m,m); RMF(L0,m,m); RMF(L1,m,m); RMF(R,m,m); RMF(M0,m,m); RMF(M1,m,m);
  RVF(F,m); RVF(F0,m); RVF(F1,m); RVF(X0,m); RVF(X1,m);
  RVF(T0,m); RVF(T1,m); RVF(r0,m); RVF(r1,m); RVF(H0,m); RVF(H1,m);
  fJ=func_del(fJ);
  return ret;
}

void rsolve_func_residual(int m, rmulti *norm, rmulti **F, rmulti **x, func_t *fF)
{
  int prec,p0,p1;
  rmulti **F0=NULL,**F1=NULL,*norm0=NULL,*norm1=NULL;
  p0=rget_prec(norm);
  p1=rvec_get_prec_max(m,F);
  prec=MAX2(p0,p1);
  RA(norm0,prec); RA(norm1,prec); RVA(F0,m,prec); RVA(F1,m,prec);
  // F=F(x)
  rvec_func_list(m,F,fF,m,x);
  // [F0,F1]=F(x)
  irvec_func_list(m,F0,F1,fF,m,x,x);
  // |F|
  rvec_max_abs(norm,m,F);
  rvec_max_abs(norm0,m,F0);
  rvec_max_abs(norm1,m,F1);
  // check
  if(rgt(norm1,norm)){ rswap(norm1,norm); rvec_swap(m,F,F1); }
  if(rgt(norm0,norm)){ rswap(norm0,norm); rvec_swap(m,F,F0); }
  // done
  RF(norm0); RF(norm1); RVF(F0,m); RVF(F1,m);
}

// Input: m,x,fF,step_max
// Output: x
int rsolve_newton(int m, rmulti **x, func_t *fF, int step_max, int debug)
{
  char *name="rsolve_newton",status=' ';
  int info,n,ret,prec;
  func_t *fJ=NULL;
  rmulti **F=NULL,**J=NULL,*eta=NULL,*eta_pre=NULL,*mu=NULL,*b=NULL,*a=NULL,**e=NULL,*emax=NULL;
  // init
  prec=rvec_get_prec_max(m,x);
  RVA(F,m,prec); RMA(J,m,m,prec); RA(eta,prec); RA(eta_pre,prec); RA(mu,prec); RA(b,prec); RA(a,prec);
  // begin
  if(step_max<0){ step_max=100; }
  fJ=func_grad(func_retain(fF),func_var1_list(m));
  if(debug>0){
    mpfr_printf("[%s] step_max=%d\n",name,step_max);
    printf("[%s] F(x)=",name); func_print(fF); printf("\n");
    printf("[%s] F'(x)=",name); func_print(fJ); printf("\n");
  }
  // compute
  if(rvec_has_nan(m,x)){ rvec_set_rand(m,x,2,-1); }
  n=0; rset_d(eta_pre,0);
  if(debug>0){ mpfr_printf("[%s] begin x[%d]\n",name,n,n); }
  do{
    // a=|x|
    rvec_max_abs(a,m,x);
    // F=F(x)
    rvec_func_list(m,F,fF,m,x);
    // b=|F|
    rvec_max_abs(b,m,F);
    // exactly F=0?
    if(req_d(b,0)){ rsolve_func_residual(m,b,F,x,fF); }
    // if F=0, then...
    if(req_d(b,0)){ status='-'; }
    else{
      // J=J(x)
      rmat_func_list2(m,m,J,m,fJ,m,x);
      // F=J\F
      rsolve(m,1,F,m,J,m,&info);
      if(info){
	// set status
	status='x';
	rset_nan(eta); rset_nan(mu);
      }else{
	// eta=|J\F|
	rvec_max_abs(eta,m,F);
	// mu=eta/eta_pre
	rdiv_r(mu,eta,eta_pre);
	// set status
	if     (n>=step_max)                   { status='/'; }
	else if(status==' ' && rle_d2(mu,1e-2)){ status='o'; }
	else if(status=='o' && rge_d2(mu,1))   { status='='; }
      }
    }
    if(debug>0){
      RVA(e,m,prec); RA(emax,prec);
      info=rsolve_krawczyk(m,e,x,fF,0);
      rvec_max_abs(emax,m,e);
      if(info){ mpfr_printf("[%s] n=%d [%c]         |x%d|=%.0Re |F(x%d)|=%.0Re |x%d-x%d|=%.0Re mu=%.0Re\n",name,n,status,n,a,n,b,n+1,n,eta,mu); }
      else    { mpfr_printf("[%s] n=%d [%c] e=%.0Re |x%d|=%.0Re |F(x%d)|=%.0Re |x%d-x%d|=%.0Re mu=%.0Re\n",name,n,status,emax,n,a,n,b,n+1,n,eta,mu); }
      RVF(e,m); RF(emax);
    }
    // update x
    if(status==' ' || status=='o'){
      // x=x-J\F
      rvec_sub_rvec(m,x,x,F);
      // next
      n++;
      rswap(eta,eta_pre);
    }
  }while(status==' ' || status=='o');
  // value
  if     (status=='/'){ ret=RSOLVE_NEWTON_SUCCESS; }
  else if(status=='-'){ ret=RSOLVE_NEWTON_SUCCESS; }
  else if(status=='='){ ret=RSOLVE_NEWTON_SUCCESS; }
  else if(status=='o'){ ret=RSOLVE_NEWTON_FAILED_MORESTEP; }
  else if(status==' '){ ret=RSOLVE_NEWTON_FAILED_STEPMAX; }
  else                { ret=RSOLVE_NEWTON_FAILED; }
  if(debug>0){ mpfr_printf("[%s] return x[%d] ret=%d\n",name,n,ret); }
  // done
  RVF(F,m); RMF(J,m,m); RF(eta); RF(eta_pre); RF(mu); RF(b); RF(a);
  fJ=func_del(fJ);
  return ret;
}

//EOF
