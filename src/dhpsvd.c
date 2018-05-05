#include"isys.h"

#define J00 (MAT(JM,0,0,LDJM))
#define J01 (MAT(JM,0,m,LDJM))
#define J10 (MAT(JM,m,0,LDJM))
#define J11 (MAT(JM,m,m,LDJM))
#define Hu  (H[0])
#define Hv  (H[m])
#define Uk  (COL(U,k,LDU))
#define Vk  (COL(V,k,LDV))


void dhpsvd_jacobi_mat(int m, int n, double *JM, int LDJM, const double *A, int LDA, const double *u, const double *v, const double *w, double sigma, double C)
{
  dmat_set_zeros(m+n,m+n,JM,LDJM);           // JM=zeros(m+n,m+n)
  dmat_diag_sub_scalar(m,&J00,LDJM,sigma);   // JM(1:m,1:m)=-sigma*Im;
  dmat_copy(m,n,&J01,LDJM,A,LDA);            // JM(1:m,(m+1):end)=A
  dmat_rank1op(m,n,&J01,LDJM,(-1.0/C),u,w);  // JM(1:m,(m+1):end)=A-(u*w')/C;
  dmat_copy_t(m,n,&J10,LDJM,A,LDA);          // JM((m+1):end,1:m)=A';
  dmat_diag_sub_scalar(n,&J11,LDJM,-sigma);  // JM((m+1):end,(m+1):end)=-sigma*In
  dmat_rank1op(n,n,&J11,LDJM,(-1.0/C),v,w);  // JM((m+1):end,(m+1):end)=-sigma_u*In-(v*wu')/Cu;
}

int dhpsvd_1pair(int m, int n, const double *A, int LDA, const double *z, double *u, double *v, double *Sigma, double *E, int *Step, int debug)
{
  int step,step_max,done,info=-1,ret=DHPSVD_NONE,LDJM;
  double *w=NULL,*H=NULL,*JM=NULL,h,e,C,eps_e,eps_r,sigma;
  // allocation
  w=dvec_allocate(n);
  H=dvec_allocate(m+n);
  LDJM=m+n; JM=dmat_allocate(LDJM,m+n);
  // constants
  step_max=DHPSVD_STEP_MAX;
  eps_e=DHPSVD_EPS_1ERROR;
  eps_r=DHPSVD_EPS_RESIDUAL;
  if(debug>0){
    printf("[dhpsvd_1pair] eps_e=%.1e log2(eps_e)=%.2f\n",eps_e,log(fabs(eps_e))/log(2));
    printf("[dhpsvd_1pair] eps_r=%.1e log2(eps_e)=%.2f\n",eps_r,log(fabs(eps_r))/log(2));
    printf("[dhpsvd_1pair] step_max=%d\n",step_max);
  }
  // constants
  dvec_lintr_t(m,n,w,A,LDA,z);                          // w=A'*z
  // compute
  done=0; h=1; e=1;
  for(step=0; !done && step<step_max; step++){
    dvec_normalize(m,u);                                // u=u/sqrt(u'*u)
    dvec_normalize(n,v);                                // v=v/sqrt(v'*v)
    C=dvec_dot(m,z,u);                                  // C=z'*u
    sigma=dvec_dot(n,w,v)/C;                            // sigma=(w'*v)/C
    dsvd_residual(m,n,H,A,LDA,u,v,sigma);               // H=[ A*v-sigma*u;  A'*u-sigma*v ]
    e=dvec_norm_max(m+n,H);                             // e=max(abs(H))
    dhpsvd_jacobi_mat(m,n,JM,LDJM,A,LDA,u,v,w,sigma,C); // Jacobi matrix
    info=dsolve(m+n,1,H,m+n,JM,LDJM);                   // H=JM\H;
    if(info<0)      { ERROR_EXIT("Error in dhpsvd_1pair(), dsolve, info=%d.\n",info); }
    else if(info>0) { done=2; ret=DHPSVD_SINGULAR; }
    else{
      h=dvec_norm_max(m+n,H);                           // h=norm_max(H)
      dvec_sub(m,u,&Hu);                                // u=u-H(1:m);    
      dvec_sub(n,v,&Hv);                                // v=v-H((m+1):end);
    }
    if(h<eps_e || e<eps_r) { done=1; ret=DHPSVD_CONVERGENT; }
    if(debug>0){ printf("[dhpsvd_1pair] done=%d ret=%d step=%02d log2(e)=%+6.2f log2(h)=%+6.2f sigma=%+.16e\n",done,ret,step,log(fabs(e))/log(2),log(fabs(h))/log(2),sigma); }
  }
  if(!done) { ret=DHPSVD_DIVERGENT; }
  dvec_normalize_sgn(m,u);                              // u=u/sqrt(u'*u)
  dvec_normalize_sgn(n,v);                              // v=v/sqrt(v'*v)
  dvec_lintr_t(m,n,w,A,LDA,u);                          // w=A'*u
  sigma=dvec_dot(n,w,v)/C;                              // sigma=(w'*v)/C
  if(sigma<0) { dvec_scale(m,u,-1); sigma=-sigma; }     // if sigma<0
  dsvd_residual(m,n,H,A,LDA,u,v,sigma);                 // H=[ A*v-sigma*u;  A'*u-sigma*v ]
  e=dvec_norm_max(m+n,H);                               // e=max(abs(H))
  if(isnan(u[0])) { ret=DHPSVD_NAN; e=1; }
  if(debug>0){ printf("[dhpsvd_1pair] done=%d ret=%d step=%02d log2(e)=%+6.2f                sigma=%+.16e\n",done,ret,step,log(fabs(e))/log(2),sigma); }
  // free
  w=dvec_free(w);
  H=dvec_free(H);
  JM=dmat_free(JM);
  // return values
  (*Sigma)=sigma;
  (*E)=e;
  (*Step)=step;
  // done
  return ret;
}

/*
int dhpsvd(int m, int n, const double *A, int LDA, double *U, int LDU, double *V, int LDV, double *Sigma, int debug)
{
  char msg;
  int t,k,itr,step,LDH,fail,flag,fail_max,info,ret=0;
  double *z=NULL,*H=NULL,*alpha=NULL,e,E,eps;
  // constants
  fail_max=n;
  eps=DHPSVD_EPS_REJECT;
  if(debug>0){
    printf("[dhpsvd] eps=%.1e log2(eps)=%.2f\n",eps,log(fabs(eps))/log(2));
  }
  // allocate
  u=dvec_allocate(m);
  v=dvec_allocate(n);
  z=dvec_allocate(m);
  LDH=m; H=dmat_allocate(m,m);
  alpha=dvec_allocate(m);
  // init
  dmat_eye(m,m,U,LDU);
  dmat_eye(n,n,V,LDV);
  dvec_zeros(m,Sigma);
  dvec_zeros(m,alpha);
  // initialize MT
  init_genrand(0);
  // loop
  t=0; done=0; itr=0;
  for(k=0; !done && k<m; k++){
    dvec_copy(m,z,&Uk); // z=U(:,k)
    if(k==0){
      dvec_set_rand(m,&Uk,2,-1);
      dvec_set_rand(n,&Vk,2,-1);
      dvec_normalize(m,Uk);
      dvec_normalize(n,Vk);
    }
    E=1; flag=1;
    for(fail=0; flag && fail<fail_max; fail++){
      info=dhpsvd_1pari(m,n,A,LDA,z,&Uk,&Vk,&Sigma[k],&E,&step,debug-1);
      t++;
      itr+=step;
      flag=0; msg='O';
      if(!flag && info==DHPSVD_CONVERGENT)                { flag=0; msg='O'; }
      else if(!flag && info!=DHPSVD_CONVERGENT && E<eps)  { flag=0; msg='*'; }
      else if(!flag && info!=DHPSVD_CONVERGENT && E>=eps) { flag=1; msg='X'; }
      if(debug>0){
	if(msg=='X') print_red(); 
	else if(msg=='*') print_yellow();
	else print_green();
	printf("[dhpeig] t=%03d fail=%03d [%c] %03d/%03d info=%d step=%02d log2(E)=%+6.2f lambda=%+.16e\n",t,fail,msg,k+1,n,info,step,log(fabs(E))/log(2),Lambda[k]);
	print_reset();
      }
      if(flag){
	dvec_rand(m,&Uk,2,-1);
	dvec_rand(n,&Vk,2,-1);
	dvec_normalize(m,&Uk);
	dvec_normalize(n,&Vk);
      }
    }
    if(flag){
      done=1;
    }
    if(!done){
      dhouseholder(0,k,Hu,LDHu,alpha_u,m,u,k,&COL(Hu,k,LDHu),&alpha_u[k]);  // u  -> Hu
      dhouseholder_right(m,m,k,&COL(Hu,k,LDHu),alpha_u[k],U,LDU);   // Hu -> U
    }

    // re-compute Sigma
    Sigma[k]=dsvd_rayleigh_quotient(type[NSVD_TYPE_UV],m,n,A,LDA,&COL(U,k,LDU),&COL(V,k,LDV));
    if(Sigma[k]<0){
      Sigma[k]=-Sigma[k];
      if(type[NSVD_TYPE_UV]=='U') dvec_scale(m,&COL(U,k,LDU),-1);
      else                        dvec_scale(n,&COL(V,k,LDV),-1);
    }

    // re-compute norm max of residual
    dsvd_residual_f(m,n,A,LDA,&COL(U,k,LDU),&COL(V,k,LDV),Sigma[k],&H[0]);
    dsvd_residual_g(m,n,A,LDA,&COL(U,k,LDU),&COL(V,k,LDV),Sigma[k],&H[m]);
    E=dvec_norm_max(m+n,H);

#ifdef SVD_SNSVD_MESSAGE
    // print
    printf("[snsvd-ht-%s] %03d/%03d step=%02d e=%.0e E=%.0e sigma=%.14e\n",type,k,m,step,e,E,Sigma[k]);
#endif

  }

  // compute v_k for k=m+1,...,n
  dsvd_compute_v_from_complement(m,n,V,LDV);

  // free
  FREE(z);
  FREE(H);
  FREE(alpha);
  // done
  return 0;
}
*/
