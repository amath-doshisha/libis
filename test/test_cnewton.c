#include<stdio.h>
#include<stdlib.h>
#include<mpfr.h>
#include<isys.h>

#define RMA(A,M,N) { A=rmat_allocate(M,N); }
#define RMF(A,M,N) { A=rmat_free(M,N,A); }
#define RVA(X,N)   { X=rvec_allocate(N); }
#define RVF(X,N)   { X=rvec_free(N,X); }
#define RA(X)      { X=rallocate(); }
#define RF(X)      { X=rmfree(X); }
#define CMA(A,M,N) { A=cmat_allocate(M,N); }
#define CMF(A,M,N) { A=cmat_free(M,N,A); }
#define CVA(X,N)   { X=cvec_allocate(N); }
#define CVF(X,N)   { X=cvec_free(N,X); }
#define CA(X)      { X=callocate(); }
#define CF(X)      { X=cfree(X); }

void usage()
{
  printf("Usage:\n");
  exit(0);
}

void eterm_show(int m, cmulti **x, func_t *fF, int bmax, int kappa)
{
  int kmax,k,db=64;
  int *b=NULL,tau;
  func_t *fJ=NULL;
  cmulti **y=NULL;
  rmulti **e=NULL,**eps=NULL,**em=NULL,**r=NULL,*rmax=NULL,**p=NULL;
  
  kmax=(bmax-64)/db; if(kmax<0){ kmax=2; }
  tau=cvec_get_exp_max(m,x);
  printf("tau=%d\n",tau);

  fJ=func_grad(func_retain(fF),func_var1_list(m));
  CVA(y,m); RVA(e,kmax); RVA(eps,kmax); RVA(em,kmax); RVA(r,kmax); RVA(p,kmax); RA(rmax);

  b=ivec_allocate(kmax);

  b[0]=53;
  for(k=1; k<kmax; k++){ b[k]=(k-1)*db+64; }

  for(k=0; k<kmax; k++){
    printf("b=%4d ",b[k]);
    set_default_prec(b[k]);
    rset_d(eps[k],1); rmul_2exp(eps[k],eps[k],-b[k]+tau);
    mpfr_printf("2^(-b+tau)=%8.1Re ",eps[k]);
    cvec_round_cvec(m,y,b[k]);
    csolve_newton_map(m,y,x,fF,fJ);
    csolve_newton_e_norm(e[k],m,y,x,fF,fJ,bmax*4);
    if(gt_rr(eps[k],e[k])){ printf("> "); }else{ printf("< "); }
    mpfr_printf("e=%8.1Re ",e[k]);
    rset_d(em[k],1); rmul_2exp(em[k],em[k],-b[k]+tau+kappa);
    if(gt_rr(e[k],em[k])){ printf("> "); }else{ printf("< "); }
    mpfr_printf("em=%8.1Re ",em[k]);

    rlog2_r(p[k],e[k]); rsub_dr(p[k],tau,p[k]); rsub_dr(p[k],b[k],p[k]);
    mpfr_printf("b-tau+log2(e)=%+6.2Rf ",p[k]);


    rlog2_r(r[k],e[k]); rsub_dr(r[k],tau,r[k]); rdiv_rd(r[k],r[k],b[k]); rsub_dr(r[k],1,r[k]);
    mpfr_printf("κ=1-(tau-log2(e))/b=%+.10Rf ",r[k]);
    printf("\n");
  }
  rmax_rvec(rmax,kmax,p);
  mpfr_printf("cancel bits=%+.2Rf\n",rmax);
  rmax_rvec(rmax,kmax,r);
  mpfr_printf("κ_max=%+.10Rf\n",rmax);
  

  // done
  fJ=func_del(fJ);
  b=ivec_free(b);
  CVF(y,m); RVF(e,kmax); RVF(eps,kmax); RVF(em,kmax); RVF(r,kmax); RVF(p,kmax); RF(rmax);
}

int main(int argc, char *argv[])
{
  char *eq[]={"{x-sin(3/2*(x+y)), y-cos(3/2*(x+y))}",
	      "{x^3+y^2+1, x^2+y^2+2}",
	      "{2*x^2+y, x+2*y+5}",
	      "{3*x^2-1,log(x)+y+3}",          //-eq 0 -eps 1e-100 --> ## 2 step 戻る ## 実数解 
	      "{x^2+y^2+x^2-1,x^2+z^2-y,x-z}",       //eq 5  ## 2step 戻る ##
	      "{-x+x^2+2*y^2+2*z^2+2*w^2, -y+2*x*y+2*y*z+2*z*w, -z+y^2+2*x*z+2*y*w, -1+x+2*y+2*z+2*w}", //eq 6 ## 2step 戻る ##
	      "{x-10000-y^2,x*y^2-100}",             // eq 8 ## 2step 戻る ##
	      "{0.3*x^2-1}",                         // eq 9 ## 2step 戻る ##
	      "{(3*x-y^2)^3-1000,log(x)^(-1)+2*y-30}",                //eq 12 ### 3step 戻る ### 
	      "{cos(0.5*x)+cos(y)-1, sin(0.5*x)-sin(y)-1}",           //eq31 // x=pi, y=0  ## 2回 2step 戻る ##	      
	      "{-x+x^2+2*y^2+2*z^2+2*w^2, -y+2*x*y+2*y*z+2*z*w, -z+y^2+2*x*z+2*y*w, -1+x+2*y+2*z+2*w}", //eq 36 4変数 x=1, y=z=w=0  ## 2step 戻る ##
	      "{x^2-y^2-1,x^4-y^4-1999}",             // 桁落ち kappa=0.06,  解 x=sqrt(1000),y=sqrt(999)
	      "{x^2+y^2-1, x^2+2*x*y+y^2-2}",         // 線形
	      "{x+y,x*y}",                            // 線形，零
	      NULL 
  };
  int debug=0,i,m,step_max0=-1,step_max=-1,prec0=53,prec=53,no=0,solve_true=0,info,seed=0,eterm=0,e_prec=2048,e_seed=-1,kappa=26,l=4;
  double mu=8;
  func_t *fF=NULL;
  cmulti **x0=NULL,**x=NULL,**x_true=NULL,**e=NULL;
  rmulti *eps=NULL,*eps_true=NULL;

  // init func_t
  func_eval(func_script("begin(x,y,z,w,v,u)"));

  // default prec
  set_default_prec(prec0);

  // allocate
  RA(eps); RA(eps_true);

  // default parameters
  rset_d(eps,1e-200);
  rset_s(eps_true,"1e-2000");

  // get options
  i=1;
  while(i<argc){
    if(STR_EQ(argv[i],"--help"))                 { usage(); }
    else if(STR_EQ(argv[i],"-v"))                { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))               { debug=2; }
    else if(STR_EQ(argv[i],"-vvv"))              { debug=3; }
    else if(STR_EQ(argv[i],"-true"))             { solve_true=1; }
    else if(STR_EQ(argv[i],"-eterm"))            { eterm=1; }
    else if(i+1<argc && STR_EQ(argv[i],"-mu"))   { mu=atof(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-l"))    { l=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-kappa")){ kappa=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-eq"))   { no=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-n"))    { step_max=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-prec0")){ prec0=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-prec")) { prec=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-e-prec")){ e_prec=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-e-seed")){ e_seed=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-eps"))  { rset_s(eps,argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-seed")) { seed=atoi(argv[++i]); }
    else                                         { usage(); }
    i++;
  }

  //
  printf("step_max=%d\n",step_max);
  printf("mu=%g\n",mu);
  printf("l=%d\n",l);
  printf("kappa=%d\n",kappa);
  fF=func_script(eq[no]);  printf("fF="); func_print(fF); printf("\n");
  m=func_asize(fF);        printf("m=%d\n",m);

  // allocate vectors and matrices
  set_default_prec(prec0);  printf("prec0=%d\n",prec0);
  CVA(x0,m); CVA(x,m); CVA(x_true,m); CVA(e,m);

  // set initial vector
  init_genrand(seed);
  cvec_set_rand(m,x0,2,-1);     cvec_print(m,x0,"x0=",'f',6);

  // solve for true solution
  printf("#-------------------\n");
  cvec_clone_cvec(m,x_true,x0);
  csolve_newton(m,x_true,fF,step_max0,debug);
    //prec=csolve_newton_adjust(m,x_true,fF,NULL,eps_true,step_max0,mu,l,kappa,debug);
  cvec_print(m,x_true,"x*=",'e',20);
  cvec_clone_cvec(m,x,x_true);

  // solve
  if(solve_true){
    printf("#-------------------\n");
    cvec_clone_cvec(m,x,x0);
    prec=csolve_newton_adjust(m,x,fF,x_true,eps,step_max,mu,l,kappa,debug);
    cvec_print(m,x,"x=",'e',20);
  }
 
  // error
  printf("#-------------------\n");
  printf("prec=%d\n",prec);
  set_default_prec(prec);
  cvec_round_cvec(m,e,prec);
  cvec_set_nan(m,e);
  info=csolve_krawczyk(m,e,x,fF,debug-2);
  if(info){ print_red(); printf("failed.\n"); } else{ print_green(); printf("succeeded.\n"); } print_reset();
  cvec_print(m,e,"e=",'e',1);

  if(eterm){
    printf("#-------------------\n");
    if(e_seed>=0){ init_genrand(e_seed); cvec_set_rand(m,x,2,-1); }
    eterm_show(m,x,fF,e_prec,kappa);
  }

  // done
  printf("#-------------------\n");
  fF=func_del(fF);
  func_clear();
  RF(eps); CVF(x,m); CVF(e,m);
  printf("func_new_del_check_sum=%d\n",func_new_del_check_sum());  
  return 0;
}

//EOF
