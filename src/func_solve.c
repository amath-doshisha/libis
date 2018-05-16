#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_cvec.h"
#include"is_rsolve.h"
#include"is_csolve.h"

#define FR(F) func_retain(F)

//////////////////////////////////////////////////////////////////////////

func_t *func_mat_solve(func_t *A, func_t *b)
{
  func_t *x=NULL,*C=NULL;
  int m,n,k,info;
  x=func_evalf(func_clone(FR(b)));
  C=func_evalf(func_clone(FR(A)));
  m=func_rows(C);
  n=func_cols(C);
  k=func_size(x);
  if(func_is(C,"rmat") && func_is(x,"rvec") && m==n && n==k){
    rsolve(n,1,func_rvec_p(x),func_rvec_size(x),func_rmat_p(C),func_rmat_ld(C),&info);
    if(info){ x=func_del(x); }
  }else if(func_is(C,"cmat") && func_is(x,"cvec") && m==n && n==k){
    csolve(n,1,func_cvec_p(x),func_cvec_size(x),func_cmat_p(C),func_cmat_ld(C),&info);
    if(info){ x=func_del(x); }
  }else if(func_is(C,"cmat") && func_is(x,"rvec")){
    x=func_rvec_get_cvec(x);
    csolve(n,1,func_cvec_p(x),func_cvec_size(x),func_cmat_p(C),func_cmat_ld(C),&info);
    if(info){ x=func_del(x); }
  }else if(func_is(C,"rmat") && func_is(x,"cvec")){
    C=func_rmat_get_cmat(C);
    csolve(n,1,func_cvec_p(x),func_cvec_size(x),func_cmat_p(C),func_cmat_ld(C),&info);
    if(info){ x=func_del(x); }
  }else{ FUNC_ERROR_ARG2("func_mat_solve",A,b); }
  C=func_del(C);
  A=func_del(A);
  b=func_del(b);
  return x;
}


//////////////////////////////////////////////////////////////////////////

func_t *func_poly_solve(func_t *f)
{
  func_t *x=NULL;
  if(func_is_poly(f) && func_poly_n_var(f)==1){
    if     (func_degree(f)==1){ x=func_poly_solve_var1(FR(f)); }
    else if(func_degree(f)==2){ x=func_poly_solve_var2(FR(f)); }
    else if(func_degree(f)>=3){ x=func_poly_solve_varn(FR(f)); }
  }
  f=func_del(f);
  return x;
}

//////////////////////////////////////////////////////////////////////////

func_t *func_poly_solve_var1(func_t *f)                   //ax+b=0 の解 
{
  func_t *x=NULL,*a=NULL,*b=NULL;
  if(func_is_mono_not_num(f)){ x=func_zero(); }       //b=0 のとき
  else if(func_is_poly(f) && func_poly_n_var(f)==1 && func_degree(f)==1){
    a=func_poly_clone_coeff_ntarm(f,1);
    b=func_number_pull(FR(f));
    x=func_poly_calculate_soltion(FR(a),FR(b));
  }
  a=func_del(a);
  b=func_del(b);
  f=func_del(f);
  return x;
}

//////////////////////////////////////////////////////////////////////////

func_t *func_poly_solve_var2(func_t *f)
{
  func_t *x=NULL;
  if(func_is_poly(f) && func_poly_n_var(f)==1 && func_degree(f)==2){
    if(func_poly_get_mono_ntarm(f,1)==NULL){ x=func_poly_solve_var2_have_ac(FR(f)); }
    else if(func_poly_get_mono_ntarm(f,0)==NULL){ x=func_poly_solve_var2_have_ab(FR(f)); }
    else{ x=func_poly_solve_var2_have_abc(FR(f)); }
    func_a_rm_null(x);
  }
  f=func_del(f);
  return x;
}

func_t *func_poly_solve_var2_have_ac(func_t *f)          //ax^2+c=0 の解 (c=0を含む) x=+-sqrt(-c/a)
{
  func_t *x=NULL,*a=NULL,*c=NULL;
  if(func_is_poly(f) && func_poly_n_var(f)==1 && func_degree(f)==2 && func_poly_get_mono_ntarm(f,1)==NULL){
    x=func_list(2);
    a=func_poly_clone_coeff_ntarm(f,2);
    c=func_number_pull(FR(f));
    if(func_is_zero(c)){ x->a[0]=func_zero(); }
    else{
      x->a[0]=func_poly_calculate_soltion(FR(a),FR(c));
      x->a[0]=func_sqrt(x->a[0]);
      x->a[1]=func_mul(func_bigint_int(-1,1),FR(x->a[0]));
    }
  }
  f=func_del(f);
  a=func_del(a);
  c=func_del(c);
  return x;
}

func_t *func_poly_solve_var2_have_ab(func_t *f)        //ax^2+bx=0 の解 x=0,-b/a
{
  func_t *x=NULL,*a=NULL,*b=NULL;
  if(func_is_poly(f) && func_poly_n_var(f)==1 && func_degree(f)==2 && func_poly_get_mono_ntarm(f,0)==NULL){
    a=func_poly_clone_coeff_ntarm(f,2);
    b=func_poly_clone_coeff_ntarm(f,1);
    x=func_list(2);
    x->a[0]=func_zero();
    x->a[1]=func_poly_calculate_soltion(FR(a),FR(b));
  }
  f=func_del(f);
  a=func_del(a);
  b=func_del(b);
  return x;
}

func_t *func_poly_solve_var2_have_abc(func_t *f)                //ax^2+bx+c=0 の解 x={-b+-sqrt(b^2-4ac)}/2a
{
  func_t *x=NULL,*a=NULL,*b=NULL,*c=NULL;
  if(func_is_poly(f) && func_poly_n_var(f)==1 && func_degree(f)==2){
    a=func_poly_clone_coeff_ntarm(f,2);
    b=func_poly_clone_coeff_ntarm(f,1);
    c=func_number_pull(FR(f));
    x=func_poly_quadratic_formula(FR(a), FR(b), FR(c));
  }
  a=func_del(a);
  b=func_del(b);
  c=func_del(c);
  f=func_del(f);
  return x;
}

//////////////////////////////////////////////////////////////////////////

func_t *func_poly_calculate_soltion(func_t *a,func_t *b)                              //x=-b/aの計算
{
  func_t *x=NULL;
  x=func_div(FR(b),FR(a));
  x=func_expand(func_mul(x,func_bigint_int(-1,1)));
  x=func_evalf(x);
  a=func_del(a);
  b=func_del(b);
  return x;
}

func_t *func_poly_quadratic_formula(func_t *a,func_t *b,func_t *c)                      //解の公式
{
  func_t *x=NULL,*x1=NULL,*x2=NULL;
  x=func_mul(func_bigint_int(4,1),func_mul(FR(a),FR(c)));              //4ac 
  x=func_sub(func_pow_n(FR(b),2),x);                          //b^2-4ac
  x=func_sqrt(x);                                                               //sqrt(b^2-4ac)
  x=func_evalf(x);
  x1=func_sub(FR(x),FR(b));
  x1=func_div(x1,func_mul(func_bigint_int(2,1),FR(a)));
  if(!func_is_zero(x)){                                                    //not重解ならば
    x2=func_mul(func_bigint_int(-1,1),func_add(FR(b),FR(x)));
    x2=func_div(x2,func_mul(func_bigint_int(2,1),FR(a)));
  }
  x=func_del(x);
  x=func_list(2);
  x->a[0]=FR(x1);
  x->a[1]=FR(x2);
  a=func_del(a);
  b=func_del(b);
  c=func_del(c);
  x1=func_del(x1);
  x2=func_del(x2);
  return x;
}

//////////////////////////////////////////////////////////////////////////

func_t *func_poly_solve_varn(func_t *f)
{
  int n=0,i,j,k=0;
  func_t *x=NULL,*g=NULL;
  cmulti **c1=NULL,**c2=NULL,**c3=NULL;
  if(func_is_mono_not_num(f)){ x=func_zero(); }       //x^n=0 のとき
  else{
    for(i=0; i<func_asize(f); i++){
      if(func_is_number(f->a[i])){
	g=FR(f);
	n=func_poly_degree_max(g);
	x=func_list(n);
    	break;	
      }else if(i==func_asize(f)-1 && !func_is_number(f->a[i])){
	g=func_poly_div_q(FR(f),FR(func_mono_get_var(f->a[i])));
	n=func_poly_degree_max(g);
	x=func_list(n+1);
	x->a[n]=func_zero();
      }
    }
    c1=cvec_allocate(n);  //係数
    c2=cvec_allocate(n);  //初期値
    c3=cvec_allocate(n);  //出力
    func_ccopy_coeff(n,c1,g);
    func_init_val(n,c2,c1);
    i=0;
    while(i!=1){
      k=k+1;
      i=func_weierstrass(n,c3,c2,c1);
      for(j=0; j<n; j++){
	cset_c(c2[j],c3[j]);
      }
      if(k>20){ printf("hanpukukaisu > 20\n"); break;}
    }
    printf("hanpuku kaisu [%d]\n",k);        //反復回数表示
    for(j=0; j<n; j++){
      x->a[j]=func_complex();
      cset_c(x->a[j]->p.cm,c3[j]);
    }
  }
  f=func_del(f);
  g=func_del(g);
  c1=cvec_free(n,c1);
  c2=cvec_free(n,c2);
  c3=cvec_free(n,c3);
  return x;
}

//////////////////////////////////////////////////////////////////////////

func_t *func_poly_list_solve(func_t *f)                          //連立方程式を解く
{
  int i=0,j,k=0,n;
  func_t *a=NULL,*b=NULL,*g=NULL,*h=NULL,*x=NULL;
  if(func_is_poly_list(f)){
    h=func_clone(FR(f));
    for(i=0; i<func_asize(h); i++){
      if(func_poly_n_var(h->a[i])==1){
	a=func_list(0);
	b=func_list(3);
	k=func_poly_var1n(h->a[i]);	
	x=func_poly_solve(h->a[i]); h->a[i]=NULL; 
	func_a_rm_null(h);
	if(func_is_number(x)){
	  b->a[0]=func_var1(k,1);
	  b->a[1]=FR(x);
	  g=func_maps(FR(h),k,FR(x));
	  for(n=0; n<func_asize(g); n++){
	    if(!func_is_zero(g->a[n]) && func_is_number(g->a[n])){
	      print_red(); printf("ERROR exist num\n"); func_print(g); printf("\n"); print_reset(); 
	      exit(0); 
	    }
	  }
	  b->a[2]=func_poly_list_solve(FR(g));
	  func_a_append(a,FR(b));
	}
	else{
	  for(j=0; j<func_asize(x); j++){
	    b->a[0]=func_var1(k,1);
	    b->a[1]=FR(x->a[j]);
	    g=func_maps(FR(h),k,FR(x->a[j]));
	    for(n=0; n<func_asize(g); n++){
	      if(!func_is_zero(g->a[n]) && func_is_number(g->a[n])){
		print_red(); printf("ERROR exist num\n"); func_print(g); printf("\n"); print_reset(); 
		exit(0);
	      }
	    }

	    b->a[2]=func_poly_list_solve(FR(g));
	    func_a_append(a,FR(b));
	    g=func_del(g);g=NULL;
	    b=func_del(b);b=func_list(3);
	  }
	}
	break;
      }
    }
  }
  f=func_del(f);
  b=func_del(b);
  g=func_del(g);
  h=func_del(h);
  x=func_del(x);
  return a;
}

func_t *func_poly_list_solve_and_arrange(func_t *f)
{
  int i,size1,size2,I[100];
  func_t *g=NULL;
  if(func_is_poly_list(f)){
    g=func_poly_list_solve(FR(f));
    if(g!=NULL){
      g=func_list_sol_convert_tree(g);
      size1=func_asize(g);
      size2=func_asize(g->a[0]->a[0]);
      for(i=0; i<size1; i++){
	func_args_sort_index(g->a[i]->a[0],I);
	func_args_reverse(g->a[i]->a[0]);
	ivec_reverse(size2,I);
	func_args_arrange(g->a[i]->a[1],I);
      }
    }
  }
  f=func_del(f);
  return g;
}

//////////////////////////////////////////////////////////////////

void func_change_cmulti(cmulti *x, func_t *f)            //(cmulti)x=(func)f
{
  func_t *g=NULL;
  g=func_evalf(FR(f));
  if(g==NULL)                 { cset_d(x,0.0); }
  else if(func_is_complex(g)) { cset_c(x,g->p.cm); }
  else if(func_is_real(g))    { cset_r(x,g->p.rm); }
  else{ FUNC_ERROR_ARG1("func_init_val_radius",f); }
  g=func_del(g);
}

func_t *func_change_cmulti_list(int n, cmulti **x)      //(list)x=(cmulti)x
{
  int i;
  func_t *f=NULL,*g=NULL,*h=NULL;
  f=func_list(0);
  g=func_list(n);
  h=func_list(n);
  for(i=0;i<n;i++){
    g->a[i]=func_var1(i,1);   
    h->a[i]=func_complex();
    cset_c(h->a[i]->p.cm,x[i]);
  }
  func_a_append(f,FR(g));
  func_a_append(f,FR(h));
  g=func_del(g);
  h=func_del(h);
  return f;
}

void func_ccopy_coeff(int n, cmulti **z, func_t *f)                  //係数
{
  int i;
  func_t *fa_0=NULL,*g=NULL;
  cmulti *a_0=NULL;
  if(func_is_poly(f)){
    a_0=callocate();
    fa_0=func_poly_clone_coeff_ntarm(f,n);           //a_0の作成
    func_change_cmulti(a_0,fa_0);
    for(i=1; i<=n; i++){
      if(n-i!=0){
	g=func_poly_clone_coeff_ntarm(f,n-i);
	func_change_cmulti(z[i-1],g);
      }else{
	g=func_number_pull(FR(f));
	func_change_cmulti(z[i-1],g);     
      }
      cdiv_cc(z[i-1],z[i-1],a_0);
      g=func_del(g); g=NULL; 
    }
  }
  fa_0=func_del(fa_0);
  g=func_del(g);
  a_0=cfree(a_0);
}

//////////////////////////////////////////////////////////////////
//Weierstrass法
//////////////////////////////////////////////////////////////////
void func_init_val_r_and_balance(rmulti *r, cmulti *g, int n, cmulti **a)  //半径r(小澤の初期値)と重心
{
  double N;
  cdiv_cd(g,a[0],(double)(n));       //重心
  cneg_c(g,g);
  N=1.0/(double)(n);
  rabs_c(r,a[n-1]);                   //半径
  rpow_rd(r,r,N);
}

//////////////////////////////////////////////////////////////////

void func_init_val(int n, cmulti **z, cmulti **a)                         //初期値
{
  int j;
  double theta;
  rmulti *r=NULL,*b=NULL;
  cmulti *g=NULL,*c=NULL;
  r=rallocate();
  b=rallocate();
  g=callocate();
  c=callocate();
  func_init_val_r_and_balance(r,g,n,a);
  for(j=1; j<=n; j++){
    theta=(2.0*M_PI*(j-1))/n+3.0/(2.0*n);                      // (2π(j-1))/n + 3/2n;
    rset_d(b,theta);
    cset_polar(c,r,b);                                         //rexp(i*theta)
    cadd_cc(z[j-1],g,c);
  }
  r=rfree(r);
  b=rfree(b);
  g=cfree(g);
  c=cfree(c);
}

//////////////////////////////////////////////////////////////////

int func_weierstrass(int n, cmulti **z1, cmulti **z0, cmulti **a)
{
  int i,j,k=0;
  double N;
  rmulti *d=NULL,*eps=NULL;
  cmulti *b=NULL,**c=NULL;
  d=rallocate();
  eps=rmepsilon(get_default_prec()/2);
  b=callocate();
  c=cvec_allocate(n);
  cvec_set_ones(n,c);
  N=(double)(n);
  for(i=0; i<n; i++){                   // f(z)
    cpow_cd(z1[i],z0[i],N);               //z^n
    cadd_cc(z1[i],z1[i],a[n-1]);             //z^n+a_n
    for(j=1; j<n; j++){
      cpow_cd(b,z0[i],N-(double)(j));     //a*z^n-1
      cmul_cc(b,b,a[j-1]);
      cadd_cc(z1[i],z1[i],b);
    }
  }
  for(i=0; i<n; i++){                   //終了条件 |f(z)|<eps
    rabs_c(d,z1[i]);
    k=lt_rr(d,eps);                       //d<eps
    if(k==0){ break; }
  }
  for(i=0; i<n-1; i++){                 //(z[i]-z[j])*...
    for(j=i+1; j<n; j++){
      csub_cc(b,z0[i],z0[j]);
      cmul_cc(c[i],c[i],b);
      cneg_c(b,b);
      cmul_cc(c[j],c[j],b);
    }
  }
  cvec_div_cvec_cvec(n,z1,z1,c);
  cvec_mul_cvec_dscalar(n,z1,z1,-1.0);             //-f()/a_0()
  cvec_add_cvec_cvec(n,z1,z1,z0);                 //z-f()/a_0()
  d=rfree(d);
  eps=rfree(eps);
  b=cfree(b);
  c=cvec_free(n,c);
  return k;
}

/////////////////////////////////////////////////////////////////////////////////


//EOF
