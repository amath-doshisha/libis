#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_rvec.h"
#include"is_cvec.h"
#include"is_irmulti.h"
#include"is_icmulti.h"

#define FR(F) func_retain(F)
#define FR(F) func_retain(F)
#define RAp(X,Y)   ((X)=rallocate_prec(rget_prec(Y)))
#define RAP(X,Y,N) ((X)=rallocate_prec(rvec_get_prec_max(N,Y)))
#define RF(X)      ((X)=rfree(X))
#define RVAp(X,Y,N){ X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)   { X=rvec_free(N,X); }
#define RA(X,Y)      ((X)=rallocate_prec(rget_prec(Y)))
#define RA2(X,Y0,Y1) ((X)=rallocate_prec(MAX2(rget_prec(Y0),rget_prec(Y1))))
#define RF(X)        ((X)=rfree(X))
#define RVA(X,Y,N)   { X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)     { X=rvec_free(N,X); }
#define RF(X)      ((X)=rfree(X))
#define CA(X,P)    ((X)=callocate_prec(P))
#define CF(X)      ((X)=cfree(X))
#define CVA(X,N,P) { X=cvec_allocate_prec(N,P); }
#define CVF(X,N)   { X=cvec_free(N,X); }
#define RVA(X,Y,N) { X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)   { X=rvec_free(N,X); }

///////////////////////////////////////////////////////////////

func_t *func_maps_var(func_t *f, int n0, func_t *x)
{
  int i,k,p;
  func_t *g=NULL,*h=NULL;
  if(!func_is_var(f)){ FUNC_ERROR_ARG1("func_maps_var",f); }
  g=func_clone(FR(f));
  if(func_is_list(x)){
    h=func_one();
    for(i=0; i<func_asize(x); i++){
      k=func_var_get_index(g,n0+i);
      if(k>=0){
	p=g->p.var->pow[k];
	g->p.var->pow[k]=0;
	h=func_mul(h,func_pow_n(FR(x->a[i]),p));
      }
    }
    g=func_var_check_mul(g);
    g=func_mul(g,FR(h));
    h=func_del(h);
  }else if(func_is(x,"rvec")){
    h=func_one();
    for(i=0; i<func_rvec_size(x); i++){
      k=func_var_get_index(g,n0+i);
      if(k>=0){
	p=g->p.var->pow[k];
	g->p.var->pow[k]=0;
	h=func_mul(h,func_pow_n(func_rvec_get(x,i),p));
      }
    }
    g=func_var_check_mul(g);
    g=func_mul(g,FR(h));
    h=func_del(h);
  }else if(func_is(x,"cvec")){
    h=func_one();
    for(i=0; i<func_cvec_size(x); i++){
      k=func_var_get_index(g,n0+i);
      if(k>=0){
	p=g->p.var->pow[k];
	g->p.var->pow[k]=0;
	h=func_mul(h,func_pow_n(func_cvec_get(x,i),p));
      }
    }
    g=func_var_check_mul(g);
    g=func_mul(g,FR(h));
    h=func_del(h);
  }else{
    k=func_var_get_index(g,n0);
    if(k>=0){
      p=g->p.var->pow[k];
      g->p.var->pow[k]=0;
      g=func_var_check_mul(g);
      g=func_mul(g,func_pow_n(FR(x),p));
    }
  }
  f=func_del(f);
  x=func_del(x);
  return g;
}

// y=f(x)
func_t *func_maps(func_t *f, int n0, func_t *x)
{
  int i;
  func_t *g=NULL;
  func_arg1_t *eval=NULL;
  if(f==NULL)            { FUNC_ERROR_ARG1("func_maps",f); }
  else if(func_is_var(f)){ g=func_maps_var(FR(f),n0,FR(x)); }
  else{
    g=func_clone(FR(f));
    for(i=0; i<func_asize(f); i++){
      func_aset(g,i,func_maps(FR(func_aget(f,i)),n0,FR(x)));
    }
    eval=func_find_eval(g);
    if(eval!=NULL){ g=eval(g); }
  }
  f=func_del(f);
  x=func_del(x);
  return g;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name rmulti型の写像に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルに関する写像 y=f(x)
 */
void rvec_func(rmulti *y, func_t *f, int n, rmulti **x)
{
  int i;
  rmulti *a=NULL,*b=NULL,*z=NULL;
  RAp(a,y); RAp(b,y); RAp(z,y);
  if(f==NULL)                { FUNC_ERROR_ARG1("rvec_func",f); }
  else if(func_is(f,"nan"))  { rset_nan(z); }
  else if(func_is(f,"inf"))  { rset_inf(z,1); }
  else if(func_is_zero(f))   { rset_i(z,0); }
  else if(func_is_one(f))    { rset_i(z,1); }
  else if(func_is_bigint(f)) { bigint_get_rmulti(z,func_bigint_p(f)); }
  else if(func_is_real(f))   { rset_r(z,func_real_p(f)); }
  else if(func_is_complex(f)){ rset_i(a,0); rset_i(b,0); rdiv_rr(z,a,b); }
  else if(func_is_var(f))    {
    rset_d(z,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
        if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
          rpow_r(a,x[func_var_num(f,i)],func_var_pow(f,i));
        }else{ rset_nan(a); }
        rmul_rr(z,z,a);
      }
    }
  }
  else if(func_is_add(f))    { rset_d(z,0); for(i=0; i<func_asize(f); i++){ rvec_func(a,func_aget(f,i),n,x); radd_rr(z,z,a); } }
  else if(func_is_mul(f))    { rset_d(z,1); for(i=0; i<func_asize(f); i++){ rvec_func(a,func_aget(f,i),n,x); rmul_rr(z,z,a); } }
  else if(func_is(f,"sqrt")) { rvec_func(a,func_aget(f,0),n,x); rsqrt_r(z,a); }
  else if(func_is(f,"exp"))  { rvec_func(a,func_aget(f,0),n,x); rexp_r(z,a); }
  else if(func_is(f,"log"))  { rvec_func(a,func_aget(f,0),n,x); rlog_r(z,a); }
  else if(func_is(f,"sin"))  { rvec_func(a,func_aget(f,0),n,x); rsin_r(z,a); }
  else if(func_is(f,"cos"))  { rvec_func(a,func_aget(f,0),n,x); rcos_r(z,a); }
  else if(func_is(f,"tan"))  { rvec_func(a,func_aget(f,0),n,x); rtan_r(z,a); }
  else if(func_is(f,"asin")) { rvec_func(a,func_aget(f,0),n,x); rasin_r(z,a); }
  else if(func_is(f,"acos")) { rvec_func(a,func_aget(f,0),n,x); racos_r(z,a); }
  else if(func_is(f,"atan")) { rvec_func(a,func_aget(f,0),n,x); ratan_r(z,a); }
  else if(func_is(f,"sinh")) { rvec_func(a,func_aget(f,0),n,x); rsinh_r(z,a); }
  else if(func_is(f,"cosh")) { rvec_func(a,func_aget(f,0),n,x); rcosh_r(z,a); }
  else if(func_is(f,"tanh")) { rvec_func(a,func_aget(f,0),n,x); rtanh_r(z,a); }
  else if(func_is(f,"asinh")){ rvec_func(a,func_aget(f,0),n,x); rasinh_r(z,a); }
  else if(func_is(f,"acosh")){ rvec_func(a,func_aget(f,0),n,x); racosh_r(z,a); }
  else if(func_is(f,"atanh")){ rvec_func(a,func_aget(f,0),n,x); ratanh_r(z,a); }
  else if(func_is(f,"pow"))  { rvec_func(a,func_aget(f,0),n,x); rvec_func(b,func_aget(f,1),n,x); rpow_rr(z,a,b); }
  else if(func_is_list(f))   { rset_nan(a); }
  else if(func_is_rvec(f))   { rset_nan(a); }
  else if(func_is_cvec(f))   { rset_nan(a); }
  else if(func_is_rmat(f))   { rset_nan(a); }
  else if(func_is_cmat(f))   { rset_nan(a); }
  else                       { rset_nan(a); }
  if(func_has_power(f))      { rpow_r(z,z,func_power(f)); }
  rset_r(y,z);
  RF(a); RF(b); RF(z);
}

/**
 @brief rmulti型のベクトルに関するベクトル写像 y=f(x)
 */
void rvec_func_list(int m, rmulti **y, func_t *f, int n, rmulti **x)
{
  int i;
  for(i=0; i<m; i++){
    if(func_is_list(f) && i<func_asize(f)){
      rvec_func(y[i],func_aget(f,i),n,x);
    }else{
      rset_nan(y[i]);
    }
  }
}

/**
 @brief rmulti型に関する写像 y=f(x0)
 */
void r1_func(rmulti *y, func_t *f, rmulti *x0)
{
  int n=1;
  rmulti **x=NULL;
  x=rvec_allocate_prec(n,rget_prec(y));
  rclone_r(x[0],x0);
  rvec_func(y,f,n,x);
  x=rvec_free(n,x);
}

/**
 @brief rmulti型に関する写像 y=f(x0,x1)
 */
void r2_func(rmulti *y, func_t *f, rmulti *x0, rmulti *x1)
{
  int n=2;
  rmulti **x=NULL;
  x=rvec_allocate_prec(n,rget_prec(y));
  rclone_r(x[0],x0);
  rclone_r(x[1],x1);
  rvec_func(y,f,n,x);
  x=rvec_free(n,x);
}

/**
 @brief rmulti型に関する写像 y=f(x0,x1,x2)
 */
void r3_func(rmulti *y, func_t *f, rmulti *x0, rmulti *x1, rmulti *x2)
{
  int n=3;
  rmulti **x=NULL;
  x=rvec_allocate_prec(n,rget_prec(y));
  rclone_r(x[0],x0);
  rclone_r(x[1],x1);
  rclone_r(x[2],x2);
  rvec_func(y,f,n,x);
  x=rvec_free(n,x);
}

/** @} */

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name rmulti型の写像に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルに関する写像 y=f(x)
 */
void cvec_func(cmulti *y, func_t *f, int n, cmulti **x)
{
  int i;
  cmulti *a=NULL,*b=NULL,*z=NULL;
  a=callocate_prec(cget_prec(y));
  b=callocate_prec(cget_prec(y));
  z=callocate_prec(cget_prec(y));
  if(f==NULL)                { FUNC_ERROR_ARG1("cvec_fun",f); }
  else if(func_is(f,"nan"))  { cset_nan(z); }
  else if(func_is(f,"inf"))  { cset_inf(z,1,1); }
  else if(func_is_zero(f))   { cset_i(z,0); }
  else if(func_is_one(f))    { cset_i(z,1); }
  else if(func_is_bigint(f)) { bigint_get_cmulti(z,func_bigint_p(f)); }
  else if(func_is_real(f))   { cset_r(z,func_real_p(f)); }
  else if(func_is_complex(f)){ cset_c(z,func_complex_p(f)); }
  else if(func_is_var(f))    {
    cset_d(z,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
        if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
          cpow_c(a,x[func_var_num(f,i)],func_var_pow(f,i));
        }else{ cset_nan(a); }
        cmul_cc(z,z,a);
      }
    }
  }
  else if(func_is_add(f))    { cset_d(z,0); for(i=0; i<func_asize(f); i++){ cvec_func(a,func_aget(f,i),n,x); cadd_cc(z,z,a); } }
  else if(func_is_mul(f))    { cset_d(z,1); for(i=0; i<func_asize(f); i++){ cvec_func(a,func_aget(f,i),n,x); cmul_cc(z,z,a); } }
  else if(func_is(f,"sqrt")) { cvec_func(a,func_aget(f,0),n,x); csqrt_c(z,a); }
  else if(func_is(f,"exp"))  { cvec_func(a,func_aget(f,0),n,x); cexp_c(z,a); }
  else if(func_is(f,"log"))  { cvec_func(a,func_aget(f,0),n,x); clog_c(z,a); }
  else if(func_is(f,"sin"))  { cvec_func(a,func_aget(f,0),n,x); csin_c(z,a); }
  else if(func_is(f,"cos"))  { cvec_func(a,func_aget(f,0),n,x); ccos_c(z,a); }
  else if(func_is(f,"tan"))  { cvec_func(a,func_aget(f,0),n,x); ctan_c(z,a); }
  else if(func_is(f,"asin")) { cvec_func(a,func_aget(f,0),n,x); casin_c(z,a); }
  else if(func_is(f,"acos")) { cvec_func(a,func_aget(f,0),n,x); cacos_c(z,a); }
  else if(func_is(f,"atan")) { cvec_func(a,func_aget(f,0),n,x); catan_c(z,a); }
  else if(func_is(f,"sinh")) { cvec_func(a,func_aget(f,0),n,x); csinh_c(z,a); }
  else if(func_is(f,"cosh")) { cvec_func(a,func_aget(f,0),n,x); ccosh_c(z,a); }
  else if(func_is(f,"tanh")) { cvec_func(a,func_aget(f,0),n,x); ctanh_c(z,a); }
  else if(func_is(f,"asinh")){ cvec_func(a,func_aget(f,0),n,x); casinh_c(z,a); }
  else if(func_is(f,"acosh")){ cvec_func(a,func_aget(f,0),n,x); cacosh_c(z,a); }
  else if(func_is(f,"atanh")){ cvec_func(a,func_aget(f,0),n,x); catanh_c(z,a); }
  else if(func_is(f,"pow"))  { cvec_func(a,func_aget(f,0),n,x); cvec_func(b,func_aget(f,1),n,x); cpow_cc(z,a,b); }
  else if(func_is_list(f))   { cset_nan(a); }
  else if(func_is_rvec(f))   { cset_nan(a); }
  else if(func_is_cvec(f))   { cset_nan(a); }
  else if(func_is_rmat(f))   { cset_nan(a); }
  else if(func_is_cmat(f))   { cset_nan(a); }
  else                       { cset_nan(a); }
  if(func_has_power(f))      { cpow_c(z,z,func_power(f)); }
  cset_c(y,z);
  a=cfree(a);
  b=cfree(b);
  z=cfree(z);
}

/**
 @brief cmulti型のベクトルに関するベクトル写像 y=f(x)
 */
void cvec_func_list(int m, cmulti **y, func_t *f, int n, cmulti **x)
{
  int i;
  if(func_is_list(f)){
    for(i=0; i<m && i<func_asize(f); i++){
      cvec_func(y[i],func_aget(f,i),n,x);
    }
  }
}

/**
 @brief cmulti型に関する写像 y=f(x0)
 */
void c1_func(cmulti *y, func_t *f, cmulti *x0)
{
  int n=1;
  cmulti **x=NULL;
  x=cvec_allocate_prec(n,cget_prec(y));
  cset_c(x[0],x0);
  cvec_func(y,f,n,x);
  x=cvec_free(n,x);
}

/**
 @brief cmulti型に関する写像 y=f(x0,x1)
 */
void c2_func(cmulti *y, func_t *f, cmulti *x0, cmulti *x1)
{
  int n=2;
  cmulti **x=NULL;
  x=cvec_allocate_prec(n,cget_prec(y));
  cset_c(x[0],x0);
  cset_c(x[1],x1);
  cvec_func(y,f,n,x);
  x=cvec_free(n,x);
}

/**
 @brief cmulti型に関する写像 y=f(x0,x1,x2)
 */
void c3_func(cmulti *y, func_t *f, cmulti *x0, cmulti *x1, cmulti *x2)
{
  int n=3;
  cmulti **x=NULL;
  x=cvec_allocate_prec(n,cget_prec(y));
  cset_c(x[0],x0);
  cset_c(x[1],x1);
  cset_c(x[2],x2);
  cvec_func(y,f,n,x);
  x=cvec_free(n,x);
}


/**
 @brief irmulti型のベクトルの写像 [y0,y1]=f([x0,x1])
 */
void irvec_func(rmulti *y0, rmulti *y1, func_t *f, int n, rmulti **x0, rmulti **x1)
{
  int i;
  rmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL,*z0=NULL,*z1=NULL;
  RA2(a0,y0,y1); RA2(a1,y0,y1); RA2(b0,y0,y1); RA2(b1,y0,y1); RA2(z0,y0,y1); RA2(z1,y0,y1);
  rset_nan(z0); rset_nan(z1);
  if(f==NULL)                { FUNC_ERROR_ARG1("irvec_func",f); }
  else if(func_is(f,"inf"))  { rset_inf(z0,1); rset_inf(z1,1); }
  else if(func_is_zero(f))   { irset_d(z0,z1,0,0); }
  else if(func_is_one(f))    { irset_d(z0,z1,1,1); }
  else if(func_is_bigint(f)) { irset_bigint(z0,z1,func_bigint_p(f)); }
  else if(func_is_real(f))   { irset_r(z0,z1,func_real_p(f),func_real_p(f)); }
  else if(func_is_var(f))    {
    irset_d(z0,z1,1,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
        if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
          irpow_r(a0,a1,x0[func_var_num(f,i)],x1[func_var_num(f,i)],func_var_pow(f,i));
        }else{ rset_nan(a0); rset_nan(a1); }
        irmul_rr(z0,z1,z0,z1,a0,a1);
      }
    }
  }
  else if(func_is_add(f))    {
    irset_d(z0,z1,0,0);
    for(i=0; i<func_asize(f); i++){
      irvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      iradd_rr(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is_mul(f))    {
    irset_d(z0,z1,1,1);
    for(i=0; i<func_asize(f); i++){
      irvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      irmul_rr(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is(f,"sqrt")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irsqrt_r(z0,z1,a0,a1); }
  else if(func_is(f,"pow"))  {
    irvec_func(a0,a1,func_aget(f,0),n,x0,x1);
    irvec_func(b0,b1,func_aget(f,1),n,x0,x1);
    rpow_rr(z0,a0,b0); ERROR_AT;
  }
  if(func_has_power(f))      { irpow_r(z0,z1,z0,z1,func_power(f)); }
  irset_r(y0,y1,z0,z1);
  RF(a0); RF(a1); RF(b0); RF(b1); RF(z0); RF(z1);
}

/**
 @brief irmulti型のベクトルのベクトル写像 [y0,y1]=f([x0,x1])
 */
void irvec_func_list(int m, rmulti **y0, rmulti **y1, func_t *f, int n, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<m; i++){
    if(func_is_list(f) && func_is_list(f) && i<func_asize(f)){
      irvec_func(y0[i],y1[i],func_aget(f,i),n,x0,x1);
    }else{
      rset_nan(y0[i]); rset_nan(y1[i]);
    }
  }
}

/* @} */


/**
 @brief 写像 [y0,y1]=f([x0,x1])
 */
void icvec_func(cmulti *y0, cmulti *y1, func_t *f, int n, cmulti **x0, cmulti **x1)
{
  int prec,i;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL,*z0=NULL,*z1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec); CA(z0,prec); CA(z1,prec);
  cset_nan(z0); cset_nan(z1);
  if(f==NULL)                { FUNC_ERROR_ARG1("icvec_func",f); }
  else if(func_is(f,"nan"))  { cset_nan(z0); cset_nan(z1); }
  else if(func_is(f,"inf"))  { cset_inf(z0,1,1); cset_inf(z1,1,1); }
  else if(func_is_zero(f))   { icset_d(z0,z1,0,0); }
  else if(func_is_one(f))    { icset_d(z0,z1,1,1); }
  else if(func_is_bigint(f)) { icset_bigint(z0,z1,func_bigint_p(f)); }
  else if(func_is_real(f))   { cset_r(z0,func_real_p(f)); cset_r(z1,func_real_p(f)); }
  else if(func_is_complex(f)){ cset_c(z0,func_complex_p(f)); cset_c(z1,func_complex_p(f)); }
  else if(func_is_var(f))    {
    icset_d(z0,z1,1,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
        if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
          icpow_c(a0,a1,x0[func_var_num(f,i)],x1[func_var_num(f,i)],func_var_pow(f,i));
        }else{ cset_nan(a0); cset_nan(a1); }
        icmul_cc(z0,z1,z0,z1,a0,a1);
      }
    }
  }
  else if(func_is_add(f))  {
    icset_d(z0,z1,0,0);
    for(i=0; i<func_asize(f); i++){
      icvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      icadd_cc(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is_mul(f))  {
    icset_d(z0,z1,1,1);
    for(i=0; i<func_asize(f); i++){
      icvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      icmul_cc(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is(f,"pow"))  {
    icvec_func(a0,a1,func_aget(f,0),n,x0,x1);
    icvec_func(b0,b1,func_aget(f,1),n,x0,x1);
    cpow_cc(z0,a0,b0); ERROR_AT;
  }
  if(func_has_power(f))      { icpow_c(z0,z1,z0,z1,func_power(f)); }
  icset_c(y0,y1,z0,z1);
  CF(a0); CF(a1); CF(b0); CF(b1); CF(z0); CF(z1);
}

/**
 @brief ベクトル写像 [y0,y1]=f([x0,x1])
 */
void icvec_func_list(int m, cmulti **y0, cmulti **y1, func_t *f, int n, cmulti **x0, cmulti **x1)
{
  int i;
  if(func_is_list(f)){
    for(i=0; i<m && i<func_asize(f); i++){
      icvec_func(y0[i],y1[i],func_aget(f,i),n,x0,x1);
    }
  }
}




//EOF
