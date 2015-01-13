#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F) func_retain(F)

/////////////////////////////////////////////////////////////////////////////////

func_t *func_diff_var(func_t *f, int var)
{
  int n,m;
  func_t *fx=NULL;
  if(!func_is_var(f)){ FUNC_ERROR_ARG1("func_diff_var",f); }
  n=ivec_first_index_if(func_var_size(f),f->p.var->num,var);
  if(n<0) { fx=func_zero(); }
  else{
    fx=func_clone(FR(f));
    m=fx->p.var->pow[n];
    fx->p.var->pow[n]--;
    fx=func_var_check_mul(fx);
    fx=func_mul(func_bigint_int(m,1),fx);
  }
  f=func_del(f);
  return fx;
}

/////////////////////////////////////////////////////////////////////////////////

func_t *func_diff_pow(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*h=NULL,*gx=NULL,*hx=NULL;
  if(!func_is(f,"pow")){ FUNC_ERROR_ARG1("func_diff_pow",f); }
  g=f->a[0];
  h=f->a[1];
  gx=func_diff(FR(g),func_var1(var,1));
  hx=func_diff(FR(h),func_var1(var,1));
  hx=func_mul(hx,func_log(FR(g)));
  gx=func_mul(gx,FR(h));
  gx=func_div(gx,FR(g));
  gx=func_add(gx,hx);
  fx=func_mul(FR(f),gx);
  f=func_del(f);
  return fx;
}

/////////////////////////////////////////////////////////////////////////////////

func_t *func_diff_pow_n(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL;
  if(!func_has_power(f) || func_power(f)==1){ FUNC_ERROR_ARG1("func_diff_pow_n",f); }
  g=func_clone(FR(f));
  func_set_power(g,1);
  g=func_diff(g,func_var1(var,1));
  fx=func_clone(FR(f));
  func_set_power(fx,func_power(fx)-1);
  fx=func_mul(func_bigint_int(func_power(f),1),fx);
  fx=func_mul(fx,g);
  f=func_del(f);
  return fx;
}

func_t *func_diff_add_pow1(func_t *f, int var)
{
  int i;
  func_t *fx=NULL;
  if(!func_is_add(f) || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_add",f); }
  fx=func_zero();
  for(i=0; i<func_asize(f); i++){
    fx=func_add(fx,func_diff(FR(f->a[i]),func_var1(var,1)));
  }
  f=func_del(f);
  return fx;
}

func_t *func_diff_mul_pow1(func_t *f, int var)
{
  int i;
  func_t *fx=NULL,*g=NULL;
  if(!func_is_mul(f) || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_mul",f); }
  fx=func_zero();
  for(i=0; i<func_asize(f); i++){
    g=func_clone(FR(f));
    g->a[i]=func_del(g->a[i]);
    g=func_mul_eval(g);
    fx=func_add(fx,func_mul(g,func_diff(FR(f->a[i]),func_var1(var,1))));
  }
  f=func_del(f);
  return fx;
}

func_t *func_diff_sqrt_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"sqrt") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_sqrt",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_div(gx,func_sqrt(FR(g))),func_bigint_int(1,2));
  f=func_del(f);
  return fx;
}

func_t *func_diff_exp_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"exp") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_exp",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_exp(FR(g)),gx);
  f=func_del(f);
  return fx;
}

func_t *func_diff_log_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"log") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_log",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_div(gx,FR(g));
  f=func_del(f);
  return fx;
}

func_t *func_diff_sin_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"sin") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_sin",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_cos(FR(g)),gx);
  f=func_del(f);
  return fx;
}

func_t *func_diff_cos_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"cos") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_cos",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_mul(func_sin(FR(g)),func_bigint_int(-1,1)),gx);
  f=func_del(f);
  return fx;
}

func_t *func_diff_tan_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"tan") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_tan",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_pow_n(func_cos(FR(g)),-2),gx);
  f=func_del(f);
  return fx;
}

func_t *func_diff_sinh_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"sinh") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_sinh",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_cosh(FR(g)),gx);
  f=func_del(f);
  return fx;
}

func_t *func_diff_cosh_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"cosh") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_cosh",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_sinh(FR(g)),gx);
  f=func_del(f);
  return fx;
}

func_t *func_diff_tanh_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL;
  if(!func_is(f,"tanh") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_tanh",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  fx=func_mul(func_pow_n(func_cosh(FR(g)),-2),gx);
  f=func_del(f);
  return fx;
}

func_t *func_diff_asin_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL,*h=NULL;
  if(!func_is(f,"asin") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_asin",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  h=func_pow_n(FR(g),2);
  h=func_sub(func_one(),h);
  h=func_sqrt(h);
  fx=func_div(gx,h);
  f=func_del(f);
  return fx;
}

func_t *func_diff_acos_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL,*h=NULL;
  if(!func_is(f,"acos") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_acosh",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  h=func_pow_n(FR(g),2);
  h=func_sub(func_one(),h);
  h=func_sqrt(h);
  fx=func_div(gx,h);
  fx=func_mul(fx,func_bigint_int(-1,1));
  f=func_del(f);
  return fx;
}

func_t *func_diff_atan_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL,*h=NULL;
  if(!func_is(f,"atan") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_atan",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  h=func_pow_n(FR(g),2);
  h=func_add(h,func_one());
  fx=func_div(gx,h);
  f=func_del(f);
  return fx;
}

func_t *func_diff_asinh_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL,*h=NULL;
  if(!func_is(f,"asinh") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_asinh",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  h=func_pow_n(FR(g),2);
  h=func_add(h,func_one());
  h=func_sqrt(h);
  fx=func_div(gx,h);
  f=func_del(f);
  return fx;
}

func_t *func_diff_acosh_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL,*h=NULL;
  if(!func_is(f,"acosh") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_acosh",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  h=func_pow_n(FR(g),2);
  h=func_add(h,func_bigint_int(-1,1));
  h=func_sqrt(h);
  fx=func_div(gx,h);
  f=func_del(f);
  return fx;
}

func_t *func_diff_atanh_pow1(func_t *f, int var)
{
  func_t *fx=NULL,*g=NULL,*gx=NULL,*h=NULL;
  if(!func_is(f,"atanh") || func_power(f)!=1){ FUNC_ERROR_ARG1("func_diff_atanh",f); }
  g=f->a[0];
  gx=func_diff(FR(g),func_var1(var,1));
  h=func_pow_n(FR(g),2);
  h=func_sub(func_one(),h);
  fx=func_div(gx,h);
  f=func_del(f);
  return fx;
}

////////////////////////////////////////////////////////////////////////////////

func_t *func_diff_list(func_t *f, int var)
{
  int i;
  func_t *fx=NULL;
  if(!func_is_list(f)){ FUNC_ERROR_ARG1("func_diff_list",f); }
  fx=func_list(func_asize(f));
  for(i=0; i<func_asize(f); i++){
    fx->a[i]=func_diff(FR(f->a[i]),func_var1(var,1));
  }
  f=func_del(f);
  return fx;
}

///////////////////////////////////////////////

func_t *func_op_diff_new(void)
{
  func_t *f=func_builtin_new("diff");
  func_builtin_p(f)->amin=2;
  func_builtin_p(f)->amax=2;
  func_builtin_p(f)->eval=func_diff_eval;
  return f;
}

func_t *func_diff(func_t *f, func_t *x)
{
  return func_diff_eval(func_arg2_new("diff",f,x));
}

func_t *func_diff_eval(func_t *g)
{
  int var;
  func_t *fx=NULL,*f=NULL,*x=NULL;
  if(g==NULL || !func_is(g,"diff") || func_asize(g)!=2){ FUNC_ERROR_ARG1("func_diff_eval",g); }
  f=func_aget(g,0);
  x=func_aget(g,1);
  if(f==NULL || x==NULL || !func_is_1var(x)){ FUNC_ERROR_ARG2("func_diff_eval",f,x); }
  var=func_var_num(x,0);
  if     (func_is_zero(f))                       { fx=func_zero(); }
  else if(func_is_one(f))                        { fx=func_zero(); }
  else if(func_is_bigint(f))                     { fx=func_zero(); }
  else if(func_is_real(f))                       { fx=func_zero(); }
  else if(func_is_complex(f))                    { fx=func_zero(); }
  else if(func_is_var(f))                        { fx=func_diff_var(FR(f),var); }
  else if(func_is(f,"pow"))                      { fx=func_diff_pow(FR(f),var); }
  else if(func_is_list(f))                       { fx=func_diff_list(FR(f),var); }
  else if(func_is(f,"rvec"))                     { fx=func_list_zeros(func_rvec_size(f)); }
  else if(func_is(f,"cvec"))                     { fx=func_list_zeros(func_cvec_size(f)); }
  else if(func_is(f,"rmat"))                     { fx=func_list_zeros2(func_rmat_rows(f),func_rmat_cols(f)); }
  else if(func_is(f,"cmat"))                     { fx=func_list_zeros2(func_cmat_rows(f),func_cmat_cols(f)); }
  else if(func_is_add(f)     && func_power(f)==1){ fx=func_diff_add_pow1(FR(f),var); }   // (add)^1
  else if(func_is_mul(f)     && func_power(f)==1){ fx=func_diff_mul_pow1(FR(f),var); }   // (mul)^1  
  else if(func_is(f,"sqrt")  && func_power(f)==1){ fx=func_diff_sqrt_pow1(FR(f),var); }  // sqrt^1
  else if(func_is(f,"exp")   && func_power(f)==1){ fx=func_diff_exp_pow1(FR(f),var); }   // exp^1
  else if(func_is(f,"log")   && func_power(f)==1){ fx=func_diff_log_pow1(FR(f),var); }   // log^1
  else if(func_is(f,"sin")   && func_power(f)==1){ fx=func_diff_sin_pow1(FR(f),var); }   // sin^1
  else if(func_is(f,"cos")   && func_power(f)==1){ fx=func_diff_cos_pow1(FR(f),var); }   // cos^1
  else if(func_is(f,"tan")   && func_power(f)==1){ fx=func_diff_tan_pow1(FR(f),var); }   // tan^1
  else if(func_is(f,"sinh")  && func_power(f)==1){ fx=func_diff_sinh_pow1(FR(f),var); }  // sinh^1
  else if(func_is(f,"cosh")  && func_power(f)==1){ fx=func_diff_cosh_pow1(FR(f),var); }  // cosh^1
  else if(func_is(f,"tanh")  && func_power(f)==1){ fx=func_diff_tanh_pow1(FR(f),var); }  // tanh^1
  else if(func_is(f,"asin")  && func_power(f)==1){ fx=func_diff_asin_pow1(FR(f),var); }  // asin^1
  else if(func_is(f,"acos")  && func_power(f)==1){ fx=func_diff_acos_pow1(FR(f),var); }  // acos^1
  else if(func_is(f,"atan")  && func_power(f)==1){ fx=func_diff_atan_pow1(FR(f),var); }  // atan^1
  else if(func_is(f,"asinh") && func_power(f)==1){ fx=func_diff_asinh_pow1(FR(f),var); } // asinh^1
  else if(func_is(f,"acosh") && func_power(f)==1){ fx=func_diff_acosh_pow1(FR(f),var); } // acosh^1
  else if(func_is(f,"atanh") && func_power(f)==1){ fx=func_diff_atanh_pow1(FR(f),var); } // atanh^1
  else if(func_has_power(f) && func_power(f)!=1) { fx=func_diff_pow_n(FR(f),var); }      // (???)^p with p!=1
  else                                           { fx=FR(g); }                           // otherwise
  g=func_del(g);
  return fx;
}

///////////////////////////////////////////////////////////////

func_t *func_op_grad_new(void)
{
  func_t *f=func_builtin_new("grad");
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=2;
  func_builtin_p(f)->eval=func_grad_eval;
  return f;
}

func_t *func_grad(func_t *f, func_t *x)
{
  if(x==NULL){ return func_grad_eval(func_arg1_new("grad",f)); }
  else       { return func_grad_eval(func_arg2_new("grad",f,x)); }
}

func_t *func_grad_eval(func_t *h)
{
  int j,n=0;
  func_t *g=NULL,*f=NULL,*x=NULL;
  if(h==NULL || !func_is(h,"grad") || func_asize(h)<1 || func_asize(h)>2){ FUNC_ERROR_ARG1("func_grad_eval",h); }  
  f=func_aget(h,0);
  if(func_asize(h)<2 && func_is_list(f)){ x=func_bigint_int(func_asize(f),1); }
  else if(func_asize(h)>1)              { x=FR(func_aget(h,1)); }
  if     (func_is_list(x))              { n=func_asize(x); }
  else if(func_bigint_is_integer(x))    { n=func_bigint_get_si(x); }
  else                                  { FUNC_ERROR_ARG2("func_grad",f,x); }
  if(func_is_list(f)){
    g=func_list(func_asize(f));
    for(j=0; j<func_asize(f); j++){
      func_aset(g,j,func_grad(FR(f->a[j]),FR(x)));
    }    
  }else{
    g=func_list(n);
    for(j=0; j<n; j++){
      if(func_is_list(x)){ func_aset(g,j,func_diff(FR(f),FR(func_aget(x,j)))); }
      else               { func_aset(g,j,func_diff(FR(f),func_var1(j,1))); }
    }
  }
  x=func_del(x);
  h=func_del(h);
  return g;
}

//EOF
