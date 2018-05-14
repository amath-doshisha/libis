#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F) (func_retain(F))

///////////////////////////////////////////////////////

func_t *func_sub(func_t *f1, func_t *f2)
{
  func_t *f=NULL;
  if     (func_in_bigint(f1)  && func_in_bigint(f2)) { f=func_bigint_sub(FR(f1),FR(f2)); }
  else if(func_in_real(f1)    && func_in_real(f2))   { f=func_real_sub(FR(f1),FR(f2)); }
  else if(func_in_complex(f1) && func_in_complex(f2)){ f=func_complex_sub(FR(f1),FR(f2)); }
  else if(func_in_vec(f1)     && func_in_vec(f2))    { f=func_vec_sub(FR(f1),FR(f2)); }
  else if(func_in_mat(f1)     && func_in_mat(f2))    { f=func_mat_sub(FR(f1),FR(f2)); }
  else                                               { f=func_add(FR(f1),func_mul(func_bigint_int(-1,1),FR(f2))); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"

func_t *func_sub_script(func_t *f)
{
  int i;
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"-",BLK0,BLK1)){
      a=strings_split(func_strings_at(f,0),"-",BLK0,BLK1,SPC);
      if(strings_size(a)>0){
	g=func_add_new(strings_size(a));
	for(i=0; i<func_asize(g); i++){
	  if(i==0 && (*func_strings_at(f,0))!='-'){
	    func_aset(g,i,func_script(strings_at(a,i)));
	  }else{
	    func_aset(g,i,func_mul(func_bigint_int(-1,1),func_script(strings_at(a,i))));
	  }
	}
	g=func_add_eval(g);
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}

///////////////////////////////////////////////////////

func_t *func_inv(func_t *f)
{
  if     (func_is_zero(f))   { f=func_del(f); return func_inf(); }
  else if(func_is_one(f))    { return f; }
  else if(func_is_bigint(f)) { return func_bigint_inv(f); }
  else if(func_is_real(f))   { return func_real_inv(f); }
  else if(func_is_complex(f)){ return func_complex_inv(f); }
  else return func_pow(f,func_bigint_int(-1,1));
}

func_t *func_div(func_t *f1, func_t *f2)
{
  func_t *f=NULL;
  if     (func_is_zero(f1)   && func_is_zero(f2))    { f=func_nan(); }
  else if(func_is_zero(f1))                          { f=func_zero(); }
  else if(func_is_zero(f2))                          { f=func_inf(); }
  else if(func_is_one(f1))                           { f=func_inv(FR(f2)); }
  else if(func_is_one(f2))                           { f=FR(f1); }
  else if(func_in_bigint(f1)  && func_in_bigint(f2)) { f=func_bigint_div(FR(f1),FR(f2)); }
  else if(func_in_real(f1)    && func_in_real(f2))   { f=func_real_div(FR(f1),FR(f2)); }
  else if(func_in_complex(f1) && func_in_complex(f2)){ f=func_complex_div(FR(f1),FR(f2)); }
  else if(func_is_var(f1)     && func_is_var(f2))    { f=func_var_div(FR(f1),FR(f2)); }
  else                                               { f=func_mul(FR(f1),func_inv(FR(f2))); }
  f1=func_del(f1);
  f2=func_del(f2); 
  return f;
}

func_t *func_div_script(func_t *f)
{
  int i;
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"/",BLK0,BLK1)){
      a=strings_split(func_strings_at(f,0),"/",BLK0,BLK1,SPC);
      if(strings_size(a)>=2){
	g=func_script(strings_at(a,0));
	for(i=1; i<strings_size(a); i++){
	  g=func_div(g,func_script(strings_at(a,i)));
	}
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}


///////////////////////////////////////////////////////

func_t *func_pow_n(func_t *f, int pow)
{
  func_t *g=NULL;
  if     (f==NULL)           { FUNC_ERROR_ARG1("func_pow_n",f); }
  else if(func_is_nan(f))    { g=FR(f); }
  else if(func_is_inf(f))    { g=FR(f); }
  else if(func_is_one(f))    { g=FR(f); }
  else if(func_is_zero(f))   { if(pow>0) { g=FR(f); } else if(pow==0){ g=func_one(); } else { g=func_inf(); } }
  else if(pow==0)            { g=func_one(); }
  else if(pow==1)            { g=FR(f); }
  else if(func_is_bigint(f)) { g=func_bigint_pow_n(FR(f),pow); }
  else if(func_is_real(f))   { g=func_real_pow_n(FR(f),pow); }
  else if(func_is_complex(f)){ g=func_complex_pow_n(FR(f),pow); }
  else if(func_is_var(f))    { g=func_var_pow_n(FR(f),pow); }
  else if(func_is_mul(f))    { g=func_mul_pow_n(FR(f),pow); }
  else if(func_is(f,"sqrt")) { g=func_sqrt_pow_n(FR(f),pow); }
  else if(func_is(f,"pow"))  { g=func_new("pow"); g->a[0]=FR(f->a[0]); g->a[1]=func_mul(FR(f->a[1]),func_bigint_int(pow,1)); }
  else if(func_has_power(f)) { g=func_clone(FR(f)); func_set_power(g,func_power(g)*pow); }
  else                       { g=func_new("pow"); g->a[0]=FR(f); g->a[1]=func_bigint_int(pow,1); }
  f=func_del(f);
  return g;
}

func_t *func_pow_script(func_t *f)
{
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"%^",BLK0,BLK1)){
      a=strings_split(func_strings_at(f,0),"^",BLK0,BLK1,SPC);
      if(strings_size(a)==2){
	g=func_pow(func_script(strings_at(a,0)),func_script(strings_at(a,1)));
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}

void func_pow_print(func_t *f)
{
  printf("(");
  func_print(f->a[0]);
  printf(")^");
  printf("(");
  func_print(f->a[1]);
  printf(")");
}

func_t *func_op_pow_new(void)
{
  func_t *f=func_builtin_new("pow");
  func_builtin_p(f)->order=FUNC_ORDER_POW;
  func_builtin_p(f)->amin=2;
  func_builtin_p(f)->amax=2;
  func_builtin_p(f)->ptype=FUNC_P_NULL;
  func_builtin_p(f)->eval=func_pow_eval;
  func_builtin_p(f)->print=func_pow_print;
  return f;
}

func_t *func_pow(func_t *a, func_t *b)
{
  return func_pow_eval(func_arg2_new("pow",a,b));
}

func_t *func_pow_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL,*b=NULL;
  if(!func_is(f,"pow") || func_asize(f)!=2){ FUNC_ERROR_ARG1("func_pow_eval",f); }
  a=func_aget(f,0);
  b=func_aget(f,1);
  if     (func_is_zero(b))                          { g=func_pow_n(FR(a),0); }
  else if(func_is_one(b))                           { g=func_pow_n(FR(a),1); }
  else if(func_bigint_is_integer(b))                { g=func_pow_n(FR(a),func_bigint_get_si(b)); }
  else if(func_is_real(a)    && func_is_real(b))    { g=func_real_pow(FR(a),FR(b)); }
  else if(func_is_real(a)    && func_in_real(b))    { g=func_real_pow(FR(a),FR(b)); }
  else if(func_in_real(a)    && func_is_real(b))    { g=func_real_pow(FR(a),FR(b)); }
  else if(func_is_complex(a) && func_is_complex(b)) { g=func_complex_pow(FR(a),FR(b)); }
  else if(func_is_complex(a) && func_in_complex(b)) { g=func_complex_pow(FR(a),FR(b)); }
  else if(func_in_complex(a) && func_is_complex(b)) { g=func_complex_pow(FR(a),FR(b)); }
  else if(func_is_real(a)    && func_is_number(b))  { g=func_pow(FR(a),func_evalf(FR(b))); }
  else if(func_is_complex(a) && func_is_number(b))  { g=func_pow(FR(a),func_evalf(FR(b))); }
  else if(func_is_number(a)  && func_is_real(b))    { g=func_pow(func_evalf(FR(a)),FR(b)); }
  else if(func_is_number(a)  && func_is_complex(b)) { g=func_pow(func_evalf(FR(a)),FR(b)); }
  else                                              { g=FR(f); }
  f=func_del(f);
  return g;  
}

/////////////////////////////////////////////////////

func_t *func_sqrt(func_t *f)
{
  return func_sqrt_eval(func_arg1_new("sqrt",f));
}

func_t *func_sqrt_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"sqrt")){ FUNC_ERROR_ARG1("func_sqrt_eval",f); }
  a=func_aget(f,0);
  if     (func_is_zero(a))                        { g=func_zero(); }
  else if(func_is_one(a))                         { g=func_one(); }
  else if(func_is_real(a) && rget_sgn(a->p.rm)>=0){ g=func_real();    rsqrt_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_real(a) && rget_sgn(a->p.rm)<0) { g=func_complex(); csqrt_r(g->p.cm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a))                     { g=func_complex(); csqrt_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); } 
  else                                            { g=FR(f); }
  f=func_del(f);
  return g;
}

func_t *func_sqrt_pow_n(func_t *f, int pow)
{
  int p;
  func_t *g=NULL;
  if(func_is(f,"sqrt")){
    p=func_power(f)*pow;
    if(p==0){ g=func_one(); }
    else if(abs(p)%2==0){ g=func_pow_n(FR(f->a[0]),p/2); }
    else if(p>0)        { g=func_mul(func_pow_n(FR(f->a[0]),p/2),func_sqrt(FR(f->a[0]))); }
    else                { g=func_mul(func_pow_n(FR(f->a[0]),(p-1)/2),func_sqrt(FR(f->a[0]))); }
  }else{ FUNC_ERROR_ARG1("func_sqrt_pow_n",f); }
  f=func_del(f);
  return g;
}

//////////////////////////////////////////////////////////////

func_t *func_exp(func_t *f)
{
  return func_exp_eval(func_arg1_new("exp",f));
}

func_t *func_exp_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"exp")){ FUNC_ERROR_ARG1("func_exp_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_one(); }
  else if(func_is_real(a))   { g=func_real();    rexp_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); cexp_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

/////////////////////////////////////////////////////////////

func_t *func_log(func_t *f)
{
  return func_log_eval(func_arg1_new("log",f));
}

func_t *func_log_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"log")){ FUNC_ERROR_ARG1("func_log_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a))                        { g=func_inf(); }
  else if(func_is_one(a))                         { g=func_zero(); }
  else if(func_is_real(a) && rget_sgn(a->p.rm)==0){ g=func_inf(); }
  else if(func_is_real(a) && rget_sgn(a->p.rm)>0) { g=func_real();    rlog_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_real(a) && rget_sgn(a->p.rm)<0) { g=func_complex(); clog_r(g->p.cm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a))                     { g=func_complex(); clog_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                                            { g=FR(f); }
  f=func_del(f);
  return g;
}

/////////////////////////////////////////////////////////////

func_t *func_sin(func_t *f)
{
  return func_sin_eval(func_arg1_new("sin",f));
}

func_t *func_sin_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"sin")){ FUNC_ERROR_ARG1("func_sin_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_zero(); }
  else if(func_is_real(a))   { g=func_real();    rsin_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); csin_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

///////////////////////////////////////////////////////////

func_t *func_cos(func_t *f)
{
  return func_cos_eval(func_arg1_new("cos",f));
}

func_t *func_cos_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"cos")){ FUNC_ERROR_ARG1("func_cos_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_one();  }
  else if(func_is_real(a))   { g=func_real();    rcos_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); ccos_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

///////////////////////////////////////////////////////////////

func_t *func_tan(func_t *f)
{
  return func_tan_eval(func_arg1_new("tan",f));
}

func_t *func_tan_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"tan")){ FUNC_ERROR_ARG1("func_tan_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_zero();  }
  else if(func_is_real(a))   { g=func_real();    rtan_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); ctan_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

////////////////////////////////////////////////////////////

func_t *func_asin(func_t *f)
{
  return func_asin_eval(func_arg1_new("asin",f));
}

func_t *func_asin_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"asin")){ FUNC_ERROR_ARG1("func_asin_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_zero(); }
  else if(func_is_real(a) && le_rd(a->p.rm,1) && ge_rd(a->p.rm,-1))
                             { g=func_real();    rasin_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_real(a))   { g=func_complex(); casin_r(g->p.cm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); casin_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

//////////////////////////////////////////////////////////////

func_t *func_acos(func_t *f)
{
  return func_acos_eval(func_arg1_new("acos",f));
}

func_t *func_acos_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"acos")){ FUNC_ERROR_ARG1("func_acos_eval",f); }
  a=f->a[0];
  if     (func_is_one(a))    { g=func_zero(); }
  else if(func_is_real(a) && le_rd(a->p.rm,1) && ge_rd(a->p.rm,-1))
                             { g=func_real();    racos_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_real(a))   { g=func_complex(); cacos_r(g->p.cm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); cacos_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

////////////////////////////////////////////////////////////////////////

func_t *func_atan(func_t *f)
{
  return func_atan_eval(func_arg1_new("atan",f));
}

func_t *func_atan_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"atan")){ FUNC_ERROR_ARG1("func_atan_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_zero(); }
  else if(func_is_real(a))   { g=func_real();    ratan_r  (g->p.rm,   a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); catan_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

/////////////////////////////////////////////////////////////////////////////

func_t *func_sinh(func_t *f)
{
  return func_sinh_eval(func_arg1_new("sinh",f));
}

func_t *func_sinh_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"sinh")){ FUNC_ERROR_ARG1("func_sinh_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_zero();  }
  else if(func_is_real(a))   { g=func_real();    rsinh_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); csinh_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

//////////////////////////////////////////////////////////////////////////////

func_t *func_cosh(func_t *f)
{
  return func_cosh_eval(func_arg1_new("cosh",f));
}

func_t *func_cosh_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"cosh")){ FUNC_ERROR_ARG1("func_cosh_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_one(); }
  else if(func_is_real(a))   { g=func_real();    rcosh_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); ccosh_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

/////////////////////////////////////////////////////////////////////////

func_t *func_tanh(func_t *f)
{
  return func_tanh_eval(func_arg1_new("tanh",f));
}

func_t *func_tanh_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"tanh")){ FUNC_ERROR_ARG1("func_tanh_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_zero(); }
  else if(func_is_real(a))   { g=func_real();    rtanh_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); ctanh_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

/////////////////////////////////////////////////////////////////////////////

func_t *func_asinh(func_t *f)
{
  return func_asinh_eval(func_arg1_new("asinh",f));
}

func_t *func_asinh_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"asinh")){ FUNC_ERROR_ARG1("func_asinh_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a)) { g=func_zero(); }
  else if(func_is_real(a))   { g=func_real();    rasinh_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); casinh_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

//////////////////////////////////////////////////////////////////////

func_t *func_acosh(func_t *f)
{
  return func_acosh_eval(func_arg1_new("acosh",f));
}

func_t *func_acosh_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"acosh")){ FUNC_ERROR_ARG1("func_acosh_eval",f); }
  a=f->a[0];
  if     (func_is_one(a))    { g=func_zero(); }
  else if(func_is_real(a) && ge_rd(a->p.rm,1))
                             { g=func_real();    racosh_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_real(a))   { g=func_complex(); cacosh_r(g->p.cm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); cacosh_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

///////////////////////////////////////////////////////////////////

func_t *func_atanh(func_t *f)
{
  return func_atanh_eval(func_arg1_new("atanh",f));
}

func_t *func_atanh_eval(func_t *f)
{
  func_t *g=NULL,*a=NULL;
  if(!func_is(f,"atanh")){ FUNC_ERROR_ARG1("func_atanh_eval",f); }
  a=f->a[0];
  if     (func_is_zero(a))   { g=func_zero(); }
  else if(func_is_real(a) && le_rd(a->p.rm,1) && ge_rd(a->p.rm,-1))
                             { g=func_real();    ratanh_r  (g->p.rm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_real(a))   { g=func_complex(); catanh_r(g->p.cm,a->p.rm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else if(func_is_complex(a)){ g=func_complex(); catanh_c(g->p.cm,a->p.cm); g=func_eval(g); g=func_pow_n(g,func_power(f)); }
  else                       { g=FR(f); }
  f=func_del(f);
  return g;
}

//EOF
