#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F) func_retain(F)

//////////////////////////////////////////////////////////

#define SPC " \t\n"

func_t *func_real_script(func_t *f)
{
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if     (str_match(func_strings_at(f,0),"^%p?%d+.%d*$",NULL,NULL))        { g=func_real_str(func_strings_at(f,0)); }
    else if(str_match(func_strings_at(f,0),"^%p?%d*.%d+$",NULL,NULL))        { g=func_real_str(func_strings_at(f,0)); }
    else if(str_match(func_strings_at(f,0),"^%p?%d+.?%d*e%p?%d+$",NULL,NULL)){ g=func_real_str(func_strings_at(f,0)); }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  return g;
}

////////////////////////////////////////////////

void func_real_print(func_t *f)
{
  if(rget_sgn(func_real_p(f))<0){
    mpfr_printf("(%.5Rg)",func_real_p(f));
  }else{
    mpfr_printf("(%.5Rg)",func_real_p(f));
  }
}

////////////////////////////////////////////////

void func_real_p_new(func_t *f)
{
  if(f->p.rm!=NULL){ func_real_p_del(f); }
  f->ptype=FUNC_P_REAL;
  f->p.rm=rallocate();
}

void func_real_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_REAL){ FUNC_ERROR_ARG1("func_real_p_del",f); }
  if(f->p.rm!=NULL){
    f->ptype=FUNC_P_NULL;
    f->p.rm=rfree(f->p.rm);
  }
}

void func_real_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_REAL){ FUNC_ERROR_ARG2("func_real_p_clone",f,g); }
  if(f->p.rm!=NULL){ func_real_p_del(f); }
  f->ptype=FUNC_P_REAL;
  f->p.rm=rallocate_clone(g->p.rm);
}

rmulti *func_real_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_REAL){ FUNC_ERROR_ARG1("func_real_p",f); }
  return f->p.rm;
}

//////////////////////////////////////////////////

func_t *func_real_eval(func_t *f)
{
  if     (!func_is_real(f))           { FUNC_ERROR_ARG1("func_real_eval",f); }
  else if(ris_nan(func_real_p(f)))    { f=func_del(f); return func_nan(); }
  else if(ris_inf(func_real_p(f)))    { f=func_del(f); return func_inf(); }
  else if(!ris_number(func_real_p(f))){ f=func_del(f); return func_nan(); }
  else if(ris_zero(func_real_p(f)))   { f=func_del(f); return func_zero(); }
  else if(req_si(func_real_p(f),1))   { f=func_del(f); return func_one(); }
  return f;
}


//////////////////////////////////////////////////

static char *__func_real="@R";

func_t *func_real(void)
{
  return func_new(__func_real);
}

func_t *func_real_rmulti(rmulti *x)
{
  func_t *f=NULL;
  f=func_new(__func_real);
  rcopy(func_real_p(f),x);
  return func_real_eval(f);
}

func_t *func_real_str(char *str)
{
  func_t *f=NULL;
  f=func_new(__func_real);
  rset_s(func_real_p(f),str);
  return func_real_eval(f);
}

////////////////////////////////////////////////

//////////////////////

func_t *func_real_get_bigint(func_t *f)
{
  func_t *g=NULL;
  if(!func_is_real(f)){ FUNC_ERROR_ARG1("func_real_get_bigint",f); }
  g=func_bigint();
  rcopy(BIGINT_NUM(g->p.bi),func_real_p(f));
  rset_si(BIGINT_DEN(g->p.bi),1);
  return func_bigint_eval(g);
}

///////////////////////////////////////////////////////////////////////

int func_real_cmp(func_t *f, func_t *g)
{
  int value;
  func_t *a=NULL,*b=NULL;
  if(func_in_bigint(f) && func_in_bigint(g)){ return func_bigint_cmp(f,g); }
  if(!func_in_real(f) || !func_in_real(g)){ FUNC_ERROR_ARG2("func_real_cmp",f,g); }
  a=func_real_cast(FR(f));
  b=func_real_cast(FR(g));
  if(!func_is_real(a) || !func_is_real(b)){ FUNC_ERROR_ARG2("func_real_cmp",f,g); }
  value=rcmp(func_real_p(a),func_real_p(b));
  a=func_del(a);
  b=func_del(b);
  return value;
}

///////////////////////////////////////////////////////////

func_t *func_real_add(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_real(f1) || !func_in_real(f2)){ FUNC_ERROR_ARG2("func_real_add",f1,f2); }
  a=func_real_cast(FR(f1));
  b=func_real_cast(FR(f2));
  if(!func_is_real(a) || !func_is_real(b)){ FUNC_ERROR_ARG2("func_real_add",f1,f2); }
  f=func_real();
  radd(func_real_p(f),func_real_p(a),func_real_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_real_eval(f);
}

func_t *func_real_sub(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_real(f1) || !func_in_real(f2)){ FUNC_ERROR_ARG2("func_real_sub",f1,f2); }
  a=func_real_cast(FR(f1));
  b=func_real_cast(FR(f2));
  if(!func_is_real(a) || !func_is_real(b)){ FUNC_ERROR_ARG2("func_real_sub",f1,f2); }
  f=func_real();
  rsub(func_real_p(f),func_real_p(a),func_real_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_real_eval(f);
}

func_t *func_real_mul(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_real(f1) || !func_in_real(f2)){ FUNC_ERROR_ARG2("func_real_mul",f1,f2); }
  a=func_real_cast(FR(f1));
  b=func_real_cast(FR(f2));
  if(!func_is_real(a) || !func_is_real(b)){ FUNC_ERROR_ARG2("func_real_mul",f1,f2); }
  f=func_real();
  rmul(func_real_p(f),func_real_p(a),func_real_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_real_eval(f);
}

func_t *func_real_div(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_real(f1) || !func_in_real(f2)){ FUNC_ERROR_ARG2("func_real_div",f1,f2); }
  a=func_real_cast(FR(f1));
  b=func_real_cast(FR(f2));
  if(!func_is_real(a) || !func_is_real(b)){ FUNC_ERROR_ARG2("func_real_div",f1,f2); }
  f=func_real();
  rdiv(func_real_p(f),func_real_p(a),func_real_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_real_eval(f);
}

func_t *func_real_pow(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_real(f1) || !func_in_real(f2)){ FUNC_ERROR_ARG2("func_real_pow",f1,f2); }
  a=func_real_cast(FR(f1));
  b=func_real_cast(FR(f2));
  if(!func_is_real(a) || !func_is_real(b)){ FUNC_ERROR_ARG2("func_real_pow",f1,f2); }
  f=func_real();
  rpow(func_real_p(f),func_real_p(a),func_real_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_real_eval(f);
}


func_t *func_real_inv(func_t *f1)
{
  func_t *f=NULL,*a=NULL;
  if(!func_in_real(f1)){ FUNC_ERROR_ARG1("func_real_inv",f1); }
  a=func_real_cast(FR(f1));
  if(!func_is_real(a)){ FUNC_ERROR_ARG1("func_real_inv",f1); }
  f=func_real();
  rinv(func_real_p(f),func_real_p(a));
  f1=func_del(f1);
  a=func_del(a);
  return func_real_eval(f);
}

func_t *func_real_pow_n(func_t *f1, int p)
{
  func_t *f=NULL,*a=NULL;
  if(!func_in_real(f1)){ FUNC_ERROR_ARG1("func_real_pow_n",f1); }
  a=func_real_cast(FR(f1));
  if(!func_is_real(a)){ FUNC_ERROR_ARG1("func_real_pow_n",f1); }
  f=func_real();
  rpow_si(func_real_p(f),func_real_p(a),p);
  f1=func_del(f1);
  a=func_del(a);
  return func_real_eval(f);
}

/////////////////////////////////////////////////////////////

func_t *func_real_cast(func_t *g)
{
  func_t *f=NULL;
  if     (func_is(g,"nan")) { f=func_real(); rset_nan(func_real_p(f)); }
  else if(func_is(g,"inf")) { f=func_real(); rset_inf(func_real_p(f),1); }
  else if(func_is_zero(g))  { f=func_real(); rset_zero(func_real_p(f)); }
  else if(func_is_one(g))   { f=func_real(); rset_one(func_real_p(f)); }
  else if(func_is_bigint(g)){ f=func_real(); bigint_get_rmulti(func_real_p(f),func_bigint_p(g)); }
  else if(func_is_real(g))  { f=FR(g); }
  g=func_del(g);
  return f;
}

/////////////////////////////////////////////////////////////

int func_is_real(func_t *f)
{
  return func_is(f,__func_real);
}

int func_in_real(func_t *f)
{
  return (f!=NULL && (func_is_real(f) || func_in_bigint(f)));
}

////////////////////////////////////////////////////////

func_t *func_op_real_new(void)
{
  func_t *f=func_builtin_new(__func_real);
  func_builtin_p(f)->order=FUNC_ORDER_REAL;
  func_builtin_p(f)->type=FUNC_SCALAR;
  func_builtin_p(f)->ptype=FUNC_P_REAL;
  func_builtin_p(f)->p_new=func_real_p_new;
  func_builtin_p(f)->p_del=func_real_p_del;
  func_builtin_p(f)->p_clone=func_real_p_clone;
  func_builtin_p(f)->p_cmp=func_real_cmp;
  func_builtin_p(f)->eval=func_real_eval;
  func_builtin_p(f)->print=func_real_print;
  return f;
}

//EOF
