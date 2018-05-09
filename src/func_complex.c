#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#define FR(F) func_retain(F)

////////////////////////////////////////////////

void func_print_complex(func_t *f)
{
  mpfr_printf("(%.5Rg%+.5Rg*I)",C_R(func_complex_p(f)),C_I(func_complex_p(f)));
}

////////////////////////////////////////////////

void func_complex_p_new(func_t *f)
{
  if(f->p.cm!=NULL){ func_complex_p_del(f); }
  f->ptype=FUNC_P_COMPLEX;
  f->p.cm=callocate();
}

void func_complex_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_COMPLEX){ FUNC_ERROR_ARG1("func_complex_p_del",f); }
  if(f->p.cm!=NULL){
    f->p.cm=cfree(f->p.cm);
    f->ptype=FUNC_P_NULL;
  }
}

void func_complex_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_COMPLEX){ FUNC_ERROR_ARG2("func_complex_p_clone",f,g); }
  if(f->p.cm!=NULL){ func_complex_p_del(f); }
  f->ptype=FUNC_P_COMPLEX;
  f->p.cm=callocate_clone(g->p.cm);
}

cmulti *func_complex_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_COMPLEX){ FUNC_ERROR_ARG1("func_complex_p",f); }
  return f->p.cm;
}

//////////////////////////////////////////////

func_t *func_complex_eval(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_real(f))      { return func_real_eval(f); }
  else if(!func_is_complex(f))  { FUNC_ERROR_ARG1("func_complex_eval",f); }
  else if(cis_nan(f->p.cm))     { f=func_del(f); return func_nan(); }
  else if(cis_inf(f->p.cm))     { f=func_del(f); return func_inf(); }
  else if(!cis_number(f->p.cm)) { f=func_del(f); return func_nan(); }
  else if(cis_zero(f->p.cm))    { f=func_del(f); return func_zero(); }
  else if(ceq_si(f->p.cm,1))    { f=func_del(f); return func_one(); }
  else if(cis_real(f->p.cm)){
    g=func_real_rmulti(C_R(f->p.cm));
    f=func_del(f);
    return g;
  }
  return f;
}

///////////////////////////////////////////

static char *__func_complex="@C";

func_t *func_complex(void)
{
  return func_new(__func_complex);
}

func_t *func_complex_cmulti(cmulti *x)
{
  func_t *f=NULL;
  f=func_new(__func_complex);
  ccopy(f->p.cm,x);
  return f;
}

///////////////////////////////////////////////////////////

func_t *func_complex_i(void)
{
  func_t *f=NULL;
  f=func_complex();
  cset_Z(func_complex_p(f),0,1);
  return f;
}

///////////////////////////////////////////////////////////

func_t *func_complex_add(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_complex(f1) || !func_in_complex(f2)){ FUNC_ERROR_ARG2("func_complex_add",f1,f2); }
  a=func_complex_cast(FR(f1));
  b=func_complex_cast(FR(f2));
  if(!func_is_complex(a) || !func_is_complex(b)){ FUNC_ERROR_ARG2("func_complex_add",f1,f2); }
  f=func_complex();
  cadd(func_complex_p(f),func_complex_p(a),func_complex_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_complex_eval(f);
}

func_t *func_complex_sub(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_complex(f1) || !func_in_complex(f2)){ FUNC_ERROR_ARG2("func_complex_sub",f1,f2); }
  a=func_complex_cast(FR(f1));
  b=func_complex_cast(FR(f2));
  if(!func_is_complex(a) || !func_is_complex(b)){ FUNC_ERROR_ARG2("func_complex_sub",f1,f2); }
  f=func_complex();
  csub(func_complex_p(f),func_complex_p(a),func_complex_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_complex_eval(f);
}

func_t *func_complex_mul(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_complex(f1) || !func_in_complex(f2)){ FUNC_ERROR_ARG2("func_complex_mul",f1,f2); }
  a=func_complex_cast(FR(f1));
  b=func_complex_cast(FR(f2));
  if(!func_is_complex(a) || !func_is_complex(b)){ FUNC_ERROR_ARG2("func_complex_mul",f1,f2); }
  f=func_complex();
  cmul(func_complex_p(f),func_complex_p(a),func_complex_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_complex_eval(f);
}

func_t *func_complex_div(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_complex(f1) || !func_in_complex(f2)){ FUNC_ERROR_ARG2("func_complex_div",f1,f2); }
  a=func_complex_cast(FR(f1));
  b=func_complex_cast(FR(f2));
  if(!func_is_complex(a) || !func_is_complex(b)){ FUNC_ERROR_ARG2("func_complex_div",f1,f2); }
  f=func_complex();
  cdiv(func_complex_p(f),func_complex_p(a),func_complex_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_complex_eval(f);
}

func_t *func_complex_pow(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_complex(f1) || !func_in_complex(f2)){ FUNC_ERROR_ARG2("func_complex_pow",f1,f2); }
  a=func_complex_cast(FR(f1));
  b=func_complex_cast(FR(f2));
  if(!func_is_complex(a) || !func_is_complex(b)){ FUNC_ERROR_ARG2("func_complex_pow",f1,f2); }
  f=func_complex();
  cpow_c(func_complex_p(f),func_complex_p(a),func_complex_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_complex_eval(f);
}

func_t *func_complex_inv(func_t *f1)
{
  func_t *f=NULL,*a=NULL;
  if(!func_in_complex(f1)){ FUNC_ERROR_ARG1("func_complex_inv",f1); }
  a=func_complex_cast(FR(f1));
  if(!func_is_complex(a)){ FUNC_ERROR_ARG1("func_complex_inv",f1); }
  f=func_complex();
  cinv(func_complex_p(f),func_complex_p(a));
  f1=func_del(f1);
  a=func_del(a);
  return func_complex_eval(f);
}

func_t *func_complex_pow_n(func_t *f1, int p)
{
  func_t *f=NULL,*a=NULL;
  if(!func_in_complex(f1)){ FUNC_ERROR_ARG1("func_complex_pow_n",f1); }
  a=func_complex_cast(FR(f1));
  if(!func_is_complex(a)){ FUNC_ERROR_ARG1("func_complex_pow_n",f1); }
  f=func_complex();
  cpow_si(func_complex_p(f),func_complex_p(a),p);
  f1=func_del(f1);
  a=func_del(a);
  return func_complex_eval(f);
}

///////////////////////////////////////////////////////////////////////

int func_complex_cmp(func_t *f, func_t *g)
{
  int value;
  func_t *a=NULL,*b=NULL;
  if     (func_in_bigint(f) && func_in_bigint(g)){ return func_bigint_cmp(f,g); }
  else if(func_in_real(f)   && func_in_real(g))  { return func_real_cmp(f,g); }
  if(!func_in_complex(f) || !func_in_complex(g)){ FUNC_ERROR_ARG2("func_complex_cmp",f,g); }
  a=func_complex_cast(FR(f));
  b=func_complex_cast(FR(g));
  value=ccmp(func_complex_p(a),func_complex_p(b));
  a=func_del(a);
  b=func_del(b);
  return value;
}

/////////////////////////////////////////////////////////////

func_t *func_complex_cast(func_t *g)
{
  func_t *f=NULL;
  if     (func_is(g,"nan"))  { f=func_complex(); cset_nan(func_complex_p(f)); }
  else if(func_is(g,"inf"))  { f=func_complex(); cset_inf(func_complex_p(f),0,0); }
  else if(func_is_zero(g))   { f=func_complex(); cset_zero(func_complex_p(f)); }
  else if(func_is_one(g))    { f=func_complex(); cset_one(func_complex_p(f)); }
  else if(func_is_bigint(g)) { f=func_complex(); cset_zero(func_complex_p(f)); bigint_get_rmulti(C_R(func_complex_p(f)),func_bigint_p(g)); }
  else if(func_is_real(g))   { f=func_complex(); ccopy_r(func_complex_p(f),func_real_p(g)); }
  else if(func_is_complex(g)){ f=FR(g); }
  g=func_del(g);
  return f;
}

///////////////////////////////////////////

int func_is_complex(func_t *f)
{
  return func_is(f,__func_complex);
}

int func_in_complex(func_t *f)
{
  return (f!=NULL && (func_is_complex(f) || func_in_real(f)));
}

///////////////////////////////////////////////////

func_t *func_op_complex_new(void)
{
  func_t *f=func_builtin_new(__func_complex);
  func_builtin_p(f)->order=FUNC_ORDER_COMPLEX;
  func_builtin_p(f)->type=FUNC_SCALAR;
  func_builtin_p(f)->ptype=FUNC_P_COMPLEX;
  func_builtin_p(f)->p_new=func_complex_p_new;
  func_builtin_p(f)->p_del=func_complex_p_del;
  func_builtin_p(f)->p_clone=func_complex_p_clone;
  func_builtin_p(f)->p_cmp=func_complex_cmp;
  func_builtin_p(f)->eval=func_complex_eval;
  func_builtin_p(f)->print=func_print_complex;
  return f;
}

//EOF
