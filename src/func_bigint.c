#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"is_bigint.h"

#define FR(F) func_retain(F)

//////////////////////////////////////////////////////////

#define SPC " \t\n"

func_t *func_bigint_script(func_t *f)
{
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"^-?%d+$",NULL,NULL)){ // Z
      g=func_bigint_str(func_strings_at(f,0),"1");
    }else if(str_match(func_strings_at(f,0),"^(-?%d+)$",NULL,NULL)){ // Z
      a=strings_split(func_strings_at(f,0),"()",NULL,NULL,SPC);
      if(strings_size(a)==1){ g=func_bigint_str(strings_at(a,0),"1"); }
    }else if(str_match(func_strings_at(f,0),"^-?%d+/%d+$",NULL,NULL)){ // Q
      a=strings_split(func_strings_at(f,0),"/",NULL,NULL,SPC);
      if(strings_size(a)==2){ g=func_bigint_str(strings_at(a,0),strings_at(a,1)); }
    }else if(str_match(func_strings_at(f,0),"^(-?%d+/%d+)$",NULL,NULL)){ // Q
      a=strings_split(func_strings_at(f,0),"/()",NULL,NULL,SPC);
      if(a->n==2){ g=func_bigint_str(strings_at(a,0),strings_at(a,1)); }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}

//////////////////////////////////////////////////////////

void func_bigint_print(func_t *f)
{
  bigint_print(func_bigint_p(f));
}

//////////////////////////////////////////////////////////

void func_bigint_p_new(func_t *f)
{
  if(f->p.bi!=NULL){ func_bigint_p_del(f); }
  f->ptype=FUNC_P_BIGINT;
  f->p.bi=bigint_allocate();
}

void func_bigint_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_BIGINT){ FUNC_ERROR_ARG1("func_bigint_p_del",f); }
  if(f->p.bi!=NULL){
    f->p.bi=bigint_free(f->p.bi);
    f->ptype=FUNC_P_NULL;
  }
}

void func_bigint_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_BIGINT){ FUNC_ERROR_ARG2("func_bigint_p_new",f,g); }
  if(f->p.bi!=NULL){ func_bigint_p_del(f); }
  f->ptype=FUNC_P_BIGINT;
  f->p.bi=bigint_allocate_clone(g->p.bi);
}

bigint *func_bigint_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BIGINT){ FUNC_ERROR_ARG1("func_bigint_p",f); }
  return f->p.bi;
}

/////////////////////////////////////////////////////

func_t *func_bigint_eval(func_t *f)
{
  if     (!func_is_bigint(f))      { FUNC_ERROR_ARG1("func_bigint_eval",f); }
  else if(bigint_is_nan(f->p.bi))  { f=func_del(f); return func_nan(); }
  else if(bigint_is_inf(f->p.bi))  { f=func_del(f); return func_inf(); }
  else if(bigint_is_zero(f->p.bi)) { f=func_del(f); return func_zero(); }
  else if(bigint_is_one(f->p.bi))  { f=func_del(f); return func_one(); }
  else                             { return f; }
}

/////////////////////////////////////////////////////

static char *__func_bigint="@Q";

func_t *func_bigint(void)
{
  return func_new(__func_bigint);
}

func_t *func_bigint_int(int num, int den)
{
  func_t *f=NULL;
  f=func_new(__func_bigint);
  bigint_set_int(f->p.bi,num,den);
  return func_bigint_eval(f);
}

func_t *func_bigint_str(char *num, char *den)
{
  func_t *f=NULL;
  f=func_new(__func_bigint);
  bigint_set_str(f->p.bi,num,den);
  return func_bigint_eval(f);
}

////////////////////////////////////////////////////////

int func_bigint_is_integer(func_t *f)
{
  return (func_is_bigint(f) && bigint_is_integer(f->p.bi));
}

int func_bigint_sgn(func_t *f)
{
  return bigint_sgn(f->p.bi);
}

int func_bigint_get_si(func_t *f)
{
  if(!func_bigint_is_integer(f)){ FUNC_ERROR_ARG1("func_bigint_int",f); }
  return bigint_get_si(f->p.bi);
}

////////////////////////////////////////////////////////////////////////////////////

int func_is_bigint(func_t *f)
{
  return func_is(f,__func_bigint);
}

int func_in_bigint(func_t *f)
{
  return (func_is_zero(f) || func_is_one(f) || func_is_bigint(f));
}

func_t *func_bigint_cast(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_zero(f))  { g=func_bigint(); bigint_set_zero(g->p.bi); }
  else if(func_is_one(f))   { g=func_bigint(); bigint_set_one(g->p.bi); }
  else if(func_is_bigint(f)){ g=FR(f); }
  f=func_del(f);
  return g;
}

////////////////////////////////////////////////////////

int func_bigint_cmp(func_t *f, func_t *g)
{
  int value;
  func_t *a=NULL,*b=NULL;
  if(!func_in_bigint(f) || !func_in_bigint(g)){ FUNC_ERROR_ARG2("func_bigint_cmp",f,g); }
  a=func_bigint_cast(FR(f));
  b=func_bigint_cast(FR(g));
  if(!func_is_bigint(a) || !func_is_bigint(b)){ FUNC_ERROR_ARG2("func_bigint_cmp",f,g); }
  value=bigint_cmp(func_bigint_p(a),func_bigint_p(b));
  a=func_del(a);
  b=func_del(b);
  return value;
}

////////////////////////////////////////////////////////

// f=f1+f2
func_t *func_bigint_add(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_bigint(f1) || !func_in_bigint(f2)){ FUNC_ERROR_ARG2("func_bigint_add",f1,f2); }
  a=func_bigint_cast(FR(f1));
  b=func_bigint_cast(FR(f2));
  if(!func_is_bigint(a) || !func_is_bigint(b)){ FUNC_ERROR_ARG2("func_bigint_add",f1,f2); }
  f=func_bigint();
  bigint_add(func_bigint_p(f),func_bigint_p(a),func_bigint_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_bigint_eval(f);
}

// f=f1-f2
func_t *func_bigint_sub(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_bigint(f1) || !func_in_bigint(f2)){ FUNC_ERROR_ARG2("func_bigint_sub",f1,f2); }
  a=func_bigint_cast(FR(f1));
  b=func_bigint_cast(FR(f2));
  if(!func_is_bigint(a) || !func_is_bigint(b)){ FUNC_ERROR_ARG2("func_bigint_sub",f1,f2); }
  f=func_bigint();
  bigint_sub(func_bigint_p(f),func_bigint_p(a),func_bigint_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_bigint_eval(f);
}

// f=f1*f2
func_t *func_bigint_mul(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_bigint(f1) || !func_in_bigint(f2)){ FUNC_ERROR_ARG2("func_bigint_mul",f1,f2); }
  a=func_bigint_cast(FR(f1));
  b=func_bigint_cast(FR(f2));
  if(!func_is_bigint(a) || !func_is_bigint(b)){ FUNC_ERROR_ARG2("func_bigint_mul",f1,f2); }
  f=func_bigint();
  bigint_mul(func_bigint_p(f),func_bigint_p(a),func_bigint_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_bigint_eval(f);
}

// f=f1/f2
func_t *func_bigint_div(func_t *f1, func_t *f2)
{
  func_t *f=NULL,*a=NULL,*b=NULL;
  if(!func_in_bigint(f1) || !func_in_bigint(f2)){ FUNC_ERROR_ARG2("func_bigint_div",f1,f2); }
  a=func_bigint_cast(FR(f1));
  b=func_bigint_cast(FR(f2));
  if(!func_is_bigint(a) || !func_is_bigint(b)){ FUNC_ERROR_ARG2("func_bigint_div",f1,f2); }
  f=func_bigint();
  bigint_div(func_bigint_p(f),func_bigint_p(a),func_bigint_p(b));
  f1=func_del(f1);
  f2=func_del(f2);
  a=func_del(a);
  b=func_del(b);
  return func_bigint_eval(f);
}

// f=1/f1
func_t *func_bigint_inv(func_t *f1)
{
  func_t *f=NULL,*a=NULL;
  if(!func_in_bigint(f1)){ FUNC_ERROR_ARG1("func_bigint_inv",f1); }
  a=func_bigint_cast(FR(f1));
  if(!func_is_bigint(a)){ FUNC_ERROR_ARG1("func_bigint_inv",f1); }
  f=func_bigint();
  bigint_inv(func_bigint_p(f),func_bigint_p(a));
  f1=func_del(f1);
  a=func_del(a);
  return func_bigint_eval(f);
}

// g=f^p
func_t *func_bigint_pow_n(func_t *f, int p)
{
  func_t *g=NULL,*a=NULL;
  if(!func_in_bigint(f)){ FUNC_ERROR_ARG1("func_bigint_pow_n",f); }
  if(p==0){ f=func_del(f); return func_one(); }
  if(p==1){ return f; }
  a=func_bigint_cast(FR(f));
  g=func_bigint();  
  bigint_pow_n(func_bigint_p(g),func_bigint_p(a),p);
  f=func_del(f);
  a=func_del(a);
  return func_bigint_eval(g);
}

////////////////////////////////////////////////

func_t *func_op_bigint_new(void)
{
  func_t *f=func_builtin_new(__func_bigint);
  func_builtin_p(f)->order=FUNC_ORDER_BIGINT;
  func_builtin_p(f)->type=FUNC_SCALAR;
  func_builtin_p(f)->ptype=FUNC_P_BIGINT;
  func_builtin_p(f)->p_new=func_bigint_p_new;
  func_builtin_p(f)->p_del=func_bigint_p_del;
  func_builtin_p(f)->p_clone=func_bigint_p_clone;
  func_builtin_p(f)->p_cmp=func_bigint_cmp;
  func_builtin_p(f)->eval=func_bigint_eval;
  func_builtin_p(f)->print=func_bigint_print;
  return f;
}

//EOF
