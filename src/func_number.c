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

#define FR(F) func_retain(F)

////////////////////////////////////////////////

static const char *__func_nan="nan";

int func_is_nan(func_t *f)
{
  return func_is(f,__func_nan);
}

func_t *func_nan()
{
  return func_new(__func_nan);
}

func_t *func_op_nan_new(void)
{
  func_t *f=func_builtin_new(__func_nan);
  func_builtin_p(f)->order=FUNC_ORDER_NAN;
  func_builtin_p(f)->print=func_nan_print;
  return f;
}

void func_nan_print(func_t *f)
{
  if(!func_is_nan(f)){ FUNC_ERROR_ARG1("func_nan_print",f); }
  printf("%s",__func_nan);
}

func_t *func_nan_script(func_t *f)
{
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1 && str_match(func_strings_at(f,0),"^nan$",NULL,NULL)){
    g=func_nan();
    f=func_del(f);
  }else{ g=f; }
  return g;
}

////////////////////////////////////////////////

static const char *__func_inf="inf";

int func_is_inf(func_t *f)
{
  return func_is(f,__func_inf);
}

func_t *func_inf()
{
  return func_new(__func_inf);
}

func_t *func_op_inf_new(void)
{
  func_t *f=func_builtin_new(__func_inf);
  func_builtin_p(f)->order=FUNC_ORDER_INF;
  func_builtin_p(f)->print=func_inf_print;
  return f;
}

void func_inf_print(func_t *f)
{
  if(!func_is_inf(f)){ FUNC_ERROR_ARG1("func_inf_print",f); }
  printf("%s",__func_inf);
}

func_t *func_inf_script(func_t *f)
{
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1 && str_match(func_strings_at(f,0),"^inf$",NULL,NULL)){
    g=func_inf();
    f=func_del(f);
  }else{ g=f; }
  return g;
}

////////////////////////////////////////////////

static const char *__func_zero="0";

func_t *func_zero(void)
{
  return func_new(__func_zero);
}

int func_is_zero(func_t *f)
{
  return func_is(f,__func_zero);
}

func_t *func_op_zero_new(void){
  func_t *f=func_builtin_new(__func_zero);
  func_builtin_p(f)->order=FUNC_ORDER_ZERO;
  func_builtin_p(f)->type=FUNC_SCALAR;
  func_builtin_p(f)->print=func_zero_print;
  return f;
}

void func_zero_print(func_t *f)
{
  if(!func_is_zero(f)){ FUNC_ERROR_ARG1("func_zero_print",f); }
  printf("%s",__func_zero);
}

func_t *func_zero_script(func_t *f)
{
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1 && str_match(func_strings_at(f,0),"^0$",NULL,NULL)){
    g=func_zero();
    f=func_del(f);
  }else{ g=f; }
  return g;
}

////////////////////////////////////////////////

static const char *__func_one="1";

func_t *func_one(void)
{
  return func_new(__func_one);
}

int func_is_one(func_t *f)
{
  return func_is(f,__func_one);
}

func_t *func_op_one_new(void){
  func_t *f=func_builtin_new(__func_one);
  func_builtin_p(f)->type=FUNC_SCALAR;
  func_builtin_p(f)->order=FUNC_ORDER_ONE;
  func_builtin_p(f)->print=func_one_print;
  return f;
}

void func_one_print(func_t *f)
{
  if(!func_is_one(f)){ FUNC_ERROR_ARG1("func_one_print",f); }
  printf("%s",__func_one);
}

func_t *func_one_script(func_t *f)
{
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1 && str_match(func_strings_at(f,0),"^1$",NULL,NULL)){
    g=func_one();
    f=func_del(f);
  }else{ g=f; }
  return g;
}

//////////////////////////////////////////////

func_t *func_evalf_eval(func_t *g)
{
  func_t *f=NULL;
  if(!func_is(g,"evalf") || func_asize(g)!=1){ FUNC_ERROR_ARG1("func_evalf_eval",g); }
  f=func_evalf(FR(func_aget(g,0)));
  g=func_del(g);
  return f;
}

//////////////////////////////////////////////

func_t *func_evalf(func_t *f)
{
  int i;
  func_t *g=NULL;
  func_arg1_t *eval=NULL;
  if(f==NULL){ FUNC_ERROR_ARG1("func_evalf",f); }
  if     (func_is_real(f))   { g=FR(f); }
  else if(func_is_complex(f)){ g=FR(f); }
  else if(func_is_var(f))    { g=FR(f); }
  else if(func_is_rvec(f))   { g=FR(f); }
  else if(func_is_cvec(f))   { g=FR(f); }
  else if(func_is_rmat(f))   { g=FR(f); }
  else if(func_is_cmat(f))   { g=FR(f); }
  if(g==NULL){ g=func_real_cast(FR(f)); }
  if(g==NULL){ g=func_complex_cast(FR(f)); }
  if(g==NULL){
    g=func_clone(FR(f));
    for(i=0; i<func_asize(f); i++){ func_aset(g,i,func_evalf(FR(f->a[i]))); }
    eval=func_find_eval(g);
    if(eval!=NULL){ g=eval(g); }
  }
  if     (func_is_zero(g)){ g=func_real_cast(g); }
  else if(func_is_one(g)) { g=func_real_cast(g); }
  f=func_del(f);
  return g;
}

int func_is_number(func_t *f)
{
  int i;
  if(f==NULL)                          { return 0; }
  if(func_in_complex(f))               { return 1; }
  if(func_find_type(f)==FUNC_SCALAR)   { return 1; }
  if(func_find_type(f)==FUNC_SCALAR_NO){ return 0; }
  if(func_find_type(f)==FUNC_COMMAND)  { return 0; }
  if(func_find_type(f)==FUNC_SCALAR_DEPEND){
    for(i=0; i<func_asize(f); i++){
      if(!func_is_number(f->a[i])){ return 0; }
    }
    return 1;
  }
  return 0;
}

int func_number_cmp(func_t *f, func_t *g)
{
  int value;
  func_t *f1=NULL,*g1=NULL;
  if     (func_in_bigint(f)  && func_in_bigint(g)) { return func_bigint_cmp(f,g); }
  else if(func_in_real(f)    && func_in_real(g))   { return func_real_cmp(f,g); }
  else if(func_in_complex(f) && func_in_complex(g)){ return func_complex_cmp(f,g); }
  else if(func_is_number(f) && func_is_number(g)){
    f1=func_evalf(FR(f));
    g1=func_evalf(FR(g));
    if     (func_is_real(f1)    && func_is_real(g1))    { value=rcmp   (func_real_p(f1),   func_real_p(g1)); }
    else if(func_is_real(f1)    && func_is_complex(g1)) { value=ccmp_r1(func_real_p(f1),   func_complex_p(g1)); }
    else if(func_is_complex(f1) && func_is_real(g1))    { value=ccmp_r2(func_complex_p(f1),func_real_p(g1)); }
    else if(func_is_complex(f1) && func_is_complex(g1)) { value=ccmp   (func_complex_p(f1),func_complex_p(g1)); }
    else{ FUNC_ERROR_ARG2("func_complex_cmp",f1,g1); }
    f1=func_del(f1);
    g1=func_del(g1);
    return value;
  }
  else{ FUNC_ERROR_ARG2("func_number_cmp",f,g); }
}

////////////////////////////////////////////////////////////////

//定数項aを引き出す x^2+x+a
func_t *func_number_pull(func_t *f)
{
  int i;
  func_t *a=func_zero();
  if(func_is_poly(f) && !func_is_mono(f)){
    for(i=0; i<func_asize(f); i++){
      if(func_is_number(f->a[i])) a=func_add(a,FR(f->a[i]));
    }
  }
  f=func_del(f);
  return a; 
}


//EOF
