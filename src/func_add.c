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
#define FGET_SP(f) (func_get_split_part(f))
#define FGET_NB(f) (func_get_number_part(f))

//////////////////////////////////////////////////////////

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"

func_t *func_add_script(func_t *f)
{
  int i;
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"%+",BLK0,BLK1)){
      a=strings_split(func_strings_at(f,0),"+",BLK0,BLK1,SPC);
      if(strings_size(a)>0){
	g=func_add_new(strings_size(a));
	for(i=0; i<func_asize(g); i++){
	  func_aset(g,i,func_script(strings_at(a,i)));
	}
	g=func_add_eval(g);
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}


/////////////////////////////////////

static char *__func_add="add";

/////////////////////////////////////

func_t *func_op_add_new(void)
{
  func_t *f=func_builtin_new(__func_add);
  func_builtin_p(f)->order=FUNC_ORDER_ADD;
  func_builtin_p(f)->amax=FUNC_NOLIMIT;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_add_eval;
  func_builtin_p(f)->print=func_add_print;
  return f;
}

/////////////////////////////////////

void func_add_print(func_t *f)
{
  int i;
  if(func_power(f)!=1) printf("(");
  for(i=0; i<func_asize(f); i++){
    if(func_is_bigint(func_aget(f,i)) && func_bigint_sgn(func_aget(f,i))<0){ } // do nothing
    else if(func_is_mul(func_aget(f,i)) && func_asize(func_aget(f,i))>0 && func_is_bigint(func_aget(f,i)->a[0]) && func_bigint_sgn(func_aget(f,i)->a[0])<0){ } // do nothing
    else if(i>=1){ printf("+"); }
    if(func_is_add(func_aget(f,i)) && func_power(func_aget(f,i))==1){ printf("("); }
    func_print(func_aget(f,i));
    if(func_is_add(func_aget(f,i)) && func_power(func_aget(f,i))==1){ printf(")"); }
  }
  if(func_power(f)!=1) printf(")");
  func_print_hat(func_power(f));
}


/////////////////////////////////////

int func_is_add(func_t *f)
{
  return func_is(f,__func_add);
}

////////////////////////////////////////////////////////////

void func_add_args_collect(func_t *f, func_is_t *fis)
{
  int i,j;
  func_t *g=NULL;
  if(!func_is_add(f)){ FUNC_ERROR_ARG1("func_add_args_collect",f); }
  // convert args to split list
  for(i=0; i<func_asize(f); i++){
    g=func_mul_split_list(FR(func_aget(f,i)),fis);
    func_aset(f,i,FR(g));
    g=func_del(g);
  }
  // add a arg to same arg
  for(i=0; i<func_asize(f); i++){
    for(j=i+1; func_aget(f,i)!=NULL && !func_is_one(func_aget(func_aget(f,i),1)) && j<func_asize(f); j++){
      if(func_aget(f,j)!=NULL && func_cmp(func_aget(func_aget(f,i),1),func_aget(func_aget(f,j),1))==0){
	g=func_list(2);
	func_aset(g,0,func_add(FR(func_aget(func_aget(f,i),0)),FR(func_aget(func_aget(f,j),0))));
	func_aset(g,1,FR(func_aget(func_aget(f,i),1)));
	func_aset(f,i,FR(g));
	func_adel(f,j);
	g=func_del(g);
      }
    }
  }
  // convert list to mul
  for(i=0; i<func_asize(f); i++){
    if(func_aget(f,i)!=NULL){
      g=func_mul(FR(func_aget(func_aget(f,i),0)),FR(func_aget(func_aget(f,i),1)));
      func_aset(f,i,FR(g));
      g=func_del(g);
    }
  }
  func_a_rm_null(f);
}

////////////////////////////////////////////////////////////////////////////

void func_add_args(func_t *f, func_is_t *fis, func_is_t *fin)
{
  int i,j;
  if(!func_is_add(f)){ FUNC_ERROR_ARG1("func_add_args",f); }
  for(i=0; i<func_asize(f); i++){
    for(j=0; fis(func_aget(f,i)) && j<func_asize(f); j++){
      if(i!=j && fin(func_aget(f,j))){
        func_aset(f,i,func_add(FR(func_aget(f,i)),FR(func_aget(f,j))));
	func_adel(f,j);
      }
    }
  }
  func_a_rm_null(f);
}

///private//////////////////////////////////////////////////

func_t *func_add_eval(func_t *f)
{
  func_t *g;
  if(!func_is_add(f)){ FUNC_ERROR_ARG1("func_add_eval",f); }
  f=func_flatten(f,__func_add);
  func_add_args(f,func_is_complex,func_is_number);
  func_add_args(f,func_is_real,func_is_number);
  func_add_args(f,func_in_bigint,func_in_bigint);
  func_add_args_collect(f,func_in_complex);
  func_add_args_collect(f,func_is_number);
  func_a_rm_op(f,func_is_zero);
  func_args_sort(f);
  func_args_reverse(f);
  if     (func_asize(f)<=0){ f=func_del(f); return func_zero(); }
  else if(func_asize(f)==1){ g=func_pow(FR(func_aget(f,0)),func_bigint_int(func_power(f),1)); f=func_del(f); return g; }
  else                     { return f; }
}

/////////////////////////////////////////////////////

func_t *func_add(func_t *f1, func_t *f2)
{
  if     (f1==NULL || f2==NULL)                      { FUNC_ERROR_ARG2("func_add",f1,f2); }
  else if(func_is_zero(f1))                          { f1=func_del(f1); return f2; }
  else if(func_is_zero(f2))                          { f2=func_del(f2); return f1; }
  else if(func_in_bigint(f1)  && func_in_bigint(f2)) { return func_bigint_add(f1,f2); }
  else if(func_in_real(f1)    && func_in_real(f2))   { return func_real_add(f1,f2); }
  else if(func_in_complex(f1) && func_in_complex(f2)){ return func_complex_add(f1,f2); }
  else if(func_is_real(f1)    && func_is_number(f2)) { return func_add(f1,func_evalf(f2)); }
  else if(func_is_complex(f1) && func_is_number(f2)) { return func_add(f1,func_evalf(f2)); }
  else if(func_is_number(f1)  && func_is_real(f2))   { return func_add(func_evalf(f1),f2); }
  else if(func_is_number(f1)  && func_is_complex(f2)){ return func_add(func_evalf(f1),f2); }
  else if(func_in_vec(f1)     && func_in_vec(f2))    { return func_vec_add(f1,f2); }
  else if(func_in_mat(f1)     && func_in_mat(f2))    { return func_mat_add(f1,f2); }
  else                                               { return func_add_eval(func_arg2_new(__func_add,f1,f2)); }
}

/////////////////////////////////////////////////////

func_t *func_add_new(int n)
{
  func_t *f=NULL;
  f=func_new(__func_add);
  func_a_resize(f,n);
  return f;
}


//EOF
